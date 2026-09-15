#!/usr/bin/env python3
"""Rename genotyped consensus FASTA headers to RdRp_VP1_cluster_node format."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def safe_token(value: str) -> str:
    value = (value or "Unknown").strip()
    value = value.replace("/", "-")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", value)
    value = re.sub(r"-+", "-", value).strip("-")
    return value or "Unknown"


def load_cluster_lookup(cluster_info: Path) -> dict[str, str]:
    lookup = {}
    with cluster_info.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            sequence_id = row.get("sequence_id", "").strip()
            cluster_id = row.get("cluster_id", "").strip()
            match = re.search(r"cluster[_-]?(\d+)", cluster_id)
            if sequence_id and match:
                lookup[sequence_id] = match.group(1)
    return lookup


def load_typing_table(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["sequence_name"]] = row
    return rows


def extract_cluster_node(sequence_name: str, cluster_lookup: dict[str, str]) -> tuple[str, str, str]:
    cluster_match = re.search(r"cluster[_-](\d+)", sequence_name)
    node_match = re.search(r"node[_-](\d+)", sequence_name)

    cluster_no = cluster_match.group(1) if cluster_match else cluster_lookup.get(sequence_name, "NA")
    node_no = node_match.group(1) if node_match else "NA"

    source = []
    source.append("name" if cluster_match else ("cluster_info" if cluster_no != "NA" else "missing"))
    source.append("name" if node_match else "missing")
    return cluster_no, node_no, "/".join(source)


def extract_accession(sequence_name: str) -> str | None:
    parts = sequence_name.split("_")
    if len(parts) < 2:
        return None
    if len(parts) >= 3 and parts[1] in {"NC", "NG", "NM", "NR", "NW", "NZ", "XM", "XR"} and parts[2].isdigit():
        return f"{parts[1]}_{parts[2]}"
    if re.match(r"^[A-Z]{1,4}\d{5,9}(?:\.\d+)?$", parts[1]):
        return parts[1]
    return None


def is_raw_accession_record(sequence_name: str) -> bool:
    if re.search(r"cluster[_-]\d+", sequence_name):
        return False
    if "consensus" in sequence_name.lower() or "enhanced" in sequence_name.lower():
        return False
    return extract_accession(sequence_name) is not None


def build_new_name(sequence_name: str, row: dict[str, str], cluster_lookup: dict[str, str]) -> dict[str, str]:
    rdrp = safe_token(row.get("rdrp_genotype", "Unknown"))
    vp1 = safe_token(row.get("vp1_genotype", "Unknown"))

    if is_raw_accession_record(sequence_name):
        accession = safe_token(extract_accession(sequence_name) or "accessionNA")
        return {
            "core_name": f"{rdrp}_{vp1}_{accession}",
            "rdrp_genotype": rdrp,
            "vp1_genotype": vp1,
            "accession_no": accession,
            "cluster_no": "",
            "node_no": "",
            "naming_mode": "raw_accession",
            "parse_source": "sequence_name/accession",
        }

    cluster_no, node_no, parse_source = extract_cluster_node(sequence_name, cluster_lookup)
    return {
        "core_name": f"{rdrp}_{vp1}_cluster{cluster_no}_node{node_no}",
        "rdrp_genotype": rdrp,
        "vp1_genotype": vp1,
        "accession_no": "",
        "cluster_no": cluster_no,
        "node_no": node_no,
        "naming_mode": "consensus_cluster_node",
        "parse_source": parse_source,
    }


def rename_fasta(fasta: Path, typing_tsv: Path, cluster_info: Path, output: Path, mapping: Path) -> None:
    typing = load_typing_table(typing_tsv)
    cluster_lookup = load_cluster_lookup(cluster_info)
    used = Counter()
    renamed_records = []
    mapping_rows = []

    for record in SeqIO.parse(fasta, "fasta"):
        original_id = record.id
        if original_id not in typing:
            raise ValueError(f"Missing typing row for FASTA record: {original_id}")
        row = typing[original_id]
        name_parts = build_new_name(original_id, row, cluster_lookup)
        core_name = name_parts["core_name"]
        used[core_name] += 1
        new_id = core_name if used[core_name] == 1 else f"{core_name}_dup{used[core_name]}"

        new_record = SeqRecord(record.seq, id=new_id, name=new_id, description=new_id)
        renamed_records.append(new_record)
        mapping_rows.append(
            {
                "old_id": original_id,
                "new_id": new_id,
                "core_name": core_name,
                "rdrp_genotype": name_parts["rdrp_genotype"],
                "vp1_genotype": name_parts["vp1_genotype"],
                "accession_no": name_parts["accession_no"],
                "cluster_no": name_parts["cluster_no"],
                "node_no": name_parts["node_no"],
                "naming_mode": name_parts["naming_mode"],
                "parse_source": name_parts["parse_source"],
                "duplicate_index": str(used[core_name]),
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(renamed_records, output, "fasta")
    with mapping.open("w", newline="") as handle:
        fieldnames = [
            "old_id",
            "new_id",
            "core_name",
            "rdrp_genotype",
            "vp1_genotype",
            "accession_no",
            "cluster_no",
            "node_no",
            "naming_mode",
            "parse_source",
            "duplicate_index",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(mapping_rows)


def maybe_backup(path: Path, backup: Path) -> None:
    if not backup.exists():
        shutil.copy2(path, backup)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-place", action="store_true", help="Overwrite the input FASTA after writing a backup.")
    parser.add_argument(
        "--backup-suffix",
        default=".before_rdrp_vp1_cluster_node_names.fasta",
        help="Suffix for the one-time backup when --in-place is used.",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        choices=("gi", "gii"),
        default=("gi", "gii"),
        help="Datasets to rename.",
    )
    args = parser.parse_args()

    config = {
        "gi": {
            "fasta": Path("results_new/genotyping_consensus/gi/gi_final_consensus_genotyped.fasta"),
            "typing": Path("results_new/genotyping_consensus/gi/gi_final_consensus_genotyping.tsv"),
            "cluster": Path("results_new/gi/02_cluster/cluster_info.csv"),
        },
        "gii": {
            "fasta": Path("results_new/genotyping_consensus/gii/gii_final_consensus_genotyped.fasta"),
            "typing": Path("results_new/genotyping_consensus/gii/gii_final_consensus_genotyping.tsv"),
            "cluster": Path("results_new/gii/02_cluster/cluster_info.csv"),
        },
    }

    for dataset in args.datasets:
        paths = config[dataset]
        fasta = paths["fasta"]
        output = fasta if args.in_place else fasta.with_suffix(".rdrp_vp1_cluster_node.fasta")
        tmp_output = fasta.with_suffix(".tmp.rdrp_vp1_cluster_node.fasta") if args.in_place else output
        mapping = fasta.with_suffix(".rdrp_vp1_cluster_node_mapping.tsv")
        backup = Path(str(fasta) + args.backup_suffix)
        source_fasta = backup if args.in_place and backup.exists() else fasta
        if args.in_place:
            maybe_backup(fasta, backup)
        rename_fasta(source_fasta, paths["typing"], paths["cluster"], tmp_output, mapping)
        if args.in_place:
            tmp_output.replace(fasta)
        print(f"{dataset}: wrote {output}")
        print(f"{dataset}: mapping {mapping}")


if __name__ == "__main__":
    main()
