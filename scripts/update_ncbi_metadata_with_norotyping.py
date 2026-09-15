#!/usr/bin/env python3
"""Merge norotyping results into downloaded NCBI metadata and rename FASTA IDs."""

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


ACCESSION_RE = re.compile(r"(?:^|_)([A-Z]{1,4}_?\d{5,9})(?:\.\d+)?(?=_|$)")


NOROTYPING_COLUMNS = [
    "norotyping_status",
    "norotyping_genogroup",
    "norotyping_sequence_length",
    "rdrp_genotype",
    "rdrp_nearest_reference",
    "rdrp_tree_distance",
    "rdrp_blast_identity_pct",
    "rdrp_blast_alignment_length",
    "rdrp_reference_coverage_pct",
    "vp1_nearest_reference",
    "vp1_tree_distance",
    "vp1_blast_identity_pct",
    "vp1_blast_alignment_length",
    "vp1_reference_coverage_pct",
    "coarse_reference",
    "coarse_identity_pct",
    "coarse_alignment_length",
    "coarse_query_coverage_pct",
    "coarse_bitscore",
]


def safe_collection_date(value):
    value = re.sub(r"[^A-Za-z0-9]+", "-", value.strip())
    return value.strip("-") or "unknown-date"


def safe_genotype(value, fallback="Unknown"):
    value = (value or "").strip()
    if not value:
        value = fallback
    value = re.sub(r"[^A-Za-z0-9.]+", "-", value)
    return value.strip("-") or fallback


def accession_from_identifier(identifier):
    match = ACCESSION_RE.search(identifier)
    if not match:
        raise ValueError(f"Cannot extract accession from identifier: {identifier}")
    return match.group(1)


def build_seqname(row, typing_row):
    rdrp = safe_genotype(typing_row.get("rdrp_genotype"), "UnknownRdRp")
    vp1 = safe_genotype(typing_row.get("vp1_genotype"), "UnknownVP1")
    accession = row["accession"].strip()
    collection_date = safe_collection_date(row.get("collection_date", ""))
    return f"{rdrp}_{vp1}_{accession}_{collection_date}"


def read_table(path):
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_table(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_genotyping(path):
    rows, _ = read_table(path)
    by_accession = {}
    for row in rows:
        accession = accession_from_identifier(row["sequence_name"])
        if accession in by_accession:
            raise ValueError(f"Duplicate genotyping row for accession: {accession}")
        by_accession[accession] = row
    return by_accession


def merged_fieldnames(original_fieldnames):
    fields = list(original_fieldnames)
    if "seqname" not in fields:
        accession_index = fields.index("accession")
        fields.insert(accession_index + 1, "seqname")
    if "previous_seqname" not in fields:
        fields.insert(fields.index("seqname") + 1, "previous_seqname")
    if "previous_vp1_genotype" not in fields:
        insert_at = fields.index("vp1_genotype") + 1 if "vp1_genotype" in fields else len(fields)
        fields.insert(insert_at, "previous_vp1_genotype")
    for column in NOROTYPING_COLUMNS:
        if column not in fields:
            fields.append(column)
    return fields


def merge_metadata(metadata_path, genotyping_path):
    rows, fieldnames = read_table(metadata_path)
    typing_by_accession = load_genotyping(genotyping_path)
    output_fields = merged_fieldnames(fieldnames)
    updated_rows = []
    excluded_rows = []
    seqname_by_accession = {}

    for row in rows:
        accession = row["accession"].strip()
        typing_row = typing_by_accession.get(accession)
        if typing_row is None:
            row["norotyping_status"] = "absent_from_input_fasta"
            excluded_rows.append(row)
            continue

        old_seqname = row.get("seqname", "")
        old_vp1 = row.get("vp1_genotype", "")
        new_seqname = build_seqname(row, typing_row)

        row["previous_seqname"] = old_seqname
        row["previous_vp1_genotype"] = old_vp1
        row["seqname"] = new_seqname
        row["vp1_genotype"] = typing_row.get("vp1_genotype", "")
        row["rdrp_genotype"] = typing_row.get("rdrp_genotype", "")
        row["norotyping_status"] = typing_row.get("status", "")
        row["norotyping_genogroup"] = typing_row.get("genogroup", "")
        row["norotyping_sequence_length"] = typing_row.get("sequence_length", "")

        for column in NOROTYPING_COLUMNS:
            if column in {
                "norotyping_status",
                "norotyping_genogroup",
                "norotyping_sequence_length",
                "rdrp_genotype",
            }:
                continue
            row[column] = typing_row.get(column, "")

        seqname_by_accession[accession] = new_seqname
        updated_rows.append(row)

    counts = Counter(seqname_by_accession.values())
    duplicates = [name for name, count in counts.items() if count > 1]
    if duplicates:
        examples = ", ".join(duplicates[:10])
        raise ValueError(f"Duplicate sequence names generated: {examples}")

    return updated_rows, excluded_rows, output_fields, seqname_by_accession


def rewrite_fasta(input_fasta, output_fasta, seqname_by_accession):
    seen = set()
    output_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        accession = accession_from_identifier(record.id)
        if accession not in seqname_by_accession:
            raise ValueError(f"FASTA accession absent from metadata: {accession}")
        seen.add(accession)
        output_records.append(
            SeqRecord(record.seq, id=seqname_by_accession[accession], description="")
        )

    missing = set(seqname_by_accession) - seen
    if missing:
        examples = ", ".join(sorted(missing)[:10])
        raise ValueError(f"Metadata accessions absent from FASTA: {examples}")

    return SeqIO.write(output_records, output_fasta, "fasta")


def write_summary(path, rows):
    counts = Counter(
        (
            row.get("rdrp_genotype", "") or "Unknown",
            row.get("vp1_genotype", "") or "Unknown",
        )
        for row in rows
    )
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["rdrp_genotype", "vp1_genotype", "count"])
        for (rdrp, vp1), count in sorted(counts.items()):
            writer.writerow([rdrp, vp1, count])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--genotyping", required=True, type=Path)
    parser.add_argument("--output-metadata", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--excluded-metadata", type=Path)
    args = parser.parse_args()

    rows, excluded_rows, fieldnames, seqname_by_accession = merge_metadata(
        args.metadata, args.genotyping
    )
    args.output_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    write_table(args.output_metadata, rows, fieldnames)
    fasta_count = rewrite_fasta(args.fasta, args.output_fasta, seqname_by_accession)
    write_summary(args.summary, rows)
    if args.excluded_metadata:
        args.excluded_metadata.parent.mkdir(parents=True, exist_ok=True)
        write_table(args.excluded_metadata, excluded_rows, fieldnames)

    print(f"Updated metadata rows: {len(rows)}")
    print(f"Excluded metadata rows absent from input FASTA: {len(excluded_rows)}")
    print(f"Renamed FASTA records: {fasta_count}")
    print(f"Wrote metadata: {args.output_metadata}")
    print(f"Wrote FASTA: {args.output_fasta}")
    print(f"Wrote summary: {args.summary}")
    if args.excluded_metadata:
        print(f"Wrote excluded metadata: {args.excluded_metadata}")


if __name__ == "__main__":
    main()
