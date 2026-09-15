#!/usr/bin/env python3
"""Norotype pre-separated GI/GII FASTA files without genome-level preclassification."""

import argparse
import csv
import importlib.util
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_norotyping_module(path):
    spec = importlib.util.spec_from_file_location("genotype_norovirus", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import norotyping module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_region_blast(module, queries, references, reference_metadata, group, work_dir, threads):
    query_fasta = work_dir / "queries.fasta"
    SeqIO.write(queries, query_fasta, "fasta")

    region_records = [
        record
        for (ref_group, region), records in references.items()
        if ref_group == group and region in {"rdrp", "vp1"}
        for record in records
    ]
    if not region_records:
        raise ValueError(f"No {group} RdRp/VP1 reference records found")

    reference_fasta = work_dir / f"{group.lower()}_region_references.fasta"
    database = work_dir / f"{group.lower()}_region_refs"
    output = work_dir / f"{group.lower()}_region_hits.tsv"
    SeqIO.write(region_records, reference_fasta, "fasta")

    module.run(["makeblastdb", "-in", reference_fasta, "-dbtype", "nucl", "-out", database])
    module.run(
        [
            "blastn",
            "-task",
            "blastn",
            "-query",
            query_fasta,
            "-db",
            database,
            "-num_threads",
            threads,
            "-max_target_seqs",
            len(region_records),
            "-max_hsps",
            "1",
            "-evalue",
            "1e-10",
            "-outfmt",
            "6 qseqid sseqid pident length qstart qend sstart send sstrand bitscore evalue",
            "-out",
            output,
        ]
    )

    hits = []
    with output.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            hits.append(
                {
                    "query": row[0],
                    "reference": row[1],
                    "identity": float(row[2]),
                    "alignment_length": int(row[3]),
                    "qstart": int(row[4]),
                    "qend": int(row[5]),
                    "sstart": int(row[6]),
                    "send": int(row[7]),
                    "sstrand": row[8],
                    "bitscore": float(row[9]),
                    "evalue": float(row[10]),
                }
            )

    return module.best_hits_by_category(hits, reference_metadata)


def format_number(module, value, decimals=3):
    return module.format_number(value, decimals)


def write_outputs(module, input_path, output_prefix, group, queries, query_metadata, tree_results, reference_metadata):
    rows = []
    annotated_records = []

    for query in queries:
        query_id = query.id
        original = query_metadata[query_id]
        row = {field: "" for field in module.OUTPUT_FIELDS}
        row.update(
            {
                "input_file": input_path.name,
                "sequence_name": original["name"],
                "sequence_length": original["length"],
                "status": "typed",
                "genogroup": group,
            }
        )

        assigned_regions = []
        for region in ("rdrp", "vp1"):
            result = tree_results.get((query_id, region))
            if result is None:
                continue
            nearest = reference_metadata[result["reference_id"]]
            blast_hit = result["blast"]
            blast_ref = reference_metadata[blast_hit["reference"]]
            row.update(
                {
                    f"{region}_nearest_reference": nearest["name"],
                    f"{region}_genotype": nearest["genotype"],
                    f"{region}_tree_distance": format_number(module, result["distance"], 6),
                    f"{region}_blast_identity_pct": format_number(module, blast_hit["identity"]),
                    f"{region}_blast_alignment_length": blast_hit["alignment_length"],
                    f"{region}_reference_coverage_pct": format_number(
                        module,
                        100.0 * module.reference_span(blast_hit) / blast_ref["length"],
                    ),
                }
            )
            assigned_regions.append(f"{region.upper()}={nearest['genotype']}")

        if not assigned_regions:
            row["status"] = "untyped"
        elif len(assigned_regions) == 1:
            row["status"] = "partial_genotype"

        annotation = f"genogroup={row['genogroup']}"
        if assigned_regions:
            annotation += " " + " ".join(assigned_regions)
        annotated_records.append(
            SeqRecord(original["sequence"], id=original["name"], description=annotation)
        )
        rows.append(row)

    table_path = Path(f"{output_prefix}_genotyping.tsv")
    fasta_path = Path(f"{output_prefix}_genotyped.fasta")
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=module.OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    SeqIO.write(annotated_records, fasta_path, "fasta")

    return table_path, fasta_path, rows


def main():
    parser = argparse.ArgumentParser(
        description="Assign RdRp and VP1 genotypes for a pre-separated GI/GII FASTA."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--group", required=True, choices=("GI", "GII"))
    parser.add_argument("--ref-dir", required=True, type=Path)
    parser.add_argument("--norotyping-script", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--min-identity", type=float, default=70.0)
    parser.add_argument("--min-reference-coverage", type=float, default=20.0)
    args = parser.parse_args()

    module = load_norotyping_module(args.norotyping_script)
    module.require_tools("mafft")
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tree_dir = Path(f"{args.output_prefix}_trees")
    tree_dir.mkdir(parents=True, exist_ok=True)

    references, reference_metadata = module.load_references(args.ref_dir)
    references = {
        key: records
        for key, records in references.items()
        if key[0] == args.group and key[1] in {"rdrp", "vp1"}
    }
    reference_metadata = {
        ref_id: meta
        for ref_id, meta in reference_metadata.items()
        if meta["group"] == args.group and meta["region"] in {"rdrp", "vp1"}
    }
    queries, query_metadata = module.load_queries(args.input)

    with tempfile.TemporaryDirectory(prefix="noro_direct_typing_") as temp:
        work_dir = Path(temp)
        category_hits = run_region_blast(
            module,
            queries,
            references,
            reference_metadata,
            args.group,
            work_dir,
            str(args.threads),
        )

        assignments = {query_id: (args.group, None) for query_id in query_metadata}
        tree_results = {}
        for region in ("rdrp", "vp1"):
            regional_results = module.build_tree_and_assign(
                args.group,
                region,
                list(query_metadata),
                query_metadata,
                references,
                reference_metadata,
                category_hits,
                work_dir,
                tree_dir,
                str(args.threads),
                args.min_identity,
                args.min_reference_coverage,
            )
            for query_id, result in regional_results.items():
                tree_results[(query_id, region)] = result

    table_path, fasta_path, rows = write_outputs(
        module,
        args.input,
        args.output_prefix,
        args.group,
        queries,
        query_metadata,
        tree_results,
        reference_metadata,
    )
    statuses = {}
    for row in rows:
        statuses[row["status"]] = statuses.get(row["status"], 0) + 1
    status_summary = ", ".join(f"{key}={value}" for key, value in sorted(statuses.items()))
    print(f"Wrote {len(rows)} sequence records to {table_path}")
    print(f"Wrote annotated FASTA to {fasta_path}")
    print(f"Wrote phylogenetic trees to {tree_dir}")
    print(f"Status summary: {status_summary}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
