#!/usr/bin/env python3
"""Assign norovirus VP1 genotypes with the IPHnano reference BLAST databases."""

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]

OUTPUT_COLUMNS = [
    "vp1_genotype",
    "vp1_assignment_status",
    "vp1_best_reference",
    "vp1_percent_identity",
    "vp1_alignment_length",
    "vp1_reference_length",
    "vp1_reference_coverage_percent",
    "vp1_query_start",
    "vp1_query_end",
    "vp1_subject_start",
    "vp1_subject_end",
    "vp1_evalue",
    "vp1_best_hsp_bitscore",
    "vp1_total_bitscore",
    "vp1_second_genotype",
    "vp1_second_genotype_bitscore",
    "vp1_bitscore_delta",
]


def genotype_from_reference(reference):
    match = re.search(r"_(GII?\.\d+)$", reference, flags=re.IGNORECASE)
    return match.group(1).upper() if match else "Unknown"


def accession_from_query(query):
    match = re.search(r"(?:^|_)([A-Z]{1,4}_?\d{5,9})(?:\.\d+)?$", query)
    return match.group(1) if match else query


def interval_union_length(intervals):
    merged = []
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def run_blast(query, database, output, threads, blastn):
    command = [
        blastn,
        "-query",
        str(query),
        "-db",
        str(database),
        "-evalue",
        "1e-20",
        "-max_target_seqs",
        "10",
        "-outfmt",
        "6 " + " ".join(BLAST_COLUMNS),
        "-num_threads",
        str(threads),
        "-out",
        str(output),
    ]
    subprocess.run(command, check=True)


def parse_blast(path):
    grouped = defaultdict(lambda: defaultdict(list))
    with open(path) as handle:
        reader = csv.DictReader(handle, fieldnames=BLAST_COLUMNS, delimiter="\t")
        for row in reader:
            for key in (
                "pident", "evalue", "bitscore",
            ):
                row[key] = float(row[key])
            for key in (
                "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "qlen", "slen",
            ):
                row[key] = int(row[key])
            grouped[row["qseqid"]][row["sseqid"]].append(row)
    return grouped


def summarize_subject(subject, hsps):
    total_alignment = sum(row["length"] for row in hsps)
    weighted_identity = (
        sum(row["pident"] * row["length"] for row in hsps) / total_alignment
    )
    subject_union = interval_union_length(
        (row["sstart"], row["send"]) for row in hsps
    )
    query_start = min(min(row["qstart"], row["qend"]) for row in hsps)
    query_end = max(max(row["qstart"], row["qend"]) for row in hsps)
    subject_start = min(min(row["sstart"], row["send"]) for row in hsps)
    subject_end = max(max(row["sstart"], row["send"]) for row in hsps)
    subject_length = hsps[0]["slen"]
    return {
        "reference": subject,
        "genotype": genotype_from_reference(subject),
        "percent_identity": weighted_identity,
        "alignment_length": total_alignment,
        "reference_length": subject_length,
        "reference_coverage_percent": 100 * subject_union / subject_length,
        "query_start": query_start,
        "query_end": query_end,
        "subject_start": subject_start,
        "subject_end": subject_end,
        "evalue": min(row["evalue"] for row in hsps),
        "best_hsp_bitscore": max(row["bitscore"] for row in hsps),
        "bitscore": sum(row["bitscore"] for row in hsps),
    }


def call_genotypes(blast_path):
    calls = {}
    for query, subjects in parse_blast(blast_path).items():
        summaries = [
            summarize_subject(subject, hsps)
            for subject, hsps in subjects.items()
        ]
        summaries.sort(
            key=lambda item: (
                item["best_hsp_bitscore"],
                item["reference_coverage_percent"],
                item["percent_identity"],
            ),
            reverse=True,
        )
        best = summaries[0]
        second = next(
            (
                item
                for item in summaries[1:]
                if item["genotype"] != best["genotype"]
            ),
            None,
        )
        calls[accession_from_query(query)] = {
            "vp1_genotype": best["genotype"],
            "vp1_assignment_status": "assigned_best_vp1_reference",
            "vp1_best_reference": best["reference"],
            "vp1_percent_identity": f"{best['percent_identity']:.3f}",
            "vp1_alignment_length": best["alignment_length"],
            "vp1_reference_length": best["reference_length"],
            "vp1_reference_coverage_percent": (
                f"{best['reference_coverage_percent']:.3f}"
            ),
            "vp1_query_start": best["query_start"],
            "vp1_query_end": best["query_end"],
            "vp1_subject_start": best["subject_start"],
            "vp1_subject_end": best["subject_end"],
            "vp1_evalue": f"{best['evalue']:.3g}",
            "vp1_best_hsp_bitscore": f"{best['best_hsp_bitscore']:.1f}",
            "vp1_total_bitscore": f"{best['bitscore']:.1f}",
            "vp1_second_genotype": second["genotype"] if second else "",
            "vp1_second_genotype_bitscore": (
                f"{second['best_hsp_bitscore']:.1f}" if second else ""
            ),
            "vp1_bitscore_delta": (
                f"{best['best_hsp_bitscore'] - second['best_hsp_bitscore']:.1f}"
                if second
                else ""
            ),
        }
    return calls


def update_metadata(metadata_path, calls, output_path):
    with open(metadata_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_columns = [
            column for column in reader.fieldnames if column not in OUTPUT_COLUMNS
        ]
        rows = list(reader)

    missing_calls = []
    for row in rows:
        accession = row["accession"]
        call = calls.get(accession)
        if call:
            row.update(call)
        else:
            row.update({column: "" for column in OUTPUT_COLUMNS})
            row["vp1_assignment_status"] = "no_significant_vp1_hit"
            missing_calls.append(accession)

    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=source_columns + OUTPUT_COLUMNS,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(rows), missing_calls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--blast-output", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--blastn", default="/usr/local/bin/blastn")
    parser.add_argument("--reuse-blast", action="store_true")
    args = parser.parse_args()

    blast_output = Path(args.blast_output)
    if not args.reuse_blast or not blast_output.exists():
        run_blast(
            args.query,
            args.database,
            blast_output,
            args.threads,
            args.blastn,
        )

    calls = call_genotypes(blast_output)
    row_count, missing = update_metadata(
        args.metadata,
        calls,
        args.output_metadata,
    )
    print(f"Metadata rows: {row_count}")
    print(f"VP1 genotype calls: {len(calls)}")
    print(f"No significant VP1 hit: {len(missing)}")
    if missing:
        print("Missing accessions: " + ",".join(missing[:20]))


if __name__ == "__main__":
    main()
