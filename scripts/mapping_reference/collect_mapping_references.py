#!/usr/bin/env python3
"""Collect per-genotype mapping references and enforce final sequence QC."""

import argparse
import csv
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-glob", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--max-n-gap-percent", type=float, default=10.0)
    args = parser.parse_args()

    files = sorted(Path(".").glob(args.input_glob))
    seen = Counter()
    kept = []
    rows = []
    for path in files:
        for record in SeqIO.parse(path, "fasta"):
            sequence = str(record.seq).upper()
            n_gap_percent = 100.0 * (sequence.count("N") + sequence.count("-")) / len(sequence)
            status = "kept" if n_gap_percent <= args.max_n_gap_percent else "removed_n_gap"
            rows.append(
                {
                    "sequence_id": record.id,
                    "source_file": path.name,
                    "length": len(sequence.replace("-", "")),
                    "n_gap_percent": f"{n_gap_percent:.6f}",
                    "status": status,
                }
            )
            if status != "kept":
                continue
            seen[record.id] += 1
            if seen[record.id] > 1:
                record.id = f"{record.id}_dup{seen[record.id]}"
                record.name = record.id
                record.description = ""
            record.seq = Seq(sequence.replace("-", ""))
            kept.append(record)

    SeqIO.write(kept, args.output, "fasta")
    with args.report.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sequence_id", "source_file", "length", "n_gap_percent", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Collected {len(files)} genotype files; retained {len(kept)} references")


if __name__ == "__main__":
    main()
