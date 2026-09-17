#!/usr/bin/env python3
"""Combine consensus FASTA files while keeping identifiers globally unique."""

import argparse
from collections import Counter
from pathlib import Path

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_glob", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()

    input_files = sorted(Path(".").glob(args.input_glob))
    seen = Counter()
    records = []

    for input_file in input_files:
        for record in SeqIO.parse(input_file, "fasta"):
            seen[record.id] += 1
            if seen[record.id] > 1:
                record.id = f"{record.id}_dup{seen[record.id]}"
                record.name = record.id
                record.description = ""
            records.append(record)

    SeqIO.write(records, args.output_file, "fasta")
    print(f"Combined {len(input_files)} files and {len(records)} sequences")


if __name__ == "__main__":
    main()
