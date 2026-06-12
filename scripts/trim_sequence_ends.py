#!/usr/bin/env python3

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


MISSING_END_CHARACTERS = "-Nn"


def trim_sequence_ends(input_fasta, output_fasta):
    with open(input_fasta) as input_handle:
        records = list(SeqIO.parse(input_handle, "fasta"))
    if not records:
        raise ValueError("FASTA file contains no sequences")

    trimmed_records = 0
    total_removed = 0
    for record in records:
        original = str(record.seq)
        trimmed = original.strip(MISSING_END_CHARACTERS)
        if not trimmed:
            raise ValueError(
                f"Sequence {record.id!r} contains only terminal missing characters"
            )
        if len(trimmed) != len(original):
            trimmed_records += 1
            total_removed += len(original) - len(trimmed)
        record.seq = Seq(trimmed)

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    return {
        "sequence_count": len(records),
        "trimmed_records": trimmed_records,
        "total_removed": total_removed,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Remove terminal N and gap characters from each FASTA sequence"
    )
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA")
    args = parser.parse_args()

    stats = trim_sequence_ends(Path(args.input), Path(args.output))
    print(
        "Sequence-end trimming completed: "
        f"{stats['trimmed_records']}/{stats['sequence_count']} records trimmed, "
        f"{stats['total_removed']} terminal characters removed"
    )


if __name__ == "__main__":
    main()
