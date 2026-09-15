#!/usr/bin/env python3

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


MISSING_END_CHARACTERS = "-Nn"


def trim_sequence_ends(input_fasta, output_fasta, remove_gaps=False):
    with open(input_fasta) as input_handle:
        records = list(SeqIO.parse(input_handle, "fasta"))
    if not records:
        raise ValueError("FASTA file contains no sequences")

    trimmed_records = 0
    total_removed = 0
    total_internal_gaps_removed = 0
    for record in records:
        original = str(record.seq)
        trimmed = original.strip(MISSING_END_CHARACTERS)
        if remove_gaps:
            gap_count = trimmed.count("-")
            trimmed = trimmed.replace("-", "")
            total_internal_gaps_removed += gap_count
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
        "total_internal_gaps_removed": total_internal_gaps_removed,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Remove terminal N/gap characters from each FASTA sequence; optionally "
            "remove all alignment gaps."
        )
    )
    parser.add_argument("--input", required=True, help="Input FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA")
    parser.add_argument(
        "--remove-gaps",
        action="store_true",
        help="Remove all '-' characters before writing final representative sequences",
    )
    args = parser.parse_args()

    stats = trim_sequence_ends(Path(args.input), Path(args.output), args.remove_gaps)
    print(
        "Sequence-end trimming completed: "
        f"{stats['trimmed_records']}/{stats['sequence_count']} records trimmed, "
        f"{stats['total_removed']} total characters removed "
        f"({stats['total_internal_gaps_removed']} alignment gaps)"
    )


if __name__ == "__main__":
    main()
