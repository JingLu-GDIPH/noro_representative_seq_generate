#!/usr/bin/env python3

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


VALID_BASES = frozenset("ACGT")


def trim_alignment_ends(input_fasta, output_fasta, min_coverage=0.5):
    with open(input_fasta) as input_handle:
        records = list(SeqIO.parse(input_handle, "fasta"))
    if not records:
        raise ValueError("Alignment contains no sequences")
    if not 0 < min_coverage <= 1:
        raise ValueError("Minimum coverage must be greater than 0 and at most 1")

    lengths = {len(record.seq) for record in records}
    if len(lengths) != 1:
        raise ValueError("Input sequences do not have equal alignment lengths")

    alignment_length = lengths.pop()
    if alignment_length == 0:
        raise ValueError("Alignment contains no columns")

    coverage = []
    for column in range(alignment_length):
        valid_count = sum(
            str(record.seq[column]).upper() in VALID_BASES for record in records
        )
        coverage.append(valid_count / len(records))

    retained = [
        column
        for column, column_coverage in enumerate(coverage)
        if column_coverage >= min_coverage
    ]
    if not retained:
        raise ValueError(
            f"No alignment column reaches {min_coverage:.1%} A/C/G/T coverage"
        )

    start = retained[0]
    end = retained[-1] + 1
    for record in records:
        record.seq = Seq(str(record.seq[start:end]))

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    return {
        "sequence_count": len(records),
        "original_length": alignment_length,
        "trimmed_length": end - start,
        "left_trimmed": start,
        "right_trimmed": alignment_length - end,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Trim low-coverage columns from both ends of a FASTA alignment. "
            "Only A/C/G/T count as covered; N and gaps do not."
        )
    )
    parser.add_argument("--input", required=True, help="Input aligned FASTA")
    parser.add_argument("--output", required=True, help="Output trimmed FASTA")
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.5,
        help="Minimum A/C/G/T sequence coverage retained at each edge (default: 0.5)",
    )
    args = parser.parse_args()

    stats = trim_alignment_ends(
        Path(args.input), Path(args.output), args.min_coverage
    )
    print(
        "End trimming completed: "
        f"{stats['sequence_count']} sequences, "
        f"{stats['original_length']} -> {stats['trimmed_length']} columns, "
        f"left removed {stats['left_trimmed']}, "
        f"right removed {stats['right_trimmed']}"
    )


if __name__ == "__main__":
    main()
