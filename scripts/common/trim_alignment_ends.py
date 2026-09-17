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

    for record in records:
        sequence = str(record.seq)
        record.seq = Seq("".join(sequence[column] for column in retained))

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    retained_set = set(retained)
    left_trimmed = 0
    while left_trimmed < alignment_length and left_trimmed not in retained_set:
        left_trimmed += 1

    right_trimmed = 0
    while (
        right_trimmed < alignment_length
        and alignment_length - right_trimmed - 1 not in retained_set
    ):
        right_trimmed += 1

    return {
        "sequence_count": len(records),
        "original_length": alignment_length,
        "trimmed_length": len(retained),
        "removed_columns": alignment_length - len(retained),
        "left_trimmed": left_trimmed,
        "right_trimmed": right_trimmed,
        "internal_removed_columns": (
            alignment_length - len(retained) - left_trimmed - right_trimmed
        ),
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Remove low-coverage columns from a FASTA alignment. "
            "Only A/C/G/T count as covered; N and gaps do not."
        )
    )
    parser.add_argument("--input", required=True, help="Input aligned FASTA")
    parser.add_argument("--output", required=True, help="Output trimmed FASTA")
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.5,
        help="Minimum A/C/G/T sequence coverage retained for each column (default: 0.5)",
    )
    args = parser.parse_args()

    stats = trim_alignment_ends(
        Path(args.input), Path(args.output), args.min_coverage
    )
    print(
        "Alignment column filtering completed: "
        f"{stats['sequence_count']} sequences, "
        f"{stats['original_length']} -> {stats['trimmed_length']} columns, "
        f"left removed {stats['left_trimmed']}, "
        f"right removed {stats['right_trimmed']}, "
        f"internal removed {stats['internal_removed_columns']}"
    )


if __name__ == "__main__":
    main()
