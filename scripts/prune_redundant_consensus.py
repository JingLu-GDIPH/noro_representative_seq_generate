#!/usr/bin/env python3
"""Remove redundant representative sequences above a similarity threshold."""

import argparse
import os
import subprocess
import tempfile

from Bio import SeqIO


def pairwise_similarity(seq1, seq2):
    """Match the similarity definition used by final validation."""
    if len(seq1) != len(seq2):
        raise ValueError("Aligned sequences must have equal lengths")

    matches = sum(
        1 for base1, base2 in zip(seq1, seq2)
        if base1 == base2 and base1 != "-"
    )
    total = sum(
        1 for base1, base2 in zip(seq1, seq2)
        if base1 != "-" or base2 != "-"
    )
    return (matches / total) * 100 if total else 0.0


def select_nonredundant_indices(sequences, threshold=95.0):
    """Retain records whose similarity to every retained record is <= threshold."""
    retained = []
    for index, sequence in enumerate(sequences):
        if all(
            pairwise_similarity(sequence, sequences[kept]) <= threshold
            for kept in retained
        ):
            retained.append(index)
    return retained


def align_records(records, threads=1):
    """Return aligned sequence strings in the same record order."""
    lengths = {len(record.seq) for record in records}
    if len(lengths) == 1:
        return [str(record.seq).upper() for record in records], "existing"

    with tempfile.TemporaryDirectory() as temp_dir:
        input_path = os.path.join(temp_dir, "input.fasta")
        output_path = os.path.join(temp_dir, "aligned.fasta")
        SeqIO.write(records, input_path, "fasta")

        with open(output_path, "w") as output_handle:
            result = subprocess.run(
                ["mafft", "--auto", "--thread", str(threads), input_path],
                stdout=output_handle,
                stderr=subprocess.PIPE,
                text=True,
            )

        if result.returncode != 0:
            raise RuntimeError(
                f"MAFFT failed (exit {result.returncode}): {result.stderr[:500]}"
            )

        aligned_records = list(SeqIO.parse(output_path, "fasta"))
        if len(aligned_records) != len(records):
            raise RuntimeError(
                f"MAFFT returned {len(aligned_records)} records, expected {len(records)}. "
                f"Stderr: {result.stderr[:300]}"
            )

        aligned_by_id = {r.id: str(r.seq).upper() for r in aligned_records}
        missing = [r.id for r in records if r.id not in aligned_by_id]
        if missing:
            raise RuntimeError(
                f"MAFFT dropped {len(missing)} sequence IDs: {missing[:5]}"
            )

        return [aligned_by_id[r.id] for r in records], "MAFFT"


def prune_records(records, threshold=95.0, threads=1):
    """Return retained original records and the alignment method."""
    if len(records) <= 1:
        return records, "not needed"

    aligned_sequences, alignment_method = align_records(records, threads)
    retained_indices = select_nonredundant_indices(
        aligned_sequences,
        threshold=threshold,
    )
    return [records[index] for index in retained_indices], alignment_method


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--threshold", type=float, default=95.0)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    records = list(SeqIO.parse(args.input, "fasta"))
    retained, alignment_method = prune_records(
        records,
        threshold=args.threshold,
        threads=args.threads,
    )
    SeqIO.write(retained, args.output, "fasta")

    print(f"Input sequences: {len(records)}")
    print(f"Retained sequences: {len(retained)}")
    print(f"Removed redundant sequences: {len(records) - len(retained)}")
    print(f"Alignment method: {alignment_method}")


if __name__ == "__main__":
    main()
