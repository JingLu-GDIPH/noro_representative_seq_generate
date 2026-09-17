#!/usr/bin/env python3
"""Normalize every sequence in a cluster to a single strand orientation.

Sequences downloaded from public databases are not always stored on the same
strand. When sequences of opposite orientation are forced into one multiple
sequence alignment, the reversed record ends up gapped across almost the whole
alignment (e.g. ``GII.P31_GII.4_KX158285``) and produces meaningless pairwise
similarities. This script detects and reverses such records **before** MAFFT so
that every cluster is aligned in a consistent orientation.

Strategy
--------
For each input cluster:

1. Pick an orientation reference: the longest sequence (A/C/G/T only). The
   longest record is the most likely to span the full genome and to represent
   the strand the majority of the cluster is stored on.
2. For every other sequence, compare how many ``k``-mers it shares with the
   reference in the forward direction versus its reverse complement. A sequence
   whose reverse complement shares clearly more ``k``-mers with the reference is
    assumed to be on the opposite strand and is reverse-complemented in place.
3. Single-record clusters are passed through unchanged (nothing to orient
   against).

Only A/C/G/T characters are considered when building ``k``-mers; ``N`` and gap
characters are skipped, so ambiguous or padded records are handled safely.
``k`` is chosen long enough (default 15) that random shared ``k``-mers between
unrelated strands are negligible, yet short enough to tolerate sequencing
divergence within a genotype.
"""

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


#: Reverse-complement translation table restricted to IUPAC DNA characters.
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

#: ``k`` for the canonical k-mer orientation test. 15-mers are specific enough
#: that two unrelated viral strands essentially never share many by chance, and
#: short enough to tolerate intra-genotype divergence.
DEFAULT_KMER_SIZE = 15

#: Reverse a sequence only when its reverse complement shares at least this many
#: more reference k-mers than its forward strand. Guarded against tiny clusters
#: where a single shared k-mer would be noise.
DEFAULT_MIN_OVERLAP = 50

#: Bases that count as "real" sequence content when picking the orientation
#: reference. N and gap characters are excluded so that padded/ambiguous records
#: cannot win the "longest sequence" comparison and become a misleading anchor.
_VALID_BASES = frozenset("ACGT")


def reverse_complement(sequence):
    """Return the reverse complement of a DNA string, preserving case."""
    return sequence.translate(_COMPLEMENT)[::-1]


def _valid_base_count(sequence):
    """Count A/C/G/T characters (case-insensitive); N and gaps are excluded."""
    return sum(1 for base in sequence.upper() if base in _VALID_BASES)


def _clean(sequence):
    """Upper-case and strip gap characters so k-mers cover real bases only."""
    return sequence.upper().replace("-", "")


def build_kmer_set(sequence, kmer_size):
    """Return the set of ``kmer_size`` substrings of ``sequence``.

    ``N`` is allowed inside a k-mer (it is still a valid, if uninformative,
    substring); gap characters are removed first.
    """
    cleaned = _clean(sequence)
    if len(cleaned) < kmer_size:
        return set()
    return {cleaned[i : i + kmer_size] for i in range(len(cleaned) - kmer_size + 1)}


def _candidate_orientation(overlap_forward, overlap_reverse, min_overlap):
    """Return ``"reverse"`` if the record should be flipped, else ``"forward"``.

    A record is reversed only when the reverse-complement overlap strictly
    exceeds the forward overlap **and** that reverse overlap clears
    ``min_overlap``, which prevents flipping on a handful of coincidental
    shared k-mers.
    """
    if overlap_reverse > overlap_forward and overlap_reverse >= min_overlap:
        return "reverse"
    return "forward"


def detect_orientation(sequence, reference_kmers, kmer_size, min_overlap):
    """Decide whether ``sequence`` matches ``reference_kmers`` forward or reverse.

    Returns a tuple ``(orientation, overlap_forward, overlap_reverse)``.
    """
    forward_kmers = build_kmer_set(sequence, kmer_size)
    reverse_kmers = build_kmer_set(reverse_complement(sequence), kmer_size)
    overlap_forward = len(forward_kmers & reference_kmers)
    overlap_reverse = len(reverse_kmers & reference_kmers)
    orientation = _candidate_orientation(
        overlap_forward, overlap_reverse, min_overlap
    )
    return orientation, overlap_forward, overlap_reverse


def normalize_orientation(input_fasta, output_fasta, kmer_size=DEFAULT_KMER_SIZE,
                          min_overlap=DEFAULT_MIN_OVERLAP):
    """Write ``output_fasta`` with every record on a single strand.

    Returns a dict with orientation statistics and the list of records that were
    reverse-complemented, so callers can report which accessions were flipped.
    """
    with open(input_fasta) as input_handle:
        records = list(SeqIO.parse(input_handle, "fasta"))
    if not records:
        raise ValueError("FASTA file contains no sequences")

    flipped_ids = []

    # Single-record clusters have nothing to orient against; copy through.
    if len(records) == 1:
        SeqIO.write(records, output_fasta, "fasta")
        return {
            "sequence_count": 1,
            "flipped_count": 0,
            "flipped_ids": [],
            "reference_id": records[0].id,
        }

    # Orientation reference: the record with the most A/C/G/T content. Using real
    # base count (not raw length) means N- or gap-padded records cannot win the
    # comparison and become a misleading anchor. The longest real genome is the
    # most likely to span the full virus and represent the cluster's strand.
    reference = max(records, key=lambda record: _valid_base_count(str(record.seq)))
    reference_kmers = build_kmer_set(str(reference.seq), kmer_size)
    if not reference_kmers:
        # Reference shorter than k (extremely fragmented input): do not flip
        # anything, we cannot make a reliable decision.
        SeqIO.write(records, output_fasta, "fasta")
        return {
            "sequence_count": len(records),
            "flipped_count": 0,
            "flipped_ids": [],
            "reference_id": reference.id,
        }

    for record in records:
        if record.id == reference.id:
            continue
        orientation, _, _ = detect_orientation(
            str(record.seq), reference_kmers, kmer_size, min_overlap
        )
        if orientation == "reverse":
            record.seq = Seq(reverse_complement(str(record.seq)))
            flipped_ids.append(record.id)

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

    return {
        "sequence_count": len(records),
        "flipped_count": len(flipped_ids),
        "flipped_ids": flipped_ids,
        "reference_id": reference.id,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reverse-complement any sequence stored opposite to the cluster "
            "majority strand, so MAFFT receives a consistently oriented cluster."
        )
    )
    parser.add_argument("--input", required=True, help="Input per-cluster FASTA")
    parser.add_argument("--output", required=True, help="Output orientation-normalized FASTA")
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=DEFAULT_KMER_SIZE,
        help=(
            "k-mer length for the forward/reverse overlap test "
            f"(default: {DEFAULT_KMER_SIZE})"
        ),
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=DEFAULT_MIN_OVERLAP,
        help=(
            "minimum reverse-complement k-mer overlap required before a record "
            f"is flipped (default: {DEFAULT_MIN_OVERLAP})"
        ),
    )
    args = parser.parse_args()

    stats = normalize_orientation(
        Path(args.input), Path(args.output), args.kmer_size, args.min_overlap
    )
    flipped = ", ".join(stats["flipped_ids"]) if stats["flipped_ids"] else "none"
    print(
        "Orientation normalization completed: "
        f"{stats['flipped_count']}/{stats['sequence_count']} records "
        f"reverse-complemented against reference '{stats['reference_id']}'; "
        f"flipped: {flipped}"
    )


if __name__ == "__main__":
    main()
