#!/usr/bin/env python3
"""
Validate internal pairwise similarity of consensus sequences within a file.
If any pairwise similarity exceeds the threshold, the file needs reclustering.
"""

import argparse
import subprocess
import tempfile
import os
import sys
import numpy as np
import re
from collections import defaultdict
from Bio import SeqIO
from Bio import AlignIO

VALID_BASES = set("ACGT")
VALID_BASE_BYTES = np.array([np.bytes_(base) for base in "ACGT"])
GENOTYPE_PAIR_PATTERN = re.compile(
    r"(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)"
)


def calculate_pairwise_similarity(seq1, seq2):
    """Calculate pairwise similarity over columns where both bases are A/C/G/T."""
    if len(seq1) != len(seq2):
        return 0.0
    
    matches = 0
    total = 0
    for a, b in zip(seq1.upper(), seq2.upper()):
        if a not in VALID_BASES or b not in VALID_BASES:
            continue
        total += 1
        if a == b:
            matches += 1
    
    if total == 0:
        return 0.0
    
    return (matches / total) * 100


def calculate_max_pairwise_similarity(alignment, chunk_size=256):
    """Calculate maximum pairwise similarity with bounded-memory NumPy chunks."""
    sequence_strings = [str(record.seq).upper() for record in alignment]
    sequence_count = len(sequence_strings)
    if sequence_count <= 1:
        return 0.0, 0
    sequence_length = len(sequence_strings[0])
    encoded = np.frombuffer(
        "".join(sequence_strings).encode("ascii"), dtype="S1"
    ).reshape(sequence_count, sequence_length)

    max_similarity = 0.0
    for i in range(sequence_count - 1):
        reference = encoded[i]
        for start in range(i + 1, sequence_count, chunk_size):
            comparison = encoded[start:min(start + chunk_size, sequence_count)]
            valid_in_both = np.isin(comparison, VALID_BASE_BYTES) & np.isin(reference, VALID_BASE_BYTES)
            matches = ((comparison == reference) & valid_in_both).sum(axis=1)
            totals = valid_in_both.sum(axis=1)
            similarities = np.divide(
                matches,
                totals,
                out=np.zeros(matches.shape, dtype=float),
                where=totals != 0,
            ) * 100
            if similarities.size:
                max_similarity = max(max_similarity, float(similarities.max()))

    pair_count = sequence_count * (sequence_count - 1) // 2
    return max_similarity, pair_count


def extract_genotype_pair(sequence_id):
    """Extract RdRp_VP1 genotype pair from a representative FASTA identifier."""
    match = GENOTYPE_PAIR_PATTERN.search(str(sequence_id))
    if match:
        return f"{match.group(1)}_{match.group(2)}"
    return "Unknown"


def validate_consensus_internal(consensus_file, output_report, similarity_threshold, group_by="none"):
    """Validate internal similarity of consensus sequences."""
    
    # Read sequences
    sequences = list(SeqIO.parse(consensus_file, "fasta"))
    
    if len(sequences) <= 1:
        with open(output_report, 'w') as f:
            f.write(f"Max similarity: 0.0\n")
            f.write(f"Needs reclustering: false\n")
            f.write(f"Reason: Single sequence file\n")
        return
    
    temp_fasta_path = None
    aligned_file = None
    try:
        sequence_lengths = {len(record.seq) for record in sequences}
        if len(sequence_lengths) == 1:
            alignment = sequences
            alignment_method = "Reused existing equal-length alignment"
        else:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_fasta:
                SeqIO.write(sequences, temp_fasta.name, "fasta")
                temp_fasta_path = temp_fasta.name

            aligned_file = temp_fasta_path + "_aligned.fasta"
            with open(aligned_file, 'w') as aligned_handle:
                subprocess.run([
                    "mafft", "--quiet", "--auto", "--maxiterate", "1000", "--localpair",
                    temp_fasta_path
                ], stdout=aligned_handle, stderr=subprocess.PIPE, check=True, text=True)

            alignment = AlignIO.read(aligned_file, "fasta")
            alignment_method = "MAFFT"
        
        # Calculate pairwise similarities
        global_max_similarity, global_pair_count = calculate_max_pairwise_similarity(alignment)

        if group_by == "genotype_pair":
            groups = defaultdict(list)
            for record in alignment:
                groups[extract_genotype_pair(record.id)].append(record)

            max_similarity = 0.0
            pair_count = 0
            failing_groups = []
            group_summaries = []
            for group_name, group_records in sorted(groups.items()):
                group_max, group_pairs = calculate_max_pairwise_similarity(group_records)
                pair_count += group_pairs
                group_summaries.append(
                    (group_name, len(group_records), group_max, group_pairs)
                )
                if group_max > max_similarity:
                    max_similarity = group_max
                if group_max >= similarity_threshold:
                    failing_groups.append(
                        (group_name, len(group_records), group_max, group_pairs)
                    )
        else:
            max_similarity = global_max_similarity
            pair_count = global_pair_count
            failing_groups = []
            group_summaries = []
        
        # Determine if reclustering is needed
        needs_reclustering = max_similarity >= similarity_threshold
        reason = f"Max pairwise similarity ({max_similarity:.4f}%) is at or above threshold ({similarity_threshold}%)" if needs_reclustering else "All pairwise similarities below threshold"
        
        # Write report
        with open(output_report, 'w') as f:
            f.write(f"Max similarity: {max_similarity:.4f}\n")
            f.write(f"Global max similarity: {global_max_similarity:.4f}\n")
            f.write(f"Needs reclustering: {str(needs_reclustering).lower()}\n")
            f.write(f"Reason: {reason}\n")
            f.write(f"Total sequences: {len(sequences)}\n")
            f.write(f"Total pairwise comparisons: {pair_count}\n")
            f.write(f"Global pairwise comparisons: {global_pair_count}\n")
            f.write(f"Alignment method: {alignment_method}\n")
            f.write(f"Validation grouping: {group_by}\n")
            f.write("Similarity rule: count only columns where both bases are A/C/G/T; N and gap are excluded from coverage\n")
            if group_summaries:
                f.write("Group summaries: genotype_pair\tsequence_count\tmax_similarity\tpairwise_comparisons\n")
                for group_name, group_count, group_max, group_pairs in group_summaries:
                    f.write(
                        f"{group_name}\t{group_count}\t{group_max:.4f}\t{group_pairs}\n"
                    )
            if failing_groups:
                f.write("Failing groups: genotype_pair\tsequence_count\tmax_similarity\tpairwise_comparisons\n")
                for group_name, group_count, group_max, group_pairs in failing_groups:
                    f.write(
                        f"{group_name}\t{group_count}\t{group_max:.4f}\t{group_pairs}\n"
                    )
        
        print(f"Validation completed: max similarity = {max_similarity:.4f}%, needs reclustering = {needs_reclustering}")
        
    except subprocess.CalledProcessError as e:
        error_message = e.stderr.strip() if e.stderr else str(e)
        with open(output_report, 'w') as f:
            f.write("Validation failed\n")
            f.write(f"Reason: MAFFT alignment failed: {error_message}\n")
        raise RuntimeError(f"MAFFT alignment failed: {error_message}") from e
    finally:
        if temp_fasta_path and os.path.exists(temp_fasta_path):
            os.unlink(temp_fasta_path)
        if aligned_file and os.path.exists(aligned_file):
            os.unlink(aligned_file)

def main():
    parser = argparse.ArgumentParser(description="Validate internal similarity of consensus sequences")
    parser.add_argument("--consensus_file", required=True, help="Input consensus FASTA file")
    parser.add_argument("--output_report", required=True, help="Output validation report file")
    parser.add_argument("--similarity_threshold", type=float, default=95.0, help="Similarity threshold (default: 95.0)")
    parser.add_argument(
        "--group_by",
        choices=["none", "genotype_pair"],
        default="none",
        help="Validate all records together or within extracted RdRp_VP1 genotype pairs",
    )
    
    args = parser.parse_args()
    
    validate_consensus_internal(
        args.consensus_file,
        args.output_report,
        args.similarity_threshold,
        args.group_by,
    )

if __name__ == "__main__":
    main() 
