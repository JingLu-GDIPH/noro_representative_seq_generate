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
from Bio import SeqIO
from Bio import AlignIO


def calculate_pairwise_similarity(seq1, seq2):
    """Calculate pairwise similarity between two sequences."""
    if len(seq1) != len(seq2):
        return 0.0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
    total = sum(1 for a, b in zip(seq1, seq2) if a != '-' or b != '-')
    
    if total == 0:
        return 0.0
    
    return (matches / total) * 100


def calculate_max_pairwise_similarity(alignment, chunk_size=256):
    """Calculate maximum pairwise similarity with bounded-memory NumPy chunks."""
    sequence_strings = [str(record.seq) for record in alignment]
    sequence_count = len(sequence_strings)
    sequence_length = len(sequence_strings[0])
    encoded = np.frombuffer(
        "".join(sequence_strings).encode("ascii"), dtype="S1"
    ).reshape(sequence_count, sequence_length)

    max_similarity = 0.0
    gap = np.bytes_("-")

    for i in range(sequence_count - 1):
        reference = encoded[i]
        for start in range(i + 1, sequence_count, chunk_size):
            comparison = encoded[start:min(start + chunk_size, sequence_count)]
            matches = ((comparison == reference) & (reference != gap)).sum(axis=1)
            totals = ((comparison != gap) | (reference != gap)).sum(axis=1)
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


def validate_consensus_internal(consensus_file, output_report, similarity_threshold):
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
                    "mafft", "--auto", "--maxiterate", "1000", "--localpair",
                    temp_fasta_path
                ], stdout=aligned_handle, stderr=subprocess.PIPE, check=True, text=True)

            alignment = AlignIO.read(aligned_file, "fasta")
            alignment_method = "MAFFT"
        
        # Calculate pairwise similarities
        max_similarity, pair_count = calculate_max_pairwise_similarity(alignment)
        
        # Determine if reclustering is needed
        needs_reclustering = max_similarity > similarity_threshold
        reason = f"Max pairwise similarity ({max_similarity:.2f}%) exceeds threshold ({similarity_threshold}%)" if needs_reclustering else "All pairwise similarities below threshold"
        
        # Write report
        with open(output_report, 'w') as f:
            f.write(f"Max similarity: {max_similarity:.2f}\n")
            f.write(f"Needs reclustering: {str(needs_reclustering).lower()}\n")
            f.write(f"Reason: {reason}\n")
            f.write(f"Total sequences: {len(sequences)}\n")
            f.write(f"Total pairwise comparisons: {pair_count}\n")
            f.write(f"Alignment method: {alignment_method}\n")
        
        print(f"Validation completed: max similarity = {max_similarity:.2f}%, needs reclustering = {needs_reclustering}")
        
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
    
    args = parser.parse_args()
    
    validate_consensus_internal(args.consensus_file, args.output_report, args.similarity_threshold)

if __name__ == "__main__":
    main() 
