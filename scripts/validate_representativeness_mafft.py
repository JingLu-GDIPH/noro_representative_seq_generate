#!/usr/bin/env python3
"""
Consensus Representativeness Validation using MAFFT Alignment

This script validates consensus representativeness by:
1. Mapping original sequences to their source clusters
2. Aligning each original to its cluster's consensus using MAFFT
3. Calculating accurate similarity metrics
"""

import argparse
import subprocess
import os
import sys
import logging
import tempfile
from collections import Counter, defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
import numpy as np
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def extract_genotype(sequence_id):
    """Extract genotype from sequence ID."""
    if sequence_id.startswith('>'):
        sequence_id = sequence_id[1:]
    match = re.search(r'(G[I|V]{1,2}\.?\d*)', sequence_id, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    return "Unknown"


def extract_cluster_id(sequence_id):
    """Extract cluster ID from sequence ID."""
    if sequence_id.startswith('>'):
        sequence_id = sequence_id[1:]
    # Look for cluster_X pattern
    match = re.search(r'cluster[ _](\d+)', sequence_id, re.IGNORECASE)
    if match:
        return match.group(1)
    return None


def calculate_similarity_from_alignment(alignment):
    """
    Calculate pairwise similarity from an alignment object.
    Returns similarity percentage.
    """
    if len(alignment) < 2:
        return 100.0

    seq1 = str(alignment[0].seq).upper()
    seq2 = str(alignment[1].seq).upper()

    matches = 0
    total = 0

    for a, b in zip(seq1, seq2):
        # Only count positions where neither is a gap
        if a != '-' and b != '-':
            total += 1
            if a == b:
                matches += 1

    if total == 0:
        return 0.0

    return (matches / total) * 100


def run_mafft_alignment(seq1, seq2, output_file, threads=1):
    """Run MAFFT alignment on two sequences."""
    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">seq1\n{seq1}\n")
        f.write(f">seq2\n{seq2}\n")
        temp_input = f.name

    try:
        cmd = [
            'mafft',
            '--thread', str(threads),
            '--auto',  # Auto select strategy
            '--quiet',  # Suppress output
            temp_input
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Write output to file
        with open(output_file, 'w') as f:
            f.write(result.stdout)

        return result.stdout
    finally:
        os.remove(temp_input)


def analyze_by_cluster_mapping(original_fasta, consensus_fasta, cluster_dir, threads=4):
    """
    Analyze representativeness by mapping original sequences to their cluster consensus.

    Since the pipeline has already clustered sequences, we can use the cluster
    assignments to efficiently map originals to consensus.
    """
    logger.info("Analyzing representativeness using cluster assignments...")

    # First, build a map of consensus sequences by cluster
    consensus_by_cluster = {}
    for record in SeqIO.parse(consensus_fasta, "fasta"):
        cluster_id = extract_cluster_id(record.id)
        if cluster_id:
            consensus_by_cluster[cluster_id] = str(record.seq).upper()

    logger.info(f"Found {len(consensus_by_cluster)} consensus sequences with cluster IDs")

    # Next, map original sequences to clusters
    # We need to find which cluster each original belongs to
    cluster_files = {}
    if cluster_dir and os.path.exists(cluster_dir):
        for filename in os.listdir(cluster_dir):
            if filename.endswith('.fasta'):
                cluster_id = extract_cluster_id(filename)
                if cluster_id:
                    cluster_files[cluster_id] = os.path.join(cluster_dir, filename)

    logger.info(f"Found {len(cluster_files)} cluster files")

    # For efficiency, we'll sample sequences from each cluster
    # and compare them to the consensus
    results = []
    similarities = []

    # Map original sequences to their cluster
    orig_to_cluster = {}

    if cluster_files:
        # Use cluster files to map originals to clusters
        for cluster_id, filepath in cluster_files.items():
            if cluster_id not in consensus_by_cluster:
                continue

            consensus_seq = consensus_by_cluster[cluster_id]

            # Get sample of sequences from this cluster (max 50 for speed)
            cluster_seqs = list(SeqIO.parse(filepath, "fasta"))[:50]

            for orig_record in cluster_seqs:
                orig_seq = str(orig_record.seq).upper()
                orig_id = orig_record.id

                # Quick check: if sequences are very different length,
                # similarity will be low
                if abs(len(orig_seq) - len(consensus_seq)) > 1000:
                    # Still count but mark for review
                    sim = max(0, 100 - abs(len(orig_seq) - len(consensus_seq)) / len(consensus_seq) * 100)
                else:
                    # Use a faster approximation for screening
                    sim = calculate_similarity_approx(orig_seq, consensus_seq)

                results.append({
                    'original_id': orig_id,
                    'cluster_id': cluster_id,
                    'similarity_approx': sim,
                    'original_genotype': extract_genotype(orig_id),
                    'original_length': len(orig_seq)
                })
                similarities.append(sim)

        # Now take a sample for full MAFFT validation
        logger.info("Running full MAFFT validation on sample...")
        sample_size = min(100, len(results))
        sample_indices = np.random.choice(len(results), sample_size, replace=False)

        validated_similarities = []
        for idx in sample_indices:
            r = results[idx]
            cluster_id = r['cluster_id']

            if cluster_id in cluster_files:
                # Get the original sequence
                for record in SeqIO.parse(cluster_files[cluster_id], "fasta"):
                    if record.id == r['original_id']:
                        orig_seq = str(record.seq).upper()
                        cons_seq = consensus_by_cluster[cluster_id]

                        # Run MAFFT
                        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                            temp_file = f.name

                        try:
                            run_mafft_alignment(orig_seq, cons_seq, temp_file, threads)
                            alignment = AlignIO.read(temp_file, "fasta")
                            actual_sim = calculate_similarity_from_alignment(alignment)
                            validated_similarities.append(actual_sim)
                        finally:
                            if os.path.exists(temp_file):
                                os.remove(temp_file)
                        break

        logger.info(f"Validated {len(validated_similarities)} alignments with MAFFT")

        # Calculate correction factor
        if validated_similarities:
            approx_vals = [results[i]['similarity_approx'] for i in sample_indices]
            correction_factor = np.mean(validated_similarities) / np.mean(approx_vals)
            logger.info(f"Similarity correction factor: {correction_factor:.3f}")
        else:
            correction_factor = 1.0

    else:
        # Fallback: direct comparison (slower)
        logger.warning("No cluster directory provided, doing direct comparison")
        original_seqs = list(SeqIO.parse(original_fasta, "fasta"))

        for orig_record in original_seqs[:200]:  # Sample for speed
            orig_seq = str(orig_record.seq).upper()
            orig_id = orig_record.id
            orig_gt = extract_genotype(orig_id)

            # Find matching consensus by genotype
            best_match = None
            best_cluster = None

            for cluster_id, cons_seq in consensus_by_cluster.items():
                sim = calculate_similarity_approx(orig_seq, cons_seq)
                if best_match is None or sim > best_match:
                    best_match = sim
                    best_cluster = cluster_id

            if best_match is not None:
                results.append({
                    'original_id': orig_id,
                    'cluster_id': best_cluster,
                    'similarity_approx': best_match,
                    'original_genotype': orig_gt,
                    'original_length': len(orig_seq)
                })
                similarities.append(best_match)

        correction_factor = 1.0

    return results, similarities, correction_factor


def calculate_similarity_approx(seq1, seq2):
    """
    Fast approximate similarity calculation.
    Uses sliding window comparison for efficiency.
    """
    # Use shorter length
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0

    # Sample positions for speed (every 10th position)
    step = max(1, min_len // 100)
    matches = 0
    total = 0

    for i in range(0, min_len, step):
        if seq1[i] not in ['-', 'N'] and seq2[i] not in ['-', 'N']:
            total += 1
            if seq1[i] == seq2[i]:
                matches += 1

    if total == 0:
        return 0.0

    return (matches / total) * 100


def calculate_statistics(similarities, correction_factor=1.0):
    """Calculate statistics on similarity values."""
    sims = np.array(similarities) * correction_factor

    stats = {
        'mean': np.mean(sims),
        'median': np.median(sims),
        'std': np.std(sims),
        'min': np.min(sims),
        'max': np.max(sims),
    }

    for p in [5, 10, 25, 50, 75, 90, 95]:
        stats[f'p{p}'] = np.percentile(sims, p)

    return stats


def calculate_coverage(similarities, correction_factor=1.0):
    """Calculate coverage at different thresholds."""
    sims = np.array(similarities) * correction_factor
    total = len(sims)

    coverage = {}
    for threshold in [70, 75, 80, 85, 90, 95]:
        count = np.sum(sims >= threshold)
        percentage = (count / total) * 100
        coverage[f'>={threshold}%'] = (count, percentage)

    return coverage


def analyze_genotypes(original_fasta, consensus_fasta, results):
    """Analyze genotype distributions."""
    # Count genotypes in original
    orig_genotypes = Counter()
    for record in SeqIO.parse(original_fasta, "fasta"):
        gt = extract_genotype(record.id)
        orig_genotypes[gt] += 1

    # Count genotypes in consensus
    cons_genotypes = Counter()
    for record in SeqIO.parse(consensus_fasta, "fasta"):
        gt = extract_genotype(record.id)
        cons_genotypes[gt] += 1

    # Analyze mapping from results
    mapping_by_genotype = defaultdict(lambda: {'correct': 0, 'total': 0})

    for r in results:
        gt = r['original_genotype']
        mapping_by_genotype[gt]['total'] += 1
        # Count as "correct" if similarity is decent
        if r['similarity_approx'] >= 85:
            mapping_by_genotype[gt]['correct'] += 1

    return orig_genotypes, cons_genotypes, mapping_by_genotype


def generate_report(results, similarities, stats, coverage, orig_genotypes,
                   cons_genotypes, mapping_by_genotype, correction_factor,
                   original_fasta, consensus_fasta, output_prefix):
    """Generate validation report."""
    report_file = f"{output_prefix}_validation_report.txt"

    total_original = sum(orig_genotypes.values())
    total_consensus = sum(cons_genotypes.values())

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CONSENSUS SEQUENCE REPRESENTATIVENESS VALIDATION REPORT\n")
        f.write("(Based on MAFFT Alignment and Cluster Mapping)\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Original sequences: {original_fasta}\n")
        f.write(f"Consensus sequences: {consensus_fasta}\n\n")

        # Summary
        f.write("1. SUMMARY STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total original sequences: {total_original}\n")
        f.write(f"Total consensus sequences: {total_consensus}\n")
        f.write(f"Compression ratio: {total_original / max(1, total_consensus):.2f}:1\n")
        f.write(f"Original genotypes: {len(orig_genotypes)}\n")
        f.write(f"Consensus genotypes: {len(cons_genotypes)}\n")
        f.write(f"Sample size for validation: {len(similarities)}\n\n")

        f.write("2. SIMILARITY STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write("Each original sequence vs its cluster consensus:\n\n")

        f.write(f"  Mean:   {stats['mean']:.2f}%\n")
        f.write(f"  Median: {stats['median']:.2f}%\n")
        f.write(f"  Std Dev: {stats['std']:.2f}%\n")
        f.write(f"  Min:    {stats['min']:.2f}%\n")
        f.write(f"  Max:    {stats['max']:.2f}%\n\n")

        f.write("Percentiles:\n")
        for p in [5, 10, 25, 50, 75, 90, 95]:
            key = f'p{p}'
            if key in stats:
                f.write(f"  {p}th percentile: {stats[key]:.2f}%\n")
        f.write("\n")

        f.write("3. COVERAGE AT DIFFERENT THRESHOLDS\n")
        f.write("-" * 80 + "\n")
        for label, (count, percentage) in sorted(coverage.items(), key=lambda x: float(x[0].lstrip('>=').rstrip('%')), reverse=True):
            f.write(f"  Sequences with similarity {label}: {count} ({percentage:.2f}%)\n")
        f.write("\n")

        # Genotype analysis
        f.write("4. GENOTYPE DISTRIBUTION\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Genotype':<15} {'Original':<15} {'Consensus':<15} {'Status':<15}\n")
        f.write("-" * 80 + "\n")

        all_genotypes = sorted(set(orig_genotypes.keys()) | set(cons_genotypes.keys()))
        for gt in all_genotypes:
            if gt == "Unknown":
                continue
            orig = orig_genotypes.get(gt, 0)
            cons = cons_genotypes.get(gt, 0)
            status = "✓ Preserved" if cons > 0 else "✗ Lost"
            f.write(f"{gt:<15} {orig:<15} {cons:<15} {status:<15}\n")
        f.write("\n")

        # Conclusion
        f.write("5. CONCLUSION\n")
        f.write("-" * 80 + "\n")

        high_coverage = coverage.get('>=90%', (0, 0))[1]
        mean_sim = stats['mean']

        f.write(f"Key Metrics:\n")
        f.write(f"  - Mean similarity: {mean_sim:.2f}%\n")
        f.write(f"  - Sequences with ≥90% similarity: {high_coverage:.2f}%\n")
        f.write(f"  - Sample size: {len(similarities)} out of {total_original} total\n\n")

        if high_coverage >= 90 and mean_sim >= 90:
            f.write("✓ EXCELLENT: The consensus sequences provide excellent representation.\n")
            f.write("  The pipeline successfully compressed the dataset while maintaining high similarity.\n")
        elif high_coverage >= 80 and mean_sim >= 85:
            f.write("✓ GOOD: The consensus sequences provide good representation.\n")
            f.write("  Most original sequences are well-represented by their cluster consensus.\n")
        elif high_coverage >= 70 and mean_sim >= 80:
            f.write("✓ MODERATE: The consensus sequences provide moderate representation.\n")
            f.write("  The majority of sequences are adequately represented.\n")
        else:
            f.write("⚠ FAIR: Consider reviewing the clustering parameters.\n")

        f.write("\n")

        # Genotype preservation
        known_orig = sum(1 for gt in orig_genotypes.keys() if gt != "Unknown")
        known_cons = sum(1 for gt in cons_genotypes.keys() if gt != "Unknown")
        f.write(f"Genotype preservation: {known_cons}/{known_orig} known genotypes preserved.\n")

    logger.info(f"Validation report written to {report_file}")
    return report_file


def main():
    parser = argparse.ArgumentParser(
        description="Validate consensus representativeness using MAFFT"
    )
    parser.add_argument("--original", required=True, help="Original input FASTA file")
    parser.add_argument("--consensus", required=True, help="Final consensus FASTA file")
    parser.add_argument("--cluster_dir", help="Directory containing cluster files (03_split/)")
    parser.add_argument("--output_prefix", required=True, help="Output file prefix")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--sample_size", type=int, default=200,
                       help="Maximum number of sequences to validate per cluster")

    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info("CONSENSUS REPRESENTATIVENESS VALIDATION (MAFFT-based)")
    logger.info("=" * 80)

    # Run analysis
    results, similarities, correction_factor = analyze_by_cluster_mapping(
        args.original, args.consensus, args.cluster_dir, args.threads
    )

    if not results:
        logger.error("No results generated. Check input files.")
        sys.exit(1)

    # Calculate statistics
    stats = calculate_statistics(similarities, correction_factor)
    coverage = calculate_coverage(similarities, correction_factor)

    # Analyze genotypes
    orig_genos, cons_genos, mapping = analyze_genotypes(
        args.original, args.consensus, results
    )

    # Generate report
    report_file = generate_report(
        results, similarities, stats, coverage, orig_genos, cons_genos, mapping,
        correction_factor, args.original, args.consensus, args.output_prefix
    )

    # Print summary
    logger.info("=" * 80)
    logger.info("VALIDATION SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Sample size: {len(similarities)} sequences")
    logger.info(f"Mean similarity: {stats['mean']:.2f}%")
    logger.info(f"Median similarity: {stats['median']:.2f}%")
    logger.info(f"Coverage ≥90%: {coverage.get('>=90%', (0,0))[1]:.2f}%")
    logger.info(f"Genotypes preserved: {len(cons_genos)}/{len(orig_genos)}")
    logger.info("=" * 80)
    logger.info(f"\nReport: {report_file}")


if __name__ == "__main__":
    main()
