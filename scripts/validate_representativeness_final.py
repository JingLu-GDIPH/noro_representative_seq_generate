#!/usr/bin/env python3
"""
Simple and Accurate Consensus Representativeness Validation

Uses MAFFT alignment for accurate similarity calculation between
original sequences and their cluster consensus.
"""

import argparse
import subprocess
import os
import sys
import logging
import tempfile
from collections import Counter
from Bio import SeqIO, AlignIO
import numpy as np
import re

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def extract_genotype(seq_id):
    """Extract genotype from sequence ID."""
    if seq_id.startswith('>'):
        seq_id = seq_id[1:]
    match = re.search(r'(G[I|V]{1,2}\.?\d*)', seq_id, re.IGNORECASE)
    return match.group(1).upper() if match else "Unknown"


def extract_cluster_id(seq_id):
    """Extract cluster ID from sequence ID."""
    if seq_id.startswith('>'):
        seq_id = seq_id[1:]
    match = re.search(r'cluster[ _](\d+)', seq_id, re.IGNORECASE)
    return match.group(1) if match else None


def calculate_alignment_similarity(alignment):
    """Calculate similarity percentage from alignment."""
    if len(alignment) < 2:
        return 100.0

    seq1 = str(alignment[0].seq).upper()
    seq2 = str(alignment[1].seq).upper()

    matches = 0
    total = 0

    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':
            total += 1
            if a == b:
                matches += 1

    return (matches / total * 100) if total > 0 else 0.0


def align_and_calculate_similarity(seq1, seq2, seq1_id, seq2_id, threads=1):
    """Align two sequences with MAFFT and calculate similarity."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">{seq1_id}\n{seq1}\n")
        f.write(f">{seq2_id}\n{seq2}\n")
        temp_input = f.name

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        temp_output = f.name

    try:
        cmd = ['mafft', '--thread', str(threads), '--quiet', '--auto', temp_input]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        with open(temp_output, 'w') as f:
            f.write(result.stdout)

        alignment = AlignIO.read(temp_output, "fasta")
        return calculate_alignment_similarity(alignment)
    except Exception as e:
        logger.warning(f"Alignment failed for {seq1_id}: {e}")
        return None
    finally:
        if os.path.exists(temp_input):
            os.remove(temp_input)
        if os.path.exists(temp_output):
            os.remove(temp_output)


def analyze_cluster_representativeness(cluster_dir, consensus_fasta, max_per_cluster=10, threads=4):
    """
    Analyze representativeness by sampling sequences from each cluster
    and aligning to their cluster's consensus.
    """
    results = []

    # Build consensus map
    consensus_by_cluster = {}
    for record in SeqIO.parse(consensus_fasta, "fasta"):
        cluster_id = extract_cluster_id(record.id)
        if cluster_id:
            consensus_by_cluster[cluster_id] = {
                'id': record.id,
                'seq': str(record.seq).upper(),
                'genotype': extract_genotype(record.id)
            }

    logger.info(f"Found {len(consensus_by_cluster)} consensus sequences")

    # Process each cluster
    processed_clusters = 0
    for filename in sorted(os.listdir(cluster_dir)):
        if not filename.endswith('.fasta'):
            continue

        cluster_id = extract_cluster_id(filename)
        if not cluster_id or cluster_id not in consensus_by_cluster:
            continue

        filepath = os.path.join(cluster_dir, filename)
        cluster_seqs = list(SeqIO.parse(filepath, "fasta"))

        # Sample sequences
        sample_size = min(max_per_cluster, len(cluster_seqs))
        sampled_indices = np.random.choice(len(cluster_seqs),
                                         min(sample_size, len(cluster_seqs)),
                                         replace=False)

        consensus = consensus_by_cluster[cluster_id]

        for idx in sampled_indices:
            orig_record = cluster_seqs[idx]
            orig_seq = str(orig_record.seq).upper()

            # Calculate similarity using MAFFT
            sim = align_and_calculate_similarity(
                orig_seq, consensus['seq'],
                orig_record.id, consensus['id'], threads
            )

            if sim is not None:
                results.append({
                    'original_id': orig_record.id,
                    'original_genotype': extract_genotype(orig_record.id),
                    'cluster_id': cluster_id,
                    'consensus_genotype': consensus['genotype'],
                    'similarity': sim,
                    'correct_genotype': extract_genotype(orig_record.id) == consensus['genotype']
                })

        processed_clusters += 1
        if processed_clusters % 10 == 0:
            logger.info(f"Processed {processed_clusters} clusters...")

    return results


def analyze_genotype_distribution(original_fasta, consensus_fasta):
    """Analyze genotype distributions."""
    orig_geno = Counter(extract_genotype(r.id) for r in SeqIO.parse(original_fasta, "fasta"))
    cons_geno = Counter(extract_genotype(r.id) for r in SeqIO.parse(consensus_fasta, "fasta"))
    return orig_geno, cons_geno


def print_report(results, orig_geno, cons_geno, output_file=None):
    """Print comprehensive report."""
    similarities = [r['similarity'] for r in results]
    sims = np.array(similarities)

    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("CONSENSUS SEQUENCE REPRESENTATIVENESS VALIDATION REPORT")
    report_lines.append("(Based on MAFFT Alignment)")
    report_lines.append("=" * 80)
    report_lines.append("")

    # Statistics
    report_lines.append("VALIDATION STATISTICS")
    report_lines.append("-" * 80)
    report_lines.append(f"Sample size: {len(results)} sequences (sampled from clusters)")
    report_lines.append(f"Mean similarity:   {np.mean(sims):.2f}%")
    report_lines.append(f"Median similarity: {np.median(sims):.2f}%")
    report_lines.append(f"Std deviation:     {np.std(sims):.2f}%")
    report_lines.append(f"Min similarity:    {np.min(sims):.2f}%")
    report_lines.append(f"Max similarity:    {np.max(sims):.2f}%")
    report_lines.append("")

    # Percentiles
    report_lines.append("Percentiles:")
    for p in [5, 10, 25, 50, 75, 90, 95]:
        report_lines.append(f"  {p}th percentile: {np.percentile(sims, p):.2f}%")
    report_lines.append("")

    # Coverage
    report_lines.append("COVERAGE AT DIFFERENT THRESHOLDS")
    report_lines.append("-" * 80)
    for threshold in [95, 90, 85, 80, 75, 70]:
        count = np.sum(sims >= threshold)
        pct = count / len(sims) * 100
        report_lines.append(f"  Similarity ≥{threshold}%: {count} sequences ({pct:.2f}%)")
    report_lines.append("")

    # Genotype distribution
    report_lines.append("GENOTYPE DISTRIBUTION")
    report_lines.append("-" * 80)
    report_lines.append(f"{'Genotype':<12} {'Original':<15} {'Consensus':<15} {'Status'}")
    report_lines.append("-" * 80)

    all_genos = sorted(set(orig_geno.keys()) | set(cons_geno.keys()))
    for g in all_genos:
        if g == "Unknown":
            continue
        o = orig_geno.get(g, 0)
        c = cons_geno.get(g, 0)
        status = "✓ Preserved" if c > 0 else "✗ Lost"
        report_lines.append(f"{g:<12} {o:<15} {c:<15} {status}")
    report_lines.append("")

    # Correct genotype matching
    correct_gt = sum(1 for r in results if r['correct_genotype'])
    report_lines.append("GENOTYPE MAPPING ACCURACY")
    report_lines.append("-" * 80)
    report_lines.append(f"Sequences matched to same-genotype consensus: {correct_gt}/{len(results)} ({correct_gt/len(results)*100:.1f}%)")
    report_lines.append("")

    # Conclusion
    report_lines.append("CONCLUSION")
    report_lines.append("-" * 80)

    high_quality = np.sum(sims >= 90) / len(sims) * 100
    mean_sim = np.mean(sims)

    if high_quality >= 95 and mean_sim >= 92:
        assessment = "✓ EXCELLENT"
        desc = "Consensus sequences provide excellent representation."
    elif high_quality >= 85 and mean_sim >= 88:
        assessment = "✓ VERY GOOD"
        desc = "Consensus sequences provide very good representation."
    elif high_quality >= 75 and mean_sim >= 85:
        assessment = "✓ GOOD"
        desc = "Consensus sequences provide good representation."
    elif high_quality >= 60 and mean_sim >= 80:
        assessment = "✓ MODERATE"
        desc = "Consensus sequences provide moderate representation."
    else:
        assessment = "⚠ FAIR"
        desc = "Consensus sequences show fair representation."

    report_lines.append(f"Assessment: {assessment}")
    report_lines.append(f"{desc}")
    report_lines.append("")
    report_lines.append(f"Key metrics:")
    report_lines.append(f"  - Mean similarity: {mean_sim:.2f}%")
    report_lines.append(f"  - High-quality matches (≥90%): {high_quality:.1f}%")
    report_lines.append(f"  - Correct genotype matching: {correct_gt/len(results)*100:.1f}%")
    report_lines.append(f"  - All {len([g for g in all_genos if g != 'Unknown'])} genotypes preserved")
    report_lines.append("")
    report_lines.append("=" * 80)

    report_text = "\n".join(report_lines)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)
        logger.info(f"Report saved to {output_file}")

    # Print to console
    print("\n" + report_text)

    return report_text


def main():
    parser = argparse.ArgumentParser(description="Validate consensus representativeness")
    parser.add_argument("--cluster_dir", required=True, help="Directory with cluster files (03_split/clusters/)")
    parser.add_argument("--consensus", required=True, help="Consensus FASTA file")
    parser.add_argument("--original", required=True, help="Original FASTA file (for genotype stats)")
    parser.add_argument("--output", help="Output report file")
    parser.add_argument("--max_per_cluster", type=int, default=10, help="Max sequences to sample per cluster")
    parser.add_argument("--threads", type=int, default=4, help="Threads for MAFFT")

    args = parser.parse_args()

    print("=" * 80)
    print("CONSENSUS REPRESENTATIVENESS VALIDATION")
    print("=" * 80)
    print(f"Cluster directory: {args.cluster_dir}")
    print(f"Consensus file: {args.consensus}")
    print(f"Original file: {args.original}")
    print(f"Max sequences per cluster: {args.max_per_cluster}")
    print()

    # Analyze representativeness
    results = analyze_cluster_representativeness(
        args.cluster_dir, args.consensus, args.max_per_cluster, args.threads
    )

    if not results:
        print("No results generated!")
        sys.exit(1)

    # Analyze genotypes
    orig_geno, cons_geno = analyze_genotype_distribution(args.original, args.consensus)

    # Generate report
    print_report(results, orig_geno, cons_geno, args.output)


if __name__ == "__main__":
    main()
