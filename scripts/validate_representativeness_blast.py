#!/usr/bin/env python3
"""
Consensus Representativeness Validation using BLAST

This script uses BLAST to accurately map each original sequence to its
most similar consensus sequence, providing proper similarity estimates.
"""

import argparse
import subprocess
import os
import sys
import logging
from collections import Counter, defaultdict
from Bio import SeqIO
import numpy as np
import re

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def extract_genotype(sequence_id):
    """Extract genotype from sequence ID."""
    # Remove leading '>'
    if sequence_id.startswith('>'):
        sequence_id = sequence_id[1:]

    # Try GII.X pattern first
    match = re.search(r'(G[I|V]{1,2}\.?\d*)', sequence_id, re.IGNORECASE)
    if match:
        return match.group(1).upper()

    return "Unknown"


def run_blast(original_fasta, consensus_fasta, output_prefix, threads=8):
    """
    Run BLASTN to map original sequences to consensus sequences.
    """
    logger.info("Running BLASTN analysis...")

    # Create BLAST database from consensus sequences
    db_name = f"{output_prefix}_consensus_db"
    cmd_makedb = [
        "makeblastdb",
        "-in", consensus_fasta,
        "-dbtype", "nucl",
        "-out", db_name,
        "-parse_seqids"
    ]

    logger.info(f"Creating BLAST database: {' '.join(cmd_makedb)}")
    result = subprocess.run(cmd_makedb, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"makeblastdb failed: {result.stderr}")
        raise RuntimeError("Failed to create BLAST database")

    # Run BLASTN
    blast_output = f"{output_prefix}_blast_results.txt"
    cmd_blast = [
        "blastn",
        "-query", original_fasta,
        "-db", db_name,
        "-out", blast_output,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-num_threads", str(int(threads)),
        "-max_target_seqs", "10",  # Keep top 10 hits per query
        "-evalue", "1e-10"
    ]

    logger.info(f"Running BLASTN: {' '.join(cmd_blast)}")
    result = subprocess.run(cmd_blast, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"blastn failed: {result.stderr}")
        raise RuntimeError("BLASTN failed")

    logger.info(f"BLAST results written to {blast_output}")

    # Clean up database files
    for ext in ['.nhr', '.nin', '.nsq']:
        f = f"{db_name}{ext}"
        if os.path.exists(f):
            os.remove(f)

    return blast_output


def parse_blast_results(blast_file):
    """
    Parse BLAST results and find the best hit for each query sequence.
    """
    logger.info("Parsing BLAST results...")

    best_hits = {}  # query_id -> best hit info
    all_hits = defaultdict(list)  # query_id -> list of all hits

    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue

            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])  # Percentage identity
            length = int(parts[3])  # Alignment length
            mismatch = int(parts[4])
            gapopen = int(parts[5])
            qstart = int(parts[6])
            qend = int(parts[7])
            sstart = int(parts[8])
            send = int(parts[9])
            evalue = float(parts[10])
            bitscore = float(parts[11])

            hit_info = {
                'sseqid': sseqid,
                'pident': pident,
                'length': length,
                'qstart': qstart,
                'qend': qend,
                'sstart': sstart,
                'send': send,
                'evalue': evalue,
                'bitscore': bitscore
            }

            all_hits[qseqid].append(hit_info)

            # Keep the best hit (highest bitscore)
            if qseqid not in best_hits or hit_info['bitscore'] > best_hits[qseqid]['bitscore']:
                best_hits[qseqid] = hit_info

    logger.info(f"Found hits for {len(best_hits)} query sequences")

    return best_hits, all_hits


def calculate_coverage_stats(original_fasta, best_hits):
    """
    Calculate coverage statistics based on BLAST results.
    """
    logger.info("Calculating coverage statistics...")

    # Get sequence lengths
    seq_lengths = {}
    for record in SeqIO.parse(original_fasta, "fasta"):
        seq_lengths[record.id] = len(record.seq)

    # Calculate statistics
    similarities = []
    coverages = []
    results = []

    for qseqid, hit in best_hits.items():
        query_len = seq_lengths.get(qseqid, 0)
        aligned_length = hit['length']
        coverage = (aligned_length / query_len * 100) if query_len > 0 else 0
        similarity = hit['pident']

        similarities.append(similarity)
        coverages.append(coverage)

        results.append({
            'query_id': qseqid,
            'subject_id': hit['sseqid'],
            'similarity': similarity,
            'coverage': coverage,
            'query_genotype': extract_genotype(qseqid),
            'subject_genotype': extract_genotype(hit['sseqid'])
        })

    # Calculate statistics
    similarities = np.array(similarities)
    coverages = np.array(coverages)

    stats = {
        'similarity': {
            'mean': np.mean(similarities),
            'median': np.median(similarities),
            'std': np.std(similarities),
            'min': np.min(similarities),
            'max': np.max(similarities),
        },
        'coverage': {
            'mean': np.mean(coverages),
            'median': np.median(coverages),
            'std': np.std(coverages),
            'min': np.min(coverages),
            'max': np.max(coverages),
        }
    }

    # Add percentiles for similarity
    for p in [5, 10, 25, 50, 75, 90, 95]:
        stats['similarity'][f'p{p}'] = np.percentile(similarities, p)

    # Count coverage at different thresholds
    coverage_counts = {}
    for threshold in [70, 75, 80, 85, 90, 95]:
        count = np.sum((similarities >= threshold) & (coverages >= 80))
        percentage = (count / len(similarities)) * 100
        coverage_counts[f'sim>={threshold}%&cov>=80%'] = (count, percentage)

    # Also count by similarity only
    for threshold in [70, 75, 80, 85, 90, 95]:
        count = np.sum(similarities >= threshold)
        percentage = (count / len(similarities)) * 100
        coverage_counts[f'sim>={threshold}%'] = (count, percentage)

    return results, stats, coverage_counts


def analyze_genotype_mapping(original_fasta, consensus_fasta, results):
    """
    Analyze genotype mapping between original and consensus sequences.
    """
    logger.info("Analyzing genotype mapping...")

    # Count genotypes in original
    original_geno_count = Counter()
    for record in SeqIO.parse(original_fasta, "fasta"):
        gt = extract_genotype(record.id)
        original_geno_count[gt] += 1

    # Count genotypes in consensus
    consensus_geno_count = Counter()
    for record in SeqIO.parse(consensus_fasta, "fasta"):
        gt = extract_genotype(record.id)
        consensus_geno_count[gt] += 1

    # Analyze mapping patterns
    mapping_patterns = defaultdict(lambda: {'correct': 0, 'cross': 0, 'total': 0})
    cross_mapping_counts = Counter()

    for r in results:
        query_gt = r['query_genotype']
        subject_gt = r['subject_genotype']

        mapping_patterns[query_gt]['total'] += 1

        if query_gt == subject_gt:
            mapping_patterns[query_gt]['correct'] += 1
        else:
            mapping_patterns[query_gt]['cross'] += 1
            cross_mapping_counts[(query_gt, subject_gt)] += 1

    return original_geno_count, consensus_geno_count, mapping_patterns, cross_mapping_counts


def generate_report(results, stats, coverage_counts, original_geno_count,
                   consensus_geno_count, mapping_patterns, cross_mapping_counts,
                   original_fasta, consensus_fasta, output_prefix):
    """
    Generate comprehensive validation report.
    """
    report_file = f"{output_prefix}_validation_report.txt"

    num_original = len(original_geno_count)
    num_consensus = len(consensus_geno_count)
    total_original_seqs = sum(original_geno_count.values())

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CONSENSUS SEQUENCE REPRESENTATIVENESS VALIDATION REPORT\n")
        f.write("(Based on BLASTN Alignment)\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Original sequences: {original_fasta}\n")
        f.write(f"Consensus sequences: {consensus_fasta}\n\n")

        # Summary statistics
        f.write("1. SUMMARY STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total original sequences: {total_original_seqs}\n")
        f.write(f"Total consensus sequences: {sum(consensus_geno_count.values())}\n")
        f.write(f"Compression ratio: {total_original_seqs / max(1, sum(consensus_geno_count.values())):.2f}:1\n")
        f.write(f"Original genotypes: {num_original}\n")
        f.write(f"Consensus genotypes: {num_consensus}\n\n")

        f.write("2. SIMILARITY & COVERAGE STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write("(Based on BLASTN alignment of each original to its best-matching consensus)\n\n")

        f.write("Similarity (BLAST % identity):\n")
        f.write(f"  Mean:   {stats['similarity']['mean']:.2f}%\n")
        f.write(f"  Median: {stats['similarity']['median']:.2f}%\n")
        f.write(f"  Std Dev: {stats['similarity']['std']:.2f}%\n")
        f.write(f"  Min:    {stats['similarity']['min']:.2f}%\n")
        f.write(f"  Max:    {stats['similarity']['max']:.2f}%\n\n")

        f.write("Similarity Percentiles:\n")
        for p in [5, 10, 25, 50, 75, 90, 95]:
            key = f'p{p}'
            if key in stats['similarity']:
                f.write(f"  {p}th percentile: {stats['similarity'][key]:.2f}%\n")
        f.write("\n")

        f.write("Coverage (alignment length / query length):\n")
        f.write(f"  Mean:   {stats['coverage']['mean']:.2f}%\n")
        f.write(f"  Median: {stats['coverage']['median']:.2f}%\n")
        f.write(f"  Min:    {stats['coverage']['min']:.2f}%\n")
        f.write(f"  Max:    {stats['coverage']['max']:.2f}%\n\n")

        f.write("Coverage at Different Thresholds:\n")
        for label, (count, percentage) in sorted(coverage_counts.items(), key=lambda x: x[1][1], reverse=True):
            f.write(f"  Sequences with {label}: {count} ({percentage:.2f}%)\n")
        f.write("\n")

        # Genotype analysis
        f.write("3. GENOTYPE DISTRIBUTION\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Genotype':<15} {'Original':<15} {'Consensus':<15} {'Status':<15}\n")
        f.write("-" * 80 + "\n")

        all_genotypes = sorted(set(original_geno_count.keys()) | set(consensus_geno_count.keys()))
        for gt in all_genotypes:
            orig = original_geno_count.get(gt, 0)
            cons = consensus_geno_count.get(gt, 0)
            status = "✓ Preserved" if cons > 0 else "✗ Lost"
            f.write(f"{gt:<15} {orig:<15} {cons:<15} {status:<15}\n")
        f.write("\n")

        # Mapping patterns
        f.write("4. GENOTYPE MAPPING ANALYSIS\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Genotype':<15} {'Total':<10} {'Same GT':<12} {'Cross GT':<12} {'Correct Rate':<15}\n")
        f.write("-" * 80 + "\n")

        for gt in sorted(mapping_patterns.keys()):
            data = mapping_patterns[gt]
            total = data['total']
            correct = data['correct']
            cross = data['cross']
            rate = (correct / total * 100) if total > 0 else 0
            f.write(f"{gt:<15} {total:<10} {correct:<12} {cross:<12} {rate:<15.1f}%\n")
        f.write("\n")

        # Cross-mapping details
        if cross_mapping_counts:
            f.write("Cross-Genotype Mappings (Top 10):\n")
            for (query, subject), count in cross_mapping_counts.most_common(10):
                f.write(f"  {query} → {subject}: {count}\n")
            f.write("\n")

        # Conclusion
        f.write("5. CONCLUSION\n")
        f.write("-" * 80 + "\n")

        high_quality_coverage = coverage_counts.get('sim>=90%&cov>=80%', (0, 0))[1]
        mean_sim = stats['similarity']['mean']
        mean_cov = stats['coverage']['mean']

        f.write(f"Key Metrics:\n")
        f.write(f"  - Mean similarity: {mean_sim:.2f}%\n")
        f.write(f"  - Mean coverage: {mean_cov:.2f}%\n")
        f.write(f"  - High-quality matches (≥90% sim, ≥80% cov): {high_quality_coverage:.1f}%\n\n")

        if high_quality_coverage >= 90 and mean_sim >= 90:
            f.write("✓ EXCELLENT: The consensus sequences provide excellent representation.\n")
        elif high_quality_coverage >= 80 and mean_sim >= 85:
            f.write("✓ GOOD: The consensus sequences provide good representation.\n")
        elif high_quality_coverage >= 70 and mean_sim >= 80:
            f.write("✓ MODERATE: The consensus sequences provide moderate representation.\n")
        else:
            f.write("⚠ FAIR: The consensus sequences show fair representation.\n")
            f.write("  Consider reviewing clustering parameters.\n")

        f.write("\n")

        # Genotype preservation
        preserved = sum(1 for gt in original_geno_count.keys() if gt in consensus_geno_count and gt != "Unknown")
        total_known = sum(1 for gt in original_geno_count.keys() if gt != "Unknown")

        if preserved == total_known:
            f.write(f"✓ All {total_known} known genotypes from original data are preserved.\n")
        else:
            f.write(f"⚠ {preserved}/{total_known} known genotypes preserved.\n")

    logger.info(f"Validation report written to {report_file}")
    return report_file


def main():
    parser = argparse.ArgumentParser(
        description="Validate consensus sequence representativeness using BLAST"
    )
    parser.add_argument("--original", required=True, help="Original input FASTA file")
    parser.add_argument("--consensus", required=True, help="Final consensus FASTA file")
    parser.add_argument("--output_prefix", required=True, help="Output file prefix")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")

    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.original):
        logger.error(f"Original file not found: {args.original}")
        sys.exit(1)

    if not os.path.exists(args.consensus):
        logger.error(f"Consensus file not found: {args.consensus}")
        sys.exit(1)

    logger.info("=" * 80)
    logger.info("CONSENSUS SEQUENCE REPRESENTATIVENESS VALIDATION (BLAST-based)")
    logger.info("=" * 80)
    logger.info(f"Original sequences: {args.original}")
    logger.info(f"Consensus sequences: {args.consensus}")

    # Run BLAST
    blast_file = run_blast(args.original, args.consensus, args.output_prefix, args.threads)

    # Parse BLAST results
    best_hits, all_hits = parse_blast_results(blast_file)

    # Calculate statistics
    results, stats, coverage_counts = calculate_coverage_stats(args.original, best_hits)

    # Analyze genotype mapping
    orig_geno, cons_geno, mapping_patterns, cross_mapping = analyze_genotype_mapping(
        args.original, args.consensus, results
    )

    # Generate report
    report_file = generate_report(
        results, stats, coverage_counts, orig_geno, cons_geno,
        mapping_patterns, cross_mapping, args.original, args.consensus, args.output_prefix
    )

    # Print summary
    logger.info("=" * 80)
    logger.info("VALIDATION SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total original sequences: {sum(orig_geno.values())}")
    logger.info(f"Total consensus sequences: {sum(cons_geno.values())}")
    logger.info(f"Mean similarity: {stats['similarity']['mean']:.2f}%")
    logger.info(f"Mean coverage: {stats['coverage']['mean']:.2f}%")
    logger.info(f"High-quality matches (≥90% sim, ≥80% cov): {coverage_counts.get('sim>=90%&cov>=80%', (0,0))[1]:.1f}%")
    logger.info(f"Genotypes preserved: {sum(1 for gt in orig_geno.keys() if gt in cons_geno and gt != 'Unknown')}/{sum(1 for gt in orig_geno.keys() if gt != 'Unknown')}")
    logger.info("=" * 80)
    logger.info(f"\nReport saved to: {report_file}")


if __name__ == "__main__":
    main()
