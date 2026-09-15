#!/usr/bin/env python3
"""Lightweight nearest-neighbour representativeness of consensus sequences.

For each RdRp_VP1 genotype, the filtered raw records and the consensus
representatives are aligned together with MAFFT; terminal low-coverage columns
are trimmed; then every raw record is assigned the p-distance to its *nearest*
consensus sequence in that genotype. Genotype-level and global summaries report
the distance distribution and the fraction of raw records covered within
2 % / 5 % / 10 % p-distance of some consensus.

This mirrors the ``nearest_consensus_alignment_metrics`` logic of
``analyze_consensus_representativeness.py`` but is self-contained: it extracts
the genotype pair straight from the FASTA IDs (the ``<RdRp>_<VP1>`` prefix is
consistent across raw and consensus IDs), so no external genotyping TSV is
required, and it skips tree inference entirely.

Inputs
------
* ``--raw``           : filtered raw FASTA (e.g. ``Rawdata/gii_merged_sequences_filtered.fa``)
* ``--consensus``     : filtered consensus FASTA (e.g. ``*_all_final_consensus_filtered.fa``)
* ``--exclude``       : optional list of consensus IDs to drop (e.g. singletons)
* ``--outdir``        : output directory
* ``--threads``       : MAFFT threads
"""

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


VALID_BASES = frozenset("ACGT")

#: Captures the leading ``<RdRp>_<VP1>`` genotype pair from any raw/consensus ID.
GENOTYPE_PAIR_PATTERN = re.compile(
    r"(G(?:I|II|IX)\.P[A-Za-z0-9]+_G(?:I|II|IX)\.[A-Za-z0-9]+)"
)

#: Genotypes with more raw records than this use MAFFT's faster ``--retree 2``.
LARGE_GENOTYPE_THRESHOLD = 300


def extract_genotype_pair(sequence_id):
    """Return the ``<RdRp>_<VP1>`` genotype label from a sequence ID, or None."""
    match = GENOTYPE_PAIR_PATTERN.search(str(sequence_id))
    return match.group(1) if match else None


def safe_label(value):
    """Filesystem-safe label for a genotype (used in directory/file names)."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(value))


def load_by_genotype(fasta_path, exclude=None):
    """Group records by genotype pair; drop IDs in ``exclude`` if provided."""
    exclude = set(exclude or [])
    groups = defaultdict(list)
    dropped = 0
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id in exclude:
            dropped += 1
            continue
        genotype = extract_genotype_pair(record.id)
        if genotype is None:
            continue
        groups[genotype].append(record)
    return groups, dropped


def trim_alignment_terminal_columns(records, min_coverage=0.5):
    """Drop alignment columns where A/C/G/T coverage is below ``min_coverage``."""
    if not records:
        return records
    alignment_length = len(records[0].seq)
    if alignment_length == 0:
        return records
    keep = []
    n = len(records)
    for column in range(alignment_length):
        covered = sum(str(record.seq[column]).upper() in VALID_BASES for record in records)
        if covered / n >= min_coverage:
            keep.append(column)
    if not keep:
        return records
    keep_set = set(keep)
    trimmed = []
    for record in records:
        sequence = "".join(base for idx, base in enumerate(str(record.seq)) if idx in keep_set)
        trimmed.append(SeqRecord(Seq(sequence), id=record.id, description=""))
    return trimmed


def run_mafft(records, threads, large=False):
    """Align records with MAFFT and return aligned SeqRecord list."""
    with tempfile.TemporaryDirectory() as tmp:
        input_path = Path(tmp) / "input.fasta"
        SeqIO.write(records, input_path, "fasta")
        cmd = ["mafft", "--quiet", "--auto", "--thread", str(threads)]
        if large:
            cmd += ["--retree", "2"]
        cmd.append(str(input_path))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr[:500]}")
        from io import StringIO
        aligned = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
    return aligned


def alignment_p_distance(seq1, seq2):
    """Pairwise p-distance over columns where both bases are A/C/G/T."""
    matches = 0
    total = 0
    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        if base1 not in VALID_BASES or base2 not in VALID_BASES:
            continue
        total += 1
        if base1 == base2:
            matches += 1
    if total == 0:
        return None
    return 1.0 - (matches / total)


def nearest_consensus_metrics(alignment, raw_ids, consensus_ids):
    """For each raw record, find its nearest consensus by p-distance.

    Returns (per_raw_rows, summary_metrics, per_consensus_load).
    """
    seqs = {record.id: str(record.seq) for record in alignment}
    rows = []
    distances = []
    consensus_load = Counter()
    for raw_id in raw_ids:
        raw_seq = seqs.get(raw_id)
        if raw_seq is None:
            continue
        best_id = None
        best_distance = math.inf
        for consensus_id in consensus_ids:
            consensus_seq = seqs.get(consensus_id)
            if consensus_seq is None:
                continue
            distance = alignment_p_distance(raw_seq, consensus_seq)
            if distance is None:
                continue
            if distance < best_distance:
                best_distance = distance
                best_id = consensus_id
        if best_id is not None:
            consensus_load[best_id] += 1
            distances.append(best_distance)
            rows.append({
                "raw_id": raw_id,
                "nearest_consensus_id": best_id,
                "nearest_p_distance": f"{best_distance:.6f}",
            })

    if not distances:
        return rows, {}, consensus_load

    array = np.array(distances, dtype=float)
    metrics = {
        "raw_count": len(distances),
        "mean_nearest_distance": float(np.mean(array)),
        "median_nearest_distance": float(np.median(array)),
        "q90_nearest_distance": float(np.quantile(array, 0.90)),
        "q95_nearest_distance": float(np.quantile(array, 0.95)),
        "max_nearest_distance": float(np.max(array)),
        "coverage_d02": float(np.mean(array <= 0.02)),
        "coverage_d05": float(np.mean(array <= 0.05)),
        "coverage_d10": float(np.mean(array <= 0.10)),
    }
    for row in rows:
        row["consensus_raw_load"] = str(consensus_load[row["nearest_consensus_id"]])
    return rows, metrics, consensus_load


def analyze(raw_path, consensus_path, outdir, exclude_ids, threads):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    per_genotype_dir = outdir / "per_genotype"
    per_genotype_dir.mkdir(parents=True, exist_ok=True)

    raw_groups, _ = load_by_genotype(raw_path)
    consensus_groups, dropped = load_by_genotype(consensus_path, exclude=exclude_ids)

    print(f"Raw genotypes: {len(raw_groups)} | Consensus genotypes: {len(consensus_groups)} "
          f"(dropped {dropped} excluded IDs)")

    genotypes = sorted(set(raw_groups) | set(consensus_groups))
    summaries = []

    for genotype in genotypes:
        raw_records = raw_groups.get(genotype, [])
        consensus_records = consensus_groups.get(genotype, [])
        raw_ids = [r.id for r in raw_records]
        consensus_ids = [r.id for r in consensus_records]

        # Skip when there is nothing to compare.
        if not raw_records or not consensus_records:
            summaries.append({
                "genotype_pair": genotype,
                "raw_count": len(raw_records),
                "consensus_count": len(consensus_records),
                "status": "skipped_no_counterpart",
                "mean_nearest_distance": "",
                "median_nearest_distance": "",
                "q90_nearest_distance": "",
                "q95_nearest_distance": "",
                "max_nearest_distance": "",
                "coverage_d02": "",
                "coverage_d05": "",
                "coverage_d10": "",
            })
            continue

        combined = raw_records + consensus_records
        large = len(combined) > LARGE_GENOTYPE_THRESHOLD
        try:
            aligned = run_mafft(combined, threads, large=large)
            aligned = trim_alignment_terminal_columns(aligned, min_coverage=0.5)
        except Exception as exc:
            print(f"  [WARN] {genotype}: alignment failed ({exc}); skipping")
            summaries.append({
                "genotype_pair": genotype,
                "raw_count": len(raw_records),
                "consensus_count": len(consensus_records),
                "status": f"alignment_failed: {exc}",
                "mean_nearest_distance": "", "median_nearest_distance": "",
                "q90_nearest_distance": "", "q95_nearest_distance": "",
                "max_nearest_distance": "",
                "coverage_d02": "", "coverage_d05": "", "coverage_d10": "",
            })
            continue

        rows, metrics, _load = nearest_consensus_metrics(aligned, raw_ids, consensus_ids)
        if not metrics:
            summaries.append({
                "genotype_pair": genotype,
                "raw_count": len(raw_records),
                "consensus_count": len(consensus_records),
                "status": "no_measurable_pairs",
                "mean_nearest_distance": "", "median_nearest_distance": "",
                "q90_nearest_distance": "", "q95_nearest_distance": "",
                "max_nearest_distance": "",
                "coverage_d02": "", "coverage_d05": "", "coverage_d10": "",
            })
            continue

        # Per-genotype TSV
        geno_tsv = per_genotype_dir / f"{safe_label(genotype)}.tsv"
        with open(geno_tsv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=[
                "raw_id", "nearest_consensus_id", "nearest_p_distance",
                "consensus_raw_load"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        summaries.append({
            "genotype_pair": genotype,
            "raw_count": len(raw_records),
            "consensus_count": len(consensus_records),
            "status": "analysed",
            "mean_nearest_distance": f"{metrics['mean_nearest_distance']:.6f}",
            "median_nearest_distance": f"{metrics['median_nearest_distance']:.6f}",
            "q90_nearest_distance": f"{metrics['q90_nearest_distance']:.6f}",
            "q95_nearest_distance": f"{metrics['q95_nearest_distance']:.6f}",
            "max_nearest_distance": f"{metrics['max_nearest_distance']:.6f}",
            "coverage_d02": f"{metrics['coverage_d02']:.4f}",
            "coverage_d05": f"{metrics['coverage_d05']:.4f}",
            "coverage_d10": f"{metrics['coverage_d10']:.4f}",
        })
        print(f"  {genotype}: raw={len(raw_records)} cons={len(consensus_records)} "
              f"median={metrics['median_nearest_distance']:.4f} "
              f"cov5%={metrics['coverage_d05']:.3f}")

    # Global summary TSV
    summary_path = outdir / "summary.tsv"
    fieldnames = ["genotype_pair", "raw_count", "consensus_count", "status",
                  "mean_nearest_distance", "median_nearest_distance",
                  "q90_nearest_distance", "q95_nearest_distance",
                  "max_nearest_distance",
                  "coverage_d02", "coverage_d05", "coverage_d10"]
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    print(f"\nSummary written: {summary_path}")

    # Global aggregate across analysed genotypes (raw-count-weighted)
    analysed = [s for s in summaries if s["status"] == "analysed"]
    if analysed:
        all_distances = []
        weighted_sum = 0.0
        total_raw = 0
        cov_counts = {0.02: 0, 0.05: 0, 0.10: 0}
        for s in analysed:
            # re-read per-genotype distances for the global distribution
            geno_tsv = per_genotype_dir / f"{safe_label(s['genotype_pair'])}.tsv"
            for row in csv.DictReader(open(geno_tsv), delimiter="\t"):
                d = float(row["nearest_p_distance"])
                all_distances.append(d)
                for thresh in cov_counts:
                    if d <= thresh:
                        cov_counts[thresh] += 1
            total_raw += int(s["raw_count"])
        arr = np.array(all_distances)
        print("\n=== GLOBAL (raw-weighted across analysed genotypes) ===")
        print(f"  raw records measured: {len(arr)}")
        print(f"  mean   nearest p-distance: {np.mean(arr):.6f}")
        print(f"  median nearest p-distance: {np.median(arr):.6f}")
        print(f"  q90    nearest p-distance: {np.quantile(arr, 0.90):.6f}")
        print(f"  q95    nearest p-distance: {np.quantile(arr, 0.95):.6f}")
        print(f"  max    nearest p-distance: {np.max(arr):.6f}")
        print(f"  coverage <=2% : {cov_counts[0.02]/len(arr):.4f}")
        print(f"  coverage <=5% : {cov_counts[0.05]/len(arr):.4f}")
        print(f"  coverage <=10%: {cov_counts[0.10]/len(arr):.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="Lightweight nearest-neighbour representativeness of consensus "
                    "sequences vs filtered raw sequences, grouped by RdRp_VP1 genotype.")
    parser.add_argument("--raw", required=True, help="Filtered raw FASTA")
    parser.add_argument("--consensus", required=True, help="Filtered consensus FASTA")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--exclude", default=None,
                        help="File with consensus IDs to drop (one per line)")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    exclude_ids = []
    if args.exclude:
        with open(args.exclude) as fh:
            exclude_ids = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
        print(f"Excluding {len(exclude_ids)} consensus IDs")

    analyze(args.raw, args.consensus, args.outdir, exclude_ids, args.threads)


if __name__ == "__main__":
    main()
