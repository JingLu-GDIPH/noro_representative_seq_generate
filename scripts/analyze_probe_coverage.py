#!/usr/bin/env python3
"""Analyze probe-set coverage and depth per genotype.

The probe-design FASTA files (``probe_design/T3339XV1_{GI,GII}.probe.fa``)
contain 100 bp sliding-window probes tiling across each representative genome.
Each probe header is ``<VP1genotype>_<rep_id>:<start>-<end> <dup> GC=.. Penalty=..``.

This script reports, per VP1 genotype:
  * number of probes (probe count = tiling depth × breadth)
  * number of distinct representative source genomes probed
  * total bp spanned by probes (union of probe intervals on each rep)
  * mean per-position probe depth (how many probes overlap each base on average)
  * breadth: fraction of the representative genome length covered by probes

Inputs
------
* ``--probe``     : probe FASTA
* ``--representatives`` : the representative genomes the probes were designed from
  (used to get each rep's full length, for breadth calculation)
* ``--outdir``    : output directory
"""

import argparse
import csv
import re
import statistics
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


#: ``>GI.1_cluster_3_node_1:0-100 ...`` -> genotype=GI.1, source=GI.1_cluster_3_node_1, start=0, end=100
PROBE_HEADER = re.compile(r"^(?P<source>[^:]+):(?P<start>\d+)-(?P<end>\d+)")

#: VP1 genotype is the first ``_`` -delimited field of the source.
def genotype_of(source_id):
    return source_id.split("_")[0]


def parse_probes(probe_fasta):
    """Return list of (genotype, source_id, start, end) per probe."""
    probes = []
    for record in SeqIO.parse(probe_fasta, "fasta"):
        header = record.description
        m = PROBE_HEADER.search(header)
        if not m:
            continue
        source = m.group("source")
        start = int(m.group("start"))
        end = int(m.group("end"))
        geno = genotype_of(source)
        probes.append((geno, source, start, end))
    return probes


def union_length(intervals):
    """Length of the union of inclusive [start, end) half-open intervals."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_end:
            cur_end = max(cur_end, e)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = s, e
    total += cur_end - cur_start
    return total


def covered_depth(intervals, genome_length):
    """Mean per-base probe depth over [0, genome_length)."""
    if genome_length <= 0:
        return 0.0
    events = []
    for s, e in intervals:
        events.append((s, 1))
        events.append((e, -1))
    events.sort()
    depth_area = 0
    cur_depth = 0
    prev_pos = 0
    for pos, delta in events:
        pos = max(0, min(pos, genome_length))
        if pos > prev_pos:
            depth_area += cur_depth * (pos - prev_pos)
        cur_depth += delta
        prev_pos = pos
    return depth_area / genome_length


def analyze(probe_fasta, representatives_fasta, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    probes = parse_probes(probe_fasta)
    # Rep length: prefer the provided representative FASTA, otherwise infer the
    # genome length from the max probe end coordinate (probes tile to the end).
    rep_lengths = {}
    if representatives_fasta and Path(representatives_fasta).exists():
        rep_lengths = {rec.id: len(rec.seq)
                       for rec in SeqIO.parse(representatives_fasta, "fasta")}
    # Infer lengths for any source not in the representative FASTA (probes may
    # have been designed from an earlier run with different IDs).
    max_end_by_source = defaultdict(int)
    for _geno, source, _start, end in probes:
        max_end_by_source[source] = max(max_end_by_source[source], end)
    for source, max_end in max_end_by_source.items():
        if source not in rep_lengths or not rep_lengths.get(source):
            rep_lengths[source] = max_end

    # group probes by (genotype, source)
    by_source = defaultdict(list)
    for geno, source, start, end in probes:
        by_source[(geno, source)].append((start, end))

    by_genotype = defaultdict(list)
    for geno, source in by_source:
        by_genotype[geno].append(source)

    rows = []
    for geno in sorted(by_genotype):
        sources = by_genotype[geno]
        geno_probes = [probe for probe in probes if probe[0] == geno]
        # per-source breadth + depth
        source_rows = []
        total_span = 0
        total_depth_area = 0
        for source in sources:
            intervals = by_source[(geno, source)]
            g_len = rep_lengths.get(source, 0)
            cov = union_length(intervals)
            depth = covered_depth(intervals, g_len) if g_len else 0.0
            breadth = cov / g_len if g_len else 0.0
            source_rows.append({
                "genotype": geno, "source_id": source,
                "probe_count": len(intervals),
                "genome_length": g_len,
                "covered_bp": cov,
                "breadth_pct": f"{breadth * 100:.2f}",
                "mean_depth": f"{depth:.2f}",
            })
            total_span += cov
            total_depth_area += depth * (g_len or 0)

        # genotype aggregate
        all_intervals_flat = [len(by_source[(geno, s)]) for s in sources]
        total_genome_bp = sum(rep_lengths.get(s, 0) for s in sources)
        agg_depth = total_depth_area / total_genome_bp if total_genome_bp else 0.0
        rows.append({
            "genotype": geno,
            "representative_count": len(sources),
            "probe_count": len(geno_probes),
            "total_genome_bp": total_genome_bp,
            "mean_probes_per_rep": f"{len(geno_probes)/len(sources):.1f}" if sources else "0",
            "mean_depth": f"{agg_depth:.2f}",
        })

    # write per-genotype summary
    geno_path = outdir / "probe_genotype_summary.tsv"
    with open(geno_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "genotype", "representative_count", "probe_count",
            "total_genome_bp", "mean_probes_per_rep", "mean_depth"],
            delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # write per-source detail
    detail_path = outdir / "probe_source_detail.tsv"
    all_source_rows = []
    for geno in sorted(by_genotype):
        for source in by_genotype[geno]:
            intervals = by_source[(geno, source)]
            g_len = rep_lengths.get(source, 0)
            cov = union_length(intervals)
            depth = covered_depth(intervals, g_len) if g_len else 0.0
            breadth = cov / g_len if g_len else 0.0
            all_source_rows.append({
                "genotype": geno, "source_id": source,
                "probe_count": len(intervals),
                "genome_length": g_len,
                "covered_bp": cov,
                "breadth_pct": f"{breadth * 100:.2f}",
                "mean_depth": f"{depth:.2f}",
            })
    with open(detail_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "genotype", "source_id", "probe_count", "genome_length",
            "covered_bp", "breadth_pct", "mean_depth"], delimiter="\t")
        w.writeheader()
        w.writerows(all_source_rows)

    # console report
    print(f"{'genotype':10} {'reps':>5} {'probes':>7} {'genome_bp':>10} "
          f"{'prob/rep':>9} {'mean_depth':>10}")
    for r in rows:
        print(f"{r['genotype']:10} {r['representative_count']:>5} {r['probe_count']:>7} "
              f"{r['total_genome_bp']:>10} {r['mean_probes_per_rep']:>9} {r['mean_depth']:>10}")
    print(f"\nPer-genotype summary: {geno_path}")
    print(f"Per-source detail:    {detail_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Report probe-set coverage breadth and depth per VP1 genotype.")
    parser.add_argument("--probe", required=True, help="Probe FASTA")
    parser.add_argument("--representatives", default=None,
                        help="Representative genomes the probes tile (optional; "
                             "lengths are inferred from probe coordinates if absent)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()
    analyze(args.probe, args.representatives, args.outdir)


if __name__ == "__main__":
    main()
