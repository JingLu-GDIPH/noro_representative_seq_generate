#!/usr/bin/env python3
"""QC analysis of raw-sequence representatives against their genotype clade root.

Motivation
----------
Some final representative sequences are raw input records that were kept as their
own representative because they differed too much from sibling sequences to be
merged into a consensus. Such a record may either be a genuinely divergent
lineage or a sequence with assembly/sequencing errors (the reason nextstrain/ncov
filters by private-mutation count and SNP clusters via Nextclade).

Without a pre-built norovirus reference tree or clade label map (the way Nextclade
uses for SARS-CoV-2), we approximate the "clade root" for each genotype by the
**majority consensus** of all sibling raw sequences aligned in that genotype's
cluster. For every raw representative we then count, relative to that clade root:

  * ``private_mutations``  — A/C/G/T positions where the representative differs
    from the clade-root consensus (gap/N positions excluded).
  * ``snp_clusters``       — number of 101 nt sliding windows that contain
    *strictly more than* 6 private substitutions (Nextclade's SNP-cluster rule).
    Multiple adjacent/overlapping windows merge into one cluster, so this is the
    cluster count (not the window count).

A representative with many private mutations and/or ``snp_clusters > 1`` is flagged
as a potential sequencing-error suspect, mirroring the nextstrain/ncov
``diagnostic.py`` logic (``snp_clusters > 1`` is an exclusion criterion there).

Inputs
------
* ``--representatives``  : the final consensus FASTA (e.g. ``all_final_consensus``).
* ``--clusters-dir``     : directory of per-cluster aligned FASTA files produced
  by the ALIGN step (filenames like
  ``aligned_cluster__cluster_NNNN__<RdRp_VP1>.fasta``).
* ``--output``           : TSV report, one row per raw representative.

Notes
-----
* Only records whose ID matches a raw accession pattern
  (``<RdRp>_<VP1>_<ACCESSION>_<date>``) are analysed; node/consensus
  representatives are skipped because they are already consensus-derived.
* The clade-root consensus is recomputed from the full cluster alignment (all
  siblings, not only representatives) so the reference is stable and independent
  of which records happened to survive pruning.
* SNP-cluster detection uses private *substitutions* only (A/C/G/T vs A/C/G/T),
  exactly like current Nextclade (deletions are not counted, per v1.10.2).
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#: A raw input record ID looks like ``GII.P4_GII.4_GQ845368_2007`` or
#: ``GII.P7_GII.6_KX158282_2015`` — i.e. the 3rd underscore field is a nucleotide
#: accession (1-2 letters then 4-6 digits) and the 4th is a date.
RAW_ID_PATTERN = re.compile(
    r"^(G(?:I|II|IX)\.P[A-Za-z0-9]+_G(?:I|II|IX)\.[A-Za-z0-9]+)_([A-Z]{1,3}\d{4,8})_"
)
GENOTYPE_PAIR_PATTERN = re.compile(
    r"(G(?:I|II|IX)\.P[A-Za-z0-9]+_G(?:I|II|IX)\.[A-Za-z0-9]+)"
)

VALID_BASES = frozenset("ACGT")

#: Nextclade SNP-cluster rule defaults (SARS-CoV-2).
SNP_CLUSTER_WINDOW = 101
SNP_CLUSTER_CUTOFF = 6  # strictly greater than -> a cluster


def is_raw_representative(sequence_id):
    """True if the ID looks like an original raw accession record."""
    return bool(RAW_ID_PATTERN.match(str(sequence_id)))


def extract_genotype_pair(sequence_id):
    match = GENOTYPE_PAIR_PATTERN.search(str(sequence_id))
    return match.group(1) if match else None


def majority_consensus(records):
    """Return the column-wise A/C/G/T majority sequence of aligned ``records``.

    Columns with no A/C/G/T base become ``N``; ties resolve to the
    lexicographically smallest base. The returned string has the alignment
    length (gaps preserved as alignment columns, but never emitted as a
    consensus call).
    """
    if not records:
        return ""
    alignment_length = len(records[0].seq)
    consensus = []
    for column in range(alignment_length):
        bases = [
            str(record.seq[column]).upper()
            for record in records
            if str(record.seq[column]).upper() in VALID_BASES
        ]
        if not bases:
            consensus.append("N")
        else:
            # most_common(1) returns a single (base, count) tuple; tie-break by base
            consensus.append(Counter(bases).most_common(1)[0][0])
    return "".join(consensus)


def count_snp_clusters(private_positions, window=SNP_CLUSTER_WINDOW,
                       cutoff=SNP_CLUSTER_CUTOFF):
    """Count Nextclade-style SNP clusters from sorted private-SNP positions.

    A cluster exists wherever some ``window``-nt span contains *more than*
    ``cutoff`` private substitutions. Adjacent/overlapping clusters merge into
    one. Returns the integer cluster count.
    """
    if not private_positions:
        return 0
    positions = sorted(private_positions)
    cluster_count = 0
    in_cluster = False
    # deque-style window: indices into positions that fall within [pos-window+1, pos]
    left = 0
    for right in range(len(positions)):
        pos = positions[right]
        while positions[left] < pos - window + 1:
            left += 1
        window_size = right - left + 1
        if window_size > cutoff:
            if not in_cluster:
                cluster_count += 1
                in_cluster = True
        else:
            in_cluster = False
    return cluster_count


def align_pairwise(reference_seq, query_seq):
    """Globally align two unaligned sequences and return them as equal-length strings."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference_seq, query_seq)[0]
    aligned_ref, aligned_query = alignment[0], alignment[1]
    return str(aligned_ref), str(aligned_query)


def compute_private_mutations(ref_aligned, query_aligned):
    """Return list of 1-based private-substitution positions in ``query`` vs ``ref``.

    A private substitution is a column where both bases are A/C/G/T and they
    differ. Positions are 1-based coordinates in the pairwise alignment.
    """
    positions = []
    for i, (ref_base, query_base) in enumerate(zip(ref_aligned, query_aligned), start=1):
        ref_base = ref_base.upper()
        query_base = query_base.upper()
        if ref_base in VALID_BASES and query_base in VALID_BASES and ref_base != query_base:
            positions.append(i)
    return positions


def build_clade_roots(clusters_dir):
    """Map genotype pair -> (majority-consensus clade root, cluster size).

    The cluster size is returned so callers can tell singleton clusters apart:
    a singleton's clade root equals the only record itself, so its private
    mutations are trivially 0 — that does NOT mean QC-passing, it means there is
    no reference to compare against.
    """
    genotype_to_root = {}
    cluster_files = sorted(Path(clusters_dir).glob("aligned_*.fasta"))
    if not cluster_files:
        # also accept unaligned per-cluster files (fallback)
        cluster_files = sorted(Path(clusters_dir).glob("cluster__*.fasta"))
    for cluster_file in cluster_files:
        records = list(SeqIO.parse(cluster_file, "fasta"))
        if not records:
            continue
        genotypes = {extract_genotype_pair(r.id) for r in records}
        genotypes.discard(None)
        if len(genotypes) != 1:
            # ambiguous mapping; skip (should not happen for VP1-stratified clusters)
            continue
        genotype = genotypes.pop()
        # Only equal-length (aligned) records give a meaningful column consensus.
        lengths = {len(r.seq) for r in records}
        if len(lengths) == 1:
            genotype_to_root[genotype] = (majority_consensus(records), len(records))
    return genotype_to_root


def analyse(representatives_path, clusters_dir, output_path):
    representatives = list(SeqIO.parse(representatives_path, "fasta"))
    if not representatives:
        raise ValueError(f"No representatives found in {representatives_path}")

    clade_roots = build_clade_roots(clusters_dir)
    print(f"Loaded {len(clade_roots)} genotype clade roots from {clusters_dir}")

    raw_reps = [r for r in representatives if is_raw_representative(r.id)]
    skipped = [r for r in representatives if not is_raw_representative(r.id)]
    print(f"Analysing {len(raw_reps)} raw representatives "
          f"(skipping {len(skipped)} node/consensus representatives)")

    rows = []
    for rep in raw_reps:
        genotype = extract_genotype_pair(rep.id)
        entry = clade_roots.get(genotype)
        if not entry:
            rows.append({
                "representative_id": rep.id,
                "genotype_pair": genotype or "NA",
                "cluster_size": "",
                "status": "no_clade_root",
                "private_mutations": "",
                "snp_clusters": "",
                "alignment_length": "",
                "flag": "NO_REFERENCE",
            })
            continue

        root, cluster_size = entry

        # Singleton cluster: the clade root IS this sequence, so private mutations
        # are trivially 0. That is NOT a QC pass — there is simply no reference to
        # compare against. Mark explicitly so it is not misread as "clean".
        if cluster_size <= 1:
            rows.append({
                "representative_id": rep.id,
                "genotype_pair": genotype,
                "cluster_size": cluster_size,
                "status": "singleton",
                "private_mutations": "",
                "snp_clusters": "",
                "alignment_length": "",
                "flag": "SINGLETON_NO_REFERENCE",
            })
            continue

        # Align the unaligned raw rep sequence to the clade-root consensus.
        rep_seq = str(rep.seq).replace("-", "").upper()
        ref_aligned, query_aligned = align_pairwise(root.replace("-", "").upper(), rep_seq)
        private_positions = compute_private_mutations(ref_aligned, query_aligned)
        snp_clusters = count_snp_clusters(private_positions)

        # Flag mirrors nextstrain/ncov diagnostic thresholds (reference defaults).
        flags = []
        if snp_clusters > 1:
            flags.append("SNP_CLUSTERS>1")
        # ncov uses contamination=5 for reversions+labeled; without labels we use a
        # generous private-mutation percentile-based heuristic only as a soft flag.
        rows.append({
            "representative_id": rep.id,
            "genotype_pair": genotype,
            "cluster_size": cluster_size,
            "status": "analysed",
            "private_mutations": len(private_positions),
            "snp_clusters": snp_clusters,
            "alignment_length": len(ref_aligned),
            "flag": ";".join(flags) if flags else "",
        })

    # Soft-flag high private-mutation counts (>P95 within its genotype, and >=10).
    by_geno = {}
    for row in rows:
        if row["status"] != "analysed":
            continue
        by_geno.setdefault(row["genotype_pair"], []).append(row)
    for geno, group_rows in by_geno.items():
        counts = sorted(int(r["private_mutations"]) for r in group_rows)
        if len(counts) >= 4:
            # 90th percentile within genotype
            p90 = counts[int(0.9 * (len(counts) - 1))]
            threshold = max(p90, 10)
            for row in group_rows:
                if int(row["private_mutations"]) >= threshold:
                    extra = f"HIGH_PRIVATE(>={threshold})"
                    row["flag"] = (row["flag"] + ";" + extra) if row["flag"] else extra

    fieldnames = [
        "representative_id",
        "genotype_pair",
        "cluster_size",
        "status",
        "private_mutations",
        "snp_clusters",
        "alignment_length",
        "flag",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    flagged = [r for r in rows if r["flag"] and r["flag"] != "SINGLETON_NO_REFERENCE"]
    singletons = [r for r in rows if r.get("status") == "singleton"]
    print(f"\nWrote {len(rows)} rows to {output_path}")
    print(f"Flagged (potential sequencing-error suspects): {len(flagged)}")
    for row in flagged:
        print(f"  {row['representative_id']}  priv={row['private_mutations']} "
              f"clusters={row['snp_clusters']}  flag={row['flag']}")
    print(f"Singleton clusters (no reference, cannot judge): {len(singletons)}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Count private mutations and SNP clusters of raw-sequence "
            "representatives vs their genotype clade-root consensus, to flag "
            "potential sequencing-error suspects (nextclade-style)."
        )
    )
    parser.add_argument("--representatives", required=True,
                        help="Final representative FASTA (e.g. all_final_consensus.fasta)")
    parser.add_argument("--clusters-dir", required=True,
                        help="Directory of per-cluster aligned FASTA files")
    parser.add_argument("--output", required=True, help="Output TSV report")
    args = parser.parse_args()

    if not Path(args.clusters_dir).is_dir():
        sys.exit(f"ERROR: clusters dir not found: {args.clusters_dir}")

    analyse(args.representatives, args.clusters_dir, args.output)


if __name__ == "__main__":
    main()
