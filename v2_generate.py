#!/usr/bin/env python3
"""Build a v2-style mapping-reference panel FROM result_4ref data.

Method (variant: by_rdrp_vp1):
  - result_4ref/{gi,gii}/04_alignments/aligned_*.fasta are ALREADY aligned per
    RdRp_VP1 cluster, with genotypes confirmed by the noro pipeline.
  - Within each cluster, greedy-dereplicate raw sequences at a similarity
    threshold (default 95%, matching the value calibrated against the original
    v2 panel for GII.17/GII.2/GII.3). Centroids prefer the most complete
    (max A/C/G/T) sequence. No MAFFT realignment needed (reuse the alignment).
  - Emit representatives in MAPREF_v2 header format, using the ACTUAL confirmed
    RdRp+VP1 (not the original v2's P0 placeholder).

Outputs to v2_test/by_rdrp_vp1/.
"""
import os, re, glob, sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REF4REF = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/result_4ref"
OUTDIR  = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/v2_test/by_rdrp_vp1"
THRESHOLD = 95.0
VALID = np.array([ord(b) for b in "ACGT"], dtype=np.uint8)

CLUSTER_RE = re.compile(r'aligned_cluster__cluster_(\d+)__(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)\.fasta$')

def valid_base_count(s):
    return sum(1 for c in s if c in "ACGT")

def greedy_derep(records, threshold):
    """records: list of SeqRecord (already aligned). Keep one per >=threshold cluster.
    Input is pre-sorted by A/C/G/T count desc so centroids are most complete."""
    if len(records) <= 1:
        return list(records)
    arr = np.array([[ord(c) for c in str(r.seq).upper()] for r in records], dtype=np.uint8)
    valid = np.isin(arr, VALID)  # (n, L) bool
    kept = [0]
    for i in range(1, len(records)):
        is_redundant = False
        for j in kept:
            v = valid[i] & valid[j]
            if v.sum() == 0:
                continue
            sim = 100.0 * ((arr[i] == arr[j]) & v).sum() / v.sum()
            if sim >= threshold:
                is_redundant = True
                break
        if not is_redundant:
            kept.append(i)
    return [records[i] for i in kept]

def parse_id(rec_id):
    """id like GII.P16_GII.4_AB123456_2020-03-18 -> (rdrp, vp1, rest)."""
    parts = rec_id.split('_', 2)
    if len(parts) < 3:
        # fallback: maybe just rdrp_vp1
        p = rec_id.split('_')
        return (p[0], p[1] if len(p)>1 else 'NA', rec_id)
    return (parts[0], parts[1], parts[2])

def build_group(group):
    aln_dir = os.path.join(REF4REF, group, "04_alignments")
    out_recs = []
    summary_rows = []  # (cluster_id, rdrp, vp1, raw_n, kept_n)
    for f in sorted(glob.glob(os.path.join(aln_dir, "aligned_*.fasta"))):
        m = CLUSTER_RE.search(os.path.basename(f))
        if not m:
            continue
        cluster_id, rdrp, vp1 = m.group(1), m.group(2), m.group(3)
        recs = list(SeqIO.parse(f, "fasta"))
        raw_n = len(recs)
        # strip alignment gaps to store RAW sequence for the reference output
        raw_recs = []
        for r in recs:
            raw = str(r.seq).replace('-', '').upper()
            raw_recs.append((r, raw))
        if raw_n == 0:
            continue
        # derep on the ALIGNED form (use alignment for correct column-wise similarity)
        # but order by raw ACGT completeness desc
        order = sorted(range(len(recs)), key=lambda i: valid_base_count(raw_recs[i][1]), reverse=True)
        ordered = [recs[i] for i in order]
        kept = greedy_derep(ordered, THRESHOLD)
        # emit kept sequences (raw, ungapped) with MAPREF_v2 header
        refn = 0
        for r in kept:
            refn += 1
            raw = str(r.seq).replace('-', '').upper()
            _, vp1_id, rest = parse_id(r.id)
            new_id = f"{rdrp}_{vp1}_MAPREF_v2_ref{refn:03d}_{vp1}_{rest}"
            out_recs.append(SeqRecord(Seq(raw), id=new_id, description=''))
        summary_rows.append((cluster_id, rdrp, vp1, raw_n, len(kept)))
    return out_recs, summary_rows

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    grand = {}
    for group in ("gi", "gii"):
        recs, rows = build_group(group)
        out_fa = os.path.join(OUTDIR, f"{group}_mapping_references_v2.fasta")
        SeqIO.write(recs, out_fa, "fasta")
        # summary tsv
        tsv = os.path.join(OUTDIR, f"{group}_summary.tsv")
        with open(tsv, 'w') as out:
            out.write("cluster_id\trdrp\tvp1\traw_n\tkept_n\n")
            for cid, rdrp, vp1, rn, kn in rows:
                out.write(f"{cid}\t{rdrp}\t{vp1}\t{rn}\t{kn}\n")
        grand[group] = (recs, rows)
        tot_raw = sum(r[3] for r in rows)
        tot_kept = sum(r[4] for r in rows)
        print(f"{group}: clusters={len(rows)}  raw={tot_raw}  kept(MAPREF)={tot_kept}  -> {out_fa}")
    # combined report
    rep = os.path.join(OUTDIR, "report.txt")
    with open(rep, 'w') as out:
        out.write("v2_test (by_rdrp_vp1, 95%% derep from result_4ref) report\n")
        out.write("="*60+"\n")
        for group in ("gi","gii"):
            recs, rows = grand[group]
            tot_raw = sum(r[3] for r in rows)
            tot_kept = sum(r[4] for r in rows)
            out.write(f"\n{group.upper()}: {len(rows)} clusters, raw={tot_raw}, kept={tot_kept}\n")
            out.write("  per VP1 kept counts:\n")
            vp1_kept = {}
            for _,_,vp1,_,kn in rows:
                vp1_kept[vp1] = vp1_kept.get(vp1,0)+kn
            for vp1 in sorted(vp1_kept):
                out.write(f"    {vp1:10s} {vp1_kept[vp1]}\n")
        out.write("\nComparison to original v2 panel (built from sewage_consensus Rawdata):\n")
        out.write("  original v2: GI=142, GII=288\n")
        out.write(f"  this v2_test: GI={len(grand['gi'][0])}, GII={len(grand['gii'][0])}\n")
        out.write("\nNote: different rawdata pool (result_4ref Rawdata vs sewage Rawdata_ncbi)\n")
        out.write("and per-RdRp_VP1 grouping (vs original v2's VP1-only P0 grouping),\n")
        out.write("so counts need not match the original 142/288.\n")
    print(f"\nreport -> {rep}")

if __name__ == "__main__":
    main()
