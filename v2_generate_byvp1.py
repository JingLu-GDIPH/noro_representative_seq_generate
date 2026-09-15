#!/usr/bin/env python3
"""v2_test variant B: faithful VP1-only grouping (like original v2).
Merge raw sequences across all RdRp clusters of the same VP1, MAFFT-realign,
greedy-derep at 95%. Output MAPREF_v2 (RdRp kept as actual cluster RdRp of the
centroid, since result_4ref genotypes are confirmed)."""
import os, re, glob, sys, subprocess
import numpy as np
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REF4REF = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/result_4ref"
OUTDIR  = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/v2_test/by_vp1_95"
THRESHOLD = 95.0
VALID = np.array([ord(b) for b in "ACGT"], dtype=np.uint8)
CLUSTER_RE = re.compile(r'aligned_cluster__cluster_(\d+)__(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)\.fasta$')

def acgt_count(s): return sum(1 for c in s if c in "ACGT")

def greedy_derep_aligned(records, threshold):
    if len(records) <= 1: return list(records)
    arr = np.array([[ord(c) for c in str(r.seq).upper()] for r in records], dtype=np.uint8)
    valid = np.isin(arr, VALID)
    kept = [0]
    for i in range(1, len(records)):
        redundant = False
        for j in kept:
            v = valid[i] & valid[j]
            if v.sum() == 0: continue
            if 100.0*((arr[i]==arr[j])&v).sum()/v.sum() >= threshold:
                redundant = True; break
        if not redundant: kept.append(i)
    return [records[i] for i in kept]

def mafft(records, tag, threads=6):
    fa = f"/tmp/v2test/vp1_{tag}.fasta"
    SeqIO.write(records, fa, "fasta")
    n = len(records)
    args = ["mafft","--quiet","--thread",str(threads)]
    if n > 300:
        args += ["--retree","2"]  # fast for large groups
    args.append(fa)
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"MAFFT failed for {tag}: {r.stderr[:300]}\n")
        return None
    return list(SeqIO.parse(StringIO(r.stdout), "fasta"))

def collect_vp1(group):
    """VP1 -> list of (raw_seq, rdrp, original_id)."""
    aln_dir = os.path.join(REF4REF, group, "04_alignments")
    by_vp1 = {}
    for f in sorted(glob.glob(os.path.join(aln_dir, "aligned_*.fasta"))):
        m = CLUSTER_RE.search(os.path.basename(f))
        if not m: continue
        rdrp, vp1 = m.group(2), m.group(3)
        for r in SeqIO.parse(f, "fasta"):
            raw = str(r.seq).replace('-','').upper()
            by_vp1.setdefault(vp1, []).append((raw, rdrp, r.id))
    return by_vp1

def build_group(group):
    by_vp1 = collect_vp1(group)
    out_recs = []
    rows = []
    for vp1 in sorted(by_vp1):
        items = by_vp1[vp1]
        # order by completeness desc
        items.sort(key=lambda t: acgt_count(t[0]), reverse=True)
        if len(items) == 1:
            kept_items = items
        else:
            recs = [SeqRecord(Seq(t[0]), id=t[2], description='') for t in items]
            aln = mafft(recs, f"{group}_{vp1.replace('.','')}")
            if aln is None:
                kept_items = items[:1]  # fallback
            else:
                kept = greedy_derep_aligned(aln, THRESHOLD)
                kept_ids = {r.id for r in kept}
                kept_items = [t for t in items if t[2] in kept_ids]
        refn = 0
        for raw, rdrp, oid in kept_items:
            refn += 1
            parts = oid.split('_', 2)
            rest = parts[2] if len(parts) > 2 else oid
            new_id = f"{rdrp}_{vp1}_MAPREF_v2_ref{refn:03d}_{vp1}_{rest}"
            out_recs.append(SeqRecord(Seq(raw), id=new_id, description=''))
        rows.append((vp1, len(items), len(kept_items)))
    return out_recs, rows

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    grand = {}
    for group in ("gi","gii"):
        recs, rows = build_group(group)
        SeqIO.write(recs, os.path.join(OUTDIR, f"{group}_mapping_references_v2.fasta"), "fasta")
        grand[group] = (recs, rows)
        print(f"{group}: VP1 groups={len(rows)}  raw={sum(r[1] for r in rows)}  kept={sum(r[2] for r in rows)}")
    rep = os.path.join(OUTDIR, "report.txt")
    with open(rep,'w') as out:
        out.write("v2_test (by_vp1_95, faithful VP1 grouping, 95%% derep) report\n"+"="*60+"\n")
        for group in ("gi","gii"):
            recs,rows = grand[group]
            out.write(f"\n{group.upper()}: kept={len(recs)}\n  per VP1 (raw->kept):\n")
            for vp1,rn,kn in rows:
                out.write(f"    {vp1:10s} {rn:>5} -> {kn}\n")
        out.write(f"\noriginal v2: GI=142, GII=288\n")
        out.write(f"this by_vp1_95: GI={len(grand['gi'][0])}, GII={len(grand['gii'][0])}\n")
    print("report ->", rep)

if __name__ == "__main__":
    main()
