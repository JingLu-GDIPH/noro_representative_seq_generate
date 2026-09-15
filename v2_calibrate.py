#!/usr/bin/env python3
"""Calibrate whether v2 panel = greedy similarity dereplication per VP1.
Reads result_4ref aligned clusters (genotype-confirmed), groups raw seqs by VP1,
MAFFT-realigns, greedy-dereps at several thresholds, compares counts to v2."""
import os, re, glob, itertools
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess, sys

REF4REF = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/result_4ref"
V2_GI = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/sewage_consensus/results/gi_mapping_references_v2.fasta"
V2_GII = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/sewage_consensus/results/gii_mapping_references_v2.fasta"
VALID = {65,67,71,84}  # A,C,G,T bytes

VP1_RE = re.compile(r'__(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)\.fasta$')

def vp1_from_header(h):
    # header like GII.P16_GII.4_AB123_2020  -> VP1 = second field
    parts = h.lstrip('>').split('_')
    if len(parts) >= 2:
        return parts[1]
    return None

def collect_vp1_raw(group, vp1):
    """Collect ungapped raw sequences for a VP1 from all aligned clusters of that VP1."""
    aln_dir = os.path.join(REF4REF, group, "04_alignments")
    out = {}
    for f in glob.glob(os.path.join(aln_dir, "aligned_*.fasta")):
        m = re.search(r'__(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)\.fasta$', os.path.basename(f))
        if not m: continue
        file_vp1 = m.group(2)
        if file_vp1 != vp1: continue
        for r in SeqIO.parse(f, "fasta"):
            raw = str(r.seq).replace('-', '').upper()
            # skip sequences too short (fragments)
            rid = r.id
            out[rid] = raw
    return out

def v2_accessions_for_vp1(v2_path, vp1):
    accs = []
    for r in SeqIO.parse(v2_path, "fasta"):
        # >GII.P0_GII.1_MAPREF_v2_ref001_GII.1_OR397767_Feb-2021
        parts = r.id.split('_')
        if len(parts) >= 4 and parts[1] == vp1:
            # accession is parts[5] typically
            accs.append(r.id)
    return accs

def mafft_align(records, tag, threads=4):
    fa = f"/tmp/v2test/{tag}.fasta"
    SeqIO.write(records, fa, "fasta")
    r = subprocess.run(["mafft","--quiet","--auto","--thread",str(threads),fa],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(r.stderr[:500])
    return list(SeqIO.parse(r.stdout.splitlines().__iter__(), "fasta")) if False else \
           _parse_text(r.stdout)

def _parse_text(text):
    from io import StringIO
    return list(SeqIO.parse(StringIO(text), "fasta"))

def greedy_derep(aln_records, threshold):
    """Greedy similarity dereplication. Keep a record if its max similarity to all
    kept records is < threshold (strictly). Returns list of kept records."""
    arr = np.array([[ord(c) for c in str(r.seq).upper()] for r in aln_records], dtype=np.uint8)
    valid = np.isin(arr, list(VALID))
    kept_idx = []
    for i in range(len(aln_records)):
        ok = True
        for j in kept_idx:
            v = valid[i] & valid[j]
            if v.sum() == 0:
                continue
            eq = ((arr[i]==arr[j]) & v).sum()
            sim = 100.0*eq/v.sum()
            if sim >= threshold:
                ok = False
                break
        if ok:
            kept_idx.append(i)
    return [aln_records[i] for i in kept_idx]

targets = [("gii", V2_GII, [("GII.1",20),("GII.17",14),("GII.2",12),("GII.3",12)]),
           ("gi",  V2_GI,  [("GI.1",9),("GI.2",7)])]

print(f"{'group':5} {'VP1':7} {'raw':>5} {'v2':>4} | " + " ".join(f"{t:>6.1f}%" for t in [99.5,99.0,98.0,97.0,95.0,90.0]) + " | best")
print("-"*90)
for group, v2_path, vp1_targets in targets:
    for vp1, v2_n in vp1_targets:
        raw = collect_vp1_raw(group, vp1)
        if len(raw) < 2:
            print(f"{group:5} {vp1:7} {len(raw):>5} {v2_n:>4} | (too few raw)")
            continue
        recs = [SeqIO.SeqRecord(Seq(s), id=k, description='') for k,s in raw.items()]
        # filter very short
        recs = [r for r in recs if len(str(r.seq))>5000]
        aln = mafft_align(recs, f"{group}_{vp1.replace('.','')}")
        if len(aln) < 2:
            print(f"{group:5} {vp1:7} {len(raw):>5} {v2_n:>4} | align failed ({len(aln)})")
            continue
        res = []
        for t in [99.5,99.0,98.0,97.0,95.0,90.0]:
            kept = greedy_derep(aln, t)
            res.append(len(kept))
        # best threshold = the one closest to v2_n
        best_t = None; best_diff=1e9
        for t,n in zip([99.5,99.0,98.0,97.0,95.0,90.0], res):
            d = abs(n-v2_n)
            if d < best_diff: best_diff=d; best_t=t
        print(f"{group:5} {vp1:7} {len(raw):>5} {v2_n:>4} | " + " ".join(f"{n:>6}" for n in res) + f" | {best_t}% (diff={best_diff})")
