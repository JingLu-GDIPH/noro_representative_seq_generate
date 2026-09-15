#!/usr/bin/env python3
"""Calibrate max-within-max and target_mapping_rate from labeled GII.4/GII.17 variants.

A) Distance window: within-variant max pairwise distance (must be BELOW threshold
   to keep a variant whole) vs adjacent-variant max cross-pair distance (must be ABOVE
   threshold to force splitting a merged clade).
B) Merged-unit target rate: build 2-variant merged consensus refs and compute member
   target_mapping_rate (direct slice comparison, 135 comparable + <=10% mismatch)
   to see whether the target criterion can discriminate variant granularity.
"""
import csv
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np
from Bio import SeqIO

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
G4 = PROJ / "v2_test/gii4_threshold_calibration/aligned_trimmed.fasta"
G17 = PROJ / "v2_test/gii17_threshold_calibration_B2000/aligned_trimmed.fasta"
G17_FULL = PROJ / "v2_test/gii17_threshold_calibration/aligned_trimmed.fasta"
G17_TSV = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/Phy_nextstrain_nf/data/GII17_genome/GII17.genome.assigned.tsv")
BASES = frozenset("ACGT")
W = 150

def load_g4():
    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(G4, "fasta")}
    groups = defaultdict(list)
    for k in aln:
        groups[k.split("|")[0]].append(k)
    return aln, groups

def load_g17():
    sc, dates = {}, {}
    for row in csv.DictReader(open(G17_TSV), delimiter="\t"):
        s, c = row["strain"].strip(), (row["clade"] or "").strip()
        if s and c:
            sc[s] = c; dates[s] = (row["date"] or "")[:4]
    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(G17, "fasta")}
    for r in SeqIO.parse(G17_FULL, "fasta"):          # add pre-2000 B members
        if sc.get(r.id) == "B" and not (dates[r.id].isdigit() and int(dates[r.id]) >= 2000):
            aln.setdefault(r.id, str(r.seq).upper())
    groups = defaultdict(list)
    for k in aln:
        groups[sc[k]].append(k)
    return aln, groups

def pdist(a, b):
    ok = [(x, y) for x, y in zip(a, b) if x in BASES and y in BASES]
    return sum(x != y for x, y in ok) / len(ok) if len(ok) > 500 else None  # 部分基因组重叠太少不算

def pairwise_within(aln, members):
    ds = [d for d in (pdist(aln[a], aln[b]) for a, b in combinations(members, 2)) if d is not None]
    return ds

def cross_stats(aln, m1, m2):
    ds = [d for d in (pdist(aln[a], aln[b]) for a in m1 for b in m2) if d is not None]
    return ds

def consensus(aln, members):
    out = []
    for col in zip(*[aln[m] for m in members]):
        b = [x for x in col if x in BASES]
        out.append(max(set(b), key=b.count) if b else "N")
    return "".join(out)

def target_rate(aln, members, ref):
    tot = tm = 0
    for m in members:
        s = aln[m]
        for st in range(0, len(s) - W + 1, 25):
            frag = s[st:st + W]
            if any(c not in BASES for c in frag):
                continue
            ts = ref[st:st + W]
            comp = [(a, b) for a, b in zip(frag, ts) if b in BASES]
            tot += 1
            if len(comp) >= 135 and sum(a != b for a, b in comp) / len(comp) <= 0.10:
                tm += 1
    return tm / tot if tot else float("nan"), tot

for tag, loader in (("GII.4", load_g4), ("GII.17", load_g17)):
    aln, groups = loader()
    print(f"\n========== {tag} ==========")
    print("--- 组内 pairwise 距离（n, 对数, max, d95）---")
    within_max = {}
    for g in sorted(groups, key=lambda g: -len(groups[g])):
        ds = pairwise_within(aln, groups[g])
        if not ds:
            print(f"  {g:14} n={len(groups[g]):>3}  (可比对不足)")
            continue
        within_max[g] = max(ds)
        print(f"  {g:14} n={len(groups[g]):>3}  pairs={len(ds):>6}  max={max(ds)*100:5.2f}%  d95={np.quantile(ds,0.95)*100:5.2f}%")
    worst_within = max(within_max.values())
    print(f"  → 组内最远对（保持变异株完整所需上限）: max = {worst_within*100:.2f}%")

    print("--- 相邻组 cross-pair 距离（合并枝被拒所需下限 = cross max）---")
    labs = sorted(groups)
    rows = []
    for a, b in combinations(labs, 2):
        ds = cross_stats(aln, groups[a], groups[b])
        if ds:
            rows.append((min(ds), max(ds), a, b))
    rows.sort(key=lambda r: r[1])   # 按 cross max 升序 = 最难分的对在前
    for mn, mx, a, b in rows[:8]:
        print(f"  {a:10} × {b:10}  cross_min={mn*100:5.2f}%  cross_max={mx*100:5.2f}%")
    tightest_cross_max = rows[0][1]
    print(f"  → 最紧相邻对的 cross_max = {tightest_cross_max*100:.2f}% ({rows[0][2]}×{rows[0][3]})")
    print(f"  → max-within-max 可行窗口: ({worst_within*100:.2f}%, {tightest_cross_max*100:.2f}%)"
          f"  {'✔ 可行' if worst_within < tightest_cross_max else '✘ 不可行（组内最远 > 相邻组间最远）'}")

# ---- B) 合并双变异株单元的 target rate ----
print("\n========== 合并双变异株单元的 target_mapping_rate ==========")
aln4, g4 = load_g4()
for a, b in [("2004", "2006a"), ("2002", "2002CN"), ("2009", "2018WI")]:
    mem = g4[a] + g4[b]
    ref = consensus(aln4, mem)
    r_all, n_all = target_rate(aln4, mem, ref)
    r_a, _ = target_rate(aln4, g4[a], ref)
    r_b, _ = target_rate(aln4, g4[b], ref)
    print(f"  GII.4 {a}+{b} 合并参考: 全体 target={r_all:.3f}  ({a}侧={r_a:.3f}, {b}侧={r_b:.3f})  windows={n_all}")
aln17, g17 = load_g17()
for a, b in [("Kawasaki_308", "Kawasaki_323"), ("Kawasaki_323", "Romania")]:
    mem = g17[a] + g17[b]
    ref = consensus(aln17, mem)
    r_all, n_all = target_rate(aln17, mem, ref)
    r_a, _ = target_rate(aln17, g17[a], ref)
    r_b, _ = target_rate(aln17, g17[b], ref)
    print(f"  GII.17 {a}+{b} 合并参考: 全体 target={r_all:.3f}  ({a}侧={r_a:.3f}, {b}侧={r_b:.3f})  windows={n_all}")
