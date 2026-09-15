#!/usr/bin/env python3
"""Forward (goal-oriented) parameter derivation. v2

Stage 1: empirical VP1 concentration factor k = d_VP1/d_genome
Stage 2: synthetic 99% sibling -> measured distinguishable rate R99 + VP1 blocks
Stage 3: sewage mixed-strain read simulation -> consensus recovery precision
         (distinguished 2-ref panel vs merged 1-ref panel), SAM-pileup based
"""
import re, subprocess, random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from Bio import SeqIO

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
G4 = PROJ / "v2_test/gii4_threshold_calibration/aligned_trimmed.fasta"
BT2 = "/Users/LuJ/Yunpan/software/Bioinf_software/bowtie2-2.3.3/bowtie2"
WORK = PROJ / "v2_test/forward_derivation"
WORK.mkdir(exist_ok=True)
BASES = frozenset("ACGT")
VP1 = (5080, 6702)
W, STEP, MIN_ALN, SMM, NMM = 150, 25, 135, 6, 2
rng = random.Random(20260910)

aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(G4, "fasta")}
groups = defaultdict(list)
for k in aln:
    groups[k.split("|")[0]].append(k)

def consensus(ms):
    out = []
    for col in zip(*[aln[m] for m in ms]):
        b = [x for x in col if x in BASES]
        out.append(Counter(b).most_common(1)[0][0] if b else "N")
    return "".join(out)

cons = {g: consensus(ms) for g, ms in groups.items()}
L = len(next(iter(cons.values())))

def region_dist(a, b, lo, hi):
    ok = [(x, y) for x, y in zip(a[lo:hi], b[lo:hi]) if x in BASES and y in BASES]
    return sum(x != y for x, y in ok) / len(ok) if len(ok) > 300 else None

# ---------- Stage 1 ----------
print("=" * 70)
print("Stage 1: 基因组距离 → VP1 距离 浓缩系数")
pairs = []
labs = sorted(cons)
PDOM = (VP1[1] - 810, VP1[1])      # P domain ≈ 3' ~810bp of VP1
for i, a in enumerate(labs):
    for b in labs[i + 1:]:
        dg = region_dist(cons[a], cons[b], 0, L)
        dv = region_dist(cons[a], cons[b], VP1[0], VP1[1])
        dp = region_dist(cons[a], cons[b], PDOM[0], PDOM[1])
        if dg and dv and dp and dg > 0.002:
            pairs.append((dg, dv, dp, a, b))
ks = sorted(dv / dg for dg, dv, _, _, _ in pairs)
kps = sorted(dp / dg for _, _, dp, _, _ in pairs)
k_med = float(np.median(ks))
kp_med = float(np.median(kps))
vp1_frac = (VP1[1] - VP1[0]) / L
pdom_frac = (PDOM[1] - PDOM[0]) / L
print(f"  n={len(pairs)} 对  k=d_VP1/d_genome: p25={np.quantile(ks,.25):.2f} median={k_med:.2f} p75={np.quantile(ks,.75):.2f}")
print(f"              k=d_P域/d_genome:  p25={np.quantile(kps,.25):.2f} median={kp_med:.2f} p75={np.quantile(kps,.75):.2f}")
for dg, dv, dp, a, b in sorted(pairs, key=lambda p: p[0])[:6]:
    print(f"    {a:8}×{b:8} 基因组 {dg*100:5.2f}%  VP1 {dv*100:5.2f}%  P域 {dp*100:5.2f}%  k_vp1={dv/dg:4.2f} k_P={dp/dg:4.2f}")
d_vp1_exp = k_med * 0.01
d_pdom_exp = kp_med * 0.01
d_rest_exp = (0.01 - vp1_frac * d_vp1_exp) / (1 - vp1_frac)
print(f"  → 1% 基因组距离 ⇒ VP1 ≈ {d_vp1_exp*100:.2f}%, P域 ≈ {d_pdom_exp*100:.2f}% (相似性 {100-d_pdom_exp*100:.1f}%), ORF1 ≈ {d_rest_exp*100:.2f}%")

# ---------- Stage 2 ----------
print("\n" + "=" * 70)
print("Stage 2: 构造 99% 相似毒株对，实测可区分率")
base = cons["2012"]
sib = list(base)
n_vp1 = round((VP1[1] - VP1[0]) * d_vp1_exp)
n_rest = round((L - (VP1[1] - VP1[0])) * d_rest_exp)
comp = {"A": "C", "C": "A", "G": "T", "T": "G"}
for p in rng.sample(range(VP1[0], VP1[1]), n_vp1):
    sib[p] = comp[base[p]]
for p in rng.sample(list(range(0, VP1[0])) + list(range(VP1[1], L)), n_rest):
    sib[p] = comp[base[p]]
sib = "".join(sib)
dg = region_dist(base, sib, 0, L)
dv = region_dist(base, sib, VP1[0], VP1[1])
print(f"  实际: 基因组 {dg*100:.2f}% (相似性 {100-dg*100:.2f}%), VP1 {dv*100:.2f}%, SNVs VP1={n_vp1} 其余={n_rest}")

def windows(seq):
    return [(st, seq[st:st + W]) for st in range(0, len(seq) - W + 1, STEP)
            if all(c in BASES for c in seq[st:st + W])]

def bowtie_rates(refs, members):
    rf = WORK / "refs.fasta"
    rf.write_text("".join(f">{k}\n{v}\n" for k, v in refs.items()))
    subprocess.run([BT2 + "-build", str(rf), str(WORK / "idx")], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    reads, meta, n = [], {}, 0
    for mk, gseq in members.items():
        for st, frag in windows(gseq):
            reads.append(f">r{n}\n{frag}\n")
            meta[f"r{n}"] = mk
            n += 1
    subprocess.run([BT2, "-x", str(WORK / "idx"), "-f", "-U", "-", "-k", "2",
                    "--very-sensitive-local", "--no-unal", "--threads", "4",
                    "-S", str(WORK / "o.sam")],
                   input="".join(reads), capture_output=True, text=True, check=True)
    hits = defaultdict(list)
    for line in open(WORK / "o.sam"):
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if f[5] == "*":
            continue
        tags = {}
        for t in f[11:]:
            p3 = t.split(":", 2)
            if len(p3) == 3:
                tags[p3[0]] = p3[2]
        hits[f[0]].append({"ref": f[2], "as": int(tags.get("AS", -999)),
                           "nm": int(tags.get("NM", 999))})
    stats = defaultdict(Counter)
    for rid, mk in meta.items():
        hh = sorted(hits.get(rid, []), key=lambda h: (-h["as"], h["nm"], h["ref"]))
        best = hh[0] if hh else None
        second = next((h for h in hh[1:] if h["ref"] != best["ref"]), None) if best else None
        sm = best["as"] - second["as"] if best and second else 999
        nm = second["nm"] - best["nm"] if best and second else 999
        own = bool(best and best["ref"] == mk)
        dist = bool(own and (sm >= SMM or nm >= NMM))
        c = stats[mk]
        c["t"] += 1
        c["d"] += int(dist)
        c["o"] += int(own)
    return {k: (c["d"] / c["t"], c["o"] / c["t"]) for k, c in stats.items()}

r2 = bowtie_rates({"A": base, "B": sib}, {"A": base, "B": sib})
print(f"  [两参考竞争]  A rate={r2['A'][0]:.3f}   B rate={r2['B'][0]:.3f}")
panel = {f"V_{g}": cons[g] for g in cons if g != "2012"}   # 排除 base 自身变异株(影子参考)
panel["A"], panel["B"] = base, sib
r22 = bowtie_rates(panel, {"A": base, "B": sib})
print(f"  [21参考竞争]  A rate={r22['A'][0]:.3f}   B rate={r22['B'][0]:.3f} (own_best A={r22['A'][1]:.3f} B={r22['B'][1]:.3f})")
R99 = min(r22["A"][0], r22["B"][0])
flags = []
for st, frag in windows(base):
    snvs = sum(1 for x, y in zip(frag, sib[st:st + W]) if x != y)
    flags.append(snvs >= 2)
run = mx = 0
for f in flags:
    run = run + 1 if f else 0
    mx = max(mx, run)
print(f"  差异块: ≥2SNV 窗 {sum(flags)}/{len(flags)}, 最长连续 {mx} 窗")
print(f"  → R99 (全panel可区分率下限) = {R99:.3f}")

# ---------- Stage 3 ----------
print("\n" + "=" * 70)
print("Stage 3: 污水混株还原模拟 (150bp, 替换错误0.3%, 总覆盖100×)")

def simulate_reads(genome, cov, tag, err=0.003):
    n = int(cov * len(genome) / W)
    out = []
    for i in range(n):
        st = rng.randrange(0, len(genome) - W + 1)
        frag = list(genome[st:st + W])
        for j in range(W):
            if rng.random() < err:
                frag[j] = rng.choice("ACGT")
        out.append((f"{tag}{i}", "".join(frag)))
    return out

def pileup_map(panel_refs, reads):
    """Map reads, assign each to its best ref, return per-ref pileups (list of Counter)."""
    rf = WORK / "refs.fasta"
    rf.write_text("".join(f">{k}\n{v}\n" for k, v in panel_refs.items()))
    subprocess.run([BT2 + "-build", str(rf), str(WORK / "idx")], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    inp = "".join(f">{n}\n{s}\n" for n, s in reads)
    subprocess.run([BT2, "-x", str(WORK / "idx"), "-f", "-U", "-", "-k", "2",
                    "--very-sensitive-local", "--no-unal", "--threads", "4",
                    "-S", str(WORK / "r.sam")],
                   input=inp, capture_output=True, text=True, check=True)
    by_id = {n: s for n, s in reads}
    best = {}
    for line in open(WORK / "r.sam"):
        if line.startswith("@"):
            continue
        f = line.split("\t")
        if f[5] == "*":
            continue
        tags = {}
        for t in f[11:]:
            p3 = t.split(":", 2)
            if len(p3) == 3:
                tags[p3[0]] = p3[2]
        qid, ref, pos, as_ = f[0], f[2], int(f[3]) - 1, int(tags.get("AS", -999))
        if qid not in best or as_ > best[qid][0]:
            best[qid] = (as_, ref, pos, f[5])
    piles = {k: [Counter() for _ in v] for k, v in panel_refs.items()}
    for qid, (as_, ref, pos, cigar) in best.items():
        seq = by_id[qid]
        qi = 0
        ri = pos
        for l, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            l = int(l)
            if op in "M=X":
                for j in range(l):
                    if ri + j < len(piles[ref]) and seq[qi + j] in BASES:
                        piles[ref][ri + j][seq[qi + j]] += 1
                qi += l; ri += l
            elif op in "IS":
                qi += l
            elif op in "DN":
                ri += l
    return piles

def call(pile, min_depth=3, majority=0.8):
    out = []
    for col in pile:
        tot = sum(col.values())
        if tot < min_depth:
            out.append("N")
        else:
            b, c = col.most_common(1)[0]
            out.append(b if c / tot >= majority else "N")
    return "".join(out)

for mix in ((0.5, 0.5), (0.9, 0.1)):
    print(f"\n  --- A:B = {mix[0]:.0%}:{mix[1]:.0%} ---")
    ra = simulate_reads(base, 100 * mix[0], "a")
    rb = simulate_reads(sib, 100 * mix[1], "b")
    snv_sites = [i for i in range(L) if base[i] in BASES and sib[i] in BASES and base[i] != sib[i]]
    piles = pileup_map({"A": base, "B": sib}, ra + rb)
    dbg = {k: sum(sum(col.values()) for col in v) for k, v in piles.items()}
    print(f"    [调试] pileup总深度: {dbg}")
    for k in ("A", "B"):
        c = call(piles[k])
        called = [i for i in range(L) if c[i] in BASES]
        if not called:
            print(f"    [分开panel] {k}: 无叫出碱基(深度不足)"); continue
        mism = sum(1 for i in called if c[i] != (base if k == "A" else sib)[i])
        ok = sum(1 for i in snv_sites if c[i] == (base if k == "A" else sib)[i])
        cross = sum(1 for i in snv_sites if c[i] == (sib if k == "A" else base)[i])
        depth_ok = sum(piles[k][i].most_common(1)[0][1] if piles[k][i] else 0 for i in range(0, L, 500)) // 15
        print(f"    [分开panel] {k}: 叫出 {len(called)}/{L} ({len(called)/L*100:.0f}%), "
              f"碱基错误率 {mism/len(called)*100:.3f}%, 判别位点正确 {ok/len(snv_sites)*100:.1f}%, "
              f"嵌合(错成对方) {cross/len(snv_sites)*100:.1f}%")
    merged_ref = "".join(base[i] if base[i] == sib[i] else
                         (Counter([base[i], sib[i]]).most_common(1)[0][0] if base[i] != sib[i] else "N")
                         for i in range(L))
    merged_ref = "".join(base[i] if base[i] == sib[i] else rng.choice([base[i], sib[i]])
                         for i in range(L))  # 共识=镶嵌(合并单元的ASR/共识等价物)
    piles2 = pileup_map({"M": merged_ref}, ra + rb)
    c = call(piles2["M"])
    called = [i for i in range(L) if c[i] in BASES]
    okA = sum(1 for i in snv_sites if c[i] == base[i])
    okB = sum(1 for i in snv_sites if c[i] == sib[i])
    mismA = sum(1 for i in called if c[i] != base[i])
    mismB = sum(1 for i in called if c[i] != sib[i])
    print(f"    [合并panel] 单一共识: 叫出 {len(called)/L*100:.0f}%, "
          f"对A错误 {mismA/len(called)*100:.2f}% / 对B错误 {mismB/len(called)*100:.2f}%, "
          f"判别位点 called-A {okA/len(snv_sites)*100:.0f}% vs called-B {okB/len(snv_sites)*100:.0f}%")
