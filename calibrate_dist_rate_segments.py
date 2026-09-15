#!/usr/bin/env python3
"""Segment-aware calibration of the distinguishable-rate threshold for GII.4 variants.

Concern: between-variant differences may be CONCENTRATED in genome segments while
whole-genome distance can even invert vs within-variant distance. Test at WINDOW level:
A) consensus-vs-consensus per-window SNV profiles for adjacent variant pairs
   (concentration, longest non-discriminating run, region location)
B) pairwise 2-ref competitive rates: adjacent between-variant pairs (keep side)
   vs within-variant pairs (merge side)
C) overlap verdict + segment-aware recommendation
"""
import re, subprocess
from collections import defaultdict, Counter
from io import StringIO
from pathlib import Path
import numpy as np
from Bio import SeqIO

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
ALN_F = PROJ / "v2_test/gii4_threshold_calibration/aligned_trimmed.fasta"
WORK = PROJ / "v2_test/segment_aware_dist_check"
BT2 = "/Users/LuJ/Yunpan/software/Bioinf_software/bowtie2-2.3.3/bowtie2"
BASES = frozenset("ACGT")
W, STEP = 150, 25
WORK.mkdir(exist_ok=True)

aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(ALN_F, "fasta")}
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

def cdist(a, b):
    ok = [(x, y) for x, y in zip(cons[a], cons[b]) if x in BASES and y in BASES]
    return sum(x != y for x, y in ok) / len(ok)

labs = sorted(groups)
pairs_by_d = sorted(((cdist(a, b), a, b) for i, a in enumerate(labs) for b in labs[i+1:]))
adjacent = [(a, b) for _, a, b in pairs_by_d[:8]]

def region(frac):
    if frac < 0.65: return "ORF1区"
    if frac < 0.90: return "VP1区"
    return "VP2/3'区"

print("===== A) 相邻变异株对的窗口级差异分布（共识 vs 共识）=====")
print(f"{'对':18} {'全基因组d':>8} {'w=0':>5} {'w=1':>5} {'w2-5':>6} {'w>5':>5} {'最长连续<2SNV窗口':>12} {'差异块所在区':>10}")
seg_profiles = {}
for a, b in adjacent:
    snvs = []
    for st in range(0, len(cons[a]) - W + 1, STEP):
        wa, wb = cons[a][st:st+W], cons[b][st:st+W]
        ok = [(x, y) for x, y in zip(wa, wb) if x in BASES and y in BASES]
        snvs.append(sum(x != y for x, y in ok))
    snvs = np.array(snvs)
    # longest run of windows with <2 SNVs (non-discriminating)
    run = best = 0
    for s in snvs:
        run = run + 1 if s < 2 else 0
        best = max(best, run)
    blocks = [i for i, s in enumerate(snvs) if s >= 2]
    if blocks:
        span = (min(blocks), max(blocks))
        regs = Counter(region((i * STEP + W/2) / len(cons[a])) for i in blocks)
        regstr = "+".join(f"{k}" for k, _ in regs.most_common())
        spanstr = f"{region((min(blocks)*STEP)/len(cons[a]))}"
    else:
        spanstr = "无"
    d = cdist(a, b)
    seg_profiles[(a, b)] = snvs
    print(f"{a+'×'+b:18} {d*100:>7.2f}% {np.mean(snvs==0)*100:>4.0f}% {np.mean(snvs==1)*100:>4.0f}% "
          f"{np.mean((snvs>=2)&(snvs<=5))*100:>5.0f}% {np.mean(snvs>5)*100:>4.0f}% "
          f"{best*STEP:>9}nt  {spanstr:>10}")

# ---- B) pairwise 2-ref competitive rates ----
def windows_of(members):
    out = []
    for m in sorted(members):
        s = aln[m]
        for st in range(0, len(s) - W + 1, STEP):
            if all(c in BASES for c in s[st:st+W]):
                out.append((m, st, s[st:st+W]))
    return out

def pairwise_rates(ref_units):
    """ref_units: {label: (members, consensus)} -> {label: rate}"""
    refs = {lab: c.replace("-", "") for lab, (_, c) in ref_units.items()}
    assign = {m: lab for lab, (ms, _) in ref_units.items() for m in ms}
    reads = windows_of([m for ms, _ in ref_units.values() for m in ms])
    rf = WORK / "refs.fasta"
    rf.write_text("".join(f">{k}\n{v}\n" for k, v in refs.items()))
    subprocess.run([BT2+"-build", str(rf), str(WORK/"idx")], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    p = subprocess.run([BT2, "-x", str(WORK/"idx"), "-f", "-U", "-", "-k", "2",
                        "--very-sensitive-local", "--no-unal", "--threads", "4", "-S", str(WORK/"o.sam")],
                       input="".join(f">r{i}\n{s}\n" for i, (_, _, s) in enumerate(reads)),
                       capture_output=True, text=True, check=True)
    hits = defaultdict(list)
    for line in p.stdout.splitlines() if False else open(WORK/"o.sam"):
        if line.startswith("@"): continue
        f = line.split("\t")
        if f[5] == "*": continue
        tags = {}
        for t in f[11:]:
            p3 = t.split(":", 2)
            if len(p3) == 3: tags[p3[0]] = p3[2]
        hits[f[0]].append({"ref": f[2], "as": int(tags.get("AS", -999)), "nm": int(tags.get("NM", 999))})
    stats = defaultdict(Counter)
    for i, (m, st, _) in enumerate(reads):
        hh = hits.get(f"r{i}", [])
        hh.sort(key=lambda h: (-h["as"], h["nm"], h["ref"]))
        best = hh[0] if hh else None
        second = next((h for h in hh[1:] if h["ref"] != best["ref"]), None) if best else None
        sm = best["as"] - second["as"] if best and second else 999
        nm = second["nm"] - best["nm"] if best and second else 999
        own = bool(best and best["ref"] == assign[m])
        ok = bool(best and (sm >= 6 or nm >= 2))
        c = stats[assign[m]]
        c["t"] += 1; c["d"] += int(own and ok)
    return {lab: c["d"]/c["t"] for lab, c in stats.items() if c["t"]}

print("\n===== B1) 变异株内成对（merge 侧，2-ref 竞争）=====")
within_rates = []
for g in sorted(groups):
    ms = groups[g]
    if len(ms) < 2:
        continue
    units = {f"{m.split('|')[1][:10]}": ([m], aln[m].replace("-", "")) for m in ms[:2]}  # 两两首两名成员
    r = pairwise_rates(units)
    vals = list(r.values())
    within_rates += vals
    print(f"  {g:8}  两成员互相竞争 rate = {', '.join(f'{v:.3f}' for v in vals)}")

print("\n===== B2) 相邻变异株对（keep 侧，2-ref 竞争）=====")
between_rates = {}
for a, b in adjacent:
    units = {a: (groups[a], cons[a]), b: (groups[b], cons[b])}
    r = pairwise_rates(units)
    between_rates[(a, b)] = r
    print(f"  {a:8}×{b:8}  rate: {a}={r.get(a, float('nan')):.3f}  {b}={r.get(b, float('nan')):.3f}")

btw = sorted(v for r in between_rates.values() for v in r.values())
wit = sorted(within_rates)
print("\n===== C) 分布对比与结论 =====")
print(f"keep 侧（株间）: min={btw[0]:.3f}  p25={np.quantile(btw,0.25):.3f}  median={np.median(btw):.3f}")
print(f"merge 侧（株内）: min={wit[0]:.3f}  median={np.median(wit):.3f}  max={wit[-1]:.3f}")
overlap = [x for x in wit if x > btw[0]]
print(f"窗口层倒挂: 株内 rate 高于株间最小值({btw[0]:.3f})的株内对有 {len(overlap)}/{len(wit)} 个: {[f'{x:.2f}' for x in overlap]}")
worst = min(between_rates.items(), key=lambda kv: min(kv[1].values()))
print(f"最难株间对: {worst[0][0]}×{worst[0][1]} (两侧最低 {min(worst[1].values()):.3f})")
