#!/usr/bin/env python3
"""Apply production-style vs rules-doc-style decision predicates (both at T=0.5)
to known GII.4 variants / GII.17 clades. Re-uses existing Run A SAMs; adds
per-member (per-raw) rates needed for the p10 checks."""
import csv, re
from collections import defaultdict, Counter
from pathlib import Path
import numpy as np
from Bio import SeqIO

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
G4 = PROJ / "v2_test/gii4_threshold_calibration"
G17 = PROJ / "v2_test/gii17_threshold_calibration_B2000"
G17_TSV = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/Phy_nextstrain_nf/data/GII17_genome/GII17.genome.assigned.tsv")
BASES = frozenset("ACGT")
W, MIN_ALN, SMM, NMM, MAXMM = 150, 135, 6, 2, 0.10
T = 0.5

def consensus(ms, aln):
    out = []
    for col in zip(*[aln[m] for m in ms]):
        b = [x for x in col if x in BASES]
        out.append(Counter(b).most_common(1)[0][0] if b else "N")
    return "".join(out)

def parse_sam(p):
    hits = defaultdict(list)
    for line in open(p):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[5] == "*":
            continue
        aln = sum(int(l) for l, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]) if op in "MI=X")
        tags = {}
        for t in f[11:]:
            p3 = t.split(":", 2)
            if len(p3) == 3:
                tags[p3[0]] = p3[2]
        hits[f[0]].append({"ref": f[2], "as": int(tags.get("AS", -999)),
                           "nm": int(tags.get("NM", 999)), "len": aln})
    return hits

def window_meta(aln, members):
    meta, n = {}, 0
    for m in sorted(members):
        s = aln[m]
        for st in range(0, len(s) - W + 1, 25):
            if any(c not in BASES for c in s[st:st + W]):
                continue
            n += 1
            meta[f"read_{n:08d}"] = (m, st)
    return meta

def per_unit_and_member(aln, groups, sam_path):
    units = {g: {"ref": consensus(ms, aln), "members": ms} for g, ms in groups.items()}
    assign = {m: g for g, u in units.items() for m in u["members"]}
    meta = window_meta(aln, list(assign))
    hits = parse_sam(sam_path)
    su = defaultdict(Counter); sm = defaultdict(Counter)
    for rid, (m, st) in meta.items():
        target = assign[m]
        ranked = sorted(hits.get(rid, []), key=lambda h: (-h["as"], h["nm"], h["ref"]))
        best = ranked[0] if ranked else None
        second = next((h for h in ranked[1:] if h["ref"] != best["ref"]), None) if best else None
        s_m = best["as"] - second["as"] if best and second else 999
        n_m = second["nm"] - best["nm"] if best and second else 999
        own = bool(best and best["ref"] == target)
        mapped = bool(best and best["len"] >= MIN_ALN)
        ts = units[target]["ref"][st:st + W]
        comp = [(a, b) for a, b in zip(aln[m][st:st + W], ts) if a in BASES and b in BASES]
        tm = bool(len(comp) >= MIN_ALN and sum(a != b for a, b in comp) / len(comp) <= MAXMM)
        dist = bool(mapped and own and (s_m >= SMM or n_m >= NMM))
        for c, key in ((su, target), (sm, m)):
            c[key]["t"] += 1; c[key]["m"] += int(mapped); c[key]["g"] += int(tm)
            c[key]["o"] += int(own); c[key]["d"] += int(dist)
    return su, sm

def p10(vals): return float(np.quantile(vals, 0.10)) if len(vals) >= 2 else (vals[0] if vals else 1.0)

def decide(r, p10_t, p10_d):
    """r: unit rates dict. Returns (prod, doc) decisions."""
    prod = ("split" if r["target"] < 0.95 else
            "merge" if r["dist"] < T else "retain")
    retain_ok = (r["dist"] >= T and r["own"] >= 0.90 and r["target"] >= 0.95
                 and p10_t >= 0.90 and p10_d >= 0.70)
    split_trig = (r["target"] < 0.95 or p10_t < 0.90 or
                  (r["target"] >= 0.95 and r["dist"] >= T and p10_d < 0.70))
    if retain_ok:
        doc = "retain"
    elif split_trig:
        doc = "split_review"
    else:
        doc = "merge_review"
    return prod, doc

def report(tag, aln_fa, groups, sam):
    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln_fa, "fasta")}
    groups = {g: [m for m in ms if m in aln] for g, ms in groups.items()}
    groups = {g: ms for g, ms in groups.items() if ms}
    su, sm = per_unit_and_member(aln, groups, sam)
    print(f"\n===== {tag} (T={T}) =====")
    print(f"{'unit':14}{'n':>3} {'target':>6} {'own':>6} {'dist':>6} {'p10T':>6} {'p10D':>6} | {'生产版':>10} {'规则文档版':>12}")
    counts = Counter()
    for g in sorted(groups, key=lambda g: su[g]["d"]/su[g]["t"]):
        c = su[g]
        r = {"target": c["g"]/c["t"], "own": c["o"]/c["t"], "dist": c["d"]/c["t"]}
        mem_d = [sm[m]["d"]/sm[m]["t"] for m in groups[g] if sm[m]["t"] > 0]
        mem_t = [sm[m]["g"]/sm[m]["t"] for m in groups[g] if sm[m]["t"] > 0]
        pt, pd = p10(mem_t), p10(mem_d)
        prod, doc = decide(r, pt, pd)
        counts[prod] += 1; counts[doc] += 1
        print(f"{g:14}{len(groups[g]):>3} {r['target']:>6.3f} {r['own']:>6.3f} {r['dist']:>6.3f} "
              f"{pt:>6.3f} {pd:>6.3f} | {prod:>10} {doc:>12}")
    return counts

# GII.4
g4g = defaultdict(list)
for r in SeqIO.parse(G4 / "aligned_trimmed.fasta", "fasta"):
    g4g[r.id.split("|")[0]].append(r.id)
c4 = report("GII.4 变异株", G4 / "aligned_trimmed.fasta", g4g,
            G4 / "runA_variant_panel/win.sam")
# GII.17 (B>=2000)
sc, dates = {}, {}
for row in csv.DictReader(open(G17_TSV), delimiter="\t"):
    s, c = row["strain"].strip(), (row["clade"] or "").strip()
    if s and c:
        sc[s] = c; dates[s] = (row["date"] or "")[:4]
drop = {s for s, c in sc.items() if c == "B" and not (dates[s].isdigit() and int(dates[s]) >= 2000)}
g17g = defaultdict(list)
for r in SeqIO.parse(G17 / "aligned_trimmed.fasta", "fasta"):
    c = sc.get(r.id)
    if c and r.id not in drop:
        g17g[c].append(r.id)
c17 = report("GII.17 clade(B≥2000)", G17 / "aligned_trimmed.fasta", g17g,
             G17 / "runA_clade_panel/win.sam")

print("\n===== 汇总（20 GII.4 变异株 + 5 GII.17 clade 的命运）=====")
for name, c in (("GII.4", c4), ("GII.17", c17)):
    print(f"{name}: 生产版 -> retain={c['retain']}, merge={c['merge']}, split={c['split']}; "
          f"规则文档版 -> retain={c['retain']}, merge_review={c['merge_review']}, split_review={c['split_review']}")
