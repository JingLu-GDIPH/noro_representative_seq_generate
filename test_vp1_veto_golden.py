#!/usr/bin/env python3
"""Gold-standard validation of the VP1 segment veto (module functions)."""
import csv
import importlib.util
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
spec = importlib.util.spec_from_file_location("bmrg", PROJ / "scripts/build_mapping_reference_group.py")
bmrg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bmrg)

G4_ALN = PROJ / "v2_test/gii4_threshold_calibration/aligned_trimmed.fasta"
G17_ALN = PROJ / "v2_test/gii17_threshold_calibration_B2000/aligned_trimmed.fasta"
REF_FA = PROJ / "reference/gii_with_genotype.fasta"
COORDS = PROJ / "reference/genome_gene_coordinates.tsv"
W, STEP, MIN_SNV, BLOCK = 150, 25, 2, 4

def consensus(aln, members):
    out = []
    for col in zip(*[aln[m] for m in members]):
        b = [x for x in col if x in bmrg.BASES]
        out.append(Counter(b).most_common(1)[0][0] if b else "N")
    return "".join(out)

def check(tag, a, b, vp1):
    veto, st = bmrg.vp1_segment_veto(a, b, vp1, W, STEP, MIN_SNV, BLOCK)
    print(f"  {tag:32} veto={str(veto):5}  VP1窗口={st['windows']:>3}  差异窗={st['divergent_windows']:>3}  最长连续={st['max_run']}")
    return veto

# ---- GII.4 ----
aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(G4_ALN, "fasta")}
groups = defaultdict(list)
for k in aln:
    groups[k.split("|")[0]].append(k)
cons_all = consensus(aln, list(aln))
vp1 = bmrg.locate_vp1_columns(cons_all, "GII.P4_GII.4", REF_FA, COORDS)
print(f"=== GII.4 校准比对上的 VP1 定位: {vp1} (比对长 {len(cons_all)}) ===")
pairs = [("2004","2006a"),("2004","2007EU"),("2002","2002CN"),("2002CN","2004"),
         ("2006a","2007EU"),("2007EU","2009"),("2002","2004"),("2004","2009")]
print("--- 相邻变异株对（应否决合并）---")
ok = 0
for a, b in pairs:
    ca, cb = consensus(aln, groups[a]), consensus(aln, groups[b])
    ok += check(f"{a} × {b}", ca, cb, vp1)
print(f"  → {ok}/{len(pairs)} 对被正确否决")

print("--- 对照（不应否决）---")
c2004 = consensus(aln, groups["2004"])
check("同一序列 vs 自身", c2004, c2004, vp1)
scattered = list(c2004)
for pos in (1000, 3000, 5500):          # 3 个散点 SNP，各占独立窗口
    if scattered[pos] in "ACGT":
        scattered[pos] = {"A":"C","C":"A","G":"T","T":"G"}[scattered[pos]]
check("自身+3个散点SNP", c2004, "".join(scattered), vp1)
m2004 = groups["2004"][0]
check("2004共识 vs 其成员之一", c2004, aln[m2004], vp1)

# ---- GII.17 ----
import csv as _csv
TSV = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/Phy_nextstrain_nf/data/GII17_genome/GII17.genome.assigned.tsv")
sc = {}
for row in _csv.DictReader(open(TSV), delimiter="\t"):
    if row["clade"]:
        sc[row["strain"]] = row["clade"]
aln17 = {r.id: str(r.seq).upper() for r in SeqIO.parse(G17_ALN, "fasta")}
g17 = defaultdict(list)
for k in aln17:
    g17[sc[k]].append(k)
cons17 = consensus(aln17, list(aln17))
vp1_17 = bmrg.locate_vp1_columns(cons17, "GII.P17_GII.17", REF_FA, COORDS)
print(f"\n=== GII.17 校准比对上的 VP1 定位: {vp1_17} (比对长 {len(cons17)}) ===")
print("--- 关键对 ---")
ck308 = consensus(aln17, g17["Kawasaki_308"])
ck323 = consensus(aln17, g17["Kawasaki_323"])
crom = consensus(aln17, g17["Romania"])
check("Kawasaki_308 × Kawasaki_323", ck308, ck323, vp1_17)
check("Kawasaki_323 × Romania", ck323, crom, vp1_17)
check("K308 共识 vs 其一成员", ck308, aln17[sorted(g17["Kawasaki_308"])[0]], vp1_17)
near1, near2 = sorted(g17["Kawasaki_308"])[:2]
check(f"K308 两个近重复成员", aln17[near1], aln17[near2], vp1_17)

# ---- 部分基因组模拟（覆盖度守卫）----
print("\n=== 部分基因组模拟（--veto-min-window-coverage 守卫）===")
c2006a = consensus(aln, groups["2006a"])
c2004b = consensus(aln, groups["2004"])
lo, hi = vp1
no_vp1 = list(c2006a)
for i in range(lo, hi + 1):                       # ORF1-only 部分序列：VP1 全 N
    no_vp1[i] = "N"
v, st = bmrg.vp1_segment_veto(c2004b, "".join(no_vp1), vp1, W, STEP, MIN_SNV, BLOCK, 0.5)
print(f"  ORF1-only部分(vs 2004)      veto={v}  evaluable={st['evaluable_windows']}/{st['windows']}  可比列占比={st['comparable_fraction']}")
no_orf1 = list(c2006a)
for i in range(0, lo):                            # 衣壳-only 部分序列：ORF1 全 N
    no_orf1[i] = "N"
v, st = bmrg.vp1_segment_veto(c2004b, "".join(no_orf1), vp1, W, STEP, MIN_SNV, BLOCK, 0.5)
print(f"  衣壳-only部分(vs 2004)      veto={v}  evaluable={st['evaluable_windows']}/{st['windows']}  差异窗={st['divergent_windows']}  max_run={st['max_run']}")
sparse = list(c2006a)
for i in range(lo, hi + 1, 3):                    # VP1 稀疏覆盖：每3列N1列(~33%缺失)
    sparse[i] = "N"
v, st = bmrg.vp1_segment_veto(c2004b, "".join(sparse), vp1, W, STEP, MIN_SNV, BLOCK, 0.5)
print(f"  VP1稀疏覆盖67%(vs 2004)     veto={v}  evaluable={st['evaluable_windows']}/{st['windows']}  差异窗={st['divergent_windows']}  max_run={st['max_run']}")
