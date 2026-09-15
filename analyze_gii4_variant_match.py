#!/usr/bin/env python3
"""Compare the converged GII.4 v3 panel against the 20 named variants.

Nearest-reference assignment in VP1-region space (variants are capsid-defined):
MAFFT labeled genomes + panel refs together, transfer VP1 columns via the
annotation machinery, then p-distance restricted to VP1.
"""
import glob
import importlib.util
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
spec = importlib.util.spec_from_file_location("bmrg", PROJ / "scripts/build_mapping_reference_group.py")
bmrg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bmrg)

G4_ALN = PROJ / "v2_test/gii4_threshold_calibration/aligned_trimmed.fasta"
REF_FA = PROJ / "reference/gii_with_genotype.fasta"
COORDS = PROJ / "reference/genome_gene_coordinates.tsv"
WORK = PROJ / "v2_test/gii4_v3_variant_match"
WORK.mkdir(exist_ok=True)
BASES = frozenset("ACGT")

# collect panel refs from the 5 converged GII.4 groups
panel = {}
for fa in sorted(glob.glob(str(PROJ / "result_4ref/gii/06_group_references/mapping_ref_*GII.4*_v4conv/final_mapping_references.fasta"))):
    for r in SeqIO.parse(fa, "fasta"):
        panel[r.id] = str(r.seq).upper()
print(f"panel GII.4 refs: {len(panel)}")
genomes = {r.id: str(r.seq).upper().replace("-", "") for r in SeqIO.parse(G4_ALN, "fasta")}
print(f"labeled genomes: {len(genomes)} (20 variants)")

recs = [SeqRecord(Seq(s), id=g, description="") for g, s in genomes.items()]
recs += [SeqRecord(Seq(s), id="REF_" + k, description="") for k, s in panel.items()]
combined = WORK / "combined.fasta"
SeqIO.write(recs, combined, "fasta")
import subprocess
from io import StringIO
proc = subprocess.run(["mafft", "--quiet", "--retree", "2", "--maxiterate", "0",
                       "--thread", "6", str(combined)], capture_output=True, text=True, check=True)
aln = {x.id: str(x.seq).upper() for x in SeqIO.parse(StringIO(proc.stdout), "fasta")}
consensus = "".join(
    max((b for b in col if b in "ACGT"), key=col.count) if any(c in "ACGT" for c in col) else "N"
    for col in zip(*[aln["REF_" + k] for k in panel])
) if panel else ""
vp1 = bmrg.locate_vp1_columns(consensus, "GII.P4_GII.4", REF_FA, COORDS)
print(f"VP1 columns on combined alignment: {vp1} (len {len(next(iter(aln.values())))})")
lo, hi = vp1

def vp1_dist(a, b):
    ok = [(x, y) for x, y in zip(a[lo:hi], b[lo:hi]) if x in BASES and y in BASES]
    return sum(x != y for x, y in ok) / len(ok) if len(ok) > 300 else None

ref_aln = {k: aln["REF_" + k] for k in panel}
tab = defaultdict(lambda: defaultdict(float))
owner = defaultdict(list)
for g in genomes:
    best, bd = None, 1e9
    for rk, rv in ref_aln.items():
        d = vp1_dist(aln[g], rv)
        if d is not None and d < bd:
            best, bd = rk, d
    label = g.split("|")[0]
    tab[label][best] += 1
    owner[best].append((label, bd))

print(f"\n{'变异株':10} {'n':>3} {'最近参考':54} {'VP1距离'}")
for lab in sorted(tab, key=lambda l: min(owner[next(iter(tab[l]))][0][1], 0)):
    row = tab[lab]
    tot = sum(row.values())
    top = max(row, key=row.get)
    dmin = min(d for l, d in owner[top] if l == lab)
    print(f"{lab:10} {tot:>3} {top.replace('GII.P','P').replace('_MAPREF','')[:54]:54} {dmin*100:5.2f}%")

ref_owner = defaultdict(set)
for ref, pairs in owner.items():
    for lab, _ in pairs:
        ref_owner[ref].add(lab)
multi = {r: l for r, l in ref_owner.items() if len(l) > 1}
solo = {r: l for r, l in ref_owner.items() if len(l) == 1}
print(f"\n独立服务单一变异株的参考: {len(solo)}/{len(panel)}")
print(f"被多变异株共用的参考: {len(multi)}")
for r, l in sorted(multi.items(), key=lambda x: -len(x[1])):
    print(f"  {r.replace('GII.P','P').replace('_MAPREF','')[:60]:60} <- {'+'.join(sorted(l))}")
variants_served = set().union(*ref_owner.values()) if ref_owner else set()
print(f"被 20 个变异株占用的参考总数: {len(ref_owner)}；未被任何标注变异株最近命中的参考: {len(panel)-len(ref_owner)}")
