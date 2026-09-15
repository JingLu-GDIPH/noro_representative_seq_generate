#!/usr/bin/env python3
"""Check whether the produced reference panel resolves named GII.4 variants
and GII.17 clades: MAFFT panel-refs + labeled genomes, nearest-ref assignment."""
import csv, subprocess, sys
from collections import defaultdict
from io import StringIO
from pathlib import Path
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
PANEL = PROJ / "result_4ref/gii_mapping_references.fasta"
G4_ALN = PROJ / "v2_test/gii4_threshold_calibration/aligned_trimmed.fasta"      # 44, VAR|ACC
G17_ALN = PROJ / "v2_test/gii17_threshold_calibration_B2000/aligned_trimmed.fasta"  # 316 (B>=2000)
G17_TSV = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/Phy_nextstrain_nf/data/GII17_genome/GII17.genome.assigned.tsv")
WORK = PROJ / "v2_test/variant_resolution_check"
BASES = frozenset("ACGT")

WORK.mkdir(parents=True, exist_ok=True)

def panel_refs(vp1):
    out = {}
    for r in SeqIO.parse(PANEL, "fasta"):
        parts = r.id.split("_")
        if len(parts) >= 2 and parts[1] == vp1:
            out[r.id] = str(r.seq).upper().replace("-", "")
    return out

def mafft(records, tag):
    fa = WORK / f"{tag}.fasta"
    SeqIO.write(records, fa, "fasta")
    n = len(records)
    args = ["mafft", "--quiet", "--thread", "6"]
    if n > 200:
        args += ["--retree", "2", "--maxiterate", "0"]
    args.append(str(fa))
    r = subprocess.run(args, capture_output=True, text=True)
    return {x.id: str(x.seq).upper() for x in SeqIO.parse(StringIO(r.stdout), "fasta")}

def pdist(a, b):
    ok = [(x, y) for x, y in zip(a, b) if x in BASES and y in BASES]
    return sum(x != y for x, y in ok) / len(ok) if ok else 1.0

def analyze(tag, genomes, refs, label_of):
    records = [SeqRecord(Seq(s), id=k, description="") for k, s in genomes.items()]
    records += [SeqRecord(Seq(s), id="REF_" + k, description="") for k, s in refs.items()]
    aln = mafft(records, tag)
    ref_aln = {k: aln["REF_" + k] for k in refs}
    assign = {}
    for gid in genomes:
        d = {rk: pdist(aln[gid], rv) for rk, rv in ref_aln.items()}
        assign[gid] = min(d, key=d.get)
    # cross-tab
    tab = defaultdict(lambda: defaultdict(int))
    for gid, rk in assign.items():
        tab[label_of(gid)][rk] += 1
    print(f"\n===== {tag}: panel={len(refs)} refs, genomes={len(genomes)} =====")
    print(f"{'标签':10} {'n':>4}  最近参考分布")
    for lab in sorted(tab, key=lambda l: -sum(tab[l].values())):
        tot = sum(tab[lab].values())
        dist = ", ".join(f"{k.split('_MAPREF')[0].split('_',2)[-1]}({v})" if False else
                         f"{k.split('MAPREF_')[-1][:28]}({v})"
                         for k, v in sorted(tab[lab].items(), key=lambda x: -x[1]))
        print(f"{lab:10} {tot:>4}  {dist}")
    # separation metrics
    ref_owner = {}
    for lab, row in tab.items():
        best = max(row, key=row.get)
        ref_owner.setdefault(best, []).append(lab)
    shared = {r: l for r, l in ref_owner.items() if len(l) > 1}
    multi = sum(1 for lab, row in tab.items() if len(row) > 1)
    print(f"\n独立参考数(被占用): {len(ref_owner)}/{len(refs)}")
    print(f"多参考共担一个标签的标签数: {multi}")
    print(f"被多个标签共用的参考: " + (", ".join(f"{r[:40]}<-{'+'.join(l)}" for r, l in shared.items()) or "无"))
    return tab

# ---- GII.4 ----
g4g = {r.id: str(r.seq).upper().replace("-", "") for r in SeqIO.parse(G4_ALN, "fasta")}
g4_label = {k: k.split("|")[0] for k in g4g}
analyze("gii4", g4g, panel_refs("GII.4"), lambda gid: g4_label[gid])

# ---- GII.17 (labels incl. pre-2000 B flagged) ----
sc, dates = {}, {}
for row in csv.DictReader(open(G17_TSV), delimiter="\t"):
    s, c = row["strain"].strip(), (row["clade"] or "").strip()
    if s and c:
        sc[s] = c; dates[s] = (row["date"] or "")[:4]
def g17_label(gid):
    c = sc.get(gid)
    if c == "B" and not (dates.get(gid, "").isdigit() and int(dates[gid]) >= 2000):
        return "B(pre2000)"
    return c or "?"
g17g = {r.id: str(r.seq).upper().replace("-", "")
        for r in SeqIO.parse(G17_ALN, "fasta")}  # B>=2000 only in this file
# add pre-2000 B members from the unfiltered set? use orig alignment file of first run:
G17_FULL = PROJ / "v2_test/gii17_threshold_calibration/aligned_trimmed.fasta"
if G17_FULL.exists():
    for r in SeqIO.parse(G17_FULL, "fasta"):
        if g17_label(r.id) == "B(pre2000)":
            g17g.setdefault(r.id, str(r.seq).upper().replace("-", ""))
analyze("gii17", g17g, panel_refs("GII.17"), g17_label)
