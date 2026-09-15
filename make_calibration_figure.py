#!/usr/bin/env python3
"""Publication figure: distinguishable_rate threshold calibration (GII.4 + GII.17). v2

A: rate distributions (keep vs merge side) + threshold lines
B: rate vs min inter-unit distance (log x, all points) + Poisson theory
C: threshold decision band
"""
import csv
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import SeqIO

PROJ = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate")
G4 = PROJ / "v2_test/gii4_threshold_calibration"
G17 = PROJ / "v2_test/gii17_threshold_calibration_B2000"
G17_TSV = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/Phy_nextstrain_nf/data/GII17_genome/GII17.genome.assigned.tsv")
OUT = PROJ / "result_paper"
BASES = frozenset("ACGT")

C_G4, C_G17 = "#0072B2", "#D55E00"
C_G4_W, C_G17_W = "#56B4E9", "#E69F00"
GREEN, YELLOW, RED = "#009E73", "#E8C547", "#E07B54"
THRESH_STYLE = {0.30: ":", 0.50: "--", 0.60: "-.", 0.80: (0, (4, 1, 1, 1, 1, 1))}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica"],
    "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

def read_rates(p):
    with open(p) as fh:
        return list(csv.DictReader(fh))

def consensus_and_min_dist(aln_fasta, groups):
    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln_fasta, "fasta")}
    cons = {}
    for lab, ms in groups.items():
        cols = list(zip(*[aln[m] for m in ms if m in aln]))
        out = []
        for col in cols:
            b = [x for x in col if x in BASES]
            out.append(Counter(b).most_common(1)[0][0] if b else "N")
        cons[lab] = "".join(out)
    labels = sorted(cons)
    mind = {}
    for i, a in enumerate(labels):
        ds = []
        for j, b in enumerate(labels):
            if i == j:
                continue
            ok = [(x, y) for x, y in zip(cons[a], cons[b]) if x in BASES and y in BASES]
            ds.append(sum(x != y for x, y in ok) / len(ok))
        mind[a] = min(ds)
    return {lab: mind[lab] for lab in labels}

# ---- data ----
g4A = read_rates(G4 / "runA_variant_rates.tsv")
g4B = read_rates(G4 / "runB_accession_rates.tsv")
g4_groups = defaultdict(list)
for r in SeqIO.parse(G4 / "aligned_trimmed.fasta", "fasta"):
    g4_groups[r.id.split("|")[0]].append(r.id)
g4_dist = consensus_and_min_dist(G4 / "aligned_trimmed.fasta", g4_groups)

g17A = read_rates(G17 / "runA_clade_rates.tsv")
g17B = read_rates(G17 / "runB_strain_rates.tsv")
strain_clade, dates = {}, {}
with open(G17_TSV) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        s, c = row["strain"].strip(), (row["clade"] or "").strip()
        if s and c:
            strain_clade[s] = c
            dates[s] = (row["date"] or "")[:4]
drop = {s for s, c in strain_clade.items()
        if c == "B" and not (dates[s].isdigit() and int(dates[s]) >= 2000)}
g17_groups = defaultdict(list)
for r in SeqIO.parse(G17 / "aligned_trimmed.fasta", "fasta"):
    c = strain_clade.get(r.id)
    if c and r.id not in drop:
        g17_groups[c].append(r.id)
g17_dist = consensus_and_min_dist(G17 / "aligned_trimmed.fasta", g17_groups)

FLOOR_G4 = min(float(r["distinguishable_rate"]) for r in g4A)
FLOOR_G17 = min(float(r["distinguishable_rate"]) for r in g17A if int(r["windows"]) > 0)

fig = plt.figure(figsize=(7.2, 2.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.22, 1.18, 1.22], wspace=0.38,
                      left=0.062, right=0.988, top=0.86, bottom=0.17)
axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1]); axC = fig.add_subplot(gs[2])

# ============ Panel A ============
cats = [
    ("GII.4\nbetween (20)", [float(r["distinguishable_rate"]) for r in g4A], C_G4),
    ("GII.4\nwithin (43)", [float(r["distinguishable_rate"]) for r in g4B], C_G4_W),
    ("GII.17\nbetween (5)", [float(r["distinguishable_rate"]) for r in g17A], C_G17),
    ("GII.17\nwithin (316)", [float(r["distinguishable_rate"]) for r in g17B], C_G17_W),
]
rng = np.random.default_rng(42)
for i, (name, vals, col) in enumerate(cats):
    v = np.array(vals)
    axA.boxplot(v, positions=[i], widths=0.55, showfliers=False, patch_artist=True,
                whis=(5, 95), medianprops=dict(color="black", lw=0.9),
                boxprops=dict(facecolor=col, alpha=0.45, edgecolor=col, lw=0.8),
                whiskerprops=dict(color=col, lw=0.8), capprops=dict(color=col, lw=0.8))
    axA.scatter(rng.uniform(i - 0.22, i + 0.22, len(v)), v, s=4, color=col,
                alpha=0.65, edgecolors="none", zorder=3)
for T, lab in ((0.30, "0.30 production"), (0.50, "0.50 recommended"),
               (0.80, "0.80 rules doc")):
    axA.axhline(T, color="0.35", lw=0.7, ls=THRESH_STYLE[T], zorder=1)
    axA.annotate(lab, xy=(1.0, T), xycoords=("axes fraction", "data"),
                 xytext=(2, 1.5), textcoords="offset points",
                 fontsize=5.4, color="0.25", ha="left", va="bottom")
axA.set_xticks(range(4))
axA.set_xticklabels([c[0] for c in cats], fontsize=6, linespacing=1.25)
axA.set_ylabel("Distinguishable rate")
axA.set_ylim(-0.03, 1.06)
axA.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axA.set_xlim(-0.55, 3.55)
axA.spines[["top", "right"]].set_visible(False)

# ============ Panel B (log x) ============
for r in g4A:
    axB.scatter(g4_dist[r["unit"]] * 100, float(r["distinguishable_rate"]),
                s=16, color=C_G4, edgecolor="white", linewidth=0.4, zorder=3)
for r in g17A:
    if int(r["windows"]) == 0:
        continue
    axB.scatter(g17_dist[r["unit"]] * 100, float(r["distinguishable_rate"]),
                s=16, color=C_G17, marker="s", edgecolor="white", linewidth=0.4, zorder=3)
dd = np.geomspace(0.018, 0.19, 250)
axB.plot(dd * 100, 1 - np.exp(-150 * dd) * (1 + 150 * dd),
         color="0.4", lw=0.9, ls="--", zorder=2)
ann = [("2006a", g4A, g4_dist, C_G4, (-14, -10), "left"),
       ("2018AG", g4A, g4_dist, C_G4, (4, -9), "left"),
       ("Kawasaki_323", g17A, g17_dist, C_G17, (0, -10), "center"),
       ("B", g17A, g17_dist, C_G17, (0, -13), "center")]
for unit, rows, dist, col, (dx, dy), ha in ann:
    row = next(r for r in rows if r["unit"] == unit)
    x, y = dist[unit] * 100, float(row["distinguishable_rate"])
    lab = unit + " (≥2000)" if unit == "B" else unit
    axB.annotate(lab, (x, y), xytext=(dx, dy), textcoords="offset points",
                 fontsize=5.6, color=col, ha=ha)
axB.set_xscale("log")
axB.set_xlim(1.7, 19)
axB.set_xticks([2, 3, 4, 6, 9, 16])
axB.set_xticklabels(["2", "3", "4", "6", "9", "16"])
axB.minorticks_off()
axB.set_ylim(0.52, 1.04)
axB.set_xlabel("Distance to nearest panel unit (%)")
axB.set_ylabel("Distinguishable rate")
axB.scatter([], [], s=16, color=C_G4, label="GII.4 variants (n=20)")
axB.scatter([], [], s=16, color=C_G17, marker="s", label="GII.17 clades (n=5)")
axB.plot([], [], color="0.4", lw=0.9, ls="--",
         label="Poisson: P(≥2 SNVs/150 bp)")
axB.legend(frameon=False, loc="lower left", handlelength=1.3, borderpad=0.2)
axB.spines[["top", "right"]].set_visible(False)

# ============ Panel C ============
axC.axvspan(0, FLOOR_G17, color=GREEN, alpha=0.16, zorder=0)
axC.axvspan(FLOOR_G17, FLOOR_G4, color=YELLOW, alpha=0.35, zorder=0)
axC.axvspan(FLOOR_G4, 1.0, color=RED, alpha=0.14, zorder=0)
for T, lab, row in ((0.30, "0.30\nproduction", 1.78), (0.50, "0.50\nrecommended", 1.52),
                    (0.60, "0.60\nmarginal", 1.78), (0.80, "0.80\nrules doc", 1.52)):
    axC.axvline(T, ymin=0, ymax=1.44 / 1.9, color="0.2", lw=0.8,
                ls=THRESH_STYLE[T], zorder=2)
    axC.text(T, row, lab, ha="center", va="bottom", fontsize=5.3, color="0.15",
             linespacing=1.2)
# verdict texts
axC.text(FLOOR_G17 / 2, 1.30, "both PASS", ha="center", fontsize=6, color="0.2")
axC.text((FLOOR_G17 + FLOOR_G4) / 2, 1.30, "GII.4 only", ha="center", fontsize=6, color="0.35")
axC.text((FLOOR_G4 + 1.0) / 2, 1.30, "both FAIL", ha="center", fontsize=6, color="0.4")
# between-label rate ticks + floors
axC.plot([float(r["distinguishable_rate"]) for r in g4A], np.full(len(g4A), 1.06),
         "|", ms=7, color=C_G4, mew=1.1, zorder=4)
axC.plot([FLOOR_G4], [1.06], "v", ms=5, color=C_G4, zorder=5)
axC.text(FLOOR_G4 + 0.017, 1.06, "GII.4 floor 0.739", fontsize=5.4, va="center", color=C_G4)
axC.plot([float(r["distinguishable_rate"]) for r in g17A if int(r["windows"]) > 0],
         np.full(sum(1 for r in g17A if int(r["windows"]) > 0), 0.78),
         "|", ms=7, color=C_G17, mew=1.1, zorder=4)
axC.plot([FLOOR_G17], [0.78], "v", ms=5, color=C_G17, zorder=5)
axC.text(FLOOR_G17 - 0.017, 0.78, "GII.17 floor 0.636", fontsize=5.4, va="center",
         ha="right", color=C_G17)
# within-label rugs
for vals, y0, col in ((g4B, 0.36, C_G4_W), (g17B, 0.05, C_G17_W)):
    v = np.array([float(r["distinguishable_rate"]) for r in vals])
    hist, edges = np.histogram(v, bins=np.linspace(0, 1, 41))
    h = hist / hist.max() * 0.20
    for hh, e0, e1 in zip(h, edges[:-1], edges[1:]):
        if hh > 0:
            axC.add_patch(plt.Rectangle((e0, y0), e1 - e0, hh,
                                        facecolor=col, edgecolor="none", alpha=0.85))
axC.text(0.99, 0.60, "GII.4 within-variant (n=43)", ha="right", fontsize=5.4, color="0.3")
axC.text(0.99, 0.285, "GII.17 within-clade (n=316)", ha="right", fontsize=5.4, color="0.3")
axC.set_xlim(0, 1.0)
axC.set_ylim(0, 1.92)
axC.set_yticks([])
axC.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
axC.set_xlabel("Distinguishable-rate threshold T")
axC.spines[["top", "right", "left"]].set_visible(False)

for ax, lab in ((axA, "A"), (axB, "B"), (axC, "C")):
    ax.text(-0.13, 1.09, lab, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top")

OUT.mkdir(exist_ok=True)
(OUT / "figure_data").mkdir(exist_ok=True)
fig.savefig(OUT / "Fig_threshold_calibration.pdf")
fig.savefig(OUT / "Fig_threshold_calibration.png", dpi=300)

with open(OUT / "figure_data" / "panel_scatter_data.tsv", "w") as fh:
    fh.write("genotype\tunit\tmin_distance_pct\tdistinguishable_rate\n")
    for r in g4A:
        fh.write(f"GII.4\t{r['unit']}\t{g4_dist[r['unit']]*100:.4f}\t{r['distinguishable_rate']}\n")
    for r in g17A:
        if int(r["windows"]) > 0:
            fh.write(f"GII.17\t{r['unit']}\t{g17_dist[r['unit']]*100:.4f}\t{r['distinguishable_rate']}\n")
with open(OUT / "figure_data" / "panel_box_data.tsv", "w") as fh:
    fh.write("category\tunit\tdistinguishable_rate\n")
    for r in g4A: fh.write(f"GII.4_between\t{r['unit']}\t{r['distinguishable_rate']}\n")
    for r in g4B: fh.write(f"GII.4_within\t{r['unit']}\t{r['distinguishable_rate']}\n")
    for r in g17A: fh.write(f"GII.17_between\t{r['unit']}\t{r['distinguishable_rate']}\n")
    for r in g17B: fh.write(f"GII.17_within\t{r['unit']}\t{r['distinguishable_rate']}\n")

print(f"floors: GII.4={FLOOR_G4:.3f} GII.17={FLOOR_G17:.3f}")
print("saved:", OUT / "Fig_threshold_calibration.png")
