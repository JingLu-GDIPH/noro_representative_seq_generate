#!/usr/bin/env python3
"""Plot per-genotype trees with consensus tips highlighted."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-noro")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from Bio import Phylo


def parse_map(path: Path) -> tuple[set[str], set[str]]:
    raw = set()
    consensus = set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["source"] == "raw":
                raw.add(row["tree_id"])
            elif row["source"] == "consensus":
                consensus.add(row["tree_id"])
    return raw, consensus


def assign_y_positions(tree):
    terminals = tree.get_terminals()
    y_positions = {terminal: index for index, terminal in enumerate(terminals)}

    def calc_y(clade):
        if clade in y_positions:
            return y_positions[clade]
        child_y = [calc_y(child) for child in clade.clades]
        y_positions[clade] = sum(child_y) / len(child_y)
        return y_positions[clade]

    calc_y(tree.root)
    return y_positions


def plot_tree(tree_path: Path, id_map_path: Path, genotype: str, output: Path) -> None:
    raw_ids, consensus_ids = parse_map(id_map_path)
    tree = Phylo.read(tree_path, "newick")
    depths = tree.depths()
    if not depths or max(depths.values()) == 0:
        depths = tree.depths(unit_branch_lengths=True)
    y_positions = assign_y_positions(tree)
    n_tips = len(tree.get_terminals())
    height = min(18, max(5, n_tips / 80))
    fig, ax = plt.subplots(figsize=(11, height))

    def draw_clade(clade):
        x_here = depths[clade]
        y_here = y_positions[clade]
        for child in clade.clades:
            x_child = depths[child]
            y_child = y_positions[child]
            ax.plot([x_here, x_child], [y_child, y_child], color="#b0b0b0", linewidth=0.35)
            ax.plot([x_here, x_here], [y_here, y_child], color="#d0d0d0", linewidth=0.35)
            draw_clade(child)

    draw_clade(tree.root)

    raw_x = []
    raw_y = []
    cons_x = []
    cons_y = []
    cons_labels = []
    for terminal in tree.get_terminals():
        if terminal.name in consensus_ids:
            cons_x.append(depths[terminal])
            cons_y.append(y_positions[terminal])
            cons_labels.append(terminal.name)
        elif terminal.name in raw_ids:
            raw_x.append(depths[terminal])
            raw_y.append(y_positions[terminal])

    ax.scatter(raw_x, raw_y, s=5, color="#778899", alpha=0.55, label="Downloaded raw sequence", zorder=3)
    ax.scatter(cons_x, cons_y, s=26, color="#d62828", edgecolors="black", linewidths=0.2, label="Consensus sequence", zorder=4)

    if len(cons_labels) <= 20:
        x_max = max(depths.values()) if depths else 1
        for label, x, y in zip(cons_labels, cons_x, cons_y):
            short = label.replace("CONSENSUS__VP1-", "").replace("__RdRp-", "|RdRp-")
            ax.text(x + x_max * 0.01, y, short[:70], fontsize=6, va="center", color="#7f0000")

    ax.set_title(f"{genotype}: raw + consensus FastTree (consensus tips highlighted)")
    ax.set_xlabel("Branch length")
    ax.set_ylabel("Tips ordered by tree traversal")
    ax.set_yticks([])
    ax.legend(loc="upper right")
    ax.invert_yaxis()
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def select_groups(summary_path: Path, requested: list[str], top_n: int) -> list[str]:
    rows = list(csv.DictReader(summary_path.open(), delimiter="\t"))
    selected = []
    for group in requested:
        if any(row["vp1_group"] == group and row["status"] == "analyzed" for row in rows):
            selected.append(group)

    analyzed = [row for row in rows if row["status"] == "analyzed" and row.get("median_nearest_distance")]
    attention = sorted(analyzed, key=lambda row: float(row["median_nearest_distance"]), reverse=True)[:top_n]
    large = [row for row in analyzed if int(row["raw_n"]) >= 200]
    for row in attention + large:
        if row["vp1_group"] not in selected:
            selected.append(row["vp1_group"])
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=Path("results_new/representativeness_analysis"))
    parser.add_argument("--groups", nargs="*", default=["GII.4", "GII.17", "GII.2", "GII.3", "GII.6", "GI.1", "GI.3"])
    parser.add_argument("--top-attention", type=int, default=5)
    args = parser.parse_args()

    summary_path = args.analysis_dir / "03_metrics" / "representativeness_summary.tsv"
    output_dir = args.analysis_dir / "04_figures" / "highlighted_trees"
    groups = select_groups(summary_path, args.groups, args.top_attention)
    for group in groups:
        group_dir = args.analysis_dir / "02_by_vp1" / group
        tree = group_dir / f"{group}_raw_plus_consensus.fasttree.nwk"
        id_map = group_dir / f"{group}_id_map.tsv"
        if tree.exists() and id_map.exists():
            print(f"Plotting {group}")
            plot_tree(tree, id_map, group, output_dir / f"{group}_combined_tree_consensus_highlighted.png")
    print(f"Tree figures: {output_dir}")


if __name__ == "__main__":
    main()
