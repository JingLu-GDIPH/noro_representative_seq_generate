#!/usr/bin/env python3
"""Assess norovirus consensus representativeness by VP1 genotype.

The workflow:
1. Rename final consensus records with VP1 and RdRp genotypes.
2. Split downloaded records and consensus records by VP1 genotype.
3. Build MAFFT alignments and, when available, FastTree trees per genotype.
4. Quantify representativeness using nearest consensus distances and topology metrics.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-noro")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo, SeqIO
from Bio.SeqRecord import SeqRecord


VALID_BASES = set("ACGT")
DEFAULT_RAW_QC_MAX_N_GAP_PERCENT = 10.0


def safe_id(value: str) -> str:
    value = value.strip() or "Unknown"
    value = value.replace("/", "-")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "-", value)
    value = re.sub(r"-+", "-", value).strip("-")
    return value or "Unknown"


def genotype_sort_key(genotype: str):
    match = re.match(r"^(G(?:I|II|IX|V|X))\.([A-Za-z]*)(\d+|NA\d*)$", genotype)
    if not match:
        return (99, 9999, genotype)
    group, prefix, number = match.groups()
    order = {"GI": 1, "GII": 2, "GIX": 9, "GV": 5, "GX": 10}.get(group, 99)
    try:
        n = int(number)
    except ValueError:
        n = 9999
    return (order, n, prefix)


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seqname = row.get("seqname", "").strip()
            accession = row.get("accession", "").strip()
            genotype = row.get("vp1_genotype", "").strip()
            if not genotype:
                continue
            if seqname:
                rows[seqname] = row
            if accession:
                rows[accession] = row
    return rows


def n_gap_percent(sequence: str) -> tuple[float, int, int]:
    """Return combined N/gap percent, ambiguous count, and ungapped/gapped length."""
    sequence = str(sequence).upper()
    length = len(sequence)
    if length == 0:
        return 100.0, 0, 0
    ambiguous = sum(1 for base in sequence if base == "N" or base == "-")
    return ambiguous / length * 100.0, ambiguous, length


def load_consensus_typing(path: Path) -> dict[str, dict[str, str]]:
    rows = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["sequence_name"]] = row
    return rows


def find_metadata(record_id: str, metadata: dict[str, dict[str, str]]) -> dict[str, str] | None:
    if record_id in metadata:
        return metadata[record_id]
    parts = record_id.split("_")
    if len(parts) >= 2 and parts[1] in metadata:
        return metadata[parts[1]]
    return None


def clone_record(record: SeqRecord, new_id: str, description: str = "") -> SeqRecord:
    new_record = SeqRecord(record.seq, id=new_id, name=new_id, description=description or new_id)
    new_record.annotations = dict(record.annotations)
    return new_record


def trim_alignment_terminal_columns(records: list[SeqRecord], min_coverage: float) -> list[SeqRecord]:
    if not records or min_coverage <= 0:
        return records
    length = len(records[0].seq)
    if length == 0:
        return records
    if any(len(record.seq) != length for record in records):
        raise ValueError("Alignment records must have equal lengths before trimming")

    coverages = []
    denom = len(records)
    for index in range(length):
        covered = 0
        for record in records:
            if str(record.seq[index]).upper() in VALID_BASES:
                covered += 1
        coverages.append(covered / denom)

    start = 0
    while start < length and coverages[start] < min_coverage:
        start += 1
    end = length
    while end > start and coverages[end - 1] < min_coverage:
        end -= 1
    if start == 0 and end == length:
        return records
    return [record[start:end] for record in records]


def run_command(command: list[str], log_path: Path, env: dict[str, str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write("Command: " + " ".join(command) + "\n\n")
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(command)}. See {log_path}")


def run_mafft(input_fasta: Path, output_fasta: Path, threads: int, log_path: Path, large: bool) -> None:
    if large:
        command = ["mafft", "--quiet", "--thread", str(threads), "--retree", "1", "--maxiterate", "0", str(input_fasta)]
    else:
        command = ["mafft", "--quiet", "--thread", str(threads), "--auto", str(input_fasta)]
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    with output_fasta.open("w") as output, log_path.open("w") as log:
        log.write("Command: " + " ".join(command) + "\n\n")
        result = subprocess.run(command, stdout=output, stderr=log, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed ({result.returncode}) for {input_fasta}. See {log_path}")


def run_fasttree(alignment: Path, tree: Path, threads: int, log_path: Path, fastest: bool = False) -> None:
    command = ["FastTree", "-quiet", "-nt", "-gtr", "-nosupport"]
    if fastest:
        command.append("-fastest")
    command.extend(["-out", str(tree), str(alignment)])
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    run_command(command, log_path, env)


def write_records(records: list[SeqRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "fasta")


def tree_splits(tree, taxon_names: set[str]) -> set[frozenset[str]]:
    splits = set()
    n_taxa = len(taxon_names)
    if n_taxa < 4:
        return splits
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            continue
        side = frozenset(term.name for term in clade.get_terminals() if term.name in taxon_names)
        if len(side) < 2 or len(side) > n_taxa - 2:
            continue
        other = frozenset(taxon_names - side)
        splits.add(side if len(side) <= len(other) else other)
    return splits


def robinson_foulds(raw_tree_path: Path, combined_tree_path: Path, raw_ids: list[str]) -> tuple[int | None, float | None]:
    if len(raw_ids) < 4:
        return None, None
    raw_taxa = set(raw_ids)
    raw_tree = Phylo.read(raw_tree_path, "newick")
    combined_tree = Phylo.read(combined_tree_path, "newick")
    raw_splits = tree_splits(raw_tree, raw_taxa)
    combined_splits = tree_splits(combined_tree, raw_taxa)
    if not raw_splits and not combined_splits:
        return 0, 0.0
    rf = len(raw_splits.symmetric_difference(combined_splits))
    denom = len(raw_splits) + len(combined_splits)
    return rf, (rf / denom if denom else 0.0)


def sampled_patristic_correlation(raw_tree_path: Path, combined_tree_path: Path, raw_ids: list[str], max_pairs: int) -> float | None:
    if len(raw_ids) < 3:
        return None
    raw_tree = Phylo.read(raw_tree_path, "newick")
    combined_tree = Phylo.read(combined_tree_path, "newick")
    pairs = []
    for i in range(len(raw_ids)):
        for j in range(i + 1, len(raw_ids)):
            pairs.append((raw_ids[i], raw_ids[j]))
    if len(pairs) > max_pairs:
        rng = np.random.default_rng(20260615)
        idx = rng.choice(len(pairs), size=max_pairs, replace=False)
        pairs = [pairs[i] for i in idx]
    raw_dist = []
    combined_dist = []
    for a, b in pairs:
        try:
            raw_dist.append(raw_tree.distance(a, b))
            combined_dist.append(combined_tree.distance(a, b))
        except Exception:
            continue
    if len(raw_dist) < 3:
        return None
    if np.std(raw_dist) == 0 or np.std(combined_dist) == 0:
        return None
    return float(np.corrcoef(raw_dist, combined_dist)[0, 1])


def nearest_consensus_metrics(tree_path: Path, raw_ids: list[str], consensus_ids: list[str]) -> tuple[list[dict[str, str]], dict[str, float | int | None]]:
    if not raw_ids or not consensus_ids:
        return [], {
            "mean_nearest_distance": None,
            "median_nearest_distance": None,
            "q90_nearest_distance": None,
            "q95_nearest_distance": None,
            "max_nearest_distance": None,
            "coverage_d02": None,
            "coverage_d05": None,
            "coverage_d10": None,
        }
    tree = Phylo.read(tree_path, "newick")
    rows = []
    distances = []
    counts = Counter()
    for raw_id in raw_ids:
        best_id = None
        best_distance = math.inf
        for consensus_id in consensus_ids:
            try:
                distance = tree.distance(raw_id, consensus_id)
            except Exception:
                continue
            if distance < best_distance:
                best_distance = distance
                best_id = consensus_id
        if best_id is not None:
            counts[best_id] += 1
            distances.append(best_distance)
            rows.append(
                {
                    "raw_tree_id": raw_id,
                    "nearest_consensus_tree_id": best_id,
                    "nearest_distance": f"{best_distance:.8f}",
                }
            )
    if not distances:
        return rows, {}
    array = np.array(distances, dtype=float)
    metrics = {
        "mean_nearest_distance": float(np.mean(array)),
        "median_nearest_distance": float(np.median(array)),
        "q90_nearest_distance": float(np.quantile(array, 0.90)),
        "q95_nearest_distance": float(np.quantile(array, 0.95)),
        "max_nearest_distance": float(np.max(array)),
        "coverage_d02": float(np.mean(array <= 0.02)),
        "coverage_d05": float(np.mean(array <= 0.05)),
        "coverage_d10": float(np.mean(array <= 0.10)),
    }
    for row in rows:
        row["nearest_consensus_raw_count"] = str(counts[row["nearest_consensus_tree_id"]])
    return rows, metrics


def alignment_p_distance(seq1: str, seq2: str) -> float | None:
    """Pairwise p-distance using positions where both records are non-gap A/C/G/T."""
    matches = 0
    total = 0
    for base1, base2 in zip(seq1.upper(), seq2.upper()):
        if base1 not in VALID_BASES or base2 not in VALID_BASES:
            continue
        total += 1
        if base1 == base2:
            matches += 1
    if total == 0:
        return None
    return 1.0 - (matches / total)


def nearest_consensus_alignment_metrics(alignment_path: Path, raw_ids: list[str], consensus_ids: list[str]) -> tuple[list[dict[str, str]], dict[str, float | int | None]]:
    """Nearest-consensus metrics from an existing MSA when tree inference is unavailable."""
    if not raw_ids or not consensus_ids:
        return [], {
            "mean_nearest_distance": None,
            "median_nearest_distance": None,
            "q90_nearest_distance": None,
            "q95_nearest_distance": None,
            "max_nearest_distance": None,
            "coverage_d02": None,
            "coverage_d05": None,
            "coverage_d10": None,
        }

    seqs = {record.id: str(record.seq) for record in SeqIO.parse(alignment_path, "fasta")}
    rows = []
    distances = []
    counts = Counter()
    for raw_id in raw_ids:
        raw_seq = seqs.get(raw_id)
        if raw_seq is None:
            continue
        best_id = None
        best_distance = math.inf
        for consensus_id in consensus_ids:
            consensus_seq = seqs.get(consensus_id)
            if consensus_seq is None:
                continue
            distance = alignment_p_distance(raw_seq, consensus_seq)
            if distance is None:
                continue
            if distance < best_distance:
                best_distance = distance
                best_id = consensus_id
        if best_id is not None:
            counts[best_id] += 1
            distances.append(best_distance)
            rows.append(
                {
                    "raw_tree_id": raw_id,
                    "nearest_consensus_tree_id": best_id,
                    "nearest_distance": f"{best_distance:.8f}",
                }
            )
    if not distances:
        return rows, {}
    array = np.array(distances, dtype=float)
    metrics = {
        "mean_nearest_distance": float(np.mean(array)),
        "median_nearest_distance": float(np.median(array)),
        "q90_nearest_distance": float(np.quantile(array, 0.90)),
        "q95_nearest_distance": float(np.quantile(array, 0.95)),
        "max_nearest_distance": float(np.max(array)),
        "coverage_d02": float(np.mean(array <= 0.02)),
        "coverage_d05": float(np.mean(array <= 0.05)),
        "coverage_d10": float(np.mean(array <= 0.10)),
    }
    for row in rows:
        row["nearest_consensus_raw_count"] = str(counts[row["nearest_consensus_tree_id"]])
    return rows, metrics


def prepare_inputs(args):
    gi_meta = read_metadata(args.gi_metadata)
    gii_meta = read_metadata(args.gii_metadata)
    metadata = {**gi_meta, **gii_meta}
    consensus_typing = load_consensus_typing(args.consensus_typing)

    out_renamed = args.outdir / "01_renamed_consensus"
    out_renamed.mkdir(parents=True, exist_ok=True)
    qc_dir = args.outdir / "00_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    excluded_raw_rows = []

    groups = defaultdict(lambda: {"raw": [], "consensus": [], "map": []})
    renamed_consensus_all = []
    renamed_consensus_by_group = defaultdict(list)

    for label, fasta in [("gi", args.gi_raw_fasta), ("gii", args.gii_raw_fasta)]:
        for record in SeqIO.parse(fasta, "fasta"):
            original_id = record.id
            percent, ambiguous, length = n_gap_percent(str(record.seq))
            if percent > args.raw_qc_max_n_gap_percent:
                excluded_raw_rows.append(
                    {
                        "group": label,
                        "sequence_id": original_id,
                        "n_gap_percent": f"{percent:.6f}",
                        "n_gap_count": str(ambiguous),
                        "length": str(length),
                        "reason": f"N_or_gap_percent_gt_{args.raw_qc_max_n_gap_percent}",
                    }
                )
                continue
            meta = find_metadata(original_id, metadata)
            if not meta:
                continue
            vp1 = meta.get("vp1_genotype", "").strip()
            if not vp1:
                continue
            tree_id = safe_id(f"RAW__{original_id}")
            new_record = clone_record(record, tree_id)
            groups[vp1]["raw"].append(new_record)
            groups[vp1]["map"].append(
                {
                    "tree_id": tree_id,
                    "source": "raw",
                    "original_id": original_id,
                    "vp1_group": vp1,
                    "vp1_genotype": vp1,
                    "rdrp_genotype": "",
                    "accession": meta.get("accession", ""),
                    "collection_date": meta.get("collection_date", ""),
                    "length": str(len(record.seq)),
                }
            )

    for label, fasta in [("gi", args.gi_consensus_fasta), ("gii", args.gii_consensus_fasta)]:
        label_records = []
        for record in SeqIO.parse(fasta, "fasta"):
            original_id = record.id
            row = consensus_typing.get(original_id)
            if not row:
                continue
            vp1 = row.get("vp1_genotype", "").strip() or "Unknown"
            rdrp = row.get("rdrp_genotype", "").strip() or "Unknown"
            name_genotype = row.get("name_genotype", "").strip()
            vp1_group = vp1
            if vp1 in {"Unknown", "NA", ""} and name_genotype:
                vp1_group = name_genotype
            new_id = safe_id(f"CONSENSUS__VP1-{vp1}__RdRp-{rdrp}__{original_id}")
            description = f"{new_id} original={original_id} vp1={vp1} rdrp={rdrp}"
            new_record = clone_record(record, new_id, description)
            groups[vp1_group]["consensus"].append(new_record)
            groups[vp1_group]["map"].append(
                {
                    "tree_id": new_id,
                    "source": "consensus",
                    "original_id": original_id,
                    "vp1_group": vp1_group,
                    "vp1_genotype": vp1,
                    "rdrp_genotype": rdrp,
                    "accession": "",
                    "collection_date": "",
                    "length": str(len(record.seq)),
                }
            )
            renamed_consensus_all.append(new_record)
            renamed_consensus_by_group[label].append(new_record)
            label_records.append(new_record)

        write_records(label_records, out_renamed / f"{label}_final_consensus.vp1_rdrp_renamed.fasta")

    write_records(renamed_consensus_all, out_renamed / "all_final_consensus.vp1_rdrp_renamed.fasta")
    excluded_path = qc_dir / "raw_sequence_qc_excluded.tsv"
    with excluded_path.open("w", newline="") as handle:
        fieldnames = ["group", "sequence_id", "n_gap_percent", "n_gap_count", "length", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded_raw_rows)
    return groups


def subset_alignment(alignment_path: Path, keep_ids: set[str], output_path: Path) -> None:
    records = [record for record in SeqIO.parse(alignment_path, "fasta") if record.id in keep_ids]
    if records:
        write_records(records, output_path)


def analyze_group(genotype: str, payload: dict, args) -> tuple[dict[str, str], list[dict[str, str]]]:
    group_dir = args.outdir / "02_by_vp1" / safe_id(genotype)
    group_dir.mkdir(parents=True, exist_ok=True)
    raw_records = payload["raw"]
    consensus_records = payload["consensus"]
    all_records = raw_records + consensus_records
    raw_ids = [record.id for record in raw_records]
    consensus_ids = [record.id for record in consensus_records]

    map_path = group_dir / f"{safe_id(genotype)}_id_map.tsv"
    with map_path.open("w", newline="") as handle:
        fieldnames = ["tree_id", "source", "original_id", "vp1_group", "vp1_genotype", "rdrp_genotype", "accession", "collection_date", "length"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(payload["map"])

    raw_fasta = group_dir / f"{safe_id(genotype)}_raw.fasta"
    consensus_fasta = group_dir / f"{safe_id(genotype)}_consensus.fasta"
    all_fasta = group_dir / f"{safe_id(genotype)}_raw_plus_consensus.fasta"
    write_records(raw_records, raw_fasta)
    write_records(consensus_records, consensus_fasta)
    write_records(all_records, all_fasta)

    summary = {
        "vp1_group": genotype,
        "raw_n": str(len(raw_records)),
        "consensus_n": str(len(consensus_records)),
        "total_n": str(len(all_records)),
        "status": "prepared",
        "group_dir": str(group_dir),
        "raw_plus_consensus_fasta": str(all_fasta),
        "alignment": "",
        "combined_tree": "",
        "raw_tree": "",
        "consensus_tree": "",
        "rf_distance": "",
        "rf_normalized": "",
        "patristic_correlation": "",
        "distance_method": "",
        "mean_nearest_distance": "",
        "median_nearest_distance": "",
        "q90_nearest_distance": "",
        "q95_nearest_distance": "",
        "max_nearest_distance": "",
        "coverage_d02": "",
        "coverage_d05": "",
        "coverage_d10": "",
    }

    if len(all_records) < 3:
        summary["status"] = "skipped_total_lt3"
        return summary, []
    if not raw_records or not consensus_records:
        summary["status"] = "skipped_missing_raw_or_consensus"
        return summary, []

    alignment = group_dir / f"{safe_id(genotype)}_raw_plus_consensus.aligned.fasta"
    trimmed_alignment = group_dir / f"{safe_id(genotype)}_raw_plus_consensus.aligned.trimmed.fasta"
    combined_tree = group_dir / f"{safe_id(genotype)}_raw_plus_consensus.fasttree.nwk"
    raw_alignment = group_dir / f"{safe_id(genotype)}_raw.aligned_from_combined.fasta"
    raw_tree = group_dir / f"{safe_id(genotype)}_raw.fasttree.nwk"
    consensus_alignment = group_dir / f"{safe_id(genotype)}_consensus.aligned_from_combined.fasta"
    consensus_tree = group_dir / f"{safe_id(genotype)}_consensus.fasttree.nwk"

    large = len(all_records) >= args.large_mafft_threshold
    fastest = len(all_records) >= args.fasttree_fastest_threshold

    use_fasttree = getattr(args, "use_fasttree", False)
    if args.skip_existing and trimmed_alignment.exists() and (not use_fasttree or combined_tree.exists()):
        pass
    else:
        run_mafft(all_fasta, alignment, args.threads, group_dir / "mafft.log", large=large)
        aligned_records = list(SeqIO.parse(alignment, "fasta"))
        trimmed_records = trim_alignment_terminal_columns(aligned_records, args.trim_end_coverage)
        write_records(trimmed_records, trimmed_alignment)
        subset_alignment(trimmed_alignment, set(raw_ids), raw_alignment)
        subset_alignment(trimmed_alignment, set(consensus_ids), consensus_alignment)
        if use_fasttree:
            run_fasttree(trimmed_alignment, combined_tree, args.threads, group_dir / "fasttree.combined.log", fastest=fastest)
            if len(raw_records) >= 3:
                run_fasttree(raw_alignment, raw_tree, args.threads, group_dir / "fasttree.raw.log", fastest=fastest)
            if len(consensus_records) >= 3:
                run_fasttree(consensus_alignment, consensus_tree, args.threads, group_dir / "fasttree.consensus.log", fastest=False)

    summary["status"] = "analyzed"
    summary["alignment"] = str(trimmed_alignment)
    summary["combined_tree"] = str(combined_tree if combined_tree.exists() else "")
    summary["raw_tree"] = str(raw_tree if raw_tree.exists() else "")
    summary["consensus_tree"] = str(consensus_tree if consensus_tree.exists() else "")

    if use_fasttree and combined_tree.exists():
        summary["distance_method"] = "tree_patristic"
        nearest_rows, nearest_metrics = nearest_consensus_metrics(combined_tree, raw_ids, consensus_ids)
    else:
        summary["distance_method"] = "alignment_p_distance"
        nearest_rows, nearest_metrics = nearest_consensus_alignment_metrics(trimmed_alignment, raw_ids, consensus_ids)
    for key, value in nearest_metrics.items():
        if value is not None:
            summary[key] = f"{value:.8f}" if isinstance(value, float) else str(value)

    if raw_tree.exists():
        rf, rf_norm = robinson_foulds(raw_tree, combined_tree, raw_ids)
        if rf is not None:
            summary["rf_distance"] = str(rf)
            summary["rf_normalized"] = f"{rf_norm:.8f}"
        corr = sampled_patristic_correlation(raw_tree, combined_tree, raw_ids, args.max_distance_pairs)
        if corr is not None:
            summary["patristic_correlation"] = f"{corr:.8f}"

    for row in nearest_rows:
        row["vp1_group"] = genotype
    return summary, nearest_rows


def write_summary_tables(args, summaries, nearest_rows):
    metrics_dir = args.outdir / "03_metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    summary_path = metrics_dir / "representativeness_summary.tsv"
    fields = [
        "vp1_group",
        "raw_n",
        "consensus_n",
        "total_n",
        "status",
        "rf_distance",
        "rf_normalized",
        "patristic_correlation",
        "distance_method",
        "mean_nearest_distance",
        "median_nearest_distance",
        "q90_nearest_distance",
        "q95_nearest_distance",
        "max_nearest_distance",
        "coverage_d02",
        "coverage_d05",
        "coverage_d10",
        "group_dir",
        "raw_plus_consensus_fasta",
        "alignment",
        "combined_tree",
        "raw_tree",
        "consensus_tree",
    ]
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)

    nearest_path = metrics_dir / "nearest_consensus_distances.tsv"
    nearest_fields = ["vp1_group", "raw_tree_id", "nearest_consensus_tree_id", "nearest_distance", "nearest_consensus_raw_count"]
    with nearest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=nearest_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(nearest_rows)
    return summary_path, nearest_path


def plot_results(args, summaries, nearest_rows):
    figures_dir = args.outdir / "04_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    analyzed = [row for row in summaries if row["status"] == "analyzed"]
    if not analyzed:
        return

    labels = [row["vp1_group"] for row in sorted(analyzed, key=lambda r: genotype_sort_key(r["vp1_group"]))]
    raw_n = [int(next(row for row in analyzed if row["vp1_group"] == label)["raw_n"]) for label in labels]
    cons_n = [int(next(row for row in analyzed if row["vp1_group"] == label)["consensus_n"]) for label in labels]
    median_dist = [
        float(next(row for row in analyzed if row["vp1_group"] == label)["median_nearest_distance"] or "nan")
        for label in labels
    ]

    x = np.arange(len(labels))
    fig, ax1 = plt.subplots(figsize=(max(10, len(labels) * 0.35), 5))
    ax1.bar(x - 0.18, raw_n, width=0.36, label="Downloaded raw sequences", color="#8fb3d9")
    ax1.bar(x + 0.18, cons_n, width=0.36, label="Consensus sequences", color="#e07a5f")
    ax1.set_yscale("log")
    ax1.set_ylabel("Sequence count (log scale)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=75, ha="right")
    ax2 = ax1.twinx()
    ax2.plot(x, median_dist, color="#2f3e46", marker="o", linewidth=1.5, label="Median nearest consensus distance")
    ax2.set_ylabel("Median distance to nearest consensus")
    lines, line_labels = ax1.get_legend_handles_labels()
    lines2, line_labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, line_labels + line_labels2, loc="upper right")
    fig.tight_layout()
    fig.savefig(figures_dir / "summary_counts_and_nearest_distance.png", dpi=220)
    plt.close(fig)

    distances_by_group = defaultdict(list)
    for row in nearest_rows:
        distances_by_group[row["vp1_group"]].append(float(row["nearest_distance"]))
    plot_labels = [label for label in labels if distances_by_group[label]]
    if plot_labels:
        data = [distances_by_group[label] for label in plot_labels]
        fig, ax = plt.subplots(figsize=(max(10, len(plot_labels) * 0.35), 5))
        ax.boxplot(data, labels=plot_labels, showfliers=False)
        ax.set_ylabel("Distance to nearest consensus")
        ax.set_title("Raw sequence coverage by nearest consensus sequence")
        ax.tick_params(axis="x", rotation=75)
        fig.tight_layout()
        fig.savefig(figures_dir / "nearest_consensus_distance_by_vp1.png", dpi=220)
        plt.close(fig)

    rf_rows = [row for row in analyzed if row.get("rf_normalized")]
    if rf_rows:
        rf_labels = [row["vp1_group"] for row in sorted(rf_rows, key=lambda r: genotype_sort_key(r["vp1_group"]))]
        rf_vals = [float(next(row for row in rf_rows if row["vp1_group"] == label)["rf_normalized"]) for label in rf_labels]
        corr_vals = [
            float(next(row for row in rf_rows if row["vp1_group"] == label)["patristic_correlation"] or "nan")
            for label in rf_labels
        ]
        x = np.arange(len(rf_labels))
        fig, ax1 = plt.subplots(figsize=(max(10, len(rf_labels) * 0.35), 5))
        ax1.bar(x, rf_vals, color="#81b29a", label="Normalized RF distance")
        ax1.set_ylabel("Normalized RF distance")
        ax1.set_ylim(0, 1.05)
        ax1.set_xticks(x)
        ax1.set_xticklabels(rf_labels, rotation=75, ha="right")
        ax2 = ax1.twinx()
        ax2.plot(x, corr_vals, color="#3d405b", marker="o", label="Patristic distance correlation")
        ax2.set_ylabel("Patristic distance correlation")
        ax2.set_ylim(-0.05, 1.05)
        lines, line_labels = ax1.get_legend_handles_labels()
        lines2, line_labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines + lines2, line_labels + line_labels2, loc="upper right")
        fig.tight_layout()
        fig.savefig(figures_dir / "tree_topology_stability_by_vp1.png", dpi=220)
        plt.close(fig)


def write_report(args, summaries, summary_path: Path, nearest_path: Path) -> None:
    report = args.outdir / "representativeness_report.md"
    analyzed = [row for row in summaries if row["status"] == "analyzed"]
    skipped = [row for row in summaries if row["status"] != "analyzed"]
    best = sorted(
        analyzed,
        key=lambda row: float(row["median_nearest_distance"] or "inf"),
    )[:5]
    worst = sorted(
        analyzed,
        key=lambda row: float(row["median_nearest_distance"] or "-inf"),
        reverse=True,
    )[:5]
    lines = [
        "# Consensus Representativeness Analysis",
        "",
        "## Inputs",
        "",
        f"- GI raw FASTA: `{args.gi_raw_fasta}`",
        f"- GII raw FASTA: `{args.gii_raw_fasta}`",
        f"- Consensus genotyping table: `{args.consensus_typing}`",
        "",
        "## Method",
        "",
        "Consensus FASTA identifiers were rewritten as `CONSENSUS__VP1-<VP1>__RdRp-<RdRp>__<old_id>`.",
        f"Downloaded raw sequences with combined N/gap content > {args.raw_qc_max_n_gap_percent}% were excluded before representativeness analysis, matching the consensus-generation QC rule.",
        "Remaining downloaded raw sequences and consensus sequences were split by VP1 genotype.",
        "For each genotype with both raw and consensus records, raw+consensus sequences were aligned with MAFFT; terminal columns with A/C/G/T coverage below the configured threshold were removed. If FastTree is available, approximate maximum-likelihood nucleotide trees are inferred; otherwise nearest-consensus distances are calculated directly from the trimmed alignment as pairwise p-distance.",
        "",
        "Representativeness was evaluated using:",
        "",
        "1. nearest-consensus distance for every downloaded raw sequence, using either combined-tree patristic distance or alignment p-distance depending on local tool availability;",
        "2. normalized Robinson-Foulds distance between the raw-only tree and the combined tree after considering raw taxa only;",
        "3. sampled pairwise patristic-distance correlation between the raw-only and combined trees.",
        "",
        "A direct RF comparison between a consensus-only tree and an all-sequence tree was not used as a primary metric because RF distance requires the same taxon set.",
        "",
        "## Existing Tools Considered",
        "",
        "- Tree distance packages such as TreeDist, phangorn, DendroPy, and ETE/ETE3 are appropriate for Robinson-Foulds or related tree-distance metrics when the compared trees share the same taxa. In this analysis, the meaningful topology comparison is therefore the raw-only tree versus the raw taxa induced by the raw+consensus tree, not consensus-only versus all-sequence directly.",
        "- Phylogenetic placement tools such as pplacer and EPA-ng are useful when placing query sequences onto a fixed reference tree. They would be a good future scaling option if the full genotype trees become too large to rebuild repeatedly.",
        "- Local software availability is checked at runtime. When FastTree is unavailable, the analysis still reports alignment-based nearest-consensus p-distance metrics, while RF and patristic-correlation fields remain blank.",
        "",
        "## Outputs",
        "",
        f"- Summary metrics: `{summary_path}`",
        f"- Per-sequence nearest consensus table: `{nearest_path}`",
        f"- Renamed consensus FASTA directory: `{args.outdir / '01_renamed_consensus'}`",
        f"- Per-genotype alignment/tree directory: `{args.outdir / '02_by_vp1'}`",
        f"- Figures: `{args.outdir / '04_figures'}`",
        f"- Consensus-highlighted tree figures: `{args.outdir / '04_figures' / 'highlighted_trees'}`",
        "",
        "## Run Summary",
        "",
        f"- Analyzed VP1 genotype groups: {len(analyzed)}",
        f"- Skipped groups: {len(skipped)}",
        f"- Threads: {args.threads}",
        f"- Terminal alignment coverage threshold: {args.trim_end_coverage}",
        "",
        "## Most Central Consensus Coverage",
        "",
        "| VP1 genotype | raw n | consensus n | median nearest distance | normalized RF | patristic corr |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in best:
        lines.append(
            f"| {row['vp1_group']} | {row['raw_n']} | {row['consensus_n']} | {row['median_nearest_distance']} | {row['rf_normalized']} | {row['patristic_correlation']} |"
        )
    lines.extend(
        [
            "",
            "## Groups Requiring More Attention",
            "",
            "| VP1 genotype | raw n | consensus n | median nearest distance | q95 nearest distance | normalized RF |",
            "| --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in worst:
        lines.append(
            f"| {row['vp1_group']} | {row['raw_n']} | {row['consensus_n']} | {row['median_nearest_distance']} | {row['q95_nearest_distance']} | {row['rf_normalized']} |"
        )
    if skipped:
        lines.extend(["", "## Skipped Groups", "", "| VP1 genotype | raw n | consensus n | reason |", "| --- | ---: | ---: | --- |"])
        for row in skipped:
            lines.append(f"| {row['vp1_group']} | {row['raw_n']} | {row['consensus_n']} | {row['status']} |")
    report.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gi-raw-fasta", type=Path, default=Path("ncbi_download_2026-06-11/ncbi_norovirus_gi_gt5800_collected_2000_onward.fasta"))
    parser.add_argument("--gii-raw-fasta", type=Path, default=Path("ncbi_download_2026-06-11/ncbi_norovirus_gii_gt5800_collected_2000_onward.fasta"))
    parser.add_argument("--gi-metadata", type=Path, default=Path("ncbi_download_2026-06-11/ncbi_norovirus_gi_gt5800_metadata_vp1_typed.tsv"))
    parser.add_argument("--gii-metadata", type=Path, default=Path("ncbi_download_2026-06-11/ncbi_norovirus_gii_gt5800_metadata_vp1_typed.tsv"))
    parser.add_argument("--gi-consensus-fasta", type=Path, default=Path("results_new/gi/07_final_results/all_final_consensus.fasta"))
    parser.add_argument("--gii-consensus-fasta", type=Path, default=Path("results_new/gii/07_final_results/all_final_consensus.fasta"))
    parser.add_argument("--consensus-typing", type=Path, default=Path("results_new/genotyping_consensus/consensus_vp1_rdrp_genotyping.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("results_new/representativeness_analysis"))
    parser.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 8) - 2))
    parser.add_argument("--trim-end-coverage", type=float, default=0.5)
    parser.add_argument("--large-mafft-threshold", type=int, default=500)
    parser.add_argument("--fasttree-fastest-threshold", type=int, default=1500)
    parser.add_argument("--max-distance-pairs", type=int, default=20000)
    parser.add_argument("--raw-qc-max-n-gap-percent", type=float, default=DEFAULT_RAW_QC_MAX_N_GAP_PERCENT)
    parser.add_argument("--skip-existing", action="store_true")
    args = parser.parse_args()

    if args.threads < 1:
        raise SystemExit("--threads must be a positive integer")
    if not shutil.which("mafft"):
        raise SystemExit("mafft not found in PATH")
    args.use_fasttree = bool(shutil.which("FastTree"))
    if not args.use_fasttree:
        print("FastTree not found in PATH; using alignment p-distance metrics and skipping tree topology metrics.", flush=True)

    args.outdir.mkdir(parents=True, exist_ok=True)
    groups = prepare_inputs(args)
    summaries = []
    nearest_rows = []
    for genotype in sorted(groups, key=genotype_sort_key):
        payload = groups[genotype]
        if not payload["raw"] and not payload["consensus"]:
            continue
        print(f"[{genotype}] raw={len(payload['raw'])} consensus={len(payload['consensus'])}", flush=True)
        summary, nearest = analyze_group(genotype, payload, args)
        summaries.append(summary)
        nearest_rows.extend(nearest)

    summary_path, nearest_path = write_summary_tables(args, summaries, nearest_rows)
    plot_results(args, summaries, nearest_rows)
    write_report(args, summaries, summary_path, nearest_path)
    print(f"Summary: {summary_path}")
    print(f"Nearest consensus table: {nearest_path}")
    print(f"Report: {args.outdir / 'representativeness_report.md'}")


if __name__ == "__main__":
    main()
