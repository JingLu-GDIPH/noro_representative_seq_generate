#!/usr/bin/env python3
"""Build wastewater mapping references for one RdRp_VP1 alignment.

The script consumes an orientation-normalized, column-masked alignment and the
matching IQ-TREE outputs. It selects supported phylogenetic boundary clades,
builds hybrid ASR/majority references with medoid fallback, and iteratively
merges references that cannot be distinguished by enough 150 nt windows.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path

import numpy as np
from Bio import AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BASES = frozenset("ACGT")
BASE_ARRAY = np.asarray(list("ACGT"))
SUPPORT_RE = re.compile(r"^(Node\d+)(?:/([0-9.]+)/([0-9.]+))?$")
RAW_ID_RE = re.compile(
    r"^(G(?:I|II|IX)\.P[A-Za-z0-9]+_G(?:I|II|IX)\.[A-Za-z0-9]+)_"
    r"([A-Z]{1,3}\d{4,8})_"
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--tree", type=Path)
    parser.add_argument("--state", type=Path)
    parser.add_argument("--mldist", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--bowtie2", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-tips", type=int, default=3)
    parser.add_argument("--min-sh-alrt", type=float, default=80.0)
    parser.add_argument("--min-ufboot", type=float, default=95.0)
    parser.add_argument("--max-within-d95", type=float, default=0.05)
    parser.add_argument(
        "--max-within-max",
        type=float,
        default=1.0,
        help=(
            "Optional hard ceiling for the maximum pair distance. The default "
            "disables this outlier-sensitive filter; d95 remains the primary "
            "within-clade criterion."
        ),
    )
    parser.add_argument("--min-parent-d95", type=float, default=0.05)
    parser.add_argument("--min-parent-delta", type=float, default=0.01)
    parser.add_argument("--asr-pp", type=float, default=0.90)
    parser.add_argument("--majority-frequency", type=float, default=0.60)
    parser.add_argument("--min-site-coverage", type=float, default=0.50)
    parser.add_argument("--max-low-pp-fraction", type=float, default=0.05)
    parser.add_argument("--max-n-fraction", type=float, default=0.005)
    parser.add_argument("--max-n-run", type=int, default=15)
    parser.add_argument("--window", type=int, default=150)
    parser.add_argument("--step", type=int, default=25)
    parser.add_argument("--min-aligned-window", type=int, default=135)
    parser.add_argument("--min-score-margin", type=int, default=6)
    parser.add_argument("--min-nm-margin", type=int, default=2)
    parser.add_argument("--min-target-mapping-rate", type=float, default=0.95)
    parser.add_argument("--min-distinguishable-rate", type=float, default=0.30)
    parser.add_argument(
        "--umbrella-own-best-floor",
        type=float,
        default=0.5,
        help=(
            "A merge-flagged unit whose own-best rate falls below this floor is "
            "treated as an umbrella (most member reads prefer other references) "
            "and split instead of merged."
        ),
    )
    parser.add_argument(
        "--umbrella-gate",
        choices=["members", "structural"],
        default="members",
        help=(
            "Umbrella split-routing gate: 'members' requires at least "
            "umbrella-min-members members; 'structural' requires the unit to "
            "have >=2 child clades each holding >= min-tips of its own members "
            "(i.e. a split can produce viable children)."
        ),
    )
    parser.add_argument(
        "--umbrella-min-members",
        type=int,
        default=25,
        help=(
            "Umbrella splitting applies only to units with at least this many "
            "members; smaller merge-flagged units are treated as cloud-diluted "
            "duplicates and routed to merge."
        ),
    )
    parser.add_argument("--merge-min-similarity", type=float, default=0.95)
    parser.add_argument("--max-iterations", type=int, default=20)
    parser.add_argument("--max-distance-pairs", type=int, default=20000)
    parser.add_argument(
        "--vp1-reference-fasta",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "reference"
        / "gii_with_genotype.fasta",
        help="Annotated reference FASTA used to locate the VP1 region.",
    )
    parser.add_argument(
        "--vp1-coords",
        type=Path,
        default=Path(__file__).resolve().parent.parent
        / "reference"
        / "genome_gene_coordinates.tsv",
        help="Gene coordinate table with VP1 start/end on the reference accessions.",
    )
    parser.add_argument(
        "--veto-block-windows",
        type=int,
        default=4,
        help="Consecutive VP1 windows (each >= veto-window-snv SNVs) that veto a merge.",
    )
    parser.add_argument(
        "--veto-window-snv",
        type=int,
        default=2,
        help="SNVs a VP1 window must carry to count as divergent for the merge veto.",
    )
    parser.add_argument(
        "--veto-min-window-coverage",
        type=float,
        default=0.5,
        help=(
            "Minimum fraction of comparable A/C/G/T columns for a veto window to be "
            "evaluable; guards against partial-genome sequences."
        ),
    )
    parser.add_argument(
        "--no-vp1-veto",
        action="store_true",
        help="Disable the VP1 segment-level merge veto.",
    )
    parser.add_argument(
        "--shadow-sweep-distance",
        type=float,
        default=0.005,
        help=(
            "After the loop stops, collapse any residual reference pairs closer "
            "than this ACGT distance (shadow references) unless the VP1 veto "
            "protects them. Set 0 to disable the final sweep."
        ),
    )
    parser.add_argument(
        "--representative",
        choices=["asr", "medoid"],
        default="asr",
        help=(
            "Reference reconstruction: 'asr' keeps the hybrid ASR/majority "
            "chain; 'medoid' emits the real observed member closest (patristic) "
            "to the clade root node, preferring QC-clean candidates."
        ),
    )
    return parser.parse_args()


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def genotype_pair(sequence_id: str) -> str:
    parts = sequence_id.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else "UnknownRdRp_UnknownVP1"


def support(node) -> tuple[str, float, float]:
    match = SUPPORT_RE.match(node.name or "")
    if not match:
        return node.name or "Root", 0.0, 0.0
    return (
        match.group(1),
        float(match.group(2) or 0.0),
        float(match.group(3) or 0.0),
    )


def load_state(path: Path | None) -> dict[str, dict[str, np.ndarray]]:
    if path is None or not path.exists():
        return {}
    states = defaultdict(list)
    posterior = defaultdict(list)
    with path.open() as handle:
        rows = (line for line in handle if not line.startswith("#"))
        for row in csv.DictReader(rows, delimiter="\t"):
            states[row["Node"]].append(row["State"].upper())
            posterior[row["Node"]].append(
                max(float(row[f"p_{base}"]) for base in "ACGT")
            )
    return {
        node: {
            "state": np.asarray(values),
            "posterior": np.asarray(posterior[node], dtype=float),
        }
        for node, values in states.items()
    }


def load_mldist(path: Path | None) -> tuple[list[str], np.ndarray] | None:
    if path is None or not path.exists():
        return None
    with path.open() as handle:
        size = int(handle.readline().strip())
        names = []
        matrix = np.zeros((size, size), dtype=np.float32)
        for row_index, line in enumerate(handle):
            fields = line.split()
            if len(fields) < size + 1:
                raise ValueError(f"Malformed IQ-TREE distance row: {line[:100]}")
            names.append(fields[0])
            matrix[row_index] = np.asarray(fields[1 : size + 1], dtype=np.float32)
    if len(names) != size:
        raise ValueError(f"Expected {size} IQ-TREE distance rows, found {len(names)}")
    return names, matrix


def distance_stats(
    member_ids: list[str],
    names_to_index: dict[str, int],
    matrix: np.ndarray,
    max_pairs: int,
) -> dict[str, float]:
    indices = np.asarray([names_to_index[name] for name in member_ids], dtype=int)
    count = len(indices)
    if count < 2:
        return {"median": 0.0, "d95": 0.0, "max": 0.0, "pairs": 0}
    row, col = np.triu_indices(count, 1)
    if len(row) > max_pairs:
        rng = np.random.default_rng(20260724)
        selected = rng.choice(len(row), max_pairs, replace=False)
        row, col = row[selected], col[selected]
    values = matrix[indices[row], indices[col]]
    return {
        "median": float(np.median(values)),
        "d95": float(np.quantile(values, 0.95)),
        "max": float(values.max()),
        "pairs": len(values),
    }


def fallback_distance_stats(
    member_ids: list[str],
    aligned_by_id: dict[str, str],
    max_pairs: int,
) -> dict[str, float]:
    pairs = [
        (left, right)
        for pos, left in enumerate(member_ids)
        for right in member_ids[pos + 1 :]
    ]
    if len(pairs) > max_pairs:
        rng = np.random.default_rng(20260724)
        pairs = [pairs[index] for index in rng.choice(len(pairs), max_pairs, replace=False)]
    values = [
        acgt_distance(aligned_by_id[left], aligned_by_id[right])
        for left, right in pairs
    ]
    values = np.asarray([value for value in values if math.isfinite(value)])
    if not len(values):
        return {"median": 0.0, "d95": 0.0, "max": 0.0, "pairs": 0}
    return {
        "median": float(np.median(values)),
        "d95": float(np.quantile(values, 0.95)),
        "max": float(values.max()),
        "pairs": len(values),
    }


def acgt_distance(left: str, right: str) -> float:
    if len(left) != len(right):
        return math.inf
    a = np.asarray(list(left.upper()))
    b = np.asarray(list(right.upper()))
    valid = np.isin(a, BASE_ARRAY) & np.isin(b, BASE_ARRAY)
    denominator = int(valid.sum())
    return float(((a != b) & valid).sum() / denominator) if denominator else math.inf


def majority_sequence(member_ids: set[str], aligned_by_id: dict[str, str]) -> str:
    matrix = np.asarray([list(aligned_by_id[name]) for name in sorted(member_ids)])
    output = []
    for column in matrix.T:
        counts = Counter(base for base in column if base in BASES)
        output.append(counts.most_common(1)[0][0] if counts else "N")
    return "".join(output)


def longest_n_run(sequence: str) -> int:
    return max((len(match.group()) for match in re.finditer(r"N+", sequence)), default=0)


def hybrid_asr(
    member_ids: set[str],
    aligned_by_id: dict[str, str],
    state: np.ndarray,
    posterior: np.ndarray,
    args: argparse.Namespace,
) -> tuple[str, dict[str, float]]:
    matrix = np.asarray([list(aligned_by_id[name]) for name in sorted(member_ids)])
    output = []
    fallback = 0
    unresolved = 0
    for position, column in enumerate(matrix.T):
        if posterior[position] >= args.asr_pp and state[position] in BASES:
            output.append(str(state[position]))
            continue
        fallback += 1
        valid = [base for base in column if base in BASES]
        counts = Counter(valid)
        if valid and len(valid) / len(column) >= args.min_site_coverage:
            base, count = counts.most_common(1)[0]
            if count / len(valid) >= args.majority_frequency:
                output.append(base)
                continue
        output.append("N")
        unresolved += 1
    sequence = "".join(output)
    length = max(1, len(sequence))
    return sequence, {
        "low_pp_fraction": fallback / length,
        "n_fraction": unresolved / length,
        "longest_n_run": longest_n_run(sequence),
    }


def count_snp_clusters(positions: list[int], window: int = 101, cutoff: int = 6) -> int:
    if not positions:
        return 0
    left = 0
    in_cluster = False
    clusters = 0
    for right, position in enumerate(positions):
        while positions[left] < position - window + 1:
            left += 1
        if right - left + 1 > cutoff:
            if not in_cluster:
                clusters += 1
                in_cluster = True
        else:
            in_cluster = False
    return clusters


def locate_vp1_columns(
    consensus: str,
    pair: str,
    reference_fasta: Path | None,
    coords_table: Path | None,
    max_candidates: int = 3,
) -> tuple[int, int] | None:
    """Locate the VP1 interval on this group's alignment columns.

    Strategy: pick up to `max_candidates` genbank accessions annotated for the
    same VP1 genotype, align each onto the group consensus with
    `mafft --add --keeplength` (consensus length == alignment length), and
    transfer the annotated 1-based VP1 [start, end] through the column
    correspondence. The best-identity candidate wins.
    """
    if reference_fasta is None or coords_table is None:
        return None
    if not reference_fasta.exists() or not coords_table.exists():
        return None
    vp1_genotype = pair.split("_")[1] if "_" in pair else ""
    group = "GI" if pair.startswith("GI.") else "GII"
    annotated: list[tuple[str, int, int]] = []
    with coords_table.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row.get("group") == group
                and row.get("protein") == "VP1"
                and row.get("genotype") == vp1_genotype
            ):
                try:
                    annotated.append(
                        (row["accession"], int(row["start"]), int(row["end"]))
                    )
                except (KeyError, ValueError):
                    continue
    if not annotated:
        return None
    reference_seqs = {}
    for record in SeqIO.parse(reference_fasta, "fasta"):
        genotype, _, accession = record.id.partition("_")
        if genotype == vp1_genotype and accession in {item[0] for item in annotated}:
            reference_seqs[accession] = str(record.seq).upper()
    candidates = [item for item in annotated if item[0] in reference_seqs][:max_candidates]
    if not candidates:
        return None

    consensus_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
    add_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
    try:
        consensus_file.write(f">consensus\n{consensus}\n")
        consensus_file.flush()
        best: tuple[float, tuple[int, int]] | None = None
        for accession, vp1_start, vp1_end in candidates:
            add_file.seek(0)
            add_file.truncate()
            add_file.write(f">{accession}\n{reference_seqs[accession]}\n")
            add_file.flush()
            proc = subprocess.run(
                [
                    "mafft",
                    "--quiet",
                    "--add",
                    add_file.name,
                    "--keeplength",
                    consensus_file.name,
                ],
                capture_output=True,
                text=True,
            )
            if proc.returncode != 0:
                continue
            added = list(SeqIO.parse(StringIO(proc.stdout), "fasta"))
            if len(added) != 2:
                continue
            ref_aligned = str(added[1].seq).upper() if added[1].id == accession else str(added[0].seq).upper()
            cons_aligned = str(added[0].seq).upper() if added[1].id == accession else str(added[1].seq).upper()
            if len(ref_aligned) != len(consensus):
                continue
            columns = [i for i, b in enumerate(ref_aligned) if b in BASES]
            identity_columns = [
                i
                for i in range(len(consensus))
                if ref_aligned[i] in BASES and cons_aligned[i] in BASES
            ]
            if not identity_columns:
                continue
            identity = sum(
                ref_aligned[i] == cons_aligned[i] for i in identity_columns
            ) / len(identity_columns)
            if vp1_end > len(columns):
                continue
            start_col = columns[vp1_start - 1]
            end_col = columns[vp1_end - 1]
            if best is None or identity > best[0]:
                best = (identity, (start_col, end_col))
        return best[1] if best else None
    finally:
        consensus_file.close()
        add_file.close()
        Path(consensus_file.name).unlink(missing_ok=True)
        Path(add_file.name).unlink(missing_ok=True)


def vp1_segment_veto(
    sequence_a: str,
    sequence_b: str,
    vp1_columns: tuple[int, int] | None,
    window: int,
    step: int,
    min_snv: int,
    block_windows: int,
    min_window_coverage: float = 0.5,
) -> tuple[bool, dict]:
    """Veto a merge when the two references carry a divergent block inside VP1.

    A window is evaluable only when at least `min_window_coverage` of its columns
    are comparable (both sequences A/C/G/T) — partial-genome sequences make some
    windows un-evaluable, which never count as divergent and interrupt runs. A
    window is divergent when it is evaluable and contains >= min_snv SNVs;
    >= block_windows consecutive divergent windows veto the merge.
    """
    if vp1_columns is None:
        return False, {
            "windows": 0,
            "evaluable_windows": 0,
            "divergent_windows": 0,
            "max_run": 0,
            "comparable_fraction": 0.0,
        }
    start_col, end_col = vp1_columns
    min_comparable = max(1, int(window * min_window_coverage))
    flags: list[bool] = []
    evaluable = 0
    comparable_total = 0
    for offset in range(start_col, end_col - window + 1, step):
        seg_a = sequence_a[offset : offset + window]
        seg_b = sequence_b[offset : offset + window]
        comparable = 0
        snvs = 0
        for a, b in zip(seg_a, seg_b):
            if a in BASES and b in BASES:
                comparable += 1
                if a != b:
                    snvs += 1
        comparable_total += comparable
        window_evaluable = comparable >= min_comparable
        evaluable += int(window_evaluable)
        flags.append(window_evaluable and snvs >= min_snv)
    max_run = run = 0
    for flag in flags:
        run = run + 1 if flag else 0
        max_run = max(max_run, run)
    span = max(1, len(flags) * window)
    return max_run >= block_windows, {
        "windows": len(flags),
        "evaluable_windows": evaluable,
        "divergent_windows": sum(flags),
        "max_run": max_run,
        "comparable_fraction": round(comparable_total / span, 4),
    }


def raw_qc(records: list[SeqRecord], aligned_by_id: dict[str, str]) -> list[dict]:
    if not records:
        return []
    if len(records) < 4:
        label = (
            "SINGLETON_NO_REFERENCE"
            if len(records) == 1
            else "SMALL_GROUP_NO_STABLE_REFERENCE"
        )
        return [
            {
                "sequence_id": record.id,
                "group_size": len(records),
                "private_mutations": "",
                "snp_clusters": "",
                "flag": label,
            }
            for record in records
        ]
    root = majority_sequence(set(aligned_by_id), aligned_by_id)
    rows = []
    for record in records:
        query = aligned_by_id[record.id]
        private = [
            position
            for position, (ref, base) in enumerate(zip(root, query), start=1)
            if ref in BASES and base in BASES and ref != base
        ]
        clusters = count_snp_clusters(private) if len(records) > 1 else 0
        rows.append(
            {
                "sequence_id": record.id,
                "group_size": len(records),
                "private_mutations": len(private),
                "snp_clusters": clusters,
                "flag": "SNP_CLUSTERS>1" if clusters > 1 else "",
            }
        )
    return rows


def choose_medoid(
    member_ids: set[str],
    aligned_by_id: dict[str, str],
    target: str,
    qc_excluded: set[str],
) -> tuple[str | None, float]:
    candidates = [name for name in sorted(member_ids) if name not in qc_excluded]
    if not candidates:
        return None, math.inf
    eligible = [
        name
        for name in candidates
        if sum(base in BASES for base in aligned_by_id[name]) / len(target) >= 0.98
    ] or candidates
    selected = min(
        eligible,
        key=lambda name: (acgt_distance(aligned_by_id[name], target), name),
    )
    return selected, acgt_distance(aligned_by_id[selected], target)


def parse_cigar(cigar: str) -> int:
    return sum(
        int(length)
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if op in {"M", "I", "=", "X"}
    )


def parse_sam(path: Path) -> dict[str, list[dict]]:
    hits = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 11 or fields[2] == "*":
                continue
            tags = {
                item.split(":", 2)[0]: item.split(":", 2)[2]
                for item in fields[11:]
                if item.count(":") >= 2
            }
            hits[fields[0]].append(
                {
                    "ref": fields[2],
                    "mapq": int(fields[4]),
                    "aligned": parse_cigar(fields[5]),
                    "as": int(tags.get("AS", "-9999")),
                    "nm": int(tags.get("NM", "9999")),
                }
            )
    return hits


def reference_similarity(left: dict, right: dict) -> float:
    distance = acgt_distance(left["aligned_sequence"], right["aligned_sequence"])
    return 1.0 - distance if math.isfinite(distance) else 0.0


def select_root_medoid(
    members: set[str],
    node,
    tree,
    aligned_by_id: dict[str, str],
    qc_excluded: set[str],
) -> tuple[str, bool] | None:
    """Real-sequence representative: member closest (patristic) to the clade root.

    QC-clean members are preferred; a QC-excluded member is only used when
    nothing else exists (tagged so outputs can flag it).
    """
    ordered = sorted(members)
    clean = [m for m in ordered if m not in qc_excluded]
    pool = clean if clean else ordered
    qc_fallback = not clean
    if tree is not None and node is not None:
        tips = {t.name: t for t in node.get_terminals()}

        def patristic(member: str) -> float:
            tip = tips.get(member)
            if tip is None:
                return float("inf")
            try:
                return tree.distance(node, tip)
            except Exception:
                return float("inf")

        best = min(pool, key=lambda m: (patristic(m), m))
        return best, qc_fallback
    # No tree: classic medoid by mean pairwise ACGT distance.
    def mean_distance(member: str) -> float:
        ds = [
            acgt_distance(aligned_by_id[member], aligned_by_id[other])
            for other in pool
            if other != member
        ]
        return sum(ds) / len(ds) if ds else 0.0

    best = min(pool, key=lambda m: (mean_distance(m), m))
    return best, qc_fallback


def build_reference(
    uid: str,
    members: set[str],
    node,
    aligned_by_id: dict[str, str],
    state_by_node: dict,
    qc_excluded: set[str],
    args: argparse.Namespace,
    tree=None,
) -> dict | None:
    if len(members) == 1:
        raw_id = next(iter(members))
        if raw_id in qc_excluded:
            return None
        return {
            "uid": uid,
            "members": set(members),
            "node": node,
            "node_name": support(node)[0] if node is not None else raw_id,
            "aligned_sequence": aligned_by_id[raw_id],
            "source_type": "observed_raw",
            "source_record": raw_id,
            "low_pp_fraction": 0.0,
            "n_fraction": aligned_by_id[raw_id].count("N") / len(aligned_by_id[raw_id]),
            "longest_n_run": longest_n_run(aligned_by_id[raw_id]),
        }

    node_name = support(node)[0] if node is not None else "Merged"
    if args.representative == "medoid":
        picked = select_root_medoid(members, node, tree, aligned_by_id, qc_excluded)
        if picked is None:
            return None
        medoid, qc_fallback = picked
        sequence = aligned_by_id[medoid]
        return {
            "uid": uid,
            "members": set(members),
            "node": node,
            "node_name": node_name,
            "aligned_sequence": sequence,
            "source_type": "medoid_root_qc_fallback" if qc_fallback else "medoid_root",
            "source_record": medoid,
            "low_pp_fraction": 0.0,
            "n_fraction": sequence.count("N") / len(sequence),
            "longest_n_run": longest_n_run(sequence),
        }
    uncertainty = {"low_pp_fraction": 1.0, "n_fraction": 0.0, "longest_n_run": 0}
    if node_name in state_by_node:
        sequence, uncertainty = hybrid_asr(
            members,
            aligned_by_id,
            state_by_node[node_name]["state"],
            state_by_node[node_name]["posterior"],
            args,
        )
        source_type = "hybrid_asr"
    else:
        sequence = majority_sequence(members, aligned_by_id)
        uncertainty = {
            "low_pp_fraction": 0.0,
            "n_fraction": sequence.count("N") / len(sequence),
            "longest_n_run": longest_n_run(sequence),
        }
        source_type = "majority_consensus"

    requires_medoid = (
        uncertainty["low_pp_fraction"] > args.max_low_pp_fraction
        or uncertainty["n_fraction"] > args.max_n_fraction
        or uncertainty["longest_n_run"] > args.max_n_run
    )
    source_record = ""
    if requires_medoid:
        medoid, _ = choose_medoid(members, aligned_by_id, sequence, qc_excluded)
        if medoid is not None:
            sequence = aligned_by_id[medoid]
            source_type = "medoid"
            source_record = medoid
        else:
            sequence = majority_sequence(members, aligned_by_id)
            group_majority = majority_sequence(set(aligned_by_id), aligned_by_id)
            sequence = "".join(
                group_base if base == "N" and group_base in BASES else base
                for base, group_base in zip(sequence, group_majority)
            )
            uncertainty["n_fraction"] = sequence.count("N") / len(sequence)
            uncertainty["longest_n_run"] = longest_n_run(sequence)
            source_type = "hierarchical_majority_no_qc_medoid"

    return {
        "uid": uid,
        "members": set(members),
        "node": node,
        "node_name": node_name,
        "aligned_sequence": sequence,
        "source_type": source_type,
        "source_record": source_record,
        **uncertainty,
    }


def make_windows(
    aligned_by_id: dict[str, str],
    active_members: set[str],
    window: int,
    step: int,
    output: Path,
) -> tuple[dict[str, tuple[str, int]], dict[str, list[str]]]:
    metadata = {}
    reads_by_raw = defaultdict(list)
    counter = 0
    with output.open("w") as handle:
        for raw_id in sorted(active_members):
            sequence = aligned_by_id[raw_id]
            for start in range(0, len(sequence) - window + 1, step):
                fragment = sequence[start : start + window]
                if any(base not in BASES for base in fragment):
                    continue
                counter += 1
                read_id = f"read_{counter:08d}"
                metadata[read_id] = (raw_id, start)
                reads_by_raw[raw_id].append(read_id)
                handle.write(f">{read_id}\n{fragment}\n")
    return metadata, reads_by_raw


def validate_units(
    units: list[dict],
    assignment: dict[str, str],
    aligned_by_id: dict[str, str],
    reads_fasta: Path,
    read_metadata: dict[str, tuple[str, int]],
    bowtie2: Path,
    threads: int,
    args: argparse.Namespace,
    iteration_dir: Path,
) -> tuple[list[dict], Counter, list[dict]]:
    ref_fasta = iteration_dir / "references.fasta"
    records = [
        SeqRecord(
            Seq(unit["aligned_sequence"].replace("-", "")),
            id=unit["uid"],
            description="",
        )
        for unit in units
    ]
    SeqIO.write(records, ref_fasta, "fasta")
    build = bowtie2.with_name("bowtie2-build")
    index = iteration_dir / "index"
    subprocess.run(
        [str(build), str(ref_fasta), str(index)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    sam = iteration_dir / "windows.sam"
    subprocess.run(
        [
            str(bowtie2),
            "-x",
            str(index),
            "-f",
            "-U",
            str(reads_fasta),
            "-k",
            "2",
            "--very-sensitive-local",
            "--no-unal",
            "--threads",
            str(threads),
            "-S",
            str(sam),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    hits = parse_sam(sam)
    unit_by_uid = {unit["uid"]: unit for unit in units}
    counts = defaultdict(Counter)
    confusions = Counter()
    details = []
    for read_id, (raw_id, start) in read_metadata.items():
        target = assignment[raw_id]
        ranked = sorted(
            hits.get(read_id, []),
            key=lambda hit: (-hit["as"], hit["nm"], hit["ref"]),
        )
        best = ranked[0] if ranked else None
        second = next(
            (hit for hit in ranked[1:] if hit["ref"] != best["ref"]),
            None,
        ) if best else None
        score_margin = best["as"] - second["as"] if best and second else 999
        nm_margin = second["nm"] - best["nm"] if best and second else 999
        own_best = bool(best and best["ref"] == target)
        mapped = bool(best and best["aligned"] >= args.min_aligned_window)
        target_window = unit_by_uid[target]["aligned_sequence"][
            start : start + args.window
        ]
        comparable = [
            (raw_base, ref_base)
            for raw_base, ref_base in zip(
                aligned_by_id[raw_id][start : start + args.window],
                target_window,
            )
            if raw_base in BASES and ref_base in BASES
        ]
        target_mapped = (
            len(comparable) >= args.min_aligned_window
            and sum(left != right for left, right in comparable) / len(comparable) <= 0.10
        )
        distinguishable = bool(
            mapped
            and own_best
            and (
                score_margin >= args.min_score_margin
                or nm_margin >= args.min_nm_margin
            )
        )
        bucket = counts[target]
        bucket["total"] += 1
        bucket["mapped"] += int(mapped)
        bucket["target_mapped"] += int(target_mapped)
        bucket["own_best"] += int(own_best)
        bucket["distinguishable"] += int(distinguishable)
        if best and best["ref"] != target:
            confusions[(target, best["ref"])] += 1
        details.append(
            {
                "read": read_id,
                "raw_sequence": raw_id,
                "target": target,
                "best": best["ref"] if best else "",
                "second": second["ref"] if second else "",
                "score_margin": score_margin if best else "",
                "nm_margin": nm_margin if best else "",
                "target_mapped": target_mapped,
                "distinguishable": distinguishable,
            }
        )

    rows = []
    for unit in units:
        bucket = counts[unit["uid"]]
        total = max(1, bucket["total"])
        target_rate = bucket["target_mapped"] / total
        distinguishable_rate = bucket["distinguishable"] / total
        own_best_rate = bucket["own_best"] / total
        if target_rate < args.min_target_mapping_rate:
            decision = "split"
        elif distinguishable_rate < args.min_distinguishable_rate:
            if own_best_rate < args.umbrella_own_best_floor and (
                splittable_children(unit, args.min_tips) >= 2
                if args.umbrella_gate == "structural"
                else len(unit["members"]) >= args.umbrella_min_members
            ):
                # Large unit whose member reads mostly best-hit other
                # references: an umbrella consensus representing nobody —
                # split it instead of merging it into a neighbour. Small
                # units with low own-best are usually cloud-diluted
                # duplicates and still route to merge.
                decision = "split"
            else:
                decision = "merge"
        else:
            decision = "retain"
        rows.append(
            {
                "reference": unit["uid"],
                "members": len(unit["members"]),
                "windows": bucket["total"],
                "mapping_rate": f"{bucket['mapped'] / total:.6f}",
                "target_mapping_rate": f"{target_rate:.6f}",
                "own_best_rate": f"{bucket['own_best'] / total:.6f}",
                "distinguishable_rate": f"{distinguishable_rate:.6f}",
                "decision": decision,
            }
        )
    return rows, confusions, details


def splittable_children(unit: dict, min_child_members: int) -> int:
    """Child clades holding at least `min_child_members` of the unit's own members."""
    node = unit.get("node")
    if node is None or node.is_terminal():
        return 0
    count = 0
    for child in node.clades:
        own = {tip.name for tip in child.get_terminals()} & unit["members"]
        if len(own) >= min_child_members:
            count += 1
    return count


def main() -> None:
    args = arguments()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not args.bowtie2.exists() or not args.bowtie2.with_name("bowtie2-build").exists():
        raise FileNotFoundError("bowtie2 and bowtie2-build must both be available")

    alignment = AlignIO.read(args.alignment, "fasta")
    records = list(alignment)
    if not records:
        raise ValueError("Alignment is empty")
    pair = genotype_pair(records[0].id)
    aligned_by_id = {record.id: str(record.seq).upper() for record in records}
    record_by_id = {record.id: record for record in records}
    if len(set(aligned_by_id.values())) == 0:
        raise ValueError("No aligned sequences found")

    qc_rows = raw_qc(records, aligned_by_id)
    qc_excluded = {
        row["sequence_id"] for row in qc_rows if row["flag"] == "SNP_CLUSTERS>1"
    }
    write_tsv(
        args.outdir / "raw_reference_qc.tsv",
        qc_rows,
        ["sequence_id", "group_size", "private_mutations", "snp_clusters", "flag"],
    )

    state_by_node = load_state(args.state)
    tree = None
    parents = {}
    node_rows = []
    boundary_nodes = set()
    metrics_cache = {}
    mldist = load_mldist(args.mldist)
    if mldist:
        distance_names, distance_matrix = mldist
        distance_index = {name: index for index, name in enumerate(distance_names)}
    else:
        distance_matrix = None
        distance_index = {}

    def metrics(node) -> dict[str, float]:
        members = tuple(terminal.name for terminal in node.get_terminals())
        if members not in metrics_cache:
            if distance_matrix is not None:
                metrics_cache[members] = distance_stats(
                    list(members),
                    distance_index,
                    distance_matrix,
                    args.max_distance_pairs,
                )
            else:
                metrics_cache[members] = fallback_distance_stats(
                    list(members),
                    aligned_by_id,
                    args.max_distance_pairs,
                )
        return metrics_cache[members]

    if len(records) >= 4:
        if args.tree is None or not args.tree.exists():
            raise FileNotFoundError("IQ-TREE tree is required for groups with >=4 sequences")
        tree = Phylo.read(args.tree, "newick")
        tree.root_at_midpoint()
        parents = {child: parent for parent in tree.find_clades() for child in parent.clades}
        for node in tree.get_nonterminals():
            name, sh_alrt, ufboot = support(node)
            current = metrics(node)
            parent = parents.get(node)
            parent_stats = metrics(parent) if parent is not None else {
                "d95": math.nan,
                "max": math.nan,
                "median": math.nan,
                "pairs": 0,
            }
            support_pass = sh_alrt >= args.min_sh_alrt and ufboot >= args.min_ufboot
            within_pass = (
                len(node.get_terminals()) >= args.min_tips
                and current["d95"] <= args.max_within_d95
                and current["max"] <= args.max_within_max
            )
            parent_pass = (
                parent is not None
                and parent_stats["d95"] > args.min_parent_d95
                and parent_stats["d95"] - current["d95"] >= args.min_parent_delta
            )
            candidate = support_pass and within_pass and parent_pass
            if candidate:
                boundary_nodes.add(node)
            node_rows.append(
                {
                    "node": name,
                    "tips": len(node.get_terminals()),
                    "sh_alrt": f"{sh_alrt:.1f}",
                    "ufboot": f"{ufboot:.1f}",
                    "within_median": f"{current['median']:.6f}",
                    "within_d95": f"{current['d95']:.6f}",
                    "within_max": f"{current['max']:.6f}",
                    "parent_d95": f"{parent_stats['d95']:.6f}",
                    "candidate": candidate,
                }
            )
    write_tsv(
        args.outdir / "node_candidate_metrics.tsv",
        node_rows,
        [
            "node",
            "tips",
            "sh_alrt",
            "ufboot",
            "within_median",
            "within_d95",
            "within_max",
            "parent_d95",
            "candidate",
        ],
    )

    selected_nodes = []
    selected_tips = []
    if tree is None:
        selected_tips = list(records)
    else:
        root_stats = metrics(tree.root)
        homogeneous_root = (
            len(records) >= args.min_tips
            and root_stats["d95"] <= args.max_within_d95
            and root_stats["max"] <= args.max_within_max
        )
        if homogeneous_root:
            selected_nodes = [tree.root]
        else:
            def select(node) -> None:
                if node in boundary_nodes:
                    selected_nodes.append(node)
                elif node.is_terminal():
                    selected_tips.append(record_by_id[node.name])
                else:
                    for child in node.clades:
                        select(child)
            select(tree.root)

    uid_counter = 0

    def next_uid() -> str:
        nonlocal uid_counter
        uid_counter += 1
        return f"U{uid_counter:04d}"

    units = []
    for node in selected_nodes:
        members = (
            {terminal.name for terminal in node.get_terminals()}
            if node is not None
            else set(aligned_by_id)
        )
        unit = build_reference(
            next_uid(),
            members,
            node,
            aligned_by_id,
            state_by_node,
            qc_excluded,
            args,
                tree=tree,
        )
        if unit:
            units.append(unit)
    for record in selected_tips:
        unit = build_reference(
            next_uid(),
            {record.id},
            next(
                (terminal for terminal in tree.get_terminals() if terminal.name == record.id),
                None,
            ) if tree else None,
            aligned_by_id,
            state_by_node,
            qc_excluded if len(records) > 1 else set(),
            args,
        )
        if unit:
            units.append(unit)
    if not units:
        fallback_node = tree.root if tree is not None else None
        fallback = build_reference(
            next_uid(),
            set(aligned_by_id),
            fallback_node,
            aligned_by_id,
            state_by_node,
            qc_excluded,
            args,
                tree=tree,
        )
        if fallback is None:
            raise RuntimeError(
                "Unable to construct a synthetic fallback after raw-reference QC"
            )
        fallback["source_type"] = f"qc_fallback_{fallback['source_type']}"
        units = [fallback]

    assignment = {
        member: unit["uid"] for unit in units for member in unit["members"]
    }
    active_members = set(assignment)
    reads_fasta = args.outdir / "sliding_windows.fasta"
    read_metadata, _ = make_windows(
        aligned_by_id,
        active_members,
        args.window,
        args.step,
        reads_fasta,
    )
    if not read_metadata:
        raise RuntimeError("No valid A/C/G/T sliding windows were generated")

    vp1_columns = None
    if not args.no_vp1_veto:
        group_consensus = majority_sequence(set(aligned_by_id), aligned_by_id)
        vp1_columns = locate_vp1_columns(
            group_consensus,
            pair,
            args.vp1_reference_fasta,
            args.vp1_coords,
        )
    vp1_veto_events: list[dict] = []
    vp1_veto_seen: set[tuple[str, str]] = set()

    iteration_rows = []
    final_validation = []
    all_details = []
    seen_states: set[tuple] = set()
    best_state: tuple | None = None
    for iteration in range(1, args.max_iterations + 1):
        iteration_dir = args.outdir / f"iteration_{iteration:02d}"
        iteration_dir.mkdir(exist_ok=True)
        validation, confusions, details = validate_units(
            units,
            assignment,
            aligned_by_id,
            reads_fasta,
            read_metadata,
            args.bowtie2,
            args.threads,
            args,
            iteration_dir,
        )
        all_details = details
        decision_by_uid = {row["reference"]: row["decision"] for row in validation}
        final_validation = validation
        flagged = sum(row["decision"] in ("merge", "split") for row in validation)
        score = (-flagged, -len(units))
        if best_state is None or score > best_state[0]:
            best_state = (score, list(units), dict(assignment), validation, details)
        iteration_rows.append(
            {
                "iteration": iteration,
                "references": len(units),
                "retain": sum(row["decision"] == "retain" for row in validation),
                "merge": sum(row["decision"] == "merge" for row in validation),
                "split": sum(row["decision"] == "split" for row in validation),
                "action": "",
            }
        )
        if all(row["decision"] == "retain" for row in validation):
            iteration_rows[-1]["action"] = "converged"
            break
        if iteration == args.max_iterations:
            iteration_rows[-1]["action"] = "max_iterations_reached"
            break
        # Split/merge actions can settle into a limit cycle (an umbrella split,
        # its children re-merged, split again...). Detect revisited states and
        # keep the best-scoring snapshot for the final output.
        signature = frozenset(
            frozenset(unit["members"]) for unit in units
        )
        if signature in seen_states:
            iteration_rows[-1]["action"] = "oscillation_detected"
            break
        seen_states.add(signature)

        unit_by_uid = {unit["uid"]: unit for unit in units}
        split_candidates = [
            unit for unit in units
            if decision_by_uid[unit["uid"]] == "split"
            and unit["node"] is not None
            and not unit["node"].is_terminal()
        ]
        if split_candidates:
            split_done = 0
            removed_by_split: set[str] = set()
            new_units: list[dict] = []
            for target in sorted(
                split_candidates, key=lambda unit: len(unit["members"]), reverse=True
            ):
                replacements = []
                for child in target["node"].clades:
                    # Intersect with the unit's OWN members: a merge-born unit's
                    # node (MRCA of the merged pair) can span members owned by
                    # other units, and pulling those in would zero out their
                    # windows.
                    members = (
                        {terminal.name for terminal in child.get_terminals()}
                        & target["members"]
                    )
                    if not members:
                        continue
                    replacement = build_reference(
                        next_uid(),
                        members,
                        child,
                        aligned_by_id,
                        state_by_node,
                        qc_excluded,
                        args,
                            tree=tree,
                    )
                    if replacement:
                        replacements.append(replacement)
                if len(replacements) >= 2:
                    removed_by_split.add(target["uid"])
                    new_units.extend(replacements)
                    for replacement in replacements:
                        for member in replacement["members"]:
                            assignment[member] = replacement["uid"]
                    split_done += 1
            if split_done:
                units = [
                    unit for unit in units if unit["uid"] not in removed_by_split
                ] + new_units
                iteration_rows[-1]["action"] = f"batch split {split_done} unit(s)"
                continue

        failing = {
            row["reference"] for row in validation if row["decision"] == "merge"
        }
        pair_candidates = []
        for left_index, left in enumerate(units):
            for right in units[left_index + 1 :]:
                if left["uid"] not in failing and right["uid"] not in failing:
                    continue
                similarity = reference_similarity(left, right)
                if similarity < args.merge_min_similarity:
                    continue
                ambiguity = (
                    confusions[(left["uid"], right["uid"])]
                    + confusions[(right["uid"], left["uid"])]
                )
                pair_candidates.append((ambiguity, similarity, left, right))
        if not pair_candidates:
            iteration_rows[-1]["action"] = "unresolved_no_high_similarity_merge"
            break
        accepted: list[tuple[float, dict, dict]] = []
        consumed: set[str] = set()
        for ambiguity, similarity, cand_left, cand_right in sorted(
            pair_candidates,
            key=lambda item: (item[0], item[1]),
            reverse=True,
        ):
            if cand_left["uid"] in consumed or cand_right["uid"] in consumed:
                continue
            vetoed, veto_stats = vp1_segment_veto(
                cand_left["aligned_sequence"],
                cand_right["aligned_sequence"],
                vp1_columns,
                args.window,
                args.step,
                args.veto_window_snv,
                args.veto_block_windows,
                args.veto_min_window_coverage,
            )
            if vetoed:
                pair_key = tuple(sorted((cand_left["uid"], cand_right["uid"])))
                if pair_key not in vp1_veto_seen:
                    vp1_veto_seen.add(pair_key)
                    vp1_veto_events.append(
                        {
                            "iteration": iteration,
                            "left": cand_left["uid"],
                            "right": cand_right["uid"],
                            "similarity": f"{similarity:.4f}",
                            "decision": "vetoed",
                            **{k: str(v) for k, v in veto_stats.items()},
                        }
                    )
                continue
            if veto_stats.get("evaluable_windows", 0) < 2 * args.veto_block_windows:
                vp1_veto_events.append(
                    {
                        "iteration": iteration,
                        "left": cand_left["uid"],
                        "right": cand_right["uid"],
                        "similarity": f"{similarity:.4f}",
                        "decision": "allowed_low_vp1_coverage",
                        **{k: str(v) for k, v in veto_stats.items()},
                    }
                )
            accepted.append((similarity, cand_left, cand_right))
            consumed.update({cand_left["uid"], cand_right["uid"]})
        if not accepted:
            iteration_rows[-1]["action"] = (
                "unresolved_vp1_veto" if vp1_veto_events else "unresolved_no_mergeable_pair"
            )
            break
        # Batch-merge every accepted (disjoint) pair, then revalidate everything.
        # Pair-only merges replace the old MRCA cascade, which could absorb many
        # unvetted units in one action.
        remaining = [unit for unit in units if unit["uid"] not in consumed]
        merged_pairs = 0
        qc_failures = 0
        for similarity, left, right in accepted:
            merge_node = None
            if tree is not None:
                merge_node = tree.common_ancestor(
                    list(left["members"] | right["members"])
                )
            merged = build_reference(
                next_uid(),
                left["members"] | right["members"],
                merge_node,
                aligned_by_id,
                state_by_node,
                qc_excluded,
                args,
                    tree=tree,
            )
            if merged is None:
                qc_failures += 1
                remaining.extend([left, right])
                continue
            remaining.append(merged)
            for member in merged["members"]:
                assignment[member] = merged["uid"]
            merged_pairs += 1
        units = remaining
        iteration_rows[-1]["action"] = (
            f"batch merge {merged_pairs} pair(s), {qc_failures} QC failure(s) "
            f"-> {len(units)} references"
        )
        if merged_pairs == 0:
            iteration_rows[-1]["action"] = "unresolved_merge_failed_qc"
            break

    if best_state is not None:
        _, best_units, best_assignment, best_validation, best_details = best_state
        if {unit["uid"] for unit in best_units} != {unit["uid"] for unit in units}:
            units = best_units
            assignment = best_assignment
            final_validation = best_validation
            all_details = best_details

    # Final shadow sweep: the split branch preempts merges every round, so
    # near-duplicate references (<= shadow-sweep-distance apart) can survive the
    # loop. Shadow references kill read distinguishability — collapse them here
    # unless the VP1 veto protects the pair.
    shadow_swept = 0
    if args.shadow_sweep_distance > 0:
        merged_again = True
        while merged_again:
            merged_again = False
            unit_by = {unit["uid"]: unit for unit in units}
            uids = sorted(unit_by)
            for i, uid_a in enumerate(uids):
                if uid_a not in unit_by:
                    continue
                for uid_b in uids[i + 1 :]:
                    if uid_b not in unit_by or uid_a not in unit_by:
                        continue
                    left, right = unit_by[uid_a], unit_by[uid_b]
                    if reference_similarity(left, right) < 1 - args.shadow_sweep_distance:
                        continue
                    vetoed, _ = vp1_segment_veto(
                        left["aligned_sequence"],
                        right["aligned_sequence"],
                        vp1_columns,
                        args.window,
                        args.step,
                        args.veto_window_snv,
                        args.veto_block_windows,
                        args.veto_min_window_coverage,
                    )
                    if vetoed:
                        continue
                    merge_node = (
                        tree.common_ancestor(list(left["members"] | right["members"]))
                        if tree is not None
                        else None
                    )
                    merged = build_reference(
                        next_uid(),
                        left["members"] | right["members"],
                        merge_node,
                        aligned_by_id,
                        state_by_node,
                        qc_excluded,
                        args,
                            tree=tree,
                    )
                    if merged is None:
                        continue
                    units = [
                        u
                        for u in units
                        if u["uid"] not in (uid_a, uid_b)
                    ] + [merged]
                    for member in merged["members"]:
                        assignment[member] = merged["uid"]
                    unit_by = {unit["uid"]: unit for unit in units}
                    shadow_swept += 1
                    merged_again = True
        if shadow_swept:
            sweep_dir = args.outdir / "final_sweep"
            sweep_dir.mkdir(exist_ok=True)
            final_validation, _, all_details = validate_units(
                units,
                assignment,
                aligned_by_id,
                reads_fasta,
                read_metadata,
                args.bowtie2,
                args.threads,
                args,
                sweep_dir,
            )

    write_tsv(
        args.outdir / "iteration_summary.tsv",
        iteration_rows,
        ["iteration", "references", "retain", "merge", "split", "action"],
    )
    write_tsv(
        args.outdir / "final_window_validation.tsv",
        final_validation,
        [
            "reference",
            "members",
            "windows",
            "mapping_rate",
            "target_mapping_rate",
            "own_best_rate",
            "distinguishable_rate",
            "decision",
        ],
    )
    write_tsv(
        args.outdir / "final_window_details.tsv",
        all_details,
        [
            "read",
            "raw_sequence",
            "target",
            "best",
            "second",
            "score_margin",
            "nm_margin",
            "target_mapped",
            "distinguishable",
        ],
    )
    membership_rows = [
        {"raw_sequence": member, "reference": unit["uid"]}
        for unit in units
        for member in sorted(unit["members"])
    ]
    write_tsv(
        args.outdir / "final_membership.tsv",
        membership_rows,
        ["raw_sequence", "reference"],
    )

    output_records = []
    source_rows = []
    for number, unit in enumerate(sorted(units, key=lambda item: item["uid"]), start=1):
        if unit["source_type"] == "observed_raw":
            output_id = f"{unit['source_record'] or next(iter(unit['members']))}_MAPREF_RAW"
        elif unit["source_type"] in ("medoid", "medoid_root"):
            output_id = f"{unit['source_record']}_MAPREF_MEDOID"
        elif unit["source_type"] == "medoid_root_qc_fallback":
            output_id = f"{unit['source_record']}_MAPREF_MEDOIDQC"
        else:
            output_id = f"{pair}_MAPREF_ref{number:03d}_{unit['node_name']}_{unit['source_type'].upper()}"
        sequence = unit["aligned_sequence"].replace("-", "")
        output_records.append(SeqRecord(Seq(sequence), id=output_id, description=""))
        source_rows.append(
            {
                "output_id": output_id,
                "internal_reference": unit["uid"],
                "genotype_pair": pair,
                "members": len(unit["members"]),
                "source_node": unit["node_name"],
                "source_type": unit["source_type"],
                "source_record": unit["source_record"],
                "low_pp_fraction": f"{unit['low_pp_fraction']:.6f}",
                "n_fraction": f"{unit['n_fraction']:.6f}",
                "longest_n_run": unit["longest_n_run"],
            }
        )
    SeqIO.write(output_records, args.outdir / "final_mapping_references.fasta", "fasta")
    write_tsv(
        args.outdir / "final_reference_sources.tsv",
        source_rows,
        [
            "output_id",
            "internal_reference",
            "genotype_pair",
            "members",
            "source_node",
            "source_type",
            "source_record",
            "low_pp_fraction",
            "n_fraction",
            "longest_n_run",
        ],
    )
    if vp1_veto_events:
        write_tsv(
            args.outdir / "vp1_veto_events.tsv",
            vp1_veto_events,
            [
                "iteration",
                "left",
                "right",
                "similarity",
                "decision",
                "windows",
                "evaluable_windows",
                "divergent_windows",
                "max_run",
                "comparable_fraction",
            ],
        )

    with (args.outdir / "summary.txt").open("w") as handle:
        handle.write(f"Genotype pair: {pair}\n")
        handle.write(f"Input sequences: {len(records)}\n")
        handle.write(
            "Raw sequences flagged by current SNP-cluster QC: "
            f"{len(qc_excluded)}\n"
        )
        handle.write(f"Initial boundary clades: {len(selected_nodes)}\n")
        handle.write(f"Valid {args.window} nt windows: {len(read_metadata)}\n")
        if args.no_vp1_veto:
            handle.write("VP1 segment veto: disabled\n")
        elif vp1_columns is None:
            handle.write("VP1 segment veto: unavailable (no annotated VP1 reference)\n")
        else:
            handle.write(
                "VP1 segment veto: columns "
                f"{vp1_columns[0] + 1}-{vp1_columns[1] + 1}, "
                f"{len(vp1_veto_events)} pair(s) vetoed\n"
            )
        handle.write(f"Final mapping references: {len(units)}\n")
        handle.write(f"Iterations: {len(iteration_rows)}\n")
        handle.write(
            "Final unresolved references: "
            f"{sum(row['decision'] != 'retain' for row in final_validation)}\n"
        )


if __name__ == "__main__":
    main()
