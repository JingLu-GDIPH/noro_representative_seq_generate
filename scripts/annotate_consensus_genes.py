#!/usr/bin/env python3
"""Annotate gene coordinates on representative consensus sequences.

Reference gene coordinates (``reference/genome_gene_coordinates.tsv``) are defined
on original database accessions, none of which are the pipeline's consensus
sequences. This script transfers those coordinates onto each consensus by:

  1. For each consensus, finding its VP1 genotype and the reference records of
     the same genotype (``reference/{gi,gii}_with_genotype.fasta``).
  2. Aligning the consensus with all same-genotype references in one MAFFT run.
  3. Picking the reference with the smallest p-distance to the consensus as the
     "nearest reference".
  4. Building a per-column position map between the nearest reference and the
     consensus (skipping gap columns) and remapping every gene interval of that
     reference onto consensus coordinates.

Output is a TSV with one row per (consensus, gene), suitable for downstream
analyses that need to know where ORF1 / p48 / NTPase / p22 / VPg / 3CLpro /
RdRp / VP1 / VP2 sit on each representative genome.

Note on naming: norovirus ORF2 encodes the major capsid protein (VP1, ~1620 bp)
and ORF3 encodes the minor capsid (VP2, ~660 bp); the reference TSV uses these
exact labels.
"""

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


VALID_BASES = frozenset("ACGT")

#: Capture ``<RdRp>_<VP1>`` from a consensus ID (e.g. ``GII.P4_GII.4``).
GENOTYPE_PAIR_PATTERN = re.compile(
    r"(G(?:I|II|IX)\.P[A-Za-z0-9]+)_(G(?:I|II|IX)\.[A-Za-z0-9]+)"
)

#: Norovirus canonical gene order (ORF1 sub-units first, then capsids).
GENE_ORDER = [
    "5'UTR", "ORF1", "p48", "NTPase", "p22", "VPg", "3CLpro", "RdRp",
    "VP1", "VP2", "3'UTR",
]

#: Core annotations required for a complete norovirus coding-coordinate set.
REQUIRED_CORE_GENES = frozenset(
    {"ORF1", "p48", "NTPase", "p22", "VPg", "3CLpro", "RdRp", "VP1", "VP2"}
)

#: Genotypes with more references than this use MAFFT's faster ``--retree 2``.
LARGE_REFERENCE_THRESHOLD = 200


def extract_vp1_genotype(sequence_id):
    """Return the VP1 genotype (2nd field, e.g. ``GII.4``) from a consensus ID."""
    match = GENOTYPE_PAIR_PATTERN.search(str(sequence_id))
    return match.group(2) if match else None


def extract_rdrp_vp1_pair(sequence_id):
    """Return ``<RdRp>_<VP1>`` (e.g. ``GII.P4_GII.4``) from a consensus ID."""
    match = GENOTYPE_PAIR_PATTERN.search(str(sequence_id))
    return match.group(0) if match else None


def reference_accession(record_id):
    """Reference FASTA IDs are ``<genotype>_<accession>`` -> return accession."""
    parts = str(record_id).split("_")
    return parts[1] if len(parts) >= 2 else str(record_id)


def reference_genotype(record_id):
    """Reference FASTA IDs are ``<genotype>_<accession>`` -> return genotype."""
    parts = str(record_id).split("_")
    return parts[0] if parts else None


def alignment_p_distance(seq1, seq2):
    """Pairwise p-distance over columns where both bases are A/C/G/T."""
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


def run_mafft(records, threads, large=False):
    """Align records with MAFFT and return aligned SeqRecord list."""
    with tempfile.TemporaryDirectory() as tmp:
        input_path = Path(tmp) / "input.fasta"
        SeqIO.write(records, input_path, "fasta")
        cmd = ["mafft", "--quiet", "--auto", "--thread", str(threads)]
        if large:
            cmd += ["--retree", "2"]
        cmd.append(str(input_path))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr[:500]}")
        aligned = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
    return aligned


def build_position_map(ref_aligned, cons_aligned):
    """Map reference positions -> consensus positions over an alignment pair.

    Returns (ref_to_cons, cons_to_ref) dicts keyed by 1-based coordinates of
    *non-gap* residues in each respective sequence. ``N`` is not valid for
    similarity calculations, but it still occupies a genomic coordinate and
    therefore increments position counters here. Reference positions follow
    the reference's own coordinate system (1-based), consensus positions follow
    the consensus's own coordinate system (1-based, ungapped).
    """
    ref_to_cons = {}
    cons_to_ref = {}
    ref_pos = 0
    cons_pos = 0
    for ref_base, cons_base in zip(ref_aligned, cons_aligned):
        ref_is_residue = ref_base != "-"
        cons_is_residue = cons_base != "-"
        if ref_is_residue:
            ref_pos += 1
        if cons_is_residue:
            cons_pos += 1
        if ref_is_residue and cons_is_residue:
            ref_to_cons[ref_pos] = cons_pos
            cons_to_ref[cons_pos] = ref_pos
    return ref_to_cons, cons_to_ref


def transfer_interval(ref_start, ref_end, ref_to_cons):
    """Transfer a 1-based inclusive reference interval to consensus coordinates.

    Returns ``(cons_start, cons_end, coverage_pct)`` where coverage_pct is the
    fraction of reference interval positions that successfully mapped to a
    consensus position. Returns ``(None, None, 0.0)`` if nothing mapped.
    """
    mapped = [ref_to_cons[p] for p in range(ref_start, ref_end + 1) if p in ref_to_cons]
    if not mapped:
        return None, None, 0.0
    coverage = len(mapped) / (ref_end - ref_start + 1)
    return min(mapped), max(mapped), coverage


def gene_sort_key(protein):
    """Sort key putting genes in canonical norovirus genome order."""
    try:
        return (GENE_ORDER.index(protein), protein)
    except ValueError:
        return (len(GENE_ORDER), protein)


def load_reference_coords(coords_path):
    """Map (group, accession) -> list of gene-interval dicts."""
    by_acc = defaultdict(list)
    with open(coords_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_acc[(row["group"], row["accession"])].append({
                "protein": row["protein"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "source": row["source"],
            })
    return by_acc


def load_reference_fasta(fasta_path, group):
    """Map (group, accession) -> SeqRecord, and genotype -> [SeqRecord, ...]."""
    by_acc = {}
    by_genotype = defaultdict(list)
    for record in SeqIO.parse(fasta_path, "fasta"):
        acc = reference_accession(record.id)
        geno = reference_genotype(record.id)
        record.description = ""
        by_acc[(group, acc)] = record
        if geno:
            by_genotype[geno].append(record)
    return by_acc, by_genotype


def infer_vp2_interval(reference, intervals):
    """Infer a missing VP2/ORF3 interval from the reference nucleotide sequence.

    Norovirus ORF3 starts close to the VP1 end and encodes a roughly 200-260 aa
    protein. Search a narrow interval around the annotated VP1 end for an ATG
    followed by an in-frame stop, retaining biologically plausible 500-1000 nt
    ORFs. Coordinates are returned as 1-based inclusive positions.
    """
    if any(interval["protein"] == "VP2" for interval in intervals):
        return None
    vp1_intervals = [
        interval for interval in intervals if interval["protein"] == "VP1"
    ]
    if not vp1_intervals:
        return None

    sequence = str(reference.seq).upper().replace("-", "")
    vp1_end = max(interval["end"] for interval in vp1_intervals)
    search_start = max(0, vp1_end - 31)
    search_end = min(len(sequence) - 2, vp1_end + 120)
    candidates = []
    for start_index in range(search_start, search_end):
        if sequence[start_index : start_index + 3] != "ATG":
            continue
        for stop_index in range(start_index + 3, len(sequence) - 2, 3):
            if sequence[stop_index : stop_index + 3] not in {"TAA", "TAG", "TGA"}:
                continue
            length = stop_index + 3 - start_index
            if 500 <= length <= 1000:
                start = start_index + 1
                end = stop_index + 3
                candidates.append((abs(start - vp1_end), -length, start, end))
            break
    if not candidates:
        return None
    _, _, start, end = min(candidates)
    return {
        "protein": "VP2",
        "start": start,
        "end": end,
        "strand": "+",
        "source": "inferred_from_orf3",
    }


def effective_reference_intervals(reference, coords, group):
    """Return reference intervals, supplementing a missing VP2 when possible."""
    accession = reference_accession(reference.id)
    intervals = [
        dict(interval) for interval in coords.get((group, accession), [])
    ]
    inferred_vp2 = infer_vp2_interval(reference, intervals)
    if inferred_vp2 is not None:
        intervals.append(inferred_vp2)
    return intervals


def annotate_group(consensus_records, references, group, vp1, coords, threads,
                   fallback_label=None):
    """Align consensus + references, transfer gene coords from nearest reference.

    Returns (output_rows, per-consensus status list). ``fallback_label`` is set
    when this call is a cross-genotype fallback so the status records it.
    """
    rows = []
    statuses = []
    status_tag = "annotated_fallback" if fallback_label else "annotated"

    complete_references = [
        reference
        for reference in references
        if REQUIRED_CORE_GENES
        <= {
            interval["protein"]
            for interval in effective_reference_intervals(reference, coords, group)
        }
    ]
    if complete_references:
        references = complete_references
    else:
        annotated_references = [
            reference
            for reference in references
            if effective_reference_intervals(reference, coords, group)
        ]
        if annotated_references:
            references = annotated_references

    large = (len(consensus_records) + len(references)) > LARGE_REFERENCE_THRESHOLD
    combined = [SeqRecord(r.seq, id=r.id, description="") for r in references] + \
               [SeqRecord(r.seq, id=r.id, description="") for r in consensus_records]
    # Ensure unique IDs (references and consensus share the `_` namespace).
    seen = set()
    safe_records = []
    for r in combined:
        rid = r.id
        counter = 1
        while rid in seen:
            rid = f"{r.id}_dup{counter}"
            counter += 1
        r.id = rid
        seen.add(rid)
        safe_records.append(r)

    try:
        aligned = run_mafft(safe_records, threads, large=large)
    except Exception as exc:
        print(f"[WARN] MAFFT failed for {vp1}: {exc}; skipping")
        for rec in consensus_records:
            statuses.append((rec.id, vp1, None, None, "mafft_failed"))
        return rows, statuses

    aligned_by_id = {r.id: str(r.seq) for r in aligned}

    for rec in consensus_records:
        cons_aln = aligned_by_id.get(rec.id)
        cons_len = len(str(rec.seq).replace("-", ""))
        if cons_aln is None:
            print(f"[WARN] consensus {rec.id} missing from {vp1} alignment")
            statuses.append((rec.id, vp1, None, None, "missing_from_alignment"))
            continue

        # Find nearest reference by p-distance.
        best_acc = None
        best_pdist = None
        best_ref_aln = None
        best_ref_record = None
        for ref_rec in references:
            ref_aln = aligned_by_id.get(ref_rec.id)
            if ref_aln is None:
                continue
            pdist = alignment_p_distance(ref_aln, cons_aln)
            if pdist is None:
                continue
            if best_pdist is None or pdist < best_pdist:
                best_pdist = pdist
                best_acc = reference_accession(ref_rec.id)
                best_ref_aln = ref_aln
                best_ref_record = ref_rec

        if best_acc is None:
            print(f"[WARN] no measurable reference for {rec.id}")
            statuses.append((rec.id, vp1, None, None, "no_measurable_ref"))
            continue

        ref_to_cons, _cons_to_ref = build_position_map(best_ref_aln, cons_aln)
        intervals = effective_reference_intervals(best_ref_record, coords, group)
        intervals = sorted(intervals, key=lambda x: gene_sort_key(x["protein"]))

        if not intervals:
            statuses.append((rec.id, vp1, best_acc, best_pdist, "ref_unannotated"))
            continue

        for iv in intervals:
            cstart, cend, cov = transfer_interval(iv["start"], iv["end"], ref_to_cons)
            rows.append({
                "consensus_id": rec.id,
                "group": group,
                "rdrp_vp1_pair": extract_rdrp_vp1_pair(rec.id) or "",
                "vp1_genotype": vp1,
                "nearest_reference_accession": best_acc,
                "nearest_reference_pdistance": f"{best_pdist:.6f}",
                "consensus_length": cons_len,
                "protein": iv["protein"],
                "ref_accession": best_acc,
                "ref_start": iv["start"],
                "ref_end": iv["end"],
                "ref_strand": iv["strand"],
                "ref_source": iv["source"],
                "consensus_start": cstart if cstart is not None else "",
                "consensus_end": cend if cend is not None else "",
                "gene_coverage_pct": f"{cov * 100:.2f}",
            })
        statuses.append((rec.id, vp1, best_acc, best_pdist, status_tag))
        print(f"  {rec.id}: nearest={best_acc} pdist={best_pdist:.4f} "
              f"genes={len(intervals)} len={cons_len}"
              f"{' [cross-genotype fallback]' if fallback_label else ''}")
    return rows, statuses


def annotate(consensus_path, reference_fasta, coords_path, group, threads, outdir,
             cross_genotype_fallback=False, output_name=None):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    coords = load_reference_coords(coords_path)
    ref_by_acc, ref_by_genotype = load_reference_fasta(reference_fasta, group)
    all_references = [rec for recs in ref_by_genotype.values() for rec in recs]

    # Group consensus by VP1 genotype so each genotype is aligned once.
    cons_by_genotype = defaultdict(list)
    for record in SeqIO.parse(consensus_path, "fasta"):
        vp1 = extract_vp1_genotype(record.id)
        if vp1 is None:
            print(f"[WARN] cannot parse VP1 genotype from {record.id}; skipping")
            continue
        cons_by_genotype[vp1].append(record)

    output_rows = []
    genotype_summary = []

    for vp1 in sorted(cons_by_genotype):
        consensus_records = cons_by_genotype[vp1]
        references = ref_by_genotype.get(vp1, [])
        if references:
            rows, statuses = annotate_group(
                consensus_records, references, group, vp1, coords, threads)
            output_rows.extend(rows)
            genotype_summary.extend(statuses)
        elif cross_genotype_fallback and all_references:
            # No same-genotype reference: fall back to nearest across ALL
            # references. Norovirus gene layout (ORF1 -> VP1 -> VP2) is conserved
            # across genotypes, so cross-genotype transfer is biologically valid,
            # though p-distance will be higher. We align each orphan consensus
            # against all references in one run (using --retree 2 for speed).
            print(f"[INFO] genotype {vp1} has no reference; using cross-genotype fallback")
            rows, statuses = annotate_group(
                consensus_records, all_references, group, vp1, coords, threads,
                fallback_label=True)
            output_rows.extend(rows)
            genotype_summary.extend(statuses)
        else:
            print(f"[WARN] no reference records for genotype {vp1}; skipping "
                  f"{len(consensus_records)} consensus "
                  f"(enable --cross-genotype-fallback to annotate via nearest other genotype)")
            for rec in consensus_records:
                genotype_summary.append((rec.id, vp1, None, None, "no_reference"))

    fieldnames = [
        "consensus_id", "group", "rdrp_vp1_pair", "vp1_genotype",
        "nearest_reference_accession", "nearest_reference_pdistance",
        "consensus_length", "protein", "ref_accession", "ref_start", "ref_end",
        "ref_strand", "ref_source", "consensus_start", "consensus_end",
        "gene_coverage_pct",
    ]
    out_path = outdir / (
        output_name or f"{group.lower()}_consensus_gene_coordinates.tsv"
    )
    with open(out_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    annotated = sum(1 for s in genotype_summary if s[4] == "annotated")
    print(f"\n{group}: annotated {annotated}/{len(genotype_summary)} consensus "
          f"-> {out_path}")
    if any(s[4] != "annotated" for s in genotype_summary):
        print("  Non-annotated:")
        for cid, vp1, acc, pd, status in genotype_summary:
            if status != "annotated":
                print(f"    {cid}: {status}")
    return out_path, genotype_summary


def main():
    parser = argparse.ArgumentParser(
        description="Transfer reference gene coordinates onto consensus "
                    "representative sequences via per-genotype MAFFT alignment "
                    "to the nearest reference.")
    parser.add_argument("--consensus", required=True, help="Filtered consensus FASTA")
    parser.add_argument("--reference-fasta", required=True,
                        help="Reference raw FASTA (gi/gii_with_genotype.fasta)")
    parser.add_argument("--coords", required=True,
                        help="reference/genome_gene_coordinates.tsv")
    parser.add_argument("--group", required=True, choices=["GI", "GII"],
                        help="Which group (GI/GII)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--cross-genotype-fallback", action="store_true",
                        help="When a consensus genotype has no same-genotype "
                             "reference, annotate it against the nearest reference "
                             "across ALL genotypes (status 'annotated_fallback').")
    parser.add_argument(
        "--output-name",
        help="Output TSV filename (default: <group>_consensus_gene_coordinates.tsv)",
    )
    args = parser.parse_args()

    annotate(args.consensus, args.reference_fasta, args.coords,
             args.group, args.threads, args.outdir,
             cross_genotype_fallback=args.cross_genotype_fallback,
             output_name=args.output_name)


if __name__ == "__main__":
    main()
