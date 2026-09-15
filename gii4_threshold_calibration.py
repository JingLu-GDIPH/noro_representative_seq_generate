#!/usr/bin/env python3
"""Calibrate distinguishable_rate threshold using known GII.4 variant genomes.

Data: 20 named variants x 2-3 accessions each (user-curated list).
Run A (keep side): units = variants (majority-consensus ref per variant, 20-ref
  competitive index) -> per-variant distinguishable_rate, production semantics.
Run B (merge side): units = individual accessions (44-ref index) -> per-accession
  rate; same-variant accessions should be indistinguishable (low rate).
Threshold should sit between the two distributions.
"""
import csv, re, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path
from io import StringIO

import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BT2 = "/Users/LuJ/Yunpan/software/Bioinf_software/bowtie2-2.3.3/bowtie2"
RAW_LOCAL = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/Rawdata/gii_merged_sequences.fasta"
FETCHED = "/tmp/gii4_missing.fasta"
VARIANT_TSV = "/tmp/gii4_variants.tsv"
WORK = Path("/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/v2_test/gii4_threshold_calibration")
WINDOW, STEP, MIN_ALIGNED = 150, 25, 135
MIN_SCORE_MARGIN, MIN_NM_MARGIN = 6, 2
TARGET_MAX_MISMATCH = 0.10
BASES = frozenset("ACGT")

def load_inputs():
    acc2var, want = {}, []
    with open(VARIANT_TSV) as fh:
        for line in fh:
            v, a = line.rstrip("\n").split("\t")
            acc2var[a] = v
            want.append(a)
    recs = {}
    # fetched (NCBI) — header ">AF145896.1 description"
    for r in SeqIO.parse(FETCHED, "fasta"):
        acc = r.id.split(".")[0]
        if acc in acc2var:
            recs[acc] = str(r.seq).upper()
    # local — header "GII.P4_GII.4_<ACC>_<date>"
    for r in SeqIO.parse(RAW_LOCAL, "fasta"):
        for acc in want:
            if acc in r.id and acc not in recs:
                recs[acc] = str(r.seq).upper().replace("-", "")
    missing = [a for a in want if a not in recs]
    if missing:
        sys.exit(f"missing accessions: {missing}")
    out = []
    for acc in want:
        out.append(SeqRecord(Seq(recs[acc]), id=f"{acc2var[acc]}|{acc}", description=""))
    return out

def mafft_and_trim(records):
    fa = WORK / "all.fasta"
    SeqIO.write(records, fa, "fasta")
    r = subprocess.run(["mafft", "--quiet", "--auto", "--thread", "6", str(fa)],
                       capture_output=True, text=True, check=True)
    aln = list(SeqIO.parse(StringIO(r.stdout), "fasta"))
    arr = np.array([[c for c in str(x.seq).upper()] for x in aln])
    keep = [j for j in range(arr.shape[1])
            if np.isin(arr[:, j], list("ACGT")).mean() >= 0.5]
    trimmed = [SeqRecord(Seq("".join(str(x.seq)[j] for j in keep)), id=x.id, description="")
               for x in aln]
    SeqIO.write(trimmed, WORK / "aligned_trimmed.fasta", "fasta")
    return trimmed, len(keep)

def ungapped(seq):
    return seq.replace("-", "")

def variant_consensus(members, aligned_by_id):
    cols = list(zip(*[aligned_by_id[m] for m in members]))
    cons = []
    for col in cols:
        bases = [b for b in col if b in BASES]
        cons.append(Counter(bases).most_common(1)[0][0] if bases else "N")
    return "".join(cons)

def make_windows(aligned_by_id, members):
    reads, meta = [], {}
    n = 0
    for m in sorted(members):
        seq = aligned_by_id[m]
        for start in range(0, len(seq) - WINDOW + 1, STEP):
            frag = seq[start:start + WINDOW]
            if any(b not in BASES for b in frag):
                continue
            n += 1
            rid = f"read_{n:08d}"
            meta[rid] = (m, start)
            reads.append(f">{rid}\n{frag}\n")
    wf = WORK / "windows.fasta"
    wf.write_text("".join(reads))
    return meta

def run_bowtie2(refs, tag):
    d = WORK / tag
    d.mkdir(exist_ok=True)
    rf = d / "refs.fasta"
    rf.write_text("".join(f">{k}\n{v}\n" for k, v in refs.items()))
    idx = d / "index"
    subprocess.run([str(BT2) + "-build", str(rf), str(idx)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    sam = d / "win.sam"
    subprocess.run([BT2, "-x", str(idx), "-f", "-U", str(WORK / "windows.fasta"),
                    "-k", "2", "--very-sensitive-local", "--no-unal",
                    "--threads", "6", "-S", str(sam)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits = defaultdict(list)
    for line in sam.open():
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        cig = f[5]
        if cig == "*":
            continue
        aln_len = sum(int(l) for l, op in re.findall(r"(\d+)([MIDNSHP=X])", cig)
                      if op in "MI=X")
        tags = {}
        for t in f[11:]:
            parts = t.split(":", 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]
        hits[f[0]].append({"ref": f[2], "as": int(tags.get("AS", -999)),
                           "nm": int(tags.get("NM", 999)), "len": aln_len})
    return hits

def evaluate(units, aligned_by_id, tag):
    """units: dict uid -> {'ref_aligned': str(gapped), 'members': [ids]}"""
    refs = {uid: ungapped(u["ref_aligned"]) for uid, u in units.items()}
    assign = {m: uid for uid, u in units.items() for m in u["members"]}
    meta = make_windows(aligned_by_id, list(assign))
    hits = run_bowtie2(refs, tag)
    stats = defaultdict(Counter)
    for rid, (m, start) in meta.items():
        target = assign[m]
        ranked = sorted(hits.get(rid, []), key=lambda h: (-h["as"], h["nm"], h["ref"]))
        best = ranked[0] if ranked else None
        second = next((h for h in ranked[1:] if h["ref"] != best["ref"]), None) if best else None
        s_margin = best["as"] - second["as"] if best and second else 999
        n_margin = second["nm"] - best["nm"] if best and second else 999
        own = bool(best and best["ref"] == target)
        mapped = bool(best and best["len"] >= MIN_ALIGNED)
        t_slice = units[target]["ref_aligned"][start:start + WINDOW]
        comp = [(a, b) for a, b in zip(aligned_by_id[m][start:start + WINDOW], t_slice)
                if a in BASES and b in BASES]
        t_mapped = bool(len(comp) >= MIN_ALIGNED and
                        sum(a != b for a, b in comp) / len(comp) <= TARGET_MAX_MISMATCH)
        dist = bool(mapped and own and (s_margin >= MIN_SCORE_MARGIN or n_margin >= MIN_NM_MARGIN))
        c = stats[target]
        c["total"] += 1; c["mapped"] += int(mapped); c["target"] += int(t_mapped)
        c["own"] += int(own); c["dist"] += int(dist)
    rows = []
    for uid in sorted(stats):
        c = stats[uid]
        rows.append({"unit": uid, "windows": c["total"],
                     "mapping_rate": c["mapped"]/c["total"],
                     "target_rate": c["target"]/c["total"],
                     "own_best_rate": c["own"]/c["total"],
                     "distinguishable_rate": c["dist"]/c["total"]})
    return rows

def main():
    WORK.mkdir(parents=True, exist_ok=True)
    records = load_inputs()
    print(f"sequences: {len(records)}, variants: {len({r.id.split('|')[0] for r in records})}")
    trimmed, aln_len = mafft_and_trim(records)
    print(f"alignment trimmed length: {aln_len}")
    aligned_by_id = {r.id: str(r.seq).upper() for r in trimmed}

    # ---- Run A: variant-level refs (keep side) ----
    by_var = defaultdict(list)
    for r in trimmed:
        by_var[r.id.split("|")[0]].append(r.id)
    unitsA = {v: {"ref_aligned": variant_consensus(ms, aligned_by_id), "members": ms}
              for v, ms in by_var.items()}
    rowsA = evaluate(unitsA, aligned_by_id, "runA_variant_panel")

    # ---- Run B: accession-level refs (merge side) ----
    unitsB = {r.id: {"ref_aligned": aligned_by_id[r.id], "members": [r.id]}
              for r in trimmed}
    rowsB = evaluate(unitsB, aligned_by_id, "runB_accession_panel")

    # variant consensus pairwise p-distances
    vids = sorted(unitsA)
    D = np.zeros((len(vids), len(vids)))
    for i in range(len(vids)):
        for j in range(i + 1, len(vids)):
            a, b = unitsA[vids[i]]["ref_aligned"], unitsA[vids[j]]["ref_aligned"]
            ok = [(x, y) for x, y in zip(a, b) if x in BASES and y in BASES]
            D[i, j] = D[j, i] = sum(x != y for x, y in ok) / len(ok)

    with open(WORK / "runA_variant_rates.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rowsA[0])); w.writeheader(); w.writerows(rowsA)
    with open(WORK / "runB_accession_rates.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rowsB[0])); w.writeheader(); w.writerows(rowsB)

    print("\n===== Run A: variant-level panel (KEEP side) =====")
    print(f"{'variant':8} {'windows':>7} {'map_rate':>8} {'target':>7} {'ownbest':>7} {'DIST':>7}")
    for r in sorted(rowsA, key=lambda r: r["distinguishable_rate"]):
        print(f"{r['unit']:8} {r['windows']:>7} {r['mapping_rate']:>8.3f} "
              f"{r['target_rate']:>7.3f} {r['own_best_rate']:>7.3f} {r['distinguishable_rate']:>7.3f}")
    keep_min = min(r["distinguishable_rate"] for r in rowsA)
    print(f"KEEP-side minimum distinguishable_rate = {keep_min:.3f}")

    print("\n===== Run B: accession-level panel (MERGE side) =====")
    same_variant = defaultdict(list)
    for r in rowsB:
        same_variant[r["unit"].split("|")[0]].append(r["distinguishable_rate"])
    print(f"{'variant':8} {'rates of its accessions':>34}")
    for v in sorted(same_variant):
        print(f"{v:8} {str([f'{x:.3f}' for x in same_variant[v]]):>34}")
    merge_max = max(max(v) for v in same_variant.values())
    print(f"MERGE-side maximum distinguishable_rate (same-variant) = {merge_max:.3f}")

    print("\n===== closest variant pairs (consensus p-distance) =====")
    pairs = sorted(((D[i, j], vids[i], vids[j]) for i in range(len(vids))
                    for j in range(i + 1, len(vids))))
    for d, a, b in pairs[:8]:
        print(f"  {a:8} vs {b:8}  d = {d:.4f}")

    print("\n===== CONCLUSION =====")
    print(f"threshold window: ({merge_max:.3f}, {keep_min:.3f})")
    for T in (0.30, 0.50, 0.80):
        ok_merge = all(max(v) < T for v in same_variant.values())
        ok_keep = keep_min >= T
        print(f"  T={T:.2f}: 合并同变异株={'PASS' if ok_merge else 'FAIL'}, 保留不同变异株={'PASS' if ok_keep else 'FAIL'}")

if __name__ == "__main__":
    main()
