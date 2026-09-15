#!/bin/bash
# Re-run SELECT_MAPPING_REFERENCES at T=0.5 reusing existing alignments+trees,
# then COLLECT. Faithful to mapping_reference.nf args except --min-distinguishable-rate.
set -u
PROJ=/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate
G=$PROJ/result_4ref/gii
PY=/Users/LuJ/mambaforge/envs/noro-consensus/bin/python3
BT2=/Users/LuJ/Yunpan/software/Bioinf_software/bowtie2-2.3.3/bowtie2
BUILD=$PROJ/scripts/build_mapping_reference_group.py
T=0.5
J=2   # concurrent groups (mirrors maxForks=2)

run_group() {
  local aln="$1"
  local key
  key=$(basename "$aln" .fasta); key=${key#aligned_}
  local out="$G/06_group_references/mapping_ref_${key}"
  [ -s "$out/final_mapping_references.fasta" ] && { echo "SKIP $key"; return 0; }
  local targs=""
  local tdir="$G/05_iqtree/iqtree_${key}"
  if [ -s "$tdir/iqtree.treefile" ]; then
    targs="--tree $tdir/iqtree.treefile --state $tdir/iqtree.state --mldist $tdir/iqtree.mldist"
  fi
  "$PY" "$BUILD" --alignment "$aln" $targs --outdir "$out" \
    --bowtie2 "$BT2" --threads 4 \
    --min-sh-alrt 80.0 --min-ufboot 95.0 \
    --max-within-d95 0.05 --max-within-max 1.0 \
    --min-parent-d95 0.05 --min-parent-delta 0.01 \
    --asr-pp 0.90 --majority-frequency 0.60 \
    --max-low-pp-fraction 0.05 --max-n-fraction 0.005 --max-n-run 15 \
    --window 150 --step 25 \
    --min-distinguishable-rate $T \
    --min-target-mapping-rate 0.95 \
    --merge-min-similarity 0.95 --max-iterations 20 \
    > "$out.log" 2>&1
  local rc=$?
  echo "DONE rc=$rc $key refs=$(grep -c '^>' "$out/final_mapping_references.fasta" 2>/dev/null || echo 0)"
  return $rc
}
export -f run_group
export G PY BT2 BUILD T

cd "$PROJ"
ls "$G"/04_alignments/aligned_*.fasta | xargs -I{} -P $J bash -c 'run_group "$@"' _ {}
echo "=== ALL GROUPS DONE ==="

# COLLECT (mirrors COLLECT_MAPPING_REFERENCES)
"$PY" "$PROJ/scripts/collect_mapping_references.py" \
  --input-glob 'mapping_ref_*/final_mapping_references.fasta' \
  --output "$G/07_final/all_mapping_references.fasta" \
  --report "$G/07_final/final_sequence_qc.tsv" \
  --max-n-gap-percent 10.0
{
  echo "Norovirus wastewater mapping-reference branch (T=0.5 rerun, trees reused)"
  echo "Input basis: gii_prev_T030_20260725 04_alignments + 05_iqtree (3697-seq run)"
  echo "Window size: 150"
  echo "Minimum distinguishable-window rate: 0.5"
  echo "Final references: $(grep -c '^>' "$G/07_final/all_mapping_references.fasta")"
} > "$G/07_final/summary.txt"
cat "$G/07_final/summary.txt"
