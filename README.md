# Norovirus Representative Reference Panel (v2.0 — strain-level medoid recipe)

> **v2.0 定案方法**：全部由**真实存在序列**（medoid，到 clade 根节点 patristic 最近成员）组成、
> 互相可被 150 bp 读段区分（支持 ≥99% 基因组相似毒株共存）的参考 panel。
> 完整步骤、参数与依据见 **[STRAIN_PANEL_METHOD.md](STRAIN_PANEL_METHOD.md)**。
>
> ```bash
> # 每组核心命令（上游比对/建树由 mapping_reference.nf 编排）
> python3 scripts/build_mapping_reference_group.py \
>   --alignment aligned.fasta --tree ... --state ... --mldist ... \
>   --outdir mapping_ref_<group> --bowtie2 <bowtie2> --threads 4 \
>   --representative medoid --min-sh-alrt 80.0 --min-ufboot 0.0 \
>   --max-within-d95 0.005 --max-within-max 0.02 \
>   --min-parent-d95 0.005 --min-parent-delta 0.001 \
>   --min-distinguishable-rate 0.5 --min-target-mapping-rate 0.95 \
>   --merge-min-similarity 0.95 --veto-block-windows 15 \
>   --umbrella-own-best-floor 0.5 --umbrella-gate structural \
>   --max-iterations 60 --shadow-sweep-distance 0.005
> ```
> 实测：GII.P16_GII.4 504→37 条（最小间距 0.64%）、GII.17 664→19+15 条
> （Kawasaki_308/323/Romania 各获专属参考，Node276 坍缩解决）。

---

# Norovirus Representative Consensus Pipeline

This repository contains a Nextflow workflow for generating representative
norovirus GI/GII consensus sequences from downloaded genome FASTA files.

The current workflow is intended for VP1-typed NCBI nucleotide records whose
FASTA IDs follow:

```text
genotype_accession_collectiondate
```

Example:

```text
GII.4_ON123456_2022-03-18
```

## Current Version

Version: `1.4.0`

Main changes in this version:

- A new `ORIENT_SEQUENCES` step puts every cluster on a single strand before
  MAFFT. Records downloaded on the reverse strand (e.g. `GII.P31_GII.4_KX158285`,
  `GII.P7_GII.6_KX158282`, `GII.P7_GII.6_KX158286`) are reverse-complemented so
  they no longer align with thousands of spurious gaps and produce misleading
  pairwise similarities.
- A representative QC step (`scripts/analyze_raw_representative_qc.py`) flags
  raw-sequence representatives that carry excessive private mutations / SNP
  clusters vs their genotype clade root (nextclade-style), to detect probable
  sequencing errors. QC-clean sets are written as
  `*_all_final_consensus_filtered.fa` (raw reps tagged `SNP_CLUSTERS>1` removed).
- VP1 genotype hard grouping before vclust clustering.
- Clustering defaults: ANI `0.90`, qcov `0.65`, rcov `0.65`, len_ratio `0.65`.
- Complete-linkage clustering to reduce chaining artifacts.
- MAFFT is run for every multi-sequence cluster after removing existing gaps.
- After MAFFT, all alignment columns with <50% A/C/G/T coverage are removed;
  `N` and gap are not counted as valid coverage.
- Final representative sequences are written without alignment gaps, so
  low-coverage insertion columns cannot be filled into consensus sequences.
- Final representatives are globally pruned so all pairwise similarities are
  below or equal to the configured threshold, default `95%`.
- Runtime defaults use `--threads 8`; vclust and final pruning use the full
  thread budget, while MAFFT/IQ-TREE tasks use up to 4 threads per task with
  Nextflow-controlled concurrency.

## Directory Layout

```text
main.nf                         Nextflow workflow
nextflow.config                 local execution and default parameter config
environment.yml                 conda/mamba environment definition
scripts/                        reusable helper scripts
tests/                          unit tests for critical helper logic
ncbi_download_2026-06-11/       downloaded FASTA and metadata tables
results_new/gi/                 final GI analysis result
results_new/gii/                final GII analysis result
```

Old test runs and Nextflow work caches have been removed from the project
workspace. Formal result directories are `results_new/gi` and `results_new/gii`.

## Run

Use the `noro-consensus` environment and Nextflow:

```bash
conda activate noro-consensus
nextflow run main.nf \
  --input_file ncbi_download_2026-06-11/ncbi_norovirus_gii_gt5800_collected_2000_onward.fasta \
  --output_dir results_new/gii \
  --threads 8
```

GI example:

```bash
nextflow run main.nf \
  --input_file ncbi_download_2026-06-11/ncbi_norovirus_gi_gt5800_collected_2000_onward.fasta \
  --output_dir results_new/gi \
  --threads 8
```

Important parameters:

```text
--ani_threshold 0.9
--query_coverage_threshold 0.65
--reference_coverage_threshold 0.65
--length_ratio_threshold 0.65
--cluster_algorithm complete
--internal_similarity_threshold 95.0
--alignment_end_min_coverage 0.5
--threads 8
```

On this workstation, `--threads 8` is the recommended default because it uses
most of the 10 logical cores while leaving two cores free for the desktop and
I/O. Use `--threads 10` only when dedicating the machine to the run.

## Main Outputs

For each output directory:

```text
01_merged/filter_report.txt
02_cluster/cluster_info.csv
02_cluster/vcluster_results/ani.same_genotype.tsv
03_split/clusters/
03b_oriented/
04_aligned/
05_consensus/
06_enhanced_consensus/
07_final_results/all_final_consensus.fasta
07_final_results/all_final_consensus_filtered.fa
07_final_results/summary_report.txt
08_validation/final_validation_report.txt
```

Downstream workflows should usually consume:

```text
results_new/gi/07_final_results/all_final_consensus_filtered.fa
results_new/gii/07_final_results/all_final_consensus_filtered.fa
results_new/gi/02_cluster/cluster_info.csv
results_new/gii/02_cluster/cluster_info.csv
```

`all_final_consensus.fasta` is the raw pipeline output; `_filtered.fa` additionally
removes raw representatives flagged `SNP_CLUSTERS>1` by the representative QC step
(see "Representative QC & Sequencing-Error Filtering" below). Use `_filtered.fa`
when you want sequencing-error-suspect raw records excluded.

## Validated Results

Latest formal results:

```text
GI/GII formal outputs should be regenerated after changing the MAFFT column
coverage mask; see `07_final_results/all_final_consensus.fasta` and
`08_validation/final_validation_report.txt` in each output directory.
```

The refreshed GII run used 3698 quality-filtered input sequences, formed 87
VP1 genotype-pure clusters, and produced 24 singleton clusters.

## Tests

Run tests with the same Python environment used by the pipeline:

```bash
/Users/LuJ/mambaforge/envs/noro-consensus/bin/python -m unittest discover -s tests -v
```

The current test suite covers genotype-pair grouping, sequence renaming,
pre-alignment strand orientation normalization, alignment column masking,
final gap removal, tree-node consensus logic, and validation helpers.

## Consensus VP1/RdRp Genotyping

Final GI/GII consensus sequences can be checked with the IPHnano dual-region
genotyping workflow:

```bash
python scripts/run_consensus_genotyping.py --threads 8 --group both
```

The helper uses the reference files in:

```text
/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_genotyping/ref_seq
```

Outputs are written under:

```text
results_new/genotyping_consensus/
```

## Consensus Representativeness Analysis

After VP1/RdRp genotyping, add genotype labels to final consensus FASTA headers
and evaluate consensus coverage within each VP1 genotype:

```bash
python scripts/analyze_consensus_representativeness.py --threads 8 --skip-existing
python scripts/plot_consensus_representativeness_trees.py
```

Main outputs:

```text
results_new/representativeness_analysis/01_renamed_consensus/
results_new/representativeness_analysis/02_by_vp1/
results_new/representativeness_analysis/03_metrics/representativeness_summary.tsv
results_new/representativeness_analysis/03_metrics/nearest_consensus_distances.tsv
results_new/representativeness_analysis/04_figures/
results_new/representativeness_analysis/representativeness_report.md
```

The tree comparison uses raw-only trees versus raw taxa induced by the
raw+consensus trees. A direct RF comparison between consensus-only and
all-sequence trees is not used because those trees do not share the same tip set.

## Representative QC & Sequencing-Error Filtering

Some final representatives are raw input records kept as their own representative
because they differed too much from sibling sequences to merge into a consensus.
Such a record may be a genuinely divergent lineage, or it may carry
assembly/sequencing errors. To tell these apart, `scripts/analyze_raw_representative_qc.py`
adapts the nextstrain/ncov (Nextclade) QC logic to norovirus:

- For each RdRp_VP1 genotype, the **clade root** is the column-wise majority
  consensus of all sibling raw sequences in that genotype's aligned cluster.
- For every raw representative it counts, relative to that clade root:
  - `private_mutations` — A/C/G/T positions where the representative differs;
  - `snp_clusters` — number of 101 nt windows containing *strictly more than* 6
    private substitutions (Nextclade's SNP-cluster rule).
- Representatives are tagged `SNP_CLUSTERS>1` (ncov's hard exclusion criterion)
  and/or `HIGH_PRIVATE(>=N)` (genotype 90th-percentile soft flag). Singleton
  clusters (only one sequence, no reference to compare against) are tagged
  `SINGLETON_NO_REFERENCE` — this is **not** a QC pass, it means undecidable.

```bash
python scripts/analyze_raw_representative_qc.py \
  --representatives results_new/gii/07_final_results/gii_all_final_consensus.fasta \
  --clusters-dir    results_new/gii/04_aligned \
  --output          results_new/qc_analysis/gii_raw_rep_qc.tsv
```

A full write-up of the method, results, and the rep-vs-sibling baseline check
lives in `results_new/qc_analysis/REPORT.md`.

To produce QC-clean representative sets, remove every sequence tagged
`SNP_CLUSTERS>1` from the final consensus:

```bash
# GI:  42 -> 29 (removed 13)   GII: 117 -> 92 (removed 25)
results_new/gi/07_final_results/gi_all_final_consensus_filtered.fa
results_new/gii/07_final_results/gii_all_final_consensus_filtered.fa
```

`*_filtered.fa` keeps all node/consensus representatives and only drops the raw
representatives flagged `SNP_CLUSTERS>1`. Singleton-cluster representatives are
retained (they are undecidable by this method, not confirmed errors).

## Notes for Downstream Calling

- The workflow assumes genotype is present in the FASTA ID before the first
  underscore.
- Untyped or cross-genotype vclust pairs are excluded before clustering.
- A singleton cluster does not always mean a sequence is globally unique; it
  means it had no same-genotype partner satisfying the configured ANI, qcov,
  rcov, and len_ratio constraints under complete linkage.
- The final FASTA may contain more sequences than the number of clusters,
  because a large cluster can be split into multiple phylogenetic boundary
  nodes before final global pruning.
