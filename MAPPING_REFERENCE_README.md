# Norovirus wastewater mapping-reference branch

This branch is independent of the existing 95%-similarity probe-reference
workflow. It builds a smaller reference set intended for competitive mapping of
150 bp wastewater sequencing reads.

## Selection rules

1. Input records are grouped by the exact `RdRp_VP1` pair. Before alignment,
   records with more than 10% combined `N` and gap content are removed and
   opposite-strand records are reverse-complemented.
2. MAFFT alignments retain only columns with at least 50% A/C/G/T coverage.
   `N` and `-` never count as covered.
3. Groups containing at least four records are analysed by IQ-TREE 3 using
   `GTR+F+R4`, 1,000 SH-aLRT replicates and 1,000 UFBoot replicates with BNNI.
4. Candidate clades require SH-aLRT >=80%, UFBoot >=95%, at least three tips,
   within-clade ML-distance d95 <=0.05, and a parent d95 that is both greater
   than 0.05 and at least 0.01 larger than the clade d95. Maximum distance is
   reported for diagnosis but is not a default hard filter because one outlier
   pair can otherwise fragment an otherwise compact clade into many raw tips.
5. IQ-TREE ASR bases with posterior probability >=0.90 are retained. Lower
   posterior sites fall back to an A/C/G/T majority call when site coverage is
   >=50% and majority frequency is >=60%; otherwise they become `N`.
6. If low-posterior sites exceed 5%, `N` exceeds 0.5%, or the longest `N` run
   exceeds 15 nt, the nearest quality-eligible observed medoid is used.
7. A/C/G/T-only 150 nt windows are generated every 25 nt from all represented
   raw sequences and mapped competitively with Bowtie2. A reference is retained
   when target mapping is >=95% and at least 30% of windows are distinguishable
   by alignment-score margin >=6 or edit-distance margin >=2.
8. References below 30% distinguishable windows are merged only when their
   pairwise A/C/G/T similarity is >=95%. The MRCA reference is reconstructed and
   the window test is repeated. Poor target mapping triggers a tree-guided split.

Groups with fewer than four records cannot provide bootstrap-supported
clade-root QC or meaningful branch support; IQ-TREE does not define bootstrap
analysis for three-tip trees. Their raw sequences are retained but
explicitly marked `SINGLETON_NO_REFERENCE` or
`SMALL_GROUP_NO_STABLE_REFERENCE`; this is an undecidable status, not a QC pass.
In larger groups, observed raw or medoid candidates are subject to the existing
raw-representative SNP-cluster filter. The final collection also reapplies the
10% `N`/gap threshold.

If every provisional reference in a group is an observed raw tip rejected by
that QC, the group is not silently lost. The workflow builds one explicitly
labelled `QC_FALLBACK` synthetic majority/ASR reference from all valid group
members and subjects it to the same 150 nt mapping validation.

If an uncertain ASR requires a medoid but every observed medoid candidate fails
raw-reference QC, the workflow emits a labelled synthetic hierarchical majority
consensus instead of retaining a high-`N` ASR or silently restoring a rejected
raw record. Clade plurality calls are used first; sites with no clade A/C/G/T
coverage fall back to the whole genotype-pair majority after the global 50%
alignment-column mask.

## Run

```bash
nextflow -C mapping_reference.config run mapping_reference.nf \
  --input_file Rawdata/gii_merged_sequences.fasta \
  --output_dir result_4ref/gii \
  --threads 8
```

Run GI with the corresponding GI input and `--output_dir result_4ref/gi`.

## Main outputs

- `01_filtered/filtered_input.fasta`: input after N/gap filtering.
- `04_alignments/`: orientation-normalized, 50%-coverage-masked alignments.
- `05_iqtree/`: IQ-TREE tree, support, distance and ASR state files.
- `06_group_references/`: per-genotype candidate metrics, iterations, window
  validation, memberships, raw-reference QC and final references.
- `07_final/all_mapping_references.fasta`: combined mapping-reference FASTA.
- `07_final/final_sequence_qc.tsv`: final N/gap filter audit.
- `pipeline_report.html`, `timeline_report.html`, `trace.txt`: run diagnostics.

The `final_reference_sources.tsv` table records whether every output is hybrid
ASR, majority consensus, observed raw sequence or medoid, together with the
source node and uncertainty metrics.
