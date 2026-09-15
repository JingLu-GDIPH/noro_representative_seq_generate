# Figure: Empirical calibration of the distinguishable-rate threshold for norovirus mapping-reference panels

**Files**: `Fig_threshold_calibration.pdf` (vector) / `Fig_threshold_calibration.png` (300 dpi)
**Plotted data**: `figure_data/panel_box_data.tsv`, `figure_data/panel_scatter_data.tsv` (tab-separated)
**Generating script**: `../make_calibration_figure.py`; underlying simulations in
`../v2_test/gii4_threshold_calibration/` and `../v2_test/gii17_threshold_calibration_B2000/`

## Experiment

150-bp sliding-window (step 25 nt) competitive Bowtie2 mapping (`-k 2 --very-sensitive-local`)
simulating the production distinguishability validation. A window is *distinguishable*
when its best hit is the reference of its own unit with alignment-score margin ≥6 or
edit-distance margin ≥2 (aligned length ≥135). GII.4: 20 named variants (44 genomes,
2–3 accessions each; consensus reference per variant, 20-reference panel).
GII.17: 5 nextstrain clades (316 genomes; clade B restricted to sequences collected ≥2000;
consensus reference per clade, 5-reference panel). Within-label rates come from
accession/strain-level panels (43 and 316 references, respectively).

## Panels

- **A** Distinguishable-rate distributions. Between-label panels (must stay separate,
  "keep side") vs within-label panels (near-duplicates, "merge side"). Candidate
  thresholds: 0.30 (production), 0.50 (recommended upper bound), 0.80 (selection-rules
  document). The two distributions overlap because variant/clade labels are coarser
  than genomic distance — the threshold separates *distances*, not labels.
- **B** Keep-side rate versus distance to the nearest panel unit (consensus p-distance,
  log scale, 2–16%). All 25 units are shown, including the ancestral GII.17 clades
  A and B (≈16%, rate ≈1.0). Empirical points sit below the Poisson expectation
  P(≥2 SNVs per 150 bp) (non-uniform variation, local-alignment trimming,
  member-vs-consensus mismatch). The hardest pairs: GII.4 2004↔2006a
  (2.66% → rate 0.739) and GII.17 Kawasaki_308↔Kawasaki_323 (1.95% → 0.880/0.772).
- **C** Threshold decision band. Green: T ≤ 0.636 keeps every GII.4 variant and every
  GII.17 clade; yellow: 0.636 < T ≤ 0.739 retains GII.4 only (GII.17 clade B would be
  merge-flagged; note B also fails target-mapping ≥0.95 and would be *split*, not merged);
  red: T > 0.739 fails both. Vertical markers: 0.30 (production, both PASS),
  0.50 (recommended upper bound, both PASS), 0.60 (marginal, 1.06× margin),
  0.80 (rules document, both FAIL — merges 2004/2006a and Kawasaki_323).

## Recommended statement

The distinguishable-rate threshold was calibrated on 20 named GII.4 variants and 5
GII.17 clades. The empirical keep-side floors were 0.739 (GII.4 2006a, closest neighbour
2004 at 2.66% genome distance) and 0.636 (GII.17 clade B). T = 0.30–0.50 passes both
calibrations with ≥1.3× margin; T = 0.80 merges named variants and was rejected.

## 中文摘要

三面板校准图：A) 组间（须保留）与组内（可合并）的 distinguishable rate 分布及候选阈值线；
B) 组间 rate 随“到最近参考单元距离”上升的散点与泊松理论曲线（2006a@2.66%→0.739、
Kawasaki_323@1.95%→0.772 为最难分离对）；C) 阈值判定带：绿色区 T≤0.636 两组均通过，
0.80（规则文档值）两组均失败。推荐 T=0.30–0.50。
