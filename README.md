# Norovirus Representative Reference Panel

从 NCBI 公共数据库的诺如病毒 GI/GII 全基因组序列生成代表性参考 panel。
仓库只保留两条生产流程；完整方法与参数依据见 **[STRAIN_PANEL_METHOD.md](STRAIN_PANEL_METHOD.md)**。

## 两条流程

| 流程 | 入口 | 用途 |
|---|---|---|
| **流程一 · 探针共识** | `main.nf` | 95% 相似性边界节点的共识/祖先代表序列，用于 tiling 探针设计 |
| **流程二 · 映射参考（株级 medoid 定案 v2.0）** | `mapping_reference.nf` + `scripts/mapping_reference/build_mapping_reference_group.py` | 150 bp 污水读段比对参考：全部真实序列（medoid），支持 ≥99% 基因组相似毒株共存 |

## 目录结构

```text
main.nf / nextflow.config            流程一（探针共识）
mapping_reference.nf / mapping_reference.config   流程二（株级 medoid 映射参考）
environment.yml                      conda 环境（noro-consensus）
scripts/common/                      两流程共用（过滤/分组/拆组/方向归一/列过滤）
scripts/probe_consensus/             流程一专用（共识生成/迭代验证/汇总）
scripts/mapping_reference/           流程二专用（核心选择+验证/汇总/审计）
reference/                           VP1 注释坐标与基因型参考（VP1 否决依赖）
STRAIN_PANEL_METHOD.md               流程二 v2.0 完整方法与参数依据
CHANGELOG.md / VERSION               版本记录
```

## 流程一 · 探针共识（main.nf）

```bash
conda activate noro-consensus
nextflow run main.nf --input_file <gii.fasta> --output_dir <dir> --threads 8
# 关键参数: similarity 95%、内部相似性 95%、列覆盖 0.5、迭代 10
```

## 流程二 · 株级 medoid 映射参考（v2.0 定案）

上游（比对/建树，每组）：

```bash
python3 scripts/common/filter_n_sequences.py --input_file IN.fasta --output_file F.fasta --threshold 10.0
python3 scripts/common/group_sequences_by_rdrp_vp1.py --input F.fasta --cluster_info ci.csv --sequences_output g.fasta --report r.tsv
python3 scripts/common/split_clusters.py --cluster_info ci.csv --sequences g.fasta --outdir clusters/
python3 scripts/common/normalize_orientation.py --input cluster.fasta --output oriented.fasta
mafft --quiet --auto --thread 4 oriented.fasta > pretrim.fasta
python3 scripts/common/trim_alignment_ends.py --input pretrim.fasta --output aligned.fasta --min-coverage 0.5
iqtree3 -s aligned.fasta -m GTR+F+R4 --alrt 1000 -B 1000 --bnni --ancestral --asr-min 0.8 -T AUTO --prefix iqtree
```

核心（每组）：

```bash
python3 scripts/mapping_reference/build_mapping_reference_group.py \
  --alignment aligned.fasta --tree iqtree.treefile --state iqtree.state --mldist iqtree.mldist \
  --outdir mapping_ref_<group> --bowtie2 <bowtie2> --threads 4 \
  --representative medoid --min-sh-alrt 80.0 --min-ufboot 0.0 \
  --max-within-d95 0.005 --max-within-max 0.02 \
  --min-parent-d95 0.005 --min-parent-delta 0.001 \
  --min-distinguishable-rate 0.5 --min-target-mapping-rate 0.95 \
  --merge-min-similarity 0.95 --veto-block-windows 15 \
  --umbrella-own-best-floor 0.5 --umbrella-gate structural \
  --max-iterations 60 --shadow-sweep-distance 0.005
```

汇总 + 审计：

```bash
cd <06_group_references 目录>
python3 scripts/mapping_reference/collect_mapping_references.py \
  --input-glob 'mapping_ref_*/final_mapping_references.fasta' \
  --output ../07_final/all_mapping_references.fasta \
  --report ../07_final/final_sequence_qc.tsv --max-n-gap-percent 10.0
python3 scripts/mapping_reference/audit_mapping_reference_qc.py --result-dir … --output …
```

或直接用 Nextflow 编排上游+核心：

```bash
nextflow -C mapping_reference.config run mapping_reference.nf \
  --input_file <gii.fasta> --output_dir <dir> --threads 8
```

## v2.0 实测

| | GII.P16_GII.4（504 序列） | GII.P17_GII.17（664 序列） |
|---|---:|---:|
| 最终参考 | 37 条 | 19 条 |
| 全部真实序列 | ✓ | ✓ |
| 参考间最小距离 | 0.64%（99.36% 相似株共存） | 0.74% |
| 影子 / 成员丢失 | 0 / 0 | 0 / 0 |

GII.17 标注 clade 对照：Kawasaki_308 / Kawasaki_323 / Romania 各获专属 medoid 参考
（K323 23/23 全中），旧版 Node276 坍缩彻底解决。

历史版本的数据准备、阈值校准与验证分析脚本已随清理移除，可在 git 历史中查阅。
