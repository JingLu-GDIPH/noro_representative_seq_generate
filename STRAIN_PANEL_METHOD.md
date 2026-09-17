# 株级代表性参考 Panel 生成方法与参数详解
（Strain-level Representative Reference Panel — Final Recipe v2.0）

本仓库从 NCBI 公共数据库的诺如病毒 GI/GII 全基因组序列出发，构建**全部由真实存在序列组成**、
**互相之间可被 150 bp 污水读段区分**（支持区分基因组相似性 ≥99% 的不同毒株）的参考 panel。

- 上游 Nextflow 流程：`mapping_reference.nf`（分组/比对/建树）
- 核心选择脚本：`scripts/mapping_reference/build_mapping_reference_group.py`（v2.0，本文档主角）
- 英文姊妹实现（standalone、任意病毒可用）：见 `generate_representative_references.py` 所在
  仓库 [viral-representative-references](https://github.com/JingLu-GDIPH/viral-representative-references)

---

## 1. 设计目标与总体思路

**目标**：panel 中 (a) 每条参考都是数据库中真实存在的序列（medoid，非 ASR/共识镶嵌体）；
(b) 任意两条参考可被 150 bp 读段区分（distinguishable rate ≥ 0.5），等效于参考间距 ≥ ~0.46%
基因组距离；(c) ≥99% 基因组相似的毒株对可共存于 panel 并在下游混合样本测序中被分别还原。

**核心结论（结果导向推导，均有实测锚点）**：
- 99% 相似毒株对的 150 bp 窗口实测可区分率 R99 = 0.805（单 SNV 窗口靠 AS margin ≥6 即可判定）
  → 验证阈值 T=0.5 有 1.6× 裕度，无需下调；
- T=0.5 的等效去重粒度 ≈ 0.46% 基因组距离 → 影子参考（更近的近重复）会被自动清除；
- 唯一瓶颈在**选枝粒度**：株级配置（d95 ≤ 0.5%）让 99% 对以独立单元进入流程后，下游机制自动守住。

## 2. 流程步骤（7 步）

```
原始 FASTA → ①质量过滤 → ②RdRp+VP1 分组 → ③拆组 → ④方向归一 + MAFFT + 50% 列过滤
           → ⑤IQ-TREE 建树 + ASR（≥4 条/组） → ⑥核心：株级选择与迭代验证 → ⑦汇总 + QC 审计
```

①–⑤、⑦ 由 `mapping_reference.nf` 编排；⑥ 为 `build_mapping_reference_group.py`（本配方）。

### ① 质量过滤（scripts/common/filter_n_sequences.py）
`(N数+gap数)/长度 > 10%` 移除。

### ② 按 RdRp+VP1 精确分组（scripts/common/group_sequences_by_rdrp_vp1.py）
ID 前两字段正则精确匹配 `G(I|II|IX)\.P*`，组号 `cluster_NNNN`。

### ③ 拆组（scripts/common/split_clusters.py）
每组一个 FASTA，文件名 `cluster__cluster_NNNN__<RdRp>_<VP1>.fasta`。

### ④ 方向归一 + 比对 + 列过滤
- `normalize_orientation.py`：15-mer / 最小重叠 50 / 参考取 A/C/G/T 最多者，反向链校正；
- `mafft --quiet --auto`（去既有 gap 后重比对）；
- `scripts/common/trim_alignment_ends.py --min-coverage 0.5`：**所有列**（端部+内部）A/C/G/T 覆盖 <50% 即删，
  N/gap 不计覆盖。

### ⑤ IQ-TREE（≥4 条/组）
```
iqtree3 -s aligned.fasta -m GTR+F+R4 --alrt 1000 -B 1000 --bnni \
    --ancestral --asr-min 0.8 -T AUTO --prefix iqtree
```
<4 条的组跳过（保留原始序列，标记 SINGLETON/SMALL_GROUP）。

### ⑥ 核心：株级选择与迭代验证（定案命令）
```bash
python3 scripts/mapping_reference/build_mapping_reference_group.py \
  --alignment aligned.fasta \
  --tree iqtree/iqtree.treefile --state iqtree/iqtree.state --mldist iqtree/iqtree.mldist \
  --outdir mapping_ref_<组> --bowtie2 <bowtie2> --threads 4 \
  --representative medoid \
  --min-sh-alrt 80.0 --min-ufboot 0.0 \
  --max-within-d95 0.005 --max-within-max 0.02 \
  --min-parent-d95 0.005 --min-parent-delta 0.001 \
  --min-distinguishable-rate 0.5 --min-target-mapping-rate 0.95 \
  --merge-min-similarity 0.95 \
  --veto-block-windows 15 --veto-window-snv 2 --veto-min-window-coverage 0.5 \
  --umbrella-own-best-floor 0.5 --umbrella-gate structural \
  --max-iterations 60 --shadow-sweep-distance 0.005
```

### ⑦ 汇总 + QC（collect_mapping_references.py / audit_mapping_reference_qc.py）
合并各组（复查 N/gap ≤10%），审计 PASS/REVIEW/UNASSESSABLE/FAIL。

---

## 3. 全部参数与阈值（含依据）

### 3.1 选枝层（初始单元产生）

| 参数 | 值 | 含义 | 依据 |
|---|---|---|---|
| `--min-tips` | 3 | 候选枝最少末端数 | 进化枝统计意义下限 |
| `--min-sh-alrt` | 80 | SH-aLRT 支持 | 分支可信度惯例下限 |
| `--min-ufboot` | **0**（定案弃用） | UFBoot 门槛 | 定案方案仅用 SH-aLRT；medoid 代表不依赖 ASR，无需双支持 |
| `--max-within-d95` | **0.005** | 枝内成对距离 95 分位 ≤0.5% | **株级标尺**：99% 对以独立单元入流程（目标反推） |
| `--max-within-max` | 0.02 | 枝内最远一对 ≤2% | 株级离群守卫 |
| `--min-parent-d95` | 0.005 | 父枝 d95 >0.5% | 株级父子对比 |
| `--min-parent-delta` | 0.001 | 父→子跃升 ≥0.1% | 边界清晰性 |

内置：`homogeneous_root` 捷径（整组满足 d95+max 时整组折一条）；未成枝末端各自成 `observed_raw`
候选单元（株级粒度的主要来源）；QC 排除 `SNP_CLUSTERS>1`（101 nt 窗 >6 私有 SNP 记 1 簇）。

### 3.2 代表序列重建（定案：medoid 真实序列）

| 参数 | 值 | 含义 |
|---|---|---|
| `--representative` | **medoid** | 多成员单元的代表 = **到 clade 根节点 patristic 距离最近的真实成员**；QC 干净者优先，全组被 QC 标记时启用最近成员并标 `MEDOIDQC` 留痕；无树小组退化为经典 medoid |

输出后缀：`_MAPREF_MEDOID` / `_MAPREF_MEDOIDQC` / `_MAPREF_RAW`（单成员 tip）——
**panel 中每一条都是真实存在的序列，无任何合成序列**。

ASR 相关参数（`--asr-pp 0.90 --majority-frequency 0.60 --min-site-coverage 0.50
--max-low-pp-fraction 0.05 --max-n-fraction 0.005 --max-n-run 15`）仅在
`--representative asr`（旧行为，hybrid ASR/多数/medoid 回退链）时生效。

### 3.3 验证层（150 bp 滑窗竞争比对）

每轮对全部当前参考做 Bowtie2 `--very-sensitive-local -k 2` 竞争比对，逐窗判定：

| 参数 | 值 | 含义 | 依据 |
|---|---|---|---|
| `--window / --step` | 150 / 25 | 窗口/步长 | 匹配污水测序读长 |
| `--min-aligned-window` | 135 | 窗口最小比对长度 | 90% of 150 |
| `--min-score-margin` | 6 | AS 分差判据 | 高质量下 1 个错配罚分=6，**单 SNV 窗即可判**（R99=0.805 的来源） |
| `--min-nm-margin` | 2 | 错配数差判据 | 稳健判据 |
| `--min-distinguishable-rate` | **0.5** | 可区分率阈值 | R99=0.805 之上 1.6×；等效去重粒度 0.46%（实测终态最小间距 0.64–0.74% 吻合） |
| `--min-target-mapping-rate` | 0.95 | 成员读段"回家率"（≤10% 错配直接比较） | 粗筛（~10% 量级），抓严重失代表 |

### 3.4 决策路由（每单元每轮）

```
target < 0.95                                    → split   （失代表：参考不像成员）
dist < 0.5 且 own_best < 0.5 且 结构性门槛通过     → split   （伞形：读段多数投别人）
dist < 0.5 其余                                   → merge   （参考互抢读段：低于读段分辨率）
其余                                              → retain
```

| 参数 | 值 | 含义 | 依据 |
|---|---|---|---|
| `--umbrella-own-best-floor` | 0.5 | 伞形判据线 | 大伞 own_best 实测 0.205 |
| `--umbrella-gate` | **structural**（定案） | 门槛=**可行性判据**：节点有 ≥2 个子枝、各含本单元自有成员 ≥ `--min-tips` | 拆分路由与可执行性一致；无数据集相关常数。备选 `members`（≥`--umbrella-min-members` 25，经验值） |

三臂对照实测（GII.P16_GII.4）：无门槛 102 条（恶化）/ 尺寸门槛 35 条 / **结构性门槛 37 条（定案）**。

### 3.5 合并/拆分执行

- **批量合并不相交对**：候选 = ≥1 方 dist<0.5 ∧ 相似度 ≥ `--merge-min-similarity 0.95`；
  按双向混淆总数降序贪心。**MRCA 级联已移除**（只合触发对，杜绝 240→2 雪崩）。
- **VP1 片段否决**：`--veto-block-windows 15 --veto-window-snv 2 --veto-min-window-coverage 0.5`
  ——两参考在 VP1 区（注释转移定位，实测落 67–89% 基因组处）存在连续 ≥15 窗、每窗 ≥2 SNV
  的差异块即**禁止合并**。依据：金标准变异株对 VP1 块 ≥29 窗、99% 对 20 窗、近重复 ≤3 窗；
  P 域浓缩系数实测 1.82（差异确实集中于 VP1）。可比列 <50% 的窗口不可判（部分基因组守卫），
  VP1 证据不足仍放行的对记 `allowed_low_vp1_coverage` 留痕。
- **批量拆分**：一轮拆掉所有 split-flagged 非末端单元（沿子枝、限本单元自有成员，要求 ≥2 替换）。

### 3.6 收敛控制与终末清扫

| 机制 | 说明 |
|---|---|
| `--max-iterations 60` | 轮数上限 |
| 振荡检测 | **成员划分签名**（与标签无关）重现即停 `oscillation_detected`（防伞拆-合回极限环） |
| 最优快照 | 全程记录 `(-flagged 数, -单元数)` 最优状态；终态输出该快照的参考集与验证表 |
| `--shadow-sweep-distance 0.005` | 停机后所有两两 ≤0.5% 且未被 VP1 否决保护的对强制合并 + `final_sweep/` 重验证——**硬保证 panel 无目标分辨率以下的影子参考**（影子会杀死读段区分度：实测 0.805→0.000） |

---

## 4. 实测验证（截至 v2.0）

| | GII.P16_GII.4（504 序列） | GII.P17_GII.17（664 序列） | GII.17 全部 7 组 |
|---|---:|---:|---:|
| 最终参考 | 37 | 19 | 34 |
| 全部真实序列 | ✓ | ✓ | ✓（10 MEDOID + 10 MEDOIDQC + 14 RAW） |
| 参考间最小距离 | 0.64% | 0.74% | — |
| 影子 <0.5% | 0 | 0 | 0 |
| 成员丢失 | 0 | 0 | 0 |
| 历史 | 240→2 雪崩已杜绝 | Node276 坍缩已解决：K308/K323/Romania 各获专属 medoid 参考（K323 23/23 全中） |

下游还原能力（150 bp 模拟，100×覆盖、0.3% 错误）：分开 panel 下 99% 相似双株
**碱基错误率 0.000%、判别位点 100% 正确、嵌合 0%**，次要株下探 10%；合并 panel 则判别位点
全 N（50:50）或全归优势株（90:10）——**还原精度由 panel 粒度决定**。

## 5. 输出文件（每组 `mapping_ref_<组>/`）

`final_mapping_references.fasta`（终 panel）、`final_reference_sources.tsv`（来源/成员/不确定性）、
`final_membership.tsv`、`final_window_validation.tsv`（4 率+决策）、`final_window_details.tsv`、
`iteration_summary.tsv`、`vp1_veto_events.tsv`、`raw_reference_qc.tsv`、`node_candidate_metrics.tsv`、
`sliding_windows.fasta`、`iteration_NN/`、`final_sweep/`（清扫后重验证）、`summary.txt`。

## 6. 复现资产

阈值校准与验证分析脚本（test_vp1_veto_golden.py、derive_params_forward.py、
gii4/gii17_threshold_calibration.py、analyze_gii4_variant_match.py 等）及出版图
（result_paper/Fig_threshold_calibration）在 2026-09-15 仓库清理中移出工作树，
完整内容见 git 历史（v2.0.0 提交）。
