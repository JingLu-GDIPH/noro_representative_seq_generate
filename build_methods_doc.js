// Generate the norovirus pipeline methods & parameters DOCX document.
// Uses the globally-installed docx library.

const docx = require('/usr/local/lib/node_modules/docx');
const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  Header, Footer, AlignmentType, LevelFormat, HeadingLevel, ImageRun,
  BorderStyle, WidthType, ShadingType, PageNumber, PageBreak, TabStopType, TabStopPosition
} = docx;

const FONT = "Arial";
const CONTENT_WIDTH = 9360; // US Letter - 2*1" margins

// ---------- helpers ----------
function h1(text) {
  return new Paragraph({ heading: HeadingLevel.HEADING_1, children: [new TextRun({ text, bold: true })] });
}
function h2(text) {
  return new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun({ text, bold: true })] });
}
function h3(text) {
  return new Paragraph({ heading: HeadingLevel.HEADING_3, children: [new TextRun({ text, bold: true })] });
}
// runs: a string, a single {text, bold?, italics?} object, or an array of them
function normRuns(runs) {
  if (typeof runs === "string") return [{ text: runs }];
  if (Array.isArray(runs)) return runs;
  return [runs];
}
function para(runs, opts = {}) {
  const children = normRuns(runs).map(r =>
    new TextRun({ text: r.text, bold: !!r.bold, italics: !!r.italics, font: FONT })
  );
  return new Paragraph({ children, spacing: { after: 120, line: 300 }, ...opts });
}
function bullet(runs, level = 0) {
  const children = normRuns(runs).map(r =>
    new TextRun({ text: r.text, bold: !!r.bold, italics: !!r.italics, font: FONT })
  );
  return new Paragraph({ numbering: { reference: "bullets", level }, children, spacing: { after: 60, line: 290 } });
}
function numbered(runs) {
  const children = normRuns(runs).map(r =>
    new TextRun({ text: r.text, bold: !!r.bold, italics: !!r.italics, font: FONT })
  );
  return new Paragraph({ numbering: { reference: "numbers", level: 0 }, children, spacing: { after: 60, line: 290 } });
}
function pageBreak() {
  return new Paragraph({ children: [new PageBreak()] });
}
function spacer() {
  return new Paragraph({ children: [new TextRun("")] });
}

// build a table from header list + rows (arrays of strings). colWidths sum to <= CONTENT_WIDTH.
function table(headers, rows, colWidths) {
  const total = colWidths.reduce((a, b) => a + b, 0);
  const border = { style: BorderStyle.SINGLE, size: 1, color: "BBBBBB" };
  const borders = { top: border, bottom: border, left: border, right: border };
  const cell = (text, w, header = false) => new TableCell({
    borders,
    width: { size: w, type: WidthType.DXA },
    shading: header ? { fill: "D5E8F0", type: ShadingType.CLEAR } : undefined,
    margins: { top: 60, bottom: 60, left: 100, right: 100 },
    children: [new Paragraph({ children: [new TextRun({ text: String(text), bold: header, font: FONT, size: 20 })] })]
  });
  const headerRow = new TableRow({ tableHeader: true, children: headers.map((t, i) => cell(t, colWidths[i], true)) });
  const bodyRows = rows.map(r => new TableRow({ children: r.map((t, i) => cell(t, colWidths[i])) }));
  return new Table({
    width: { size: total, type: WidthType.DXA },
    columnWidths: colWidths,
    rows: [headerRow, ...bodyRows]
  });
}

// ====================================================================
// DOCUMENT CONTENT
// ====================================================================
const children = [];

// ---- Cover ----
children.push(new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { before: 2400, after: 240 },
  children: [new TextRun({ text: "诺如病毒代表性序列生成流水线", bold: true, size: 44, font: FONT })]
}));
children.push(new Paragraph({
  alignment: AlignmentType.CENTER,
  spacing: { after: 480 },
  children: [new TextRun({ text: "分析方法与参数设置原理说明文档", bold: true, size: 36, font: FONT })]
}));
children.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 120 },
  children: [new TextRun({ text: "Norovirus GI / GII Representative Sequence Generation", italics: true, size: 24, font: FONT, color: "666666" })]
}));
children.push(new Paragraph({
  alignment: AlignmentType.CENTER, spacing: { after: 2400 },
  children: [new TextRun({ text: "IPHnano · 分子流行病学与基因组学分析流程", size: 22, font: FONT, color: "666666" })]
}));

const coverInfo = table(
  ["项目", "说明"],
  [
    ["流程代码版本", "main.nf v1.4.0；mapping_reference.nf + build_mapping_reference_group.py v2.0.0（株级 medoid 定案配方）"],
    ["文档生成日期", "2026-09-15"],
    ["工作目录", "noro_representative_seq_generate/"],
    ["运行环境", "conda 环境 noro-consensus（Python 3.8+、Biopython、MAFFT、IQ-TREE / IQ-TREE 3、Bowtie2、FastTree）"],
    ["并发与资源", "默认 --threads 8（10 核工作站，预留 2 核）"],
    ["数据来源", "NCBI GenBank 诺如病毒 GI/GII 核苷酸记录（序列长度 >5800 bp，采集日期 ≥2000 年）"],
  ],
  [2400, 6960]
);
children.push(coverInfo);
children.push(pageBreak());

// ====================================================================
// 1. 项目概述
// ====================================================================
children.push(h1("1  项目概述"));

children.push(h2("1.1  研究目的"));
children.push(para("本仓库提供两条独立的自动化生物信息学流水线，用于从大规模诺如病毒（Norovirus）GI/GII 基因组序列中生成高质量的代表性参考序列。两条流水线共享相同的序列预处理与基因型分组逻辑，但在“代表性序列”的定义、建树方法与去冗余策略上存在本质差异，以服务两类不同的下游应用。"));

children.push(h2("1.2  两条流程的定位与关系"));
children.push(table(
  ["流程", "入口文件", "代表性序列定义", "目标下游应用", "结果目录"],
  [
    ["流程一 · 探针共识流程", "main.nf (v1.4.0)", "95% 相似性边界节点对应的祖先/共识序列", "tiling 探针设计、分子流行病学参考集", "results_4probe/"],
    ["流程二 · 映射参考流程（株级 medoid 定案）", "mapping_reference.nf + build_mapping_reference_group.py v2.0", "全部由真实存在序列（medoid）组成，支持 ≥99% 基因组相似毒株共存，经 150 bp 读段可区分性验证", "150 bp 污水测序读段的竞争性比对参考", "result_4ref/"],
  ],
  [1900, 1700, 2300, 2060, 1400]
));
children.push(para("两条流程相对独立：流程一强调序列间的遗传多样性边界（成对相似性 ≤ 95%），流程二强调参考序列在短读段比对下的“可区分性”（differentiability）。它们共用序列质量过滤、基因型分组、方向归一化、MAFFT 比对与覆盖度列过滤等前处理步骤。"));

children.push(h2("1.3  实际运行结果总览"));
children.push(table(
  ["指标", "GI", "GII"],
  [
    ["输入序列（质量过滤后）", "374", "3697"],
    ["RdRp+VP1 基因型分组数", "17", "61"],
    ["映射参考 panel（v2.0 定案，全真实序列）", "—（GI 待跑）", "GII.P16_GII.4 组 37 条 / GII.P17_GII.17 组 19 条 / GII.17 七组共 34 条；参考间最小间距 0.64–0.74%（≥99% 相似株共存），影子 0"],
    ["探针共识序列（results_4probe，原始/过滤后）", "42 / 29", "117 / 92"],
    ["最终验证最大组内相似性", "94.45%", "94.97%"],
    ["探针覆盖基因型数", "9 / 9（全覆盖）", "17 / 28（缺 11 个型）"],
  ],
  [4200, 2580, 2580]
));
children.push(para([
  { text: "说明：", bold: true },
  { text: "“过滤后”指进一步移除被 nextclade 风格原始代表序列 QC 标记为 SNP_CLUSTERS>1 的疑似测序错误序列及单例（singleton）代表序列后的结果。“组内最大相似性”为 08_validation 步骤计算的所有成对相似性最大值，均低于 95% 阈值，满足多样性约束。" },
]));

children.push(pageBreak());

// ====================================================================
// 2. 软件环境
// ====================================================================
children.push(h1("2  软件工具与运行环境"));

children.push(h2("2.1  工作流管理"));
children.push(para("两条流程均以 Nextflow（≥20.0.0）实现。Nextflow 基于数据流编程模型，将各处理步骤（process）以 DAG 形式连接，支持本地并行、断点续跑（-resume）与可重复执行。流程通过 nextflow.config / mapping_reference.config 统一配置本地执行器（local executor）、线程预算、每个 process 的 CPU/内存/时间上限与 conda 环境（/Users/LuJ/mambaforge/envs/noro-consensus）。"));

children.push(h2("2.2  核心生物信息学软件"));
children.push(table(
  ["软件", "版本/调用", "在本流程中的作用"],
  [
    ["MAFFT", "--auto（多序列簇）；--retree 1/2（大簇）", "多序列比对；去 gap 后重新比对以避免误把 padding 当作同源位点"],
    ["IQ-TREE", "-m JC -asr（流程一）", "流程一：最大似然建树 + 祖先序列重建（ASR），用于相似性边界节点定位"],
    ["IQ-TREE 3", "-m GTR+F+R4 --alrt 1000 -B 1000 --bnni --ancestral --asr-min 0.8（流程二）", "流程二：带分支支持（SH-aLRT + UFBoot + BNNI）的 ML 建树与 ASR"],
    ["Bowtie2", "--very-sensitive-local -k 2", "流程二：150 bp 滑动窗口读段的竞争性比对，验证参考可区分性"],
    ["FastTree", "-nt -gtr -nosupport", "代表性评估：按 VP1 基因型快速建树，计算最近共识距离等指标"],
    ["BLAST+", "makeblastdb / blastn -task blastn -evalue 1e-10", "基因型分型：RdRp/VP1 分区参考序列比对"],
    ["vclust", "（早期版本用于 ANI 聚类，当前流程已弃用）", "历史聚类工具；当前流程一改用 RdRp+VP1 精确分组"],
  ],
  [1700, 3300, 4360]
));

children.push(h2("2.3  编程语言与库"));
children.push(para("自定义脚本使用 Python 3.8+ 编写，主要依赖：Biopython（FASTA/比对/系统发育树读写与操作）、Pandas（表格数据处理）、NumPy（向量化成对相似性计算）。流程二的滑窗验证使用 Bowtie2 命令行而非 Python 比对库，以保证与实际读段映射软件一致。"));

children.push(pageBreak());

// ====================================================================
// 3. 公共预处理
// ====================================================================
children.push(h1("3  公共预处理步骤（两条流程共用）"));
children.push(para("以下步骤在两条流程中实现一致，是后续建树与代表性序列生成的基础。"));

// 3.1 quality filter
children.push(h2("3.1  序列质量过滤（filter_n_sequences.py）"));
children.push(para("对输入 FASTA 逐条扫描，计算未知碱基（N）与缺口（gap，“-”）的合并占比，移除超过阈值的低质量序列。"));
children.push(para([
  { text: "质量评分公式：", bold: true },
  { text: "  Quality = (N_count + Gap_count) / L_sequence × 100%。当 Quality > 10% 时移除该序列。N 与 gap 的合并统计可同时惩罚测序模糊碱基与比对/拼接缺口。" },
]));
children.push(table(
  ["参数", "设定值", "说明"],
  [
    ["--threshold", "10.0", "N 与 gap 合并百分比上限（严大于判定）"],
  ],
  [3000, 2000, 4360]
));

// 3.2 genotyping chain
children.push(h2("3.2  RdRp + VP1 双区基因型分型链"));
children.push(para("代表性序列的分组以 RdRp（RNA 依赖 RNA 聚合酶，ORF1）与 VP1（主要衣壳蛋白，ORF2）双区基因型为硬分层基础。分型由三步链式脚本完成："));
children.push(numbered([
  { text: "direct_group_norotyping.py（分型）：", bold: true },
  { text: "对 GI/GII 分别构建 RdRp+VP1 参考 BLAST 库（makeblastdb -dbtype nucl；blastn -task blastn -max_hsps 1 -evalue 1e-10）。对每个区分别执行 BLAST + 系统发育置树（MAFFT 比对后建树，将 query 放到参考树最近位置），输出 rdrp/vp1 最近参考、基因型、树距离、BLAST 一致性与参考覆盖度。状态分为 typed（双区均定出）/ partial_genotype / untyped。" },
]));
children.push(numbered([
  { text: "update_ncbi_metadata_with_norotyping.py（落库 + 重命名）：", bold: true },
  { text: "将分型结果合并入 NCBI 下载元数据表，并把 FASTA 记录重命名为规范 seqname：<RdRp>_<VP1>_<accession>_<collection_date>。该校验序列无重复 seqname，并保证元数据与 FASTA 一一对应。" },
]));
children.push(numbered([
  { text: "group_sequences_by_rdrp_vp1.py（分组）：", bold: true },
  { text: "按规范 ID 前两字段（RdRp、VP1）精确字符串匹配分组，输出 cluster_info.csv 与分组 FASTA。" },
]));
children.push(para([
  { text: "重要更正：", bold: true },
  { text: "早期版本（及旧版方法文档）使用 vclust 的 ANI 层次聚类（ANI 90%、完全连锁、qcov/rcov/len_ratio 65%）。当前两条流程均改为按 RdRp+VP1 基因型精确分组，不再计算 ANI，亦不再进行层次聚类。其生物学依据是：诺如病毒 RdRp-VP1 重组常见，按双区基因型对硬分层可从源头避免不同基因型或重组型被局部高相似比对错误合并；同一基因型内的进一步分化交由下游系统发育建树处理。" },
]));

// 3.3 orientation
children.push(h2("3.3  序列方向归一化（normalize_orientation.py）"));
children.push(para("公共数据库下载的序列并非都以正义链存储（例如 GII.P31_GII.4_KX158285、GII.P7_GII.6_KX158282 / KX158286 以负链下载）。若直接送入 MAFFT，反向链会被强行对齐到错误方向，产生约 1500 个虚假 gap 与无意义的相似性。本步骤在比对前统一链方向。"));
children.push(para("方法：对每个分组，选取 A/C/G/T 有效碱基数最多（非原始长度）的序列作为方向参考；对其余每条序列分别构建正义链与反向互补链的 k-mer 集合，统计与参考 k-mer 集合的重叠数。仅当反向互补链重叠严格大于正义链、且达到最低重叠阈值时，对该序列做反向互补。单序列分组与参考过短（短于 k）时直接透传。"));
children.push(table(
  ["参数", "默认值", "说明"],
  [
    ["--kmer-size", "15", "15-mer 足够特异，两条无关病毒链几乎不可能偶然共享大量 15-mer"],
    ["--min-overlap", "50", "触发反向互补所需的最低 RC k-mer 重叠数，避免短偶然重叠误判"],
    ["方向参考选择", "max(A/C/G/T count)", "以有效碱基数最多者为参考，而非最长序列，避免 N/padding 干扰"],
  ],
  [2500, 2200, 4660]
));

// 3.4 alignment
children.push(h2("3.4  多序列比对与列覆盖度过滤"));
children.push(para([
  { text: "重新比对：", bold: true },
  { text: "等长序列可能仅是 padding 结果，并不代表同源位点已对齐。流程对每个多序列簇先移除已有 gap，再用 MAFFT --auto（大簇改用 --retree 1/2 加速）重新比对；单序列簇直接保留。" },
]));
children.push(para([
  { text: "列覆盖度过滤（trim_alignment_ends.py）：", bold: true },
  { text: "比对完成后，删除所有 A/C/G/T 覆盖度低于阈值的 alignment 列。N 与 gap 不计为有效覆盖，因此低覆盖的内部插入列与端部低覆盖区域同时被屏蔽，避免人工拉长共识序列。" },
]));
children.push(table(
  ["参数", "设定值", "说明"],
  [
    ["--min-coverage", "0.5", "列内 A/C/G/T 占比下限；仅 A/C/G/T 计有效覆盖"],
    ["有效碱基集合", "{A,C,G,T}", "N 与 “-” 均不计覆盖，也不参与相似性分母（贯穿全流程）"],
  ],
  [3000, 2000, 4360]
));

children.push(pageBreak());

// ====================================================================
// 4. Pipeline 1
// ====================================================================
children.push(h1("4  流程一：探针代表性共识序列生成（main.nf v1.4.0）"));
children.push(para("本流程以“95% 成对相似性”为多样性边界，输出用于探针设计与流行病学参考的代表序列。包含 9 个 Nextflow process。"));

children.push(h2("4.1  流程步骤总览"));
children.push(table(
  ["步骤", "Process", "脚本", "说明"],
  [
    ["1", "MERGE_FASTA", "filter_n_sequences.py", "复制输入并按 10% N+gap 阈值过滤"],
    ["2", "CLUSTER", "group_sequences_by_rdrp_vp1.py", "按 RdRp+VP1 基因型精确分组"],
    ["3", "SPLIT_CLUSTERS", "split_clusters.py", "按簇拆分为独立 FASTA（文件名含基因型）"],
    ["4", "ORIENT_SEQUENCES", "normalize_orientation.py", "比对前统一链方向（15-mer）"],
    ["5", "ALIGN", "MAFFT + trim_alignment_ends.py", "去 gap 重比对 + 50% 列覆盖度过滤"],
    ["6", "CONSENSUS", "generate_consensus.py", "JC 模型建树 + ASR + 相似性边界节点共识"],
    ["7", "ENHANCED_CONSENSUS", "enhanced_consensus_with_validation.py + trim_sequence_ends.py", "迭代验证 + 树引导合并 + 最终去 gap"],
    ["8", "COLLECT_RESULTS", "combine_consensus_fastas.py", "合并各簇共识序列，保证 ID 全局唯一"],
    ["9", "FINAL_VALIDATION", "validate_consensus_internal.py", "全局成对相似性最终验证并出具报告"],
  ],
  [600, 2200, 3200, 3360]
));

children.push(h2("4.2  基于系统发育的共识序列生成（generate_consensus.py）"));
children.push(h3("4.2.1  建树与祖先序列重建"));
children.push(para("对 ≥3 条序列的簇调用 IQ-TREE，使用 Jukes-Cantor（JC）核酸替换模型并开启祖先序列重建（-asr）。JC 模型假设四种核苷酸间替换概率相等，是核酸序列最基础模型，对诺如病毒这类近缘、高度保守的基因组比较足够准确且计算高效。"));
children.push(para([
  { text: "中点定根：", bold: true },
  { text: "采用 Biopython 的 tree.root_at_midpoint() 实现中点定根——找到树中最远的两个末端，在其路径中点设根，使树具备有根的层次结构以支持“父-子”节点关系判定。" },
]));
children.push(para([
  { text: "ASR 后验概率说明：", bold: true },
  { text: "本首遍脚本仅使用 -asr（未显式设 --asr-min），直接读取 IQ-TREE .state 文件中各节点的 ML 状态序列；显式的 80% 后验概率阈值（--asr-min 0.8）在步骤 7 的增强迭代脚本中施加（见 4.3）。" },
]));

children.push(h3("4.2.2  最优（相似性边界）节点选择"));
children.push(para("对中点定根后的树遍历所有节点（内部 + 末端）。对每个节点 v，定义其内部相似性 S_min(v) 为其所有后代末端序列两两相似性的最小值。节点 v 被选为“最优节点”当且仅当："));
children.push(para([
  { text: "条件一：", bold: true },
  { text: " S_min(v) ≥ θ（95%）——该分支内部已足够相似；" },
]));
children.push(para([
  { text: "条件二：", bold: true },
  { text: " S_min(parent(v)) < θ（95%）——其父分支已跨越相似性边界。" },
]));
children.push(para("即从高相似区域向低相似区域转变的“边界节点”。对满足条件的内部节点输出其祖先序列；对不满足的末端则保留原始序列为代表（含低于阈值的末端回退与根回退兜底）。该策略无需预先指定分组数，自适应定位遗传多样性的自然分界。"));

children.push(h3("4.2.3  成对相似性公式"));
children.push(para("Similarity = M / N_valid × 100%，其中 M 为两序列在相同位置均为同一 A/C/G/T 的位点数，N_valid 为两序列相同位置均属 {A,C,G,T} 的位点数。N 与 gap 同时排除在分子与分母之外，避免低覆盖区域虚增相似性。该公式以 NumPy 向量化实现（VALID_BASE_BYTES 比对，分块 chunk_size=256）。"));

children.push(h3("4.2.4  少序列簇的简化处理"));
children.push(bullet("单序列：直接输出，ID 形如 <genotype>_cluster_<n>_node_1。"));
children.push(bullet("双序列：若成对相似性 ≥ 95% 则取首条作为唯一共识；否则两条均作为代表（node_1/node_2）。"));
children.push(bullet("多序列建树失败：异常时回退到逐位点多数投票（仅统计 A/C/G/T）。"));

children.push(h2("4.3  增强迭代验证（enhanced_consensus_with_validation.py）"));
children.push(para("对步骤 6 输出的共识序列进行迭代多样性保证。核心循环："));
children.push(numbered("计算当前共识集合的最大成对相似性 S_max。"));
children.push(numbered("若 S_max < 95%（严格小于），判定收敛，输出结果；否则进入下一轮。"));
children.push(numbered("对超阈序列重新执行 IQ-TREE 建树（JC 模型）+ 祖先序列重建，施加 --asr-min 0.8（80% 后验概率阈值）。"));
children.push(numbered("采用“树引导合并”策略：对相似性 ≥ 阈值的序列对，不直接删除，而是合并为新的多数共识序列（按系统发育距离排序成对处理，每条仅消费一次）。"));
children.push(numbered("重复直到满足条件或达到最大迭代次数。"));
children.push(table(
  ["参数", "设定值", "说明"],
  [
    ["--similarity_threshold", "95.0", "迭代收敛阈值（最大成对相似性须低于此值）"],
    ["--max_iterations", "10", "最大迭代次数；实际通常 3–5 次收敛"],
    ["IQ-TREE 模型", "JC (-m JC)", "与首遍一致"],
    ["--asr-min", "0.8", "祖先状态后验概率下限 80%（仅此步骤施加）"],
    ["--threads", "4（每任务）", "multiprocessing 池上限 min(cpu_count, 8)"],
  ],
  [2800, 1800, 4760]
));

children.push(h2("4.4  最终修剪与去冗余"));
children.push(para([
  { text: "去 gap 写出（trim_sequence_ends.py --remove-gaps）：", bold: true },
  { text: "最终代表序列写出前，先剥离两端 N/n/-，再移除所有由比对引入的 gap，避免低覆盖插入列被填入共识序列。" },
]));
children.push(para([
  { text: "合并（combine_consensus_fastas.py）：", bold: true },
  { text: "合并各簇 FASTA，对重复 ID 追加 _dupN 后缀保证全局唯一。当前版本不再进行额外的贪心跨簇删除——冗余已由迭代树引导共识处理。" },
]));
children.push(para([
  { text: "最终验证（validate_consensus_internal.py）：", bold: true },
  { text: "若序列不等长，使用更严格的 MAFFT --maxiterate 1000 --localpair 重比对；随后计算所有成对相似性（同样仅 {A,C,G,T} 计有效位点）。当最大成对相似性 ≥ 阈值时标记 Needs reclustering=true，并支持按 genotype_pair 分组报告。验证采用更高质量比对模式以保证下游判定稳健。" },
]));

children.push(pageBreak());

// ====================================================================
// 5. Pipeline 2
// ====================================================================
children.push(h1("5  流程二：株级 medoid 映射参考生成（mapping_reference.nf + build_mapping_reference_group.py v2.0 定案）"));
children.push(para("本流程面向 150 bp 污水读段的竞争性比对参考。v2.0 定案配方的三个核心特征：(1) 代表序列全部为真实存在的序列（medoid）——取到 clade 根节点系统发育距离最近的成员，杜绝 ASR/共识镶嵌体偏离真实毒株的问题；(2) 株级选枝（枝内 d95≤0.5%）使 ≥99% 基因组相似的毒株以独立单元进入流程；(3) 150 bp 滑窗竞争性验证保证任意两条参考可被读段区分（等效距离 ≥约 0.46%），并由 VP1 片段否决保护真实谱系差异。"));

children.push(h2("5.1  流程步骤总览"));
children.push(table(
  ["步骤", "Process/脚本", "说明"],
  [
    ["1", "FILTER_INPUT / filter_n_sequences.py", "10% N+gap 过滤"],
    ["2", "GROUP_BY_GENOTYPE / group_sequences_by_rdrp_vp1.py", "RdRp+VP1 基因型精确分组"],
    ["3", "SPLIT_GROUPS / split_clusters.py", "拆分为独立 FASTA"],
    ["4", "ORIENT_AND_ALIGN / normalize_orientation + MAFFT + trim_alignment_ends", "方向归一（15-mer/重叠50）+ 重比对 + 50% 列覆盖过滤"],
    ["5", "IQTREE_ASR / iqtree3", "≥4 序列/组：GTR+F+R4、SH-aLRT/UFBoot 1000、BNNI、ASR"],
    ["6", "SELECT_MAPPING_REFERENCES / build_mapping_reference_group.py", "株级选枝 + medoid 代表 + 150bp 验证 + 批量合并/拆分 + VP1 否决 + 影子清扫"],
    ["7", "COLLECT_MAPPING_REFERENCES / collect_mapping_references.py", "汇总 + 二次 10% N/gap 复查"],
    ["8", "audit_mapping_reference_qc.py（线下）", "PASS/REVIEW/UNASSESSABLE/FAIL 审计"],
  ],
  [600, 2600, 6160]
));

children.push(h2("5.2  IQ-TREE 建树（上游，不变）"));
children.push(para("仅对 ≥4 条序列的组执行：iqtree3 -m GTR+F+R4 --alrt 1000 -B 1000 --bnni --ancestral --asr-min 0.8。GTR+F+R4（广义时间可逆+经验频率+4类自由速率异质性）配 SH-aLRT+UFBoot+BNNI 提供可靠分支支持与 ML 距离矩阵（.mldist）。<4 条的组保留原始序列（标记 SINGLETON/SMALL_GROUP，不可判定状态）。"));

children.push(h2("5.3  株级选枝（v2.0 关键改动）"));
children.push(table(
  ["判据", "参数", "定案值"],
  [
    ["分支支持 SH-aLRT", "--min-sh-alrt", "≥ 80（唯一支持门槛）"],
    ["分支支持 UFBoot", "--min-ufboot", "0（定案弃用：medoid 代表不依赖 ASR，无需双支持）"],
    ["枝内 d95（成对 ML 距离 95 分位）", "--max-within-d95", "≤ 0.005（株级标尺）"],
    ["枝内 max（最远一对）", "--max-within-max", "≤ 0.02（株级离群守卫）"],
    ["父枝 d95", "--min-parent-d95", "> 0.005"],
    ["父-子 d95 跃升", "--min-parent-delta", "≥ 0.001"],
    ["最少末端数", "--min-tips", "≥ 3"],
  ],
  [3400, 2400, 3560]
));
children.push(para("株级标尺的含义：99% 基因组相似的毒株对（距离约 1%）不再被装进同一个边界枝——它们以独立单元（候选枝或末端 observed_raw）进入验证层。实测 GII.P16_GII.4 初始产生 220 个株级单元。未成枝的末端各自成为 observed_raw 候选；QC 排除 SNP_CLUSTERS>1（101 nt 窗 >6 私有 SNP 记 1 簇）的序列。"));

children.push(h2("5.4  medoid 代表序列（v2.0 核心改动）"));
children.push(para([
  { text: "--representative medoid：", bold: true },
  { text: "多成员单元的代表 = 到该 clade 根节点 patristic 距离最近的真实成员（用 IQ-TREE 树直接计算）。QC 干净成员优先；全组被 SNP 簇 QC 标记时启用距离最近的被标记成员并以后缀 MEDOIDQC 显式留痕（保证谱系零丢失）。无树小组退化为经典 medoid（平均成对距离最小者）。单成员单元即该序列本身（observed_raw）。" }
]));
children.push(para([
  { text: "输出 ID：", bold: true },
  { text: "<RdRp>_<VP1>_MAPREF_refNNN_<节点>_<MEDOID/MEDOIDQC>，原始序列型为 <序列ID>_MAPREF_RAW——panel 中每一条都是数据库中真实存在的序列，无任何合成序列。旧 ASR 重建链（hybrid ASR：后验≥0.90 取祖先碱基→覆盖≥50%且多数≥60% 取多数→medoid 回退）完整保留为 --representative asr 选项，供对照。" }
]));

children.push(h2("5.5  150 bp 滑窗可区分性验证"));
children.push(para("从每条活跃成员生成 150 bp 窗口（步长 25，全 A/C/G/T 窗才入分母），对全部当前参考做 Bowtie2 -k 2 --very-sensitive-local 竞争比对。逐窗计算三个率：target_mapping_rate（读段与自家参考同坐标切片直接比较，≥135 可比列且错配≤10%，不经比对软件）、own_best_rate（最佳命中为自家参考）、distinguishable_rate（最佳命中为自家且对次佳参考有 AS 分差≥6 或错配差≥2 的决定性优势）。单 SNV 窗口靠 AS margin≥6 即可判定——这是 99% 相似对可分率高达 0.805 的机理。"));

children.push(h2("5.6  决策路由与批量合并/拆分"));
children.push(para("每单元每轮按四条路径判定："));
children.push(bullet([{ text: "target < 0.95 → split：", bold: true }, { text: "参考不像其成员（约 10% 量级失代表），沿子枝拆分；" }]));
children.push(bullet([{ text: "dist < 0.5 且 own_best < 0.5 且通过结构性门槛 → split（伞形）：", bold: true }, { text: "成员读段多数最佳命中落在其他参考上——伞罩不住成员。结构性门槛（--umbrella-gate structural，定案）= 节点有 ≥2 个子枝、各含本单元自有成员 ≥ min_tips——即拆得出合格孩子，无数据集相关常数。备选 members 门槛（≥25 成员）为经验值；" }]));
children.push(bullet([{ text: "dist < 0.5 其余 → merge：", bold: true }, { text: "参考互抢读段（低于读段分辨率 ≈0.46%），批量合并不相交对（≥1 方失败 ∧ 相似度 ≥0.95 ∧ VP1 否决放行），按双向混淆总数降序贪心。MRCA 级联合并已移除（只合触发对，杜绝 240→2 雪崩）；" }]));
children.push(bullet([{ text: "其余 → retain。", bold: true }]));
children.push(para([
  { text: "VP1 片段级合并否决（--veto-block-windows 15 --veto-window-snv 2 --veto-min-window-coverage 0.5）：", bold: true },
  { text: "两条参考在 VP1 区（注释转移定位：对每组共识用 mafft --add --keeplength 加入同 VP1 型 genbank 注释参考，转移其 VP1 坐标，实测落在基因组 67–89% 处）存在连续 ≥15 窗、每窗 ≥2 SNV 的差异块即禁止合并——保护整体相似但衣壳有真实块差的谱系。可比列 <50% 的窗口不可判（部分基因组守卫）；VP1 证据不足仍放行的对记 allowed_low_vp1_coverage 留痕。" }
]));

children.push(h2("5.7  收敛控制与终末影子清扫"));
children.push(bullet([{ text: "振荡检测：", bold: true }, { text: "成员划分签名（frozenset，与单元标签无关）重现即停——防伞拆-合回极限环（实测 3 轮环）；" }]));
children.push(bullet([{ text: "最优快照：", bold: true }, { text: "全程记录 (-flagged 数, -单元数) 最优的状态，终态输出该快照的参考集与验证表（split 分支每轮优先于 merge，最后一轮未必最优）；" }]));
children.push(bullet([{ text: "影子清扫（--shadow-sweep-distance 0.005）：", bold: true }, { text: "停机后所有两两 ≤0.5% 且未被 VP1 否决保护的对强制合并并在 final_sweep/ 重验证——硬保证 panel 无目标分辨率以下的影子参考（实测影子参考会把 99% 对的可分率从 0.805 杀到 0.000）。" }]));

children.push(h2("5.8  汇总与输出"));
children.push(para("collect_mapping_references.py 合并各组（二次 10% N/gap 复查、去重 ID）。每组输出：final_mapping_references.fasta、final_reference_sources.tsv（来源/成员/节点）、final_membership.tsv、final_window_validation.tsv（四率+决策）、vp1_veto_events.tsv、iteration_summary.tsv、raw_reference_qc.tsv、node_candidate_metrics.tsv、final_sweep/ 等。"));

children.push(pageBreak());

// ====================================================================
// 6. QC & analysis modules
// ====================================================================
children.push(h1("6  质量控制模块（线下）"));
children.push(para([
  { text: "audit_mapping_reference_qc.py（scripts/mapping_reference/）：", bold: true },
  { text: "对最终映射参考重跑流程自有 QC 规则——重建每条参考的比对形态，相对 local root（自身成员多数）与 group root（全组多数）重算私有 SNP 与 SNP 簇（101 nt 窗 >6 记 1 簇），核对 IUPAC 歧义码、N/gap、来源记录的原始 QC 标记与窗口验证决策。审计状态：PASS / REVIEW（local_SNP_clusters>1）/ UNASSESSABLE（singleton/small-group 来源）/ FAIL（非 ACGTN 碱基、N/gap>10%、来源被 SNP_CLUSTERS>1 标记）。" }
]));
children.push(para("2026-09-15 仓库清理后，流程一的原始代表 QC、代表性评估、探针覆盖度与基因注释等一次性分析模块已移出工作树（完整内容见 git 历史 v2.0.0 提交）；生产流程仅保留两条管线各自调用的脚本（scripts/common、scripts/probe_consensus、scripts/mapping_reference）。"));

children.push(pageBreak());

// ====================================================================
// 7. Parameter rationale
// ====================================================================
children.push(h1("7  关键参数设置原理（结果导向推导）"));
children.push(para("v2.0 的参数不是从机制正向设定，而是从下游目标反推：最终 panel 须能区分并还原基因组相似性 99% 的毒株。推导链的每一环都有实测锚点（校准脚本见 SCRIPT 版本记录，出版图见图 7-1）。"));
children.push(numbered([{ text: "浓缩系数实测（190 对标注变异株）：", bold: true }, { text: "1% 基因组距离 ⇒ 整段 VP1 约 1.0%（k≈1.00）、P 结构域约 1.8%（k_P 中位 1.82，极端 2.93）、ORF1 约 1.0%——衣壳判别信号确实集中于 VP1（尤其 P 域），为 VP1 否决锚定 VP1 区提供依据；" }]));
children.push(numbered([{ text: "R99 实测（构造真实 99.01% 姊妹株）：", bold: true }, { text: "150 bp 窗口可区分率 0.805（两参考竞争与 21 参考竞争一致）——单 SNV 窗靠 AS margin≥6 即可判定。⇒ T 必须低于 0.805；" }]));
children.push(numbered([{ text: "去重下限推导：", bold: true }, { text: "可分率=1−e^(−150d)，T=0.5 ⇔ 等效距离 0.46%——低于此距离的近重复（影子参考）自动合并。实测影子参考会把 99% 对可分率从 0.805 杀到 0.000，故影子必须清除；⇒ T=0.5 恰落在 (0.46% 去重线, 0.805 保留线) 安全带中点（1.6×/保守裕度）；" }]));
children.push(numbered([{ text: "VP1 否决块阈值：", bold: true }, { text: "金标准变异株对 VP1 连续差异块 ≥29 窗、99% 对 20 窗、近重复 ≤3 窗 ⇒ --veto-block-windows 15 同时覆盖 99% 对与变异株对（2×/5× 保护带）；" }]));
children.push(numbered([{ text: "选枝粒度：", bold: true }, { text: "变异株级 d95=5% 会把 99% 对装进同一枝（target≈1.0 永不触发拆分，结构信息不可恢复）⇒ 株级 d95=0.005 使 99% 对以独立单元入流程；" }]));
children.push(numbered([{ text: "伞形路由门槛三臂对照（GII.P16_GII.4）：", bold: true }, { text: "无门槛 102 条（恶化，50% retain）/ 尺寸门槛≥25：35 条（74%）/ 结构性门槛：37 条（62%，最小间距 0.64%）——结构性门槛以可行性判据（拆得出 ≥2×3 成员的子枝）取得同等稳定性且无经验常数，定为默认。" }]));
children.push(para([{ text: "图 7-1  阈值校准三面板：", bold: true }, { text: "A 组间/组内可分率分布与候选阈值；B 组间 rate 随最近参考距离上升（含泊松理论线，2006a@2.66%→0.739 为最难对）；C 阈值判定带（绿区 T≤0.636 两组均通过，0.80 在两基因型上均被证伪）。" }]));
children.push(new Paragraph({
  alignment: AlignmentType.CENTER,
  children: [new ImageRun({ type: "png",
    data: fs.readFileSync("result_paper/Fig_threshold_calibration.png"),
    transformation: { width: 624, height: 225 },
    altText: { title: "阈值校准", description: "Threshold calibration three panels", name: "calib" } })]
}));
children.push(para([{ text: "下游还原能力（混株模拟，150 bp、100×覆盖、0.3% 替换错误）：", bold: true }, { text: "panel 中两参考分开存在时，99% 相似双株可被同时还原——碱基错误率 0.000%、75 个株判别位点 100% 正确、嵌合率 0%，次要株比例下探 10%；若 panel 将两株合并为一条参考：50:50 时判别位点全部 N、90:10 时次要株 100% 被优势株遮蔽——还原精度由 panel 粒度决定，不由测序决定。" }]));

children.push(pageBreak());

// ====================================================================
// 8. Results summary
// ====================================================================
children.push(h1("8  结果汇总（v2.0 定案配方）"));
children.push(h2("8.1  GII 实测"));
children.push(table(
  ["指标", "GII.P16_GII.4（504 序列）", "GII.P17_GII.17（664 序列）", "GII.17 全部 7 组"],
  [
    ["最终参考数", "37", "19", "34"],
    ["全部真实序列", "✓（MEDOID/MEDOIDQC/RAW）", "✓", "✓ 10 MEDOID + 10 MEDOIDQC + 14 RAW"],
    ["参考间最小距离", "0.64%（99.36% 相似株共存）", "0.74%", "—"],
    ["影子 <0.5%", "0", "0", "0"],
    ["成员丢失", "0", "0", "0"],
    ["retain 比例", "62%", "63%", "小组全部收敛"],
    ["历史问题", "240→2 级联雪崩已杜绝", "Node276 坍缩已解决", "—"],
  ],
  [2600, 2300, 2300, 2160]
));
children.push(h2("8.2  标注变异株/clade 对照"));
children.push(para("GII.4（VP1 区最近参考）：20 个命名变异株中 14 个获得专属参考（VP1 距离 0.31–4.6%），未能专属的 6 个全部可归因于源数据缺失（其 accession 不在 rawdata 中）而非方法缺陷；63 条参考服务标签之外的株级多样性。GII.17（332 条标注基因组映射）：Kawasaki_308（191/209 集中）、Kawasaki_323（23/23 全中专属 medoid）、Romania（25/69 主参考+7 次级）各获独立代表——对照旧版全部挤在 Node276 一条参考上的坍缩。"));
children.push(h2("8.3  与旧版本对照"));
children.push(table(
  ["版本", "GII.P16_GII.4", "GII.P17_GII.17", "代表类型"],
  [
    ["旧生产版（T=0.3，ASR，级联合并）", "2 条（240→2 雪崩）", "3 条（Node276 坍缩）", "ASR 镶嵌体"],
    ["v1（T=0.5 + 批量 + 否决）", "10 条（伞形残留）", "3 条", "ASR"],
    ["v2.0（定案：medoid + 株级 + 结构门槛）", "37 条，最小间距 0.64%", "19 条，clade 分离", "全部真实序列"],
  ],
  [2800, 2200, 2200, 2160]
));

children.push(pageBreak());

// ====================================================================
// 9. Limitations
// ====================================================================
children.push(h1("9  局限性与注意事项"));
children.push(bullet("MEDOIDQC 参考需人工复核：其所在 clade 全部成员被 SNP 簇 QC 标记时，代表取最近被标记成员（保证谱系不丢），如 GII.P17_GII.17 主组 19 条中 10 条为 MEDOIDQC——该组 266/504 成员带 SNP 簇标记，QC 标准本身或需按组校准。"));
children.push(bullet("收尾振荡依赖兜底机制：伞形拆分与合并可在 21↔26 条间打转（3 轮极限环），由振荡检测（成员划分签名）+ 最优快照（-flagged 数最少）止损；终态仍可带少量 flagged 单元（本质是 VP1 有真实差异但 150 bp margin 不足的固有分辨率冲突）。"));
children.push(bullet("部分基因组序列：列覆盖 <50% 的窗口不可判（守卫生效），但 ORF1-only 部分序列的 VP1 否决盲区以 allowed_low_vp1_coverage 留痕放行；极短序列（如 422 bp 有效碱基）不产生任何窗口。"));
children.push(bullet("株级 panel 规模是变异株级的 2–3 倍（37 vs 16 条/组）：99% 分辨的固有价格，下游 EM 丰度估计的参考竞争随之增加；如仅需变异株级分辨率，可回退 d95=0.05 配方（其余参数不变）。"));
children.push(bullet("模拟理想化声明：还原精度 0.000% 基于无 indel、覆盖均匀、无质量波动的模拟；实测口径通常比理论低 20–30%（局部裁剪、成员-共识偏差），实际错误率会更高，但分开可还原/合并必丢失的定性结论不受影响。"));
children.push(bullet("10% N/gap 过滤无法捕获 IUPAC 歧义码（如 Y），审计中曾发现一例——汇总前建议补加字符集校验。"));

// ====================================================================
// 10. References
// ====================================================================
children.push(h1("10  主要参考文献与软件"));
children.push(bullet("Di Tommaso P, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology, 2017, 35(4):316–319."));
children.push(bullet("Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7. Molecular Biology and Evolution, 2013, 30(4):772–780."));
children.push(bullet("Nguyen LT, et al. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. MBE, 2015, 32(1):268–274."));
children.push(bullet("Minh BQ, et al. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. MBE, 2020, 37(5):1530–1534."));
children.push(bullet("Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nature Methods, 2012, 9(4):357–359."));
children.push(bullet("Price MN, et al. FastTree 2 — approximately maximum-likelihood trees for large alignments. PLoS ONE, 2010, 5(3):e9490."));
children.push(bullet("Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 2009, 25(11):1422–1423."));
children.push(bullet("Harris CR, et al. Array programming with NumPy. Nature, 2020, 585(7825):357–362."));
children.push(bullet("Nextstrain / Nextclade QC 逻辑（private mutations、SNP clusters）：nextstrain/ncov bioinformatics pipelines."));
children.push(spacer());
children.push(para({ text: "— 文档完 —", bold: false }, { alignment: AlignmentType.CENTER }));

// ====================================================================
// BUILD DOCUMENT
// ====================================================================
const doc = new Document({
  creator: "IPHnano noro pipeline",
  title: "诺如病毒代表性序列生成流水线 — 分析方法与参数设置说明",
  styles: {
    default: { document: { run: { font: FONT, size: 22 } } }, // 11pt
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 32, bold: true, font: FONT, color: "1F3864" },
        paragraph: { spacing: { before: 320, after: 160 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: FONT, color: "2E5496" },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 23, bold: true, font: FONT, color: "2E5496" },
        paragraph: { spacing: { before: 180, after: 100 }, outlineLevel: 2 } },
    ]
  },
  numbering: {
    config: [
      { reference: "bullets", levels: [
        { level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 600, hanging: 300 } } } },
        { level: 1, format: LevelFormat.BULLET, text: "◦", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 1080, hanging: 300 } } } },
      ] },
      { reference: "numbers", levels: [
        { level: 0, format: LevelFormat.DECIMAL, text: "%1.", alignment: AlignmentType.LEFT,
          style: { paragraph: { indent: { left: 600, hanging: 300 } } } },
      ] },
    ]
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 }
      }
    },
    headers: {
      default: new Header({ children: [new Paragraph({
        alignment: AlignmentType.RIGHT,
        children: [new TextRun({ text: "诺如病毒代表性序列生成流水线 · 分析方法说明", font: FONT, size: 18, color: "888888" })]
      })] })
    },
    footers: {
      default: new Footer({ children: [new Paragraph({
        alignment: AlignmentType.CENTER,
        children: [
          new TextRun({ text: "第 ", font: FONT, size: 18, color: "888888" }),
          new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 18, color: "888888" }),
          new TextRun({ text: " 页 / 共 ", font: FONT, size: 18, color: "888888" }),
          new TextRun({ children: [PageNumber.TOTAL_PAGES], font: FONT, size: 18, color: "888888" }),
          new TextRun({ text: " 页", font: FONT, size: 18, color: "888888" }),
        ]
      })] })
    },
    children
  }]
});

const OUT = "/Users/LuJ/Yunpan/NGS/pipline/IPHnano/noro_representative_seq_generate/诺如病毒代表性序列生成流水线_分析方法说明.docx";
Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(OUT, buf);
  console.log("WROTE " + OUT + "  (" + buf.length + " bytes)");
});
