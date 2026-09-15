# 病毒共识序列生成流水线：方法与技术路线

## 1. 研究背景与目的

本研究开发了一套自动化生物信息学流水线，用于从大规模病毒基因组序列数据中生成高质量的代表性共识序列。该方法整合了序列质量控制、聚类分析、多序列比对、系统发育分析和祖先序列重建等技术，通过迭代验证策略确保最终输出的共识序列之间具有足够的遗传多样性（序列相似性≤95%）。

## 2. 材料与方法

### 2.1 软件工具

本研究使用的生物信息学软件及其功能描述如下：

#### 2.1.1 工作流管理系统

**Nextflow** 是一种基于数据流编程模型的工作流管理系统，能够实现生物信息学流水线的可重复性执行和并行计算。该工具支持跨平台部署，可在本地服务器和云计算环境中运行。

#### 2.1.2 多序列比对工具

**MAFFT**（Multiple Alignment using Fast Fourier Transform）是一种高效的多序列比对软件。本流水线对每个多序列 cluster 去除既有 gap 后重新执行 MAFFT `--auto` 比对，避免将下载序列中已有的 padding 误认为同源位点比对结果。比对后删除所有 A/C/G/T 覆盖度低于 50% 的 alignment column；N 和 gap 不计为有效覆盖，因此低覆盖内部插入列和端部低覆盖列都会被屏蔽。

#### 2.1.3 系统发育分析工具

**IQ-TREE** 是一种基于最大似然法（Maximum Likelihood, ML）的系统发育树构建软件。本研究使用IQ-TREE进行两项主要任务：（1）构建系统发育树以揭示序列间的进化关系；（2）进行祖先序列重建（Ancestral Sequence Reconstruction, ASR）以推断内部节点的序列状态。

系统发育分析采用Jukes-Cantor（JC）核酸替换模型，该模型假设四种核苷酸之间的替换概率相等，是核酸序列分析中最基础的替换模型，适用于近缘序列的比较分析。祖先序列重建过程中，设置后验概率阈值为0.8（80%），即只有当某位点某碱基的后验概率达到80%以上时，才将该碱基确定为祖先状态，否则标记为不确定位点。

#### 2.1.4 序列聚类工具

**vclust** 是一种专门针对病毒序列设计的聚类工具，基于平均核苷酸一致性（Average Nucleotide Identity, ANI）进行序列相似性计算和聚类分析。该工具采用三步策略：首先通过预过滤快速筛选候选序列对，然后计算序列对之间的ANI值，最后基于ANI阈值进行层次聚类。

#### 2.1.5 编程语言与生物信息学库

流水线的自定义脚本使用**Python**编程语言编写，主要依赖以下库：**Biopython**用于序列文件解析、比对结果读取和系统发育树操作；**Pandas**用于表格数据处理；**NumPy**用于数值计算和矩阵运算。

### 2.2 分析流程

本流水线包含九个主要处理步骤，各步骤之间通过数据流自动连接。

#### 2.2.1 步骤一：序列质量过滤

**目的**：移除低质量序列，确保后续分析的可靠性。

**方法**：对输入的FASTA格式序列文件进行逐条扫描，计算每条序列中未知碱基（N）和比对缺口（gap，以"-"表示）的比例。设定质量阈值为10%，即当序列中N碱基和gap的总比例超过10%时，该序列被判定为低质量序列并从数据集中移除。

**计算公式**：

$$Quality\ Score = \frac{N_{count} + Gap_{count}}{L_{sequence}} \times 100\%$$

其中，$N_{count}$为序列中N碱基的数量，$Gap_{count}$为缺口数量，$L_{sequence}$为序列总长度。当$Quality\ Score > 10\%$时，移除该序列。

#### 2.2.2 步骤二：序列聚类分析

**目的**：将遗传相似性较高的序列归入同一聚类，为后续分组分析奠定基础。

**方法**：采用基于ANI的层次聚类方法。首先计算所有序列对之间的ANI值，ANI定义为两条序列之间相同核苷酸位点占可比对位点总数的比例。随后仅保留 VP1 genotype 相同的序列对（如 GI.3-GI.3、GII.4-GII.4、GIX.1-GIX.1），去除跨 genotype 或无法解析 genotype 的配对。最后采用完全连锁法（Complete Linkage）进行层次聚类，该方法以聚类间最大距离作为聚类间距离度量，降低链式聚类导致的跨型别合并风险。

**参数设置**：
- ANI聚类阈值：0.90（90%）
- query coverage阈值：0.65
- reference coverage阈值：0.65
- 序列长度比例阈值：0.65
- 聚类算法：完全连锁法

**生物学意义**：90%的ANI阈值能够在不同长度基因组片段之间保持适度容忍度；qcov、rcov和长度比例阈值用于避免短局部高相似区域导致错误合并；VP1 genotype 硬分层用于防止不同衣壳型别在低覆盖局部比对下被并入同一 cluster。

#### 2.2.3 步骤三：聚类分割

**目的**：将聚类结果分割为独立的序列文件，便于后续并行处理。

**方法**：根据聚类分析结果，将属于同一聚类的序列提取并保存为独立的FASTA文件。同时从序列标识符中提取基因型信息用于文件命名，便于结果追溯和管理。

#### 2.2.4 步骤四：序列方向统一

**目的**：在多序列比对前，将同一聚类内所有序列统一到同一条链方向，避免反向互补链存储的序列被强行比对到错误方向，产生大量 gap 和错误的相似性。

**方法**：公共数据库下载的序列并非都以正义链存储（例如 `GII.P31_GII.4_KX158285_2015`、`GII.P7_GII.6_KX158282_2015` 等以负链形式下载）。对每个聚类，选取 A/C/G/T 有效碱基数最多的序列作为方向参考，对其余每条序列分别计算与参考在正义链和反向互补链上的 15-mer 重叠数；当反向互补链重叠显著高于正义链且达到最低重叠阈值时，将该序列反向互补后再用于比对。仅含单条序列的聚类无需处理。该步骤仅改变链方向，不改变正向存储序列的内容。

#### 2.2.5 步骤五：多序列比对

**目的**：对每个聚类内的序列进行多序列比对，为系统发育分析提供输入数据。

**方法**：对于多序列聚类，先移除输入序列中已有 gap，再使用 MAFFT `--auto` 进行多序列比对。对于仅含单条序列的聚类，直接保留原序列，无需进行比对操作。比对完成后，对所有 alignment column 进行覆盖度过滤：若该列中 A/C/G/T 的覆盖比例低于 50%，则删除该列；N 和 gap 不作为有效覆盖，内部低覆盖插入列不进入后续建树和 consensus。

**算法特点**：
- 对每个多序列 cluster 进行真实 MAFFT 重比对
- 通过 `--auto` 自动选择适合该 cluster 的比对策略
- 对全 alignment 低覆盖列进行过滤，降低内部插入列导致的人工拉长 consensus 风险

#### 2.2.6 步骤六：共识序列生成

**目的**：基于系统发育分析生成代表性共识序列。

**方法**：本步骤采用创新的基于系统发育树节点分析的共识序列生成策略，具体包括以下子步骤：

**（1）系统发育树构建**

使用最大似然法构建系统发育树。采用Jukes-Cantor核酸替换模型，该模型的替换速率矩阵假设所有核苷酸之间的替换概率相等，替换速率为$\mu$，模型可表示为：

$$Q = \begin{pmatrix} -3\mu & \mu & \mu & \mu \\ \mu & -3\mu & \mu & \mu \\ \mu & \mu & -3\mu & \mu \\ \mu & \mu & \mu & -3\mu \end{pmatrix}$$

**（2）祖先序列重建**

利用最大似然祖先序列重建方法推断系统发育树内部节点的序列状态。对于每个内部节点的每个位点，计算四种核苷酸的后验概率，选择后验概率最高的核苷酸作为该位点的祖先状态。设置后验概率阈值为80%，当最高后验概率低于该阈值时，该位点标记为不确定。

**（3）树根定位**

采用中点定根法（Midpoint Rooting）确定树的根节点位置。该方法首先找到树中距离最远的两个末端节点，然后在连接这两个节点的路径中点处设置根节点。

**（4）最优节点选择**

遍历系统发育树的所有节点（包括内部节点和末端节点），对每个节点计算其包含的所有后代序列之间的最低成对相似性。基于以下两个条件选择最优节点：

- 条件一：节点内部序列的最低成对相似性 ≥ 相似性阈值（95%）
- 条件二：该节点的父节点内部序列最低相似性 < 相似性阈值（95%）

满足上述条件的节点即为"相似性边界节点"，代表了在系统发育树上从高相似性区域向低相似性区域的转变点。

**（5）共识序列输出**

对于满足条件的内部节点，输出其重建的祖先序列作为该分支的代表性共识序列；对于末端节点，保留原始序列作为输出。

**序列相似性计算公式**：

$$Similarity(\%) = \frac{M}{N_{valid}} \times 100$$

其中，$M$为两条序列在相同位置具有相同核苷酸（均非gap）的位点数，$N_{valid}$为两条序列在相同位置均不是gap的位点总数。

#### 2.2.7 步骤七：增强共识序列验证

**目的**：通过迭代验证确保最终输出的共识序列满足多样性要求。

**方法**：对步骤五输出的共识序列进行成对相似性验证。计算所有序列对之间的相似性，找出最高相似性值。若最高相似性超过95%阈值，则判定当前结果不满足多样性要求，需要进行进一步处理。

**迭代优化流程**：

1. 计算当前共识序列集合的最高成对相似性
2. 判断是否满足阈值条件：
   - 若最高相似性 ≤ 95%，验证通过，输出结果
   - 若最高相似性 > 95%，进入迭代优化
3. 对超过阈值的序列重新进行系统发育分析
4. 应用节点分析策略生成新的共识序列
5. 返回步骤1，重复验证
6. 当满足条件或达到最大迭代次数（10次）时终止

**终止条件**：
- 条件一：所有序列对的相似性均 ≤ 95%
- 条件二：迭代次数达到上限（10次）

#### 2.2.8 步骤八：结果收集

**目的**：合并所有聚类产生的最终共识序列。

**方法**：将各聚类经过验证的最终共识序列合并为单一FASTA文件，同时生成流水线运行摘要报告，包括处理的序列数量、生成的共识序列数量等统计信息。

#### 2.2.9 步骤九：最终质量验证

**目的**：对合并后的完整结果进行最终质量检验。

**方法**：对所有最终共识序列进行多序列比对，计算所有序列对的相似性，验证是否所有成对相似性均不超过95%阈值。生成最终验证报告，包括最高相似性值、总序列数、成对比较总数等信息。

### 2.3 关键参数说明

#### 2.3.1 核心参数

| 参数名称 | 设定值 | 生物学意义 |
|----------|--------|------------|
| 序列质量阈值 | 10% | 控制输入数据质量，移除含有过多未知碱基或缺口的低质量序列 |
| ANI聚类阈值 | 90% | 区分不同进化分支的边界，同一聚类内的序列应具有较高的遗传相似性 |
| query/reference coverage阈值 | 65% / 65% | 避免短局部高相似比对导致跨型别或远缘序列错误合并 |
| 序列长度比例阈值 | 65% | 允许5000-6000 bp片段与约7200 bp近完整基因组共同分析，同时排除过短局部重叠 |
| 聚类算法 | 完全连锁法 | 确保聚类内任意两序列相似性均达到阈值，避免链式聚类问题 |
| VP1 genotype过滤 | 相同型别内聚类 | 避免不同VP1型别因局部ANI较高被合并为mixed genotype cluster |
| 相似性阈值 | 95% | 判定序列是否足够相似可合并的边界，同时确保最终结果的多样性 |
| 后验概率阈值 | 80% | 祖先序列重建的置信度要求，保证重建结果的可靠性 |
| 最大迭代次数 | 10次 | 防止算法陷入无限循环，实际应用中通常3-5次即可收敛 |

#### 2.3.2 参数选择依据

**ANI阈值90%的选择依据**：病毒种内变异通常在85-95%的ANI范围内，选择90%作为聚类阈值可以有效区分不同的进化分支，同时将亲缘关系较近的变异株归为同一组进行联合分析。

**相似性阈值95%的选择依据**：95%的序列相似性通常被认为是同一病毒株系或血清型内部变异的上限。设置此阈值可确保最终输出的共识序列代表不同的病毒变异类型，具有足够的遗传多样性。

**后验概率阈值80%的选择依据**：在分子进化分析中，80%的后验概率被广泛接受为置信度的合理下限。低于此阈值的祖先状态推断存在较大不确定性，不宜直接采用。

**JC模型的选择依据**：Jukes-Cantor模型是最简单的核酸替换模型，假设条件最为宽松。对于病毒基因组这类高度保守的近缘序列比较，简单模型通常能够提供足够准确的结果，同时计算效率较高。如需更精确的分析，可考虑采用更复杂的模型如GTR（广义时间可逆模型）。

## 3. 算法原理

### 3.1 基于系统发育的共识序列生成策略

传统的共识序列生成方法通常采用简单多数投票或加权投票策略，即对比对后的每个位点统计各核苷酸的出现频率，选择频率最高者作为共识碱基。然而，这种方法忽略了序列之间的进化关系，可能导致生成的共识序列不能准确代表序列集合的遗传特征。

本研究提出的方法创新性地将系统发育信息整合到共识序列生成过程中。核心思想是：在系统发育树上寻找"相似性边界"，即在树的层次结构中找到子节点满足高相似性条件而父节点不满足的位置。这些边界节点代表了遗传多样性的自然分界点，其对应的祖先序列能够较好地代表该分支的整体遗传特征。

### 3.2 最优节点选择的数学描述

设系统发育树为$T$，包含节点集合$V$。对于任意节点$v \in V$，定义其后代序列集合为$D(v)$。节点$v$的内部相似性定义为：

$$S_{min}(v) = \min_{i,j \in D(v), i \neq j} Similarity(seq_i, seq_j)$$

即节点$v$包含的所有后代序列之间的最低成对相似性。

对于非根节点$v$，设其父节点为$parent(v)$。节点$v$被选为最优节点的条件为：

$$S_{min}(v) \geq \theta \quad \text{且} \quad S_{min}(parent(v)) < \theta$$

其中$\theta$为相似性阈值（本研究设为95%）。

### 3.3 迭代验证的收敛性

迭代验证过程旨在确保最终输出满足相似性约束。设第$k$次迭代后的共识序列集合为$C_k$，最高成对相似性为$S_{max}(C_k)$。迭代过程可形式化为：

$$C_{k+1} = f(C_k) \quad \text{当} \quad S_{max}(C_k) > \theta$$

其中$f$为系统发育分析和节点选择函数。

由于每次迭代都会对高相似性序列进行进一步分割，序列数量单调不增或序列间差异单调增加，因此算法具有收敛性。实践中，通常在3-5次迭代内即可达到收敛。

## 4. 结果输出

### 4.1 输出文件说明

流水线产生的主要输出文件包括：

| 输出内容 | 说明 |
|----------|------|
| 过滤后序列 | 通过质量控制的序列集合 |
| 聚类信息 | 序列到聚类的映射关系及统计信息 |
| 比对结果 | 各聚类的多序列比对结果 |
| 初始共识序列 | 系统发育分析产生的初步共识序列 |
| 最终共识序列 | 经迭代验证后的最终代表性序列 |
| 验证报告 | 最终结果的质量验证报告 |

### 4.2 输出序列命名规则

最终共识序列采用统一的命名格式，包含以下信息：基因型标识、聚类编号、处理轮次和节点编号。这种命名方式便于追溯序列来源和处理过程。

## 5. 方法学特点

### 5.1 创新点

1. **进化信息整合**：不同于传统的基于频率统计的共识序列生成方法，本方法充分利用系统发育信息，通过祖先序列重建获得更具生物学意义的代表性序列。

2. **自适应边界检测**：通过在系统发育树上自动识别"相似性边界"，实现了对序列多样性的自适应分割，无需预先指定分组数量。

3. **迭代质量保证**：通过迭代验证机制确保输出结果满足预设的多样性要求，避免产生过于相似的冗余序列。

4. **模块化设计**：流水线各步骤相对独立，参数可灵活调整，便于根据具体研究需求进行定制。

### 5.2 适用范围

本方法适用于以下应用场景：
- 病毒基因组代表性序列库的构建
- 疫苗候选株的筛选
- 分子流行病学研究中的序列去冗余
- 病毒进化分析的参考序列集准备

### 5.3 局限性

1. 对于高度分化的序列集（如不同物种），90%的ANI阈值可能需要调整。
2. JC模型假设较为简单，对于存在明显碱基组成偏好的序列可能不够准确。
3. 计算资源需求随序列数量增加而显著增长，大规模数据集需要较长处理时间。

## 6. 参考文献引用格式

### 6.1 软件工具引用

**工作流管理系统**
- Di Tommaso P, Chatzou M, Floden EW, et al. Nextflow enables reproducible computational workflows. Nature Biotechnology, 2017, 35(4): 316-319.

**多序列比对**
- Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution, 2013, 30(4): 772-780.

**系统发育分析**
- Nguyen LT, Schmidt HA, von Haeseler A, et al. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 2015, 32(1): 268-274.
- Minh BQ, Schmidt HA, Chernomor O, et al. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 2020, 37(5): 1530-1534.

**序列聚类**
- vclust: 病毒序列聚类工具. https://github.com/refresh-bio/vclust

**编程库**
- Cock PJ, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 2009, 25(11): 1422-1423.
- Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature, 2020, 585(7825): 357-362.

### 6.2 方法描述示例（中文）

> 本研究采用自主开发的自动化流水线生成病毒代表性共识序列。首先对输入序列进行质量过滤，移除未知碱基和缺口比例超过10%的低质量序列。随后按相同RdRp genotype和VP1 genotype对序列进行分组，以避免不同基因型或重组型序列被错误合并。在比对前先统一序列方向：对每个分组以A/C/G/T有效碱基数最多的序列为方向参考，比较其余序列正义链与反向互补链的15-mer重叠，将反向存储的序列反向互补到与参考一致的链上。对各分组移除既有gap后使用MAFFT `--auto` 重新进行多序列比对；比对后删除所有A/C/G/T覆盖度低于50%的列，其中N和gap均不计为有效覆盖，从而去除低覆盖内部插入列和端部低覆盖区域。最终代表序列写出前移除所有由比对引入的gap，避免低覆盖插入列被填入共识序列。系统发育树构建和祖先序列重建使用IQ-TREE软件完成，选用Jukes-Cantor核酸替换模型，祖先状态后验概率阈值设为80%。
>
> 共识序列的选择基于系统发育树节点分析策略：对中点定根后的系统发育树遍历所有节点，计算各节点内部序列的最低成对相似性，选择满足以下条件的节点——节点内部最低相似性≥95%且父节点内部最低相似性<95%——作为最优节点，输出其祖先序列作为该分支的代表。通过迭代验证确保最终所有共识序列之间的相似性均不超过95%。

### 6.3 方法描述示例（英文）

> Representative consensus sequences were generated using an automated bioinformatics pipeline. Input sequences were first filtered to remove low-quality sequences with more than 10% ambiguous bases and gaps combined. Sequences were grouped by matched RdRp and VP1 genotypes to avoid erroneous merging of distinct genotypes or recombinant lineages. Within each group, sequences were oriented to a single strand before alignment: every record was compared, by 15-mer overlap in the forward versus reverse-complement direction, against the group's longest A/C/G/T sequence, and any record stored on the opposite strand was reverse-complemented. For each multi-sequence group, pre-existing gaps were removed and sequences were realigned using MAFFT `--auto`; all alignment columns with <50% A/C/G/T coverage were removed, with both `N` and gap characters excluded from valid coverage. This masks low-occupancy internal insertion columns as well as low-coverage terminal regions. Before final representative sequences were written, all alignment-induced gaps were removed to prevent low-coverage insertion columns from being incorporated into consensus sequences. Phylogenetic tree construction and ancestral sequence reconstruction were performed using IQ-TREE with the Jukes-Cantor nucleotide substitution model and a posterior probability threshold of 80% for ancestral state assignment.
>
> Consensus sequence selection was based on a phylogenetic node analysis strategy. After midpoint rooting, all nodes in the phylogenetic tree were traversed to calculate the minimum pairwise similarity among descendant sequences. Optimal nodes were selected based on two criteria: (1) minimum internal similarity ≥ 95%, and (2) parent node minimum internal similarity < 95%. Ancestral sequences of optimal internal nodes were output as representative consensus sequences for each branch. Iterative validation ensured that all final consensus sequences had pairwise similarities not exceeding 95%.

### 6.4 关键术语中英对照

| 中文术语 | 英文术语 | 缩写 |
|----------|----------|------|
| 平均核苷酸一致性 | Average Nucleotide Identity | ANI |
| 多序列比对 | Multiple Sequence Alignment | MSA |
| 最大似然法 | Maximum Likelihood | ML |
| 祖先序列重建 | Ancestral Sequence Reconstruction | ASR |
| 系统发育树 | Phylogenetic Tree | - |
| 中点定根 | Midpoint Rooting | - |
| 层次聚类 | Hierarchical Clustering | - |
| 完全连锁法 | Complete Linkage | - |
| 后验概率 | Posterior Probability | PP |
| 核酸替换模型 | Nucleotide Substitution Model | - |
| 共识序列 | Consensus Sequence | - |
| 成对相似性 | Pairwise Similarity | - |

---

*本文档用于辅助科研论文方法部分撰写。文中描述的参数值为本研究的默认设置，可根据具体研究对象和目的进行调整。*
