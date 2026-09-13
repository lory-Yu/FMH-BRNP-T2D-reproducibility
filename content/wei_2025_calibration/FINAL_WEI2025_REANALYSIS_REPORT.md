# Wei 2025 Rb1 高/低转化组公开数据再分析：最终报告

## 结论先说

这项工作已经解决了“公开数据能否做**分组级复现**”的问题：答案是**可以**。50 名健康供者的 D1–D50、SRA、年龄、性别和 Rb1 High/Middle/Low 标签已经逐一对应；100 个 FASTQ 完整；DADA2 到 SILVA 注释和统计均已跑通。

它仍未解决“供者级**连续数值校准**”：公开材料没有每人的 Rb1 转化百分比、Compound K 生成量或 0/24 h 原始浓度。因此不能计算连续相关、回归校准误差，也不能声称预测了真实 CK 产量。

## 数据与流程审计

- 研究：Wei et al., 2025，50 名 18–31 岁健康供者；并非 T2D 队列。
- 分组：High 10、Middle 30、Low 10；主比较固定为 High vs Low。
- SRA：50 个供者对应 50 个唯一 Run、50 个唯一 BioSample，全部为 paired-end。
- FASTQ：100 个文件，共 3,473,618 对 reads；ENA MD5、gzip、双端 read 数及抽样 read ID 全部通过。
- DADA2：输入 3,473,618 对；过滤后 2,746,501；合并 2,567,893；非嵌合 1,916,623；最终 4,141 个细菌/古菌 ASV。
- 所有 50 个样本通过预设 QC；过滤保留率中位数 79.1%，过滤后合并率中位数 93.9%，输入到非嵌合保留率中位数 55.2%。
- 环境：R 4.5.3、DADA2 1.38.0、vegan 2.7.5；SILVA NR99 138.2，Zenodo DOI 10.5281/zenodo.14169026。
- 因桌面沙箱禁止 BiocParallel 本地端口，采用可复现的 SerialParam；这影响速度，不改变统计定义。

## High vs Low 主要结果

### α 多样性

- Observed_ASV：High 中位数 284.000，Low 224.000；Wilcoxon P=0.045；年龄/性别调整 P=0.012，q=0.028。
- Chao1：High 中位数 287.401，Low 231.724；Wilcoxon P=0.038；年龄/性别调整 P=0.014，q=0.028。
- Shannon：High 中位数 3.548，Low 3.246；Wilcoxon P=0.104；年龄/性别调整 P=0.053，q=0.070。
- Simpson：High 中位数 0.940，Low 0.904；Wilcoxon P=0.076；年龄/性别调整 P=0.264，q=0.264。

解释：高转化组的**丰富度**（Observed ASV、Chao1）较高并在调整后通过 BH；Shannon 和 Simpson 方向一致但证据不足。30,000 reads 稀释敏感性分析得到相同结论。

### β 多样性

- 年龄/性别调整 PERMANOVA：R²=0.066，P=0.131。
- 与原文方法一致的 ANOSIM：R=0.058，P=0.160。
- PERMDISP：P=0.092。

解释：本次冻结的 DADA2/SILVA 流程**没有复现原文显著的整体群落分离**。这不是流程失败，而是独立处理流程下的部分复现结果；原文使用 Majorbio 流程与 OTU/ANOSIM，方法差异可能贡献不一致。

### 原文预设候选菌

8 个原文候选中，SILVA 138.2 在相近分类层级可检测 6 个，6 个的点估计均为 High 较高。Eubacterium hallii group 与 UBA1819 未在冻结注释层级检出，不做替代归类。

- Anaerostipes：High−Low CLR=2.031；调整 P=0.052，候选集 q=0.070；Wilcoxon P=0.031，候选集 q=0.037。
- norank_f_Eubacterium_coprostanoligenes_group：High−Low CLR=3.139；调整 P=0.009，候选集 q=0.053；Wilcoxon P=0.006，候选集 q=0.032。
- Coprococcus：High−Low CLR=2.262；调整 P=0.026，候选集 q=0.070；Wilcoxon P=0.018，候选集 q=0.032。
- Barnesiella：High−Low CLR=2.241；调整 P=0.039，候选集 q=0.070；Wilcoxon P=0.019，候选集 q=0.032。
- norank_f_Oscillospiraceae：High−Low CLR=0.387；调整 P=0.661，候选集 q=0.661；Wilcoxon P=0.165，候选集 q=0.165。
- Family_XIII_AD3011_group：High−Low CLR=1.644；调整 P=0.058，候选集 q=0.070；Wilcoxon P=0.022，候选集 q=0.032。

解释：按原文方法的相对丰度 Wilcoxon，在 6 个可检测候选内 BH 后有 5 个 q<0.05；但年龄/性别调整 CLR 模型的候选集 BH 后为 0 个，且对全部 115 个属级特征校正后没有任何 q<0.05。最稳妥措辞是“**方向和部分候选得到复现，但协变量调整及全特征多重校正后证据减弱**”。

### 三等级趋势

把补充表中的“/”解释为论文所述 Middle 组后，完成了 50 人 Low→Middle→High 探索性趋势分析；共 115 个特征，无一通过全特征 BH-FDR。由于没有连续转化率，三等级不能替代连续剂量—反应校准。

## 与三个 T2D shotgun 队列的探索性桥接

4 个可直接映射的候选属分别在 Qin、Karlsson、MetaCardis 内部建模后合并。只有 Coprococcus 在三个队列均为 T2D 偏低，但随机效应结果为 β=-1.196（95% CI -4.753 至 2.361），P=0.285，BH q=0.285，I²=91.1%；证据不显著且异质性很高。

因此只能写：**Coprococcus 是连接健康人 Rb1 高转化表型与 T2D 菌群改变的候选线索。**不能写它已被证明负责 Rb1→CK，也不能写 T2D 导致 Rb1 转化下降。

## 解决了什么、没解决什么

### 已解决

1. 公开表型与 SRA 是否能逐样本对齐：能，50/50。
2. 能否从原始 16S 独立重跑：能，完整 DADA2/SILVA 流程已成功。
3. 能否复现 High/Low 菌群差异：能做严格检验；丰富度与部分预设候选得到支持，β 多样性和全属 FDR 未复现。
4. 能否建立可审计的分组级外部证据：能；所有输入、参数、版本、日志和校验值均留存。

### 仍未解决

1. 供者级连续 Rb1 转化率和 CK 生成量缺失，故不能完成连续校准。
2. 16S 只能给分类组成，不能直接确认 β-葡萄糖苷酶基因或 Rb1→Rd→F2→CK 反应通量。
3. Wei 队列为健康年轻人，不能验证 T2D 或二甲双胍对转化的影响。
4. 公共 SRA 未发现提取/PCR 阴性对照 FASTQ，污染控制只能列为局限。

## 论文可用的一句话

> Independent reprocessing of public 16S data from 50 healthy donors partially reproduced the microbiome stratification of high versus low Rb1 converters: richness and several prespecified taxa were higher in efficient converters, whereas overall beta-diversity separation and global genus-level FDR signals were not reproduced. These data support inter-individual microbiome-associated biotransformation phenotypes, but do not establish continuous CK production or T2D-specific impairment.

## 来源

- 论文：https://link.springer.com/article/10.1186/s13020-025-01190-2
- SRA：https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1268742
- SILVA DADA2 138.2：https://zenodo.org/records/14169026
