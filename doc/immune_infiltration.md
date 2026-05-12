# 免疫浸润分析

## 概述

免疫浸润分析旨在从 bulk 转录组表达矩阵中推断 22 种免疫细胞类型的相对丰度或富集分数。本工具基于 [TumorDecon](https://github.com/shawnstat/TumorDecon) 库，提供四种反卷积算法，无需安装 R 环境，完全基于 Python 生态。

## 依赖

| 依赖 | 用途 |
|------|------|
| `TumorDecon` | 核心反卷积引擎，内嵌 LM22 特征矩阵 |
| `scipy` | DeconRNASeq 的二次规划求解、统计检验 |
| `scikit-learn` | CIBERSORT 的 NuSVR 回归 |
| `pySingScore` | SingScore 秩和评分 |
| `pandas / numpy` | 数据处理 |
| `matplotlib / seaborn` | 可视化 |

## 输入数据要求

### 表达矩阵

- 行 = 基因（Hugo_Symbol），列 = 样本
- 支持 microarray 和 RNA-seq 数据
- 基因名列（`Gene`、`Symbol`、`SYMBOL` 等）会被自动识别并设为索引

### 自动预处理

免疫分析模式下，系统会自动执行以下预处理：

1. **基因名列识别** → 设为 DataFrame 索引
2. **TPM 归一化** → 每列除以列总和 × 10⁶（模拟 TPM 尺度）
3. **跳过 log2 转换** → 免疫模式不执行 log2(x+1)
4. **基因交集** → 与 LM22 的 547 个基因取交集

## 四种反卷积方法

通过 `config.yaml` 中的 `immune_method` 控制，默认为 `"DeconRNASeq"`。

### DeconRNASeq（推荐）

| 项目 | 说明 |
|------|------|
| 原理 | 二次规划解 `min ||S·F - M||²`，约束 `F ≥ 0, ΣF = 1` |
| 输出 | 细胞比例（0~1），所有细胞类型之和 = 1 |
| 优势 | 输出直观，适合堆叠柱状图 |
| 劣势 | 对基因覆盖敏感，LM22 基因缺失较多时解不稳定 |

### CIBERSORT

| 项目 | 说明 |
|------|------|
| 原理 | Nu-支持向量回归 + 网格搜索最优参数 |
| 输出 | 细胞比例（0~1），所有细胞类型之和 = 1 |
| 优势 | 业界广泛使用，文献接受度高 |
| 劣势 | 每个样本约需 10~30 秒，大样本集较慢 |

### ssGSEA

| 项目 | 说明 |
|------|------|
| 原理 | 对每种细胞类型的上调基因集在全部基因中算富集分数（ES） |
| 输出 | 富集分数（无上界），各细胞类型之间不互斥 |
| 优势 | 对平台差异（microarray vs RNA-seq）鲁棒，基因缺失容忍度高 |
| 劣势 | 输出不是比例，堆叠柱状图 Y 轴语义不准确（系统会给出 warning） |

### SingScore

| 项目 | 说明 |
|------|------|
| 原理 | 双向/单向秩和评分 |
| 输出 | 富集分数（无上下界） |
| 优势 | 计算快，适合大规模数据 |
| 劣势 | 与 ssGSEA 类似，非比例输出 |

### 方法选择建议

| 场景 | 推荐方法 |
|------|----------|
| 数据质量好、基因覆盖全 | DeconRNASeq |
| microarray 与 RNA-seq 混合分析 | ssGSEA |
| 跨数据集横向比较 | ssGSEA |
| 关注某种细胞类型的绝对变化 | ssGSEA 或 SingScore |
| 需要直观的比例堆叠图 | DeconRNASeq 或 CIBERSORT |
| 文献投稿（需要 Cite CIBERSORT 原文） | CIBERSORT |

## 配置项

```yaml
# config.yaml 中免疫相关的配置
immune_method: "DeconRNASeq"     # 反卷积算法：DeconRNASeq / CIBERSORT / ssGSEA / SingScore
plot_data_warning: true          # 画图数据合理性 warning（ssGSEA/SingScore 时建议保持开启）
```

## 输出文件

### 数据文件

| 文件 | 路径 | 说明 |
|------|------|------|
| 免疫浸润结果 | `data/{GSE_ID}/pkl/{GSE_ID}_immune_summary.pkl` | 样本 × 22 种免疫细胞的完整数据 |
| 免疫浸润 CSV | `res/csv/{GSE_ID}_immune_summary.csv` | 同上，CSV 格式 |
| 相关性数据 | `res/csv/{GSE_ID}_{gene}_immune_correlation.csv` | 基因表达与免疫细胞丰度的 Spearman R 值 |

### 图片输出

| 图 | 路径 | 说明 |
|------|------|------|
| 堆叠柱状图 | `res/figures/immune/{GSE_ID}_{gene}_immune_stacked_bar.png` | 每个样本一根柱子，22 种细胞按颜色堆叠 |
| 箱线图（全量） | `res/figures/immune/{GSE_ID}_{gene}_immune_boxplot.png` | 每种细胞类型对照 vs 实验组，Mann-Whitney U 检验 + Bonferroni 校正 |
| 箱线图（过滤） | `res/figures/immune/{GSE_ID}_{gene}_immune_boxplot_filtered.png` | 同上，剔除两组中位数均 < 0.001 的低丰度细胞 |
| 相关性热图 | `res/figures/immune/{GSE_ID}_{gene}_immune_heatmap.png` | 目标基因表达 × 免疫细胞丰度的 Spearman 相关性 |

## 图表说明

### 堆叠柱状图

- 纵轴为 `Cell Fraction (%)`，假设数据为比例（行和 ≈ 1）
- 使用 Tab20 / Tab20b 配色区分 22 种细胞类型
- 按 Control → Experiment 分组排列，组间虚线分隔
- 若使用 ssGSEA 或 SingScore（非比例数据），系统会检测到行和偏离 1 并输出 warning，图表仍会生成但 Y 轴百分比含义不准确

### 箱线图

- 每种免疫细胞类型一个子图，5 列网格排列
- Control vs Experiment 两两对比
- Mann-Whitney U 检验 + Bonferroni 多重检验校正
- `*` p < 0.05、`**` p < 0.01、`***` p < 0.001
- 细胞按 adjusted p-value 升序排列（显著的排在前面）
- 过滤版会剔除两组中位数均 < 0.001 的细胞类型

### 相关性热图

- 横轴 = 22 种免疫细胞（按平均 |R| 降序排列）
- 纵轴 = 目标基因（支持 `tar_gene` 单基因或 `multi_gene` 多基因）
- Spearman 秩相关，中心为 0，RdBu_r 配色
- 格内标注 R 值 + 显著性星号

## 22 种免疫细胞类型

| 类别 | 细胞类型 |
|------|----------|
| B 细胞 | B cells naive, B cells memory, Plasma cells |
| T 细胞 | T cells CD8, T cells CD4 naive, T cells CD4 memory resting, T cells CD4 memory activated, T cells follicular helper, T cells regulatory (Tregs), T cells gamma delta |
| NK 细胞 | NK cells resting, NK cells activated |
| 髓系 | Monocytes, Macrophages M0, Macrophages M1, Macrophages M2 |
| 树突状细胞 | Dendritic cells resting, Dendritic cells activated |
| 肥大细胞 | Mast cells resting, Mast cells activated |
| 粒细胞 | Eosinophils, Neutrophils |

## 常见问题

### Q: 使用 ssGSEA 或 SingScore 时，堆叠柱状图 warning 是什么意思？

ssGSEA 和 SingScore 输出的是富集分数（ES），不是细胞比例。每个样本所有细胞类型分数之和通常不等于 1，可能远超 1。堆叠柱状图的 Y 轴标注为 `Cell Fraction (%)`，此时这个标注在数学上不准确，但图表仍可展示各组间的趋势差异。如果你不希望看到该 warning，将 `config.yaml` 中 `plot_data_warning` 设为 `false`。

### Q: 如何手动输入目标基因做相关性热图？

免疫浸润本身不需要基因输入（只依赖表达矩阵）。相关性热图需要的基因通过以下方式指定：

```yaml
tar_gene: "APEX1"                              # 单个基因
multi_gene: "APEX1,TP53,EGFR"                  # 多个基因（逗号分隔）
multi_gene: "/path/to/gene_list.txt"           # 从文件读取（每行一个基因名）
```

`tar_gene` 与 `multi_gene` 互斥，不能同时存在。

### Q: 不同 GSE_ID 的免疫浸润结果能直接比较吗？

不建议直接横向比较绝对值。不同数据集样本人群、组织类型、平台不同，绝对分数不可比。每个数据集内部做 Control vs Experiment 对比是有意义的。ssGSEA 对跨数据集的相对趋势对比稍好于比例方法，但仍需谨慎。

### Q: 转录组数据做免疫浸润有什么局限？

LM22 特征矩阵基于 mRNA 表达构建，无法捕获翻译后修饰（磷酸化、乙酰化等）引起的变化。如果免疫功能变化主要由蛋白水平的调控驱动，可能出现假阴性。详情见[该项目的局限性提示]或相关文献讨论。

### Q: DeconRNASeq 速度很慢怎么办？

DeconRNASeq 对每个样本解二次规划，样本多时较慢。可以改用 SingScore（最快）或 ssGSEA（较快）。
