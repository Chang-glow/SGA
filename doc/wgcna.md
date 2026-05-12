# WGCNA 数据准备

## 概述

WGCNA 模式（`analysis_mode: "wgcna"`）为 R 语言的 WGCNA（Weighted Gene Co-expression Network Analysis）分析准备输入数据，并调用 [PyWGCNA](https://github.com/shawnstat/PyWGCNA) 自动执行完整的共表达网络分析流水线。

整个流程包括：基因筛选 → 表达矩阵输出 → 性状矩阵输出 → PyWGCNA 模块检测 → 自动绘图。

## 输入数据

- **表达矩阵**：行 = 基因，列 = 样本（来自阶段 1 数据包）
- **元数据**：`meta` / `meta_full`，用于样本分组
- **不做 log2 转换**（WGCNA 内部自行处理标准化）

## 基因筛选（目标 ~3000 个）

基因数硬编码为 3000，选取策略分两级：

### 优先级 1：diff / hilo 缓存

若 `debug: false` 且存在以下缓存文件之一：
- `data/{GSE_ID}/pkl/{GSE_ID}_differential_summary.pkl`
- `data/{GSE_ID}/pkl/{GSE_ID}_highlow_summary.pkl`

则读取缓存，按 `padj` 升序排序，去重后取 top 3000。

### 优先级 2：表达方差

若无缓存可用，按基因行方差（`df.var(axis=1)`）降序取 top 3000。

### 兜底

若目标基因 `tar_gene` 被过滤掉，会强制加回。

## 输出文件

### 核心输出

| 文件 | 路径 | 说明 |
|------|------|------|
| 表达矩阵 | `res/csv/{GSE_ID}_{GENE}_datExpr.csv` | 样本 × 基因（筛选后），带 Sample 索引 |
| 性状矩阵 | `res/csv/{GSE_ID}_{GENE}_datTraits.csv` | 样本 × 性状（Group + 目标基因表达），带 Sample 索引 |

### 辅助输出

| 文件 | 说明 |
|------|------|
| `wgcna_module_genes.csv` | 每个基因的模块分配（moduleColors, moduleLabels, signed kME） |
| `wgcna_hub_genes.csv` | 每个模块的 top 10 hub 基因（按 kME） |

### datTraits.csv 列说明

| 列 | 说明 |
|----|------|
| `Group` | Control 或 Experiment（按 `exp_type` 命名） |
| `{GENE}_exp` | 目标基因在各样本中的表达值 |

样本按 Control 优先排列。

## PyWGCNA 流水线

```python
wgcna = WGCNA(
    name=f"{gse_id}_{tar_gene}",
    species=organism,
    geneExpPath=expr_csv_path,
    sampleInfo=datTraits,
    save=True,
    outputPath=res/figures/wgcna/,
)
```

自动执行以下步骤：

| 步骤 | 方法 | 说明 |
|------|------|------|
| 1 | `runWGCNA()` | 预处理 → 软阈值筛选 → TOM → 层次聚类 → 模块检测 |
| 2 | 自定义热图 | 模块 eigengene × 性状 Pearson 相关热图（含显著性星星） |
| 3 | `CalculateSignedKME()` | 计算每个基因在各模块中的 signed kME |
| 4 | `plotModuleEigenGene()` | 每个模块的 eigengene 条形图 |
| 5 | `functional_enrichment_analysis()` | 每个模块的 GO 富集点图 |
| 6 | `saveWGCNA()` | 持久化所有结果 |

所有图表在**阶段 2 内由 PyWGCNA 完成**，阶段 3 自动跳过（无需单独绘图步骤）。

## 配置项

| 配置 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `analysis_mode` | string | diff | 设为 `"wgcna"` 启动 |
| `tar_gene` | string | — | 目标基因（用于 datTraits 和文件名） |
| `organism` | string | human | 物种（影响 GO 富集） |
| `exp_type` | string | Fibrosis | datTraits 中实验组的标签名 |
| `debug` | bool | false | true 时跳过 diff/hilo 缓存，只用方差筛选 |

## 图表

所有图表输出至 `res/figures/wgcna/figures/`：

| 图 | 说明 |
|----|------|
| 模块-性状热图 | 各模块 eigengene 与 Group / {GENE}_exp 的 Pearson 相关性，含 `*/**/***` 标注 |
| 模块 eigengene 条形图 | 每个模块一张，展示 eigengene 在各样本中的值 |
| GO 富集点图 | 每个模块一张，气泡图展示 GO 通路富集 |

## 常见问题

### Q: WGCNA 模式可以直接在 R 中运行吗？

WGCNA 模式生成两个 CSV 文件（`datExpr.csv`、`datTraits.csv`）可以直接导入 R 的 WGCNA 包。如果不需要 PyWGCNA 自动分析，取这两个文件后即可。

### Q: 为什么基因数固定为 3000？

WGCNA 的计算复杂度与基因数的平方相关，3000 是实践中平衡覆盖度和计算时间的常用值。需要调整可修改 `modules/strategies/wgcna.py` 中的 `n_genes` 变量。

### Q: 优先使用 diff/hilo padj 还是方差筛选更好？

- **优先 padj**：差异显著的基因更有生物学意义，构建的模块更可能与性状相关
- **方差筛选**：无偏筛选，适合没有明确差异分组的数据集

若 diff/hilo 结果存在且 `debug: false`，系统默认使用 padj 优先。设为 `debug: true` 可强制使用方差筛选。

### Q: WGCNA 运行很慢怎么办？

- 减少基因数（修改源码中的 `n_genes`）
- 确保已运行过 diff/hilo（padj 筛选通常比方差筛选得到的基因更有信息量）
- PyWGCNA 的软阈值计算和 TOM 矩阵构建是主要耗时步骤，无法通过配置跳过

## 示例

以下示例基于 `GSE143318`，目标基因 `APEX1`，`exp_type: "Fibrosis"`。

### 输入

表达矩阵：27 个样本 × 基因。基因数筛选目标 3000，优先使用 diff 结果（`differential_summary.pkl`），按 padj 升序取 top。

### 输出：性状矩阵（`res/csv/GSE143318_APEX1_datTraits.csv`）

```csv
Sample,Group,APEX1_exp
GSM4257053,Fibrosis,1444.0
GSM4257054,Fibrosis,1363.0
GSM4257055,Fibrosis,1422.0
GSM4257056,Fibrosis,1581.0
GSM4257057,Fibrosis,1422.0
GSM4257058,Control,1191.0
GSM4257059,Control,369.0
...
```

> 表达矩阵 `datExpr.csv` 较大（样本 × 3000 基因），此处不展示完整内容。文件路径：`res/csv/GSE143318_APEX1_datExpr.csv`。

### 输出：PyWGCNA 图表

PyWGCNA 在阶段 2 内自动生成以下图表至 `res/figures/wgcna/figures/`：

| 图 | 文件名 |
|----|--------|
| 模块-性状热图 | `{name}_module-traitRelationships.pdf` |
| 样本聚类清洗 | `{name}_sample_clustering_cleaning.pdf` |
| 软阈值筛选 | `{name}_summary_power.pdf` |
| 模块 eigengene 图 | `{name}_module_heatmap_eigengene_{module}.pdf` |
| GO 富集点图 | `{name}_GO/` 目录下各模块的点图 |
