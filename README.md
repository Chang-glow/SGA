# **SGA: Simple GEO Analyzer (v1.1)**

**SGA** 是一个为生信初学者和课题组日常科研设计的轻量级、自动化 GEO 数据处理工具。它将 GEO 数据下载、清洗、分析与可视化流程集成在一起，帮助你快速从 GEO 原始数据得到差异分析或相关性结果。

项目最初在生物科学本科学习阶段开发，旨在减少课题组在转录组数据处理、Polb/肝纤维化相关分析中的重复劳动。

---

## **🌟 核心特性**

- **⚙️ 自动流水线**：只需提供 GSE ID，自动检索 GEO 数据并提取可用矩阵文件。  
- **🔄 缓存机制**：自动检测本地 `data/{GSE_ID}/pkl/` 缓存，跳过重复下载与重复清洗。  
- **📊 多分析模式**：支持差异分析（`diff`）和相关性分析（`corr`），可灵活切换。  
- **📈 可视化输出**：差异分析生成分组箱线图，相关性分析生成散点图并自动添加拟合线。  
- **🧠 智能分组映射**：支持元数据自动分组，可根据原始表达矩阵列名进行样本映射。  
- **🔧 Hydra 配置**：通过配置文件或命令行覆盖参数，无需修改源码即可切换分析目标。

---

## **📂 项目结构**

```
SGA/
├── main.py                      # 项目入口，调度数据下载、分析与可视化流程
├── conf/                        # Hydra 配置文件存放目录
│   ├── config.yaml.template     # 配置模板文件
│   └── config.yaml              # 实际运行配置（首次运行时自动生成）
├── modules/                     # 核心功能模块
│   ├── __init__.py
│   ├── data_loader.py           # GEO 数据下载、矩阵加载、分组与缓存打包
│   ├── calculater.py            # 差异分析与相关性分析逻辑
│   └── fig_plotter.py           # 差异箱线图与相关性散点图绘制
├── utils/                       # 工具模块
│   ├── __init__.py
│   ├── config_manager.py        # Hydra 配置映射与数据传递对象
│   ├── loggers.py               # 日志配置与记录
│   ├── parse_user_input.py      # 交互式输入解析
│   └── paths.py                 # 项目路径管理与目录初始化
├── data/                        # 下载的原始数据与缓存数据包
├── error_logs/                  # 运行时错误日志
├── res/                         # 可视化输出目录
│   └── figures/                 # 绘制图像存放位置
├── environment.yml              # Conda 环境配置
├── requirements.txt             # Python 依赖清单
└── README.md                    # 项目说明文档
```

---

## **🚀 快速上手**

### 1. 环境准备

```bash
git clone https://github.com/YourUsername/SGA.git
cd SGA

conda env create -f environment.yml
conda activate sga

# 或者
pip install -r requirements.txt
```

### 2. config.yaml.template 用法说明

`conf/config.yaml.template` 是项目默认配置模板，用于定义分析参数与运行行为。

- 首次运行 `python main.py` 时，如果 `conf/config.yaml` 不存在，程序会自动将 `conf/config.yaml.template` 复制为 `conf/config.yaml`。
- 建议编辑 `conf/config.yaml` 来设置你的分析目标；不要直接修改模板文件，模板文件用于初始化配置和恢复默认参数。
- 如果需要重置配置，删除 `conf/config.yaml` 后重新运行 `python main.py` 即可重新生成默认配置。

模板中的常用配置项：

- `tar_gene`：目标基因名称（例如 `APEX1`）。  
- `gse_id`：GEO 系列编号（例如 `GSE143318`）。  
- `data_dir`：项目数据根目录，默认 `data`。  
- `group_select_col`：用于自动分组的元数据列名。  
- `control_label` / `fib_label`：识别 Control / Fibrosis 组的关键词列表。  
- `analysis_mode`：分析模式，支持 `diff`（差异分析）和 `corr`（相关性分析）。  
- `storage`：是否保存缓存包与分析结果。  
- `strict_mode`：是否严格对齐表达矩阵与元数据样本。  
- `debug`：开启调试模式，打印更多检查信息。  
- `p_threshold`：相关性分析显著性阈值。  
- `log_threshold`：判断是否已 log 转换的阈值。  
- `signs`：相关性方向过滤选项，支持 `positive` / `negative`。

### 3. 运行分析

```bash
python main.py
```

命令行覆盖参数：

```bash
python main.py tar_gene="Acta2" gse_id="GSE123456" analysis_mode="diff"
```

---

## **📁 输出目录说明**

- `data/{GSE_ID}/pkl/`：保存清洗后的数据包和分析结果。  
- `data/{GSE_ID}/csv/`：保存结果表格。  
- `res/figures/`：保存生成的差异箱线图或相关性散点图。

---

## **💡 使用提示**

- 如果目标数据集包含多个补充矩阵，程序会提示你选择需要加载的文件。  
- 如果自动分组失败，会进入交互式分组流程，帮助你手动选择目标组。  
- 表达矩阵若未 log 转换，程序会自动执行 `log2(x + 1)`。

---

## **📝 未来改进方向**

- [ ] 增加多基因热图与层次聚类展示。  
- [ ] 支持更多数据源，如 ArrayExpress、TCGA。  
- [ ] 扩展统计分析方法，如斯皮尔曼相关、更多差异分析模型。

---

## **🤝 贡献与反馈**

如果你在处理特定 GEO 数据集时遇到问题，欢迎提交 Issue 或 Pull Request。非常欢迎关于分析逻辑、分组映射或可视化效果的改进建议。

*📅 Update at: 2026-04-28*
