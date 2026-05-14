# **SGA: Simple GEO Analyzer** ![version](https://img.shields.io/github/v/tag/Chang-glow/SGA)

**sga** 是一个为生信初学者日常科研设计的轻量级、自动化 GEO 数据处理工具。它将 GEO 数据下载、清洗、分析与可视化流程集成在一起，可以快速从 GEO 原始数据得到差异分析或相关性结果。

---

## **🌟 核心特性**

- **自动流水线**：GSE ID → 下载 → 清洗 → 分析 → 画图
- **缓存机制**：跳过重复下载/计算
- **六种分析模式**：分组差异分析`diff` / 相关性分析`corr` / 高低表达分析`hilo` / 通路富集`enrich` / 免疫浸润`immune` / 加权基因共表达网络`wgcna`
- **可视化**：箱线图、聚类热图、散点拟合线、堆叠柱状图、相关性热图、火山图、气泡图
- **免疫浸润**：四种去卷积算法（DeconRNASeq、CIBERSORT、ssGSEA、SingScore），纯 Python 无 R 依赖
- **智能分组 + 分组记忆**
- **Hydra 配置驱动**，命令行覆盖
- **bash 预解析**：`version` / `help` / `config` 子命令秒出结果（零 Python 启动开销）

---

## **📂 项目结构**

```
SGA/
├── main.py                      # 项目入口，调度数据下载、分析与可视化流程
├── conf/                         # Hydra 配置文件存放目录
│   ├── config.yaml.template      # 配置模板文件
│   ├── config.yaml               # 实际运行配置（由 setup.sh 生成）
│   └── help.yaml                 # 自定义帮助信息
├── doc/                          # 文档
├── modules/                      # 核心功能模块
│   ├── __init__.py
│   ├── data_loader.py            # GEO 数据下载、矩阵加载、分组与缓存打包
│   ├── data_packer.py            # 数据打包/解包
│   ├── calculater.py             # 差异分析、相关性分析、免疫浸润等计算逻辑
│   ├── fig_plotter.py            # 差异箱线图、相关性散点图、免疫浸润图等绘制
│   └── strategies/               # 分析策略（按模式分文件）
│       ├── correlation.py
│       ├── difference.py
│       ├── highlow.py
│       ├── enrichment.py
│       ├── immune.py
│       └── wgcna.py
├── utils/                        # 工具模块
│   ├── __init__.py
│   ├── config_manager.py         # Hydra 配置映射与数据传递对象
│   ├── loggers.py                # 日志配置与记录
│   ├── parse_user_input.py       # 交互式输入解析
│   └── paths.py                  # 项目路径管理与目录初始化
├── data/                         # 下载的原始数据与缓存数据包
├── error_logs/                   # 运行时错误日志
├── res/                          # 可视化输出目录
│   ├── csv/                      # CSV 结果表格
│   └── figures/                  # 绘制图像存放位置（按分析模式分子目录）
├── environment.yml               # Conda 环境配置
├── setup.sh                      # 一键安装脚本
├── requirements.txt              # Python 依赖清单
└── README.md                     # 项目说明文档
```

---

## **🚀 快速上手**

```bash
# Linux / macOS
git clone https://github.com/YourUsername/SGA.git
cd SGA
bash setup.sh && source ~/.bashrc   # 创建环境 + 生成配置 + 注册命令
# 编辑 conf/config.yaml 设置 tar_gene、gse_id 等参数
sga help                            # 查看所有配置项参考
```

```bash
# 或者手动安装（所有平台）
conda env create -f environment.yml            # 通过文件创建环境
conda activate sga                             # 启动环境
cp conf/config.yaml.template conf/config.yaml  # 手动生成配置
python main.py                                 # 直接运行（该安装方式下无 sga 命令）
```

### 查看当前配置

```bash
sga config               # 核心配置摘要（含版本号）
sga config --all         # 全部配置
sga config immune_method # 查询单个配置项
sga version              # 查看版本号（来自 git tag）
```

### 运行分析

```bash
sga tar_gene=APEX1 gse_id=GSE143318 analysis_mode=immune
```

命令行覆盖参数：

```bash
sga tar_gene=Acta2 gse_id=GSE123456 analysis_mode=diff
sga tar_gene=APEX1 gse_id=GSE143318 analysis_mode=immune immune_method=ssGSEA
```

**多基因输入**（enrich 模式可直接输入基因列表，无需表达数据）：

```bash
sga analysis_mode=enrich multi_gene=TP53_EGFR_KRAS
sga analysis_mode=enrich multi_gene=./my_genes.txt
```

执行 `sga help` 查看所有配置项的完整参考。

---

## **📁 输出目录说明**

- `data/{GSE_ID}/pkl/`：保存清洗后的数据包和分析结果。
- `res/csv/`：保存结果表格。
- `res/figures/{analysis_mode}/`：保存生成的差异箱线图或相关性散点图。

---

## **💡 使用提示**

- 如果目标数据集包含多个补充矩阵，程序会提示选择需要加载的文件。
- 如果自动分组失败，会进入交互式分组流程，手动选择目标组。
- 表达矩阵若未 log 转换，程序会自动执行 `log2(x + 1)`。
- 免疫浸润分析等需要原始数据的策略不会进行自动 log 转换。
- 免疫浸润使用 ssGSEA 或 SingScore 时，堆叠柱状图会提示数据非比例（可通过 `plot_data_warning: false` 关闭提示）。详见 `doc/immune_infiltration.md`。
- `tar_gene` 与 `multi_gene` 互斥，同时设置时 `multi_gene` 优先。`multi_gene` 在 CLI 中使用 `_` 分隔基因名，YAML 内支持 `,` `;` `+` `_`。

---

## **🤝 贡献与反馈**

如果你在处理特定 GEO 数据集时遇到问题，欢迎提交 Issue 或 Pull Request。非常欢迎关于分析逻辑、分组映射或可视化效果的改进建议。

---

### **注：本项目由 Claude code 参与辅助开发**

*📅 Update at: 2026-05-14*
