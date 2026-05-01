# **SGA: Simple GEO Analyzer (v1.3)**

**SGA** 是一个为生信初学者和课题组日常科研设计的轻量级、自动化 GEO 数据处理工具。它将 GEO 数据下载、清洗、分析与可视化流程集成在一起，帮助你快速从 GEO 原始数据得到差异分析或相关性结果。

项目最初在生物科学本科学习阶段开发，旨在减少课题组在转录组数据处理、Polb/肝纤维化相关分析中的重复劳动。

---

## **🌟 核心特性**

- **自动流水线**：GSE ID → 下载 → 清洗 → 分析 → 画图
- **缓存机制**：跳过重复下载/计算
- **五种分析模式**：`diff` / `corr` / `hilo` / `enrich` / `immune`
- **可视化**：箱线图、聚类热图、散点拟合线、堆叠柱状图、相关性热图、火山图、气泡图
- **智能分组 + 分组记忆**
- **Hydra 配置驱动**，命令行覆盖

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
│   ├── csv/                     # CSV 结果表格
│   └── figures/                 # 绘制图像存放位置（按分析模式分子目录）
├── environment.yml              # Conda 环境配置
├── setup.sh                     # 一键安装脚本
├── requirements.txt             # Python 依赖清单
└── README.md                    # 项目说明文档
```

---

## **🚀 快速上手**

```bash
# Linux / macOS
git clone https://github.com/YourUsername/SGA.git
cd SGA
bash setup.sh && source ~/.bashrc   # 创建环境 + 注册 SGA 命令
SGA help                            # 查看所有配置项参考
```

```bash
# 或者手动安装（所有平台）
conda env create -f environment.yml
conda activate sga
pip install -r requirements.txt
python main.py help
```

### 运行分析

```bash
SGA tar_gene=APEX1 gse_id=GSE143318 analysis_mode=immune
```

命令行覆盖参数：

```bash
SGA tar_gene=Acta2 gse_id=GSE123456 analysis_mode=diff
```

执行 `SGA help` 查看所有配置项的完整参考。

---

## **📁 输出目录说明**

- `data/{GSE_ID}/pkl/`：保存清洗后的数据包和分析结果。
- `res/csv/`：保存结果表格。
- `res/figures/{analysis_mode}/`：保存生成的差异箱线图或相关性散点图。

---

## **💡 使用提示**

- 如果目标数据集包含多个补充矩阵，程序会提示你选择需要加载的文件。
- 如果自动分组失败，会进入交互式分组流程，帮助你手动选择目标组。
- 表达矩阵若未 log 转换，程序会自动执行 `log2(x + 1)`。

---

## **🤝 贡献与反馈**

如果你在处理特定 GEO 数据集时遇到问题，欢迎提交 Issue 或 Pull Request。非常欢迎关于分析逻辑、分组映射或可视化效果的改进建议。

*📅 Update at: 2026-05-01*
