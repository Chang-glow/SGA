import os
import re
import pandas as pd
from dataclasses import dataclass, field
from typing import Optional, List, Dict
from hydra.core.config_store import ConfigStore

from utils.paths import DATA_DIR


@dataclass
class Config:
    tar_gene: str = ""
    gse_id: str = ""
    version: str = "0.3.3"
    group_select_col: Optional[str] = "source_name_ch1"
    control_label: Optional[List[str]] = field(default_factory=lambda: ["control", "ND"])
    exp_label: Optional[List[str]] = field(default_factory=lambda: ["CCl4"])
    exp_type: Optional[str] = "Fibrosis"
    analysis_mode: str = "diff"
    data_dir: str = DATA_DIR
    storage: bool = True
    strict_mode: bool = False
    debug: bool = False
    log_threshold: int = 50
    p_threshold: float = 0.05
    signs: List[str] = field(default_factory=lambda: ["positive", "negative"])
    enrichment_source_mode: str = "diff"
    enrichment_gene_sets: List[str] = field(default_factory=lambda: [
        "KEGG_2026",
        "GO_Biological_Process_2025",
        "GO_Molecular_Function_2025",
        "GO_Cellular_Component_2025",
    ])
    organism: str = "human"
    enrichment_min_genes_fallback: int = 10  # 显著基因不足时回退到 |log2FC| 排序
    strict_filter: bool = True
    log2fc_threshold: float = 0.0
    max_input_genes: int = 500
    max_output_genes: int = 0
    min_samples_per_group: int = 3   # t-test / 免疫检验最小样本数
    gene_blacklist: List[str] = field(default_factory=list)
    group_memory_enabled: bool = False
    group_memory_use: bool = False
    tar_tuple: str = ""
    custom_marker_dict: Dict[str, List[str]] = field(default_factory=dict)  # 相关性分析自定义标志物基因集
    immune_method: str = "DeconRNASeq"
    immune_low_abundance_threshold: float = 0.001  # 免疫低丰度细胞过滤阈值
    multi_gene: str = ""           # 多基因：逗号分隔 OR 文件路径（每行一个基因）
    overwrite_figures: bool = False  # 默认不覆盖，迭代 (1)(2)...
    figure_dpi: int = 300            # 图片分辨率
    volcano_top_n_labels: int = 15   # 火山图标注基因数
    heatmap_top_n_genes: int = 20    # 热图展示 top DEG 数
    immune_boxplot_n_cols: int = 5   # 免疫箱线图每行子图数
    enrich_plot_top_terms: int = 15  # 富集图展示 top term 数
    plot_data_warning: bool = True   # 画图数据合理性 warning 开关（通用）
    process: str = "123"
    wgcna_top_n_genes: int = 3000    # WGCNA 筛选基因数
    _batch_suffix: str = field(default="", init=False)  # 批处理文件后缀（运行时注入，非 Hydra 配置项）

    def __post_init__(self):
        # 互斥检查：multi_gene 优先，自动清 tar_gene
        has_tar = bool(self.tar_gene)
        has_multi = bool(self.multi_gene)
        if has_tar and has_multi:
            object.__setattr__(self, "tar_gene", "")
            has_tar = False
        if not has_tar and not has_multi:
            raise ValueError("tar_gene 和 multi_gene 必须存在一个")


def parse_tar_genes(tar_gene: str, multi_gene: str) -> List[str]:
    """从 tar_gene / multi_gene 解析目标基因列表（供 FigurePlotter 等使用）"""
    if multi_gene:
        if os.path.isfile(multi_gene):
            with open(multi_gene) as f:
                return [l.strip() for l in f if l.strip() and not l.startswith("#")]
        return [g.strip() for g in re.split(r"[,;+_]", multi_gene) if g.strip()]
    return [tar_gene.strip()] if tar_gene else []


@dataclass
class DataHandler:
    meta_matrix_pack: Optional[dict] = None
    gene_corr_table: Optional[pd.DataFrame] = None
    gene_diff_table: Optional[pd.DataFrame] = None
    gene_enrich_table: Optional[pd.DataFrame] = None
    gene_immune_table: Optional[pd.DataFrame] = None
    gene_wgcna_table: Optional[pd.DataFrame] = None


cs = ConfigStore.instance()
cs.store(name="Config", node=Config)
