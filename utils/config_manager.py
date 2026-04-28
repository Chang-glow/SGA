import pandas as pd
from dataclasses import dataclass, field
from typing import Optional, List
from hydra.core.config_store import ConfigStore

from utils.paths import DATA_DIR


@dataclass
class Config:
    tar_gene: str
    gse_id: str
    version: str = "1.0"
    group_select_col: Optional[str] = "source_name_ch1"
    control_label: Optional[List[str]] = field(default_factory=lambda: ["control"])
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


@dataclass
class DataHandler:
    meta_matrix_pack: Optional[dict] = None
    gene_corr_table: Optional[pd.DataFrame] = None
    gene_diff_table: Optional[pd.DataFrame] = None


cs = ConfigStore.instance()
cs.store(name="Config", node=Config)
