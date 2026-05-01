import os
from typing import Optional

import numpy as np
import pandas as pd

from .difference import DiffStrategy


class HighLowStrategy:
    """实现 HighLow 分析策略，比较 High 组和 Low 组的基因表达差异"""
    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.cfg = analyzer.cfg
        self._logger = analyzer._logger

    def calculate(self) -> Optional[pd.DataFrame]:
        # 首先检查是否存在缓存结果
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        gene_hilo_path = os.path.join(data_dir, "pkl", f"{self.cfg.gse_id}_highlow_summary.pkl")

        diff_strategy = DiffStrategy(self.analyzer)
        cached_result = diff_strategy._read_cached_summary(gene_hilo_path)
        if cached_result is not None:
            return cached_result

        if not self.analyzer._meta_matrix_pack:
            self._logger.info("正在读取数据...")
            self.analyzer.roaming_data()

        # 从数据包中提取表达矩阵和分组信息
        pack = self.analyzer._meta_matrix_pack
        meta = pack.get("meta")

        if meta is None or "group" not in meta.columns or meta["group"].isnull().all():
            self._logger.error("未找到分组信息，无法进行 hilo 分析")
            raise ValueError("未找到分组信息，无法进行 hilo 分析")

        low_samples = meta[meta["group"] == "Low"].index.tolist()
        high_samples = meta[meta["group"] == "High"].index.tolist()
        if not low_samples or not high_samples:
            self._logger.error("hilo 分组信息不完整，缺失 Low 或 High 样本")
            raise ValueError("hilo 分组信息不完整，缺失 Low 或 High 样本")

        expr_keys = [k for k in pack.keys() if k not in {"meta", "meta_full", "group_info"}]
        if not expr_keys:
            raise KeyError("No expression matrix found in data pack")

        expr_df = pack[expr_keys[0]]
        if isinstance(expr_df, dict) and "matrix_aligned" in expr_df:
            expr_df = expr_df["matrix_aligned"]

        diff_strategy = DiffStrategy(self.analyzer)
        expr_df = diff_strategy._align_expression_matrix(expr_df, low_samples + high_samples, pack)
        expr_df = diff_strategy._prefilter_expression_matrix(expr_df)
        result = diff_strategy._calculate_diff(expr_df, low_samples, high_samples)

        return result
