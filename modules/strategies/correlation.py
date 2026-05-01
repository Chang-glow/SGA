import pandas as pd
from typing import Optional

from modules.calculater import fetch_gene_vector


class CorrelationStrategy:
    """实现 Correlation 分析策略，计算目标基因与常见标识物基因的相关性"""
    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.cfg = analyzer.cfg
        self._logger = analyzer._logger

    def calculate(self) -> Optional[pd.DataFrame]:
        if not self.analyzer._meta_matrix_pack:
            self._logger.info("正在读取数据...")
            self.analyzer.roaming_data()

        meta_matrix_pack = self.analyzer._meta_matrix_pack
        if meta_matrix_pack:
            self._logger.info("读取成功，将继续分析")

        tar_gene = self.cfg.tar_gene
        results_list = []

        for name, df in meta_matrix_pack.items():
            if name == "meta":
                continue

            self._logger.info(f"--- 当前处理数据{name} ---")
            self._logger.info(f"提取目标基因{tar_gene}的数据中...")
            target_vec = fetch_gene_vector(df, tar_gene)

            self._logger.info("将以常见标识物分类进行计算并储存")
            for category, gene_list in self.analyzer.hfm_dict.items():
                for gene in gene_list:
                    self._logger.debug(f"提取标识物基因{gene}的数据中...")
                    marker_vec = fetch_gene_vector(df, gene)

                    r, p = None, None
                    if not marker_vec.empty:
                        self._logger.debug("数据提取完成，计算相关性中...")
                        r, p = self.analyzer.pearson_analyze(target_vec, marker_vec)

                    if r is not None:
                        self._logger.debug("相关性计算完成！")
                        results_list.append({
                            "Matrix": name,
                            "Category": category,
                            "Gene": gene,
                            "R": r,
                            "P_value": p
                        })

        if results_list:
            self._logger.info("相关性计算完成！")
            return pd.DataFrame(results_list)

        return None
