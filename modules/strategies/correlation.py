import pandas as pd
from typing import Optional

from modules.calculater import fetch_gene_vector, detect_gene_case_convention, normalize_gene_symbol


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
            if name in ("meta", "meta_full") or not isinstance(df, pd.DataFrame):
                continue

            self._logger.info(f"--- 当前处理数据{name} ---")
            self._logger.info(f"提取目标基因{tar_gene}的数据中...")
            target_vec = fetch_gene_vector(df, tar_gene)

            self._logger.info("将以常见标识物分类进行计算并储存")
            convention = detect_gene_case_convention(df.index)
            homolog_map = getattr(self.cfg, "homolog_map", {}) or {}
            for category, gene_list in self.analyzer._get_hfm_dict().items():
                for gene in gene_list:
                    marker_vec = self._find_marker(df, gene, convention, homolog_map)
                    if marker_vec is None:
                        continue

                    r, p = None, None
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

    def _find_marker(self, df: pd.DataFrame, gene: str, convention: str,
                     homolog_map: dict) -> pd.Series | None:
        """在表达矩阵中查找标志物基因。依次尝试：原始名 → 规范化名 → 同源映射。"""
        normalized = normalize_gene_symbol(gene, convention)
        candidates = [gene]
        if normalized != gene:
            candidates.append(normalized)
        if homolog_map:
            m = homolog_map.get(gene) or homolog_map.get(normalized)
            if m:
                candidates.append(m)
            rev = {v: k for k, v in homolog_map.items()}
            rm = rev.get(gene) or rev.get(normalized)
            if rm:
                candidates.append(rm)
        # 先查 index（避免 fetch_gene_vector 对中间候选打无意义的 WARNING）
        index_upper = df.index.astype(str).str.upper()
        for name in candidates:
            if name.upper() in index_upper.values:
                return fetch_gene_vector(df, name)
        # 回退：列查找（兼容非标准格式）
        for name in candidates:
            vec = fetch_gene_vector(df, name)
            if not vec.empty:
                return vec
        self._logger.warning(f"未能在矩阵中找到基因: {gene}（已尝试: {', '.join(candidates)}）")
        return None
