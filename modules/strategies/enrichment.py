import os
import re
from typing import Optional, List

import numpy as np
import pandas as pd

import gseapy

from utils.paths import RESULT_DIR


_PROBE_PREFIXES = (
    "ILMN_", "AFFY-", "A_", "CUST_", "GE_", "TC0", "TC1",
    "Hs.", "Mm.", "Rn.", "DDB_", "ENS", "FBgn_",
)
_MIN_GENES_FALLBACK = 10


class EnrichStrategy:
    """富集分析策略，基于差异分析结果做 KEGG / GO 通路富集。

    读取 diff 或 hilo 分析结果（PKL 优先，CSV 回退），筛选显著差异
    基因列表，通过 Enrichr Web API 进行富集分析，返回合并结果。
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.cfg = analyzer.cfg
        self._logger = analyzer._logger

    def calculate(self) -> Optional[pd.DataFrame]:
        source_mode = getattr(self.cfg, "enrichment_source_mode", "diff")
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        source_mapping = {"diff": "differential", "hilo": "highlow"}
        source_key = source_mapping.get(source_mode, "differential")

        csv_path = os.path.join(RESULT_DIR, "csv", f"{self.cfg.gse_id}_{source_key}_summary.csv")
        pkl_path = os.path.join(data_dir, "pkl", f"{self.cfg.gse_id}_{source_key}_summary.pkl")

        source_df = self._read_source_results(csv_path, pkl_path)
        if source_df is None or source_df.empty:
            raise FileNotFoundError(
                f"未找到 {source_mode} 分析结果。请先运行 diff 或 hilo 分析。"
            )

        gene_list = self._filter_significant_genes(source_df)
        if len(gene_list) == 0:
            self._logger.warning("没有基因通过筛选阈值。")
            return None
        self._logger.info(f"{len(gene_list)} 个显著基因入选富集分析。")

        gene_sets = getattr(self.cfg, "enrichment_gene_sets", ["KEGG_2016"])
        organism = getattr(self.cfg, "organism", "human")
        all_results = []

        for gs in gene_sets:
            self._logger.info(f"正在对基因集 {gs} 进行富集分析...")
            try:
                enr = gseapy.enrichr(
                    gene_list=gene_list,
                    gene_sets=gs,
                    organism=organism,
                    cutoff=self.cfg.p_threshold,
                    outdir=None,
                    no_plot=True,
                    verbose=False,
                )
                if enr.res2d is not None and not enr.res2d.empty:
                    df = enr.res2d.copy()
                    df["Gene_set"] = gs
                    all_results.append(df)
                    self._logger.info(f"{gs} 富集完成，共 {len(df)} 个条目。")
            except Exception as e:
                self._logger.error(f"{gs} 富集失败: {e}")
                continue

        if not all_results:
            self._logger.error("所有富集分析均返回空结果。")
            return None

        combined = pd.concat(all_results, ignore_index=True)
        before = len(combined)
        combined = self._deduplicate_enrichment(combined)
        after = len(combined)
        if before != after:
            self._logger.info(f"去重移除 {before - after} 个重复条目。")
        self._logger.info(f"富集分析完成: {len(combined)} 个条目，{len(all_results)} 个基因集。")
        return combined

    def _deduplicate_enrichment(self, df: pd.DataFrame) -> pd.DataFrame:
        """三方面去重：指纹精确匹配、数值接近、基因高重叠。

        同一 Gene_set 内，对 Term 名有子串关系的两个条目，只要满足
        任一条件即视为同一通路，保留较短的 Term。
        """
        if df.empty or "Gene_set" not in df.columns:
            return df

        results = []
        for gs, grp in df.groupby("Gene_set", sort=False):
            rows = grp.sort_values(
                by="Term", key=lambda s: s.str.len(), na_position="last"
            ).reset_index(drop=True)
            n = len(rows)
            removed = set()

            for i in range(n):
                if i in removed:
                    continue
                for j in range(i + 1, n):
                    if j in removed:
                        continue
                    # strip GO ID 后缀避免不同 GO ID 阻止子串匹配
                    _go_re = re.compile(r"\s*\(GO:\d+\)")
                    term_i = _go_re.sub("", str(rows.iloc[i]["Term"]))
                    term_j = _go_re.sub("", str(rows.iloc[j]["Term"]))
                    if term_i not in term_j and term_j not in term_i:
                        continue
                    if (self._is_fingerprint_match(rows, i, j)
                            or self._is_score_close(rows, i, j)
                            or self._is_gene_overlap_high(rows, i, j)):
                        removed.add(j)

            for i in range(n):
                if i not in removed:
                    results.append(rows.iloc[[i]])

        return pd.concat(results, ignore_index=True)

    @staticmethod
    def _is_fingerprint_match(rows: pd.DataFrame, i: int, j: int) -> bool:
        """指纹匹配：关键数值列 round 6 位后全等，且 Overlap/Genes 字符串一致。"""
        for col in ["Adjusted P-value", "Odds Ratio", "Combined Score"]:
            if col in rows.columns:
                a, b = rows.iloc[i][col], rows.iloc[j][col]
                try:
                    if round(float(a), 6) != round(float(b), 6):
                        return False
                except (TypeError, ValueError):
                    if str(a) != str(b):
                        return False
        for col in ["Overlap", "Genes"]:
            if col in rows.columns:
                if str(rows.iloc[i][col]).strip() != str(rows.iloc[j][col]).strip():
                    return False
        return True

    @staticmethod
    def _is_score_close(rows: pd.DataFrame, i: int, j: int) -> bool:
        """数值接近：Combined Score 保留 2 位小数一致，且 Overlap 相同。"""
        if "Combined Score" not in rows.columns:
            return False
        a, b = rows.iloc[i]["Combined Score"], rows.iloc[j]["Combined Score"]
        try:
            if round(float(a), 2) != round(float(b), 2):
                return False
        except (TypeError, ValueError):
            if str(a) != str(b):
                return False
        if "Overlap" in rows.columns:
            if str(rows.iloc[i]["Overlap"]).strip() != str(rows.iloc[j]["Overlap"]).strip():
                return False
        return True

    @staticmethod
    def _is_gene_overlap_high(rows: pd.DataFrame, i: int, j: int,
                              threshold: float = 0.8) -> bool:
        """基因重叠高：Genes 列交集占较小基因集的比例 ≥ 阈值。"""
        if "Genes" not in rows.columns:
            return False
        g_i = {g.strip() for g in str(rows.iloc[i]["Genes"]).split(";") if g.strip()}
        g_j = {g.strip() for g in str(rows.iloc[j]["Genes"]).split(";") if g.strip()}
        if not g_i or not g_j:
            return False
        overlap = len(g_i & g_j)
        smaller = min(len(g_i), len(g_j))
        return smaller > 0 and overlap / smaller >= threshold

    def _read_source_results(self, csv_path: str, pkl_path: str) -> Optional[pd.DataFrame]:
        if os.path.exists(pkl_path):
            self._logger.info(f"从 PKL 读取源结果: {pkl_path}")
            return pd.read_pickle(pkl_path)
        elif os.path.exists(csv_path):
            self._logger.info(f"从 CSV 读取源结果: {csv_path}")
            return pd.read_csv(csv_path)
        return None

    def _filter_significant_genes(self, df: pd.DataFrame) -> List[str]:
        p_thr = self.cfg.p_threshold
        log2fc_thr = getattr(self.cfg, "log2fc_threshold", 0.0)
        max_genes = getattr(self.cfg, "max_input_genes", 500)

        if "padj" not in df.columns:
            raise ValueError("源结果缺少 padj 列，无法筛选。")
        if "Gene" not in df.columns:
            raise ValueError("源结果缺少 Gene 列，无法筛选。")

        df = df.copy()
        df["Gene"] = df["Gene"].astype(str).str.strip()
        original_count = len(df)

        # 清洗无效基因名
        df = df[~df["Gene"].isin(["nan", "None", ""])]

        # 过滤只含数字的（Entrez ID 或未知 ID）
        df = df[~df["Gene"].str.fullmatch(r"\d+")]

        # 过滤探针 ID
        pattern = "^(?:" + "|".join(_PROBE_PREFIXES) + ")"
        df = df[~df["Gene"].str.match(pattern)]

        cleaned_count = original_count - len(df)
        if cleaned_count:
            self._logger.info(f"清洗掉 {cleaned_count} 个无效/探针 ID 基因名。")

        # 按 padj 筛选
        mask = df["padj"].notna() & (df["padj"] < p_thr)
        if "log2FC" in df.columns and log2fc_thr > 0:
            mask = mask & (df["log2FC"].abs() >= log2fc_thr)

        filtered = df.loc[mask].drop_duplicates(subset="Gene").copy()
        strictly_significant = len(filtered)

        # 如果显著基因太少，回退到按 |log2FC| 排名取 top N
        if strictly_significant < _MIN_GENES_FALLBACK:
            self._logger.warning(
                f"padj < {p_thr} 仅筛出 {strictly_significant} 个基因，"
                f"不足 {_MIN_GENES_FALLBACK}，将回退为按 |log2FC| 排名取前 {max_genes} 个。"
            )
            fallback = df.drop_duplicates(subset="Gene").copy()
            if "log2FC" in fallback.columns:
                fallback = fallback.sort_values(
                    by="log2FC", key=abs, ascending=False, na_position="last"
                )
            else:
                fallback = fallback.sort_values(by="padj", ascending=True, na_position="last")
            filtered = fallback.head(max_genes)
        else:
            sort_col = "log2FC" if "log2FC" in filtered.columns else "padj"
            filtered = filtered.sort_values(
                by=sort_col, key=abs, ascending=False, na_position="last"
            )
            if len(filtered) > max_genes:
                filtered = filtered.head(max_genes)

        genes = filtered["Gene"].tolist()

        self._logger.info(
            f"筛选得到 {len(genes)} 个基因 "
            f"(padj < {p_thr}"
            f"{', |log2FC| >= ' + str(log2fc_thr) if log2fc_thr > 0 else ''}"
            f")"
        )
        return genes
