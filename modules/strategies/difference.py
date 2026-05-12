import os
import re
import scipy
import numpy as np
import pandas as pd
from typing import Optional
from statsmodels.stats.multitest import fdrcorrection

_PROBE_RE = re.compile(
    r"^(ILMN_|AFFY-|A_\d|CUST_|GI_|NM_|NR_|XM_|XR_"
    r"|TC\d|Hs\.|Mm\.|Rn\.|DDB_|ENS|FBgn_|agc_|GE_)"
)
# 非编码 RNA + 假基因过滤规则
_NON_CODING_RES = [
    re.compile(r"^MIR", re.IGNORECASE),                # miRNA (含 MIRLET 家族)
    re.compile(r"^SNOR[A-Z]\d", re.IGNORECASE),          # snoRNA
    re.compile(r"^SCARNA\d", re.IGNORECASE),             # scaRNA
    re.compile(r"^RNU\d", re.IGNORECASE),                # snRNA
    re.compile(r"^RNVU\d", re.IGNORECASE),               # variant U snRNA
    re.compile(r"^LINC\d", re.IGNORECASE),               # lncRNA (lincRNA)
    re.compile(r"^CASC\d", re.IGNORECASE),               # cancer susceptibility lncRNA
    re.compile(r"^FLJ\d", re.IGNORECASE),                # FLJ clones (mostly non-coding)
    re.compile(r"^RNA5", re.IGNORECASE),                 # rRNA / RNA5 family
    re.compile(r"^RN7SKP\d", re.IGNORECASE),             # RN7SK pseudogenes
    re.compile(r"-AS\d", re.IGNORECASE),                 # 反义RNA
    re.compile(r"-IT\d+$", re.IGNORECASE),               # intronic lncRNA (FGF14-IT1)
    re.compile(r"-DT\d*$", re.IGNORECASE),               # divergent lncRNA (APCDD1L-DT)
    re.compile(r"^LOC\d", re.IGNORECASE),                # 未表征位点
    re.compile(r"^RP[LS]\d*[A-Z]*P\d", re.IGNORECASE),   # 核糖体蛋白假基因
    re.compile(r"[A-Z][A-Z0-9]*[A-Z]P\d+$", re.IGNORECASE),  # 通用假基因 (≥2字母+P+数字)
    re.compile(r"\d+P\d"),                                # digit-P-digit 假基因
    re.compile(r"\d+P$"),                                 # RNA假基因 (末尾 digit+P)
]


class DiffStrategy:
    """实现 Difference 分析策略，比较 Control 组和 Experiment 组的基因表达差异"""
    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.cfg = analyzer.cfg
        self._logger = analyzer._logger

    def calculate(self) -> Optional[pd.DataFrame]:
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        gene_diff_path = os.path.join(data_dir, "pkl", f"{self.cfg.gse_id}_differential_summary.pkl")

        cached_result = self._read_cached_summary(gene_diff_path)
        if cached_result is not None:
            return cached_result

        if not self.analyzer._meta_matrix_pack:
            self._logger.info("正在读取数据...")
            self.analyzer.roaming_data()

        pack = self.analyzer._meta_matrix_pack
        meta = pack.get("meta")

        if 'group' in meta.columns and not meta['group'].isnull().all():
            self._logger.info("读取成功，将继续分析")
        else:
            self._logger.error("未找到分组信息，无法进行差异分析")
            raise ValueError("未找到分组信息，无法进行差异分析")

        expr_keys = [k for k, v in pack.items() if k not in {"meta", "meta_full", "group_info"}]
        if not expr_keys:
            raise KeyError("No expression matrix found in data pack")

        expr_df = pack[expr_keys[0]]
        if isinstance(expr_df, dict) and 'matrix_aligned' in expr_df:
            expr_df = expr_df['matrix_aligned']

        exp_type = self.cfg.exp_type if self.cfg.exp_type else 'Experiment'
        control_samples = meta[meta['group'] == 'Control'].index.tolist()
        exp_samples = meta[meta['group'] == exp_type].index.tolist()
        if not control_samples or not exp_samples:
            self._logger.error(f"分组信息不完整,无法找到Control或{exp_type}组的样本")
            raise ValueError(f"分组信息不完整,无法找到Control或{exp_type}组的样本")

        expr_df = self._align_expression_matrix(expr_df, control_samples + exp_samples, pack)
        expr_df = self._prefilter_expression_matrix(expr_df)
        return self._calculate_diff(expr_df, control_samples, exp_samples)

    def _align_expression_matrix(self, expr_df: pd.DataFrame, common_samples: list, pack: dict) -> pd.DataFrame:
        if set(common_samples).issubset(expr_df.index):
            expr_df = expr_df.T

        missing_samples = [s for s in common_samples if s not in expr_df.columns]
        if missing_samples:
            self._logger.warning(f"样本名未在表达矩阵列中找到: {missing_samples[:10]}{'...' if len(missing_samples) > 10 else ''}")
            if 'meta_full' in pack:
                expr_df = self.analyzer._rename_expr_columns_by_meta_order(expr_df, pack['meta_full'])
                missing_samples = [s for s in common_samples if s not in expr_df.columns]
                if not missing_samples:
                    sample_columns = [c for c in common_samples if c in expr_df.columns]
                    other_columns = [c for c in expr_df.columns if c not in sample_columns]
                    expr_df = expr_df.loc[:, other_columns + sample_columns]
                    self._logger.info("已按原始元数据顺序将表达矩阵列映射为 GSM 样本名")
                else:
                    self._logger.error("样本名与原始元数据顺序无法完全对应,无法重新映射表达矩阵")
                    raise KeyError("样本名与表达矩阵列不匹配，请检查数据和元数据")
            else:
                self._logger.error("无法从表达矩阵列名中推断Control/Experiment分组，请检查数据")
                raise KeyError("样本名与表达矩阵列不匹配，请检查数据和元数据")
        else:
            sample_columns = [c for c in common_samples if c in expr_df.columns]
            other_columns = [c for c in expr_df.columns if c not in sample_columns]
            expr_df = expr_df.loc[:, other_columns + sample_columns]
        return expr_df

    def _prefilter_expression_matrix(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """剔除探针行并按基因符号去重（保留表达量最高的行），减少假阳性"""
        label_columns = ['SYMBOL', 'GENE', 'GENENAME', 'ENSEMBL', 'ENTREZID', 'ID_REF', 'TARGETID']
        gene_col = None
        for col in label_columns:
            if col in expr_df.columns:
                gene_col = col
                break

        if gene_col is None:
            self._logger.warning("未找到基因标注列，跳过探针过滤与去重")
            return expr_df

        before = len(expr_df)
        gene_labels = expr_df[gene_col].astype(str).str.strip()

        keep = (
            gene_labels.notna()
            & (gene_labels != "")
            & (gene_labels != "nan")
            & (gene_labels != "None")
            & ~gene_labels.str.fullmatch(r"\d+")
        )
        if not getattr(self.cfg, "tar_tuple", ""):
            keep = keep & ~gene_labels.str.match(_PROBE_RE)
            for ncrna_re in _NON_CODING_RES:
                keep = keep & ~gene_labels.str.contains(ncrna_re, regex=True, na=False)

        blacklist = getattr(self.cfg, "gene_blacklist", [])
        if blacklist:
            black_upper = [g.upper() for g in blacklist]
            keep = keep & ~gene_labels.str.upper().isin(black_upper)

        expr_df = expr_df.loc[keep]
        gene_labels = gene_labels[keep]

        sample_cols = self.analyzer._get_sample_columns(expr_df)
        expr_df = expr_df.copy()
        expr_df["_mean_expr_"] = expr_df[sample_cols].mean(axis=1)
        expr_df["_gene_label_"] = gene_labels.values
        expr_df = expr_df.sort_values("_mean_expr_", ascending=False).drop_duplicates(
            subset="_gene_label_"
        )
        expr_df = expr_df.drop(columns=["_mean_expr_", "_gene_label_"])

        removed = before - len(expr_df)
        self._logger.info(
            f"表达矩阵预过滤完成：剔除 {removed} 行（探针+重复基因），"
            f"保留 {len(expr_df)} 个唯一基因，有效样本列 {len(sample_cols)} 个"
        )
        return expr_df

    def _calculate_diff(self, expr_df: pd.DataFrame, control_samples: list, exp_samples: list) -> pd.DataFrame:
        control_values = expr_df[control_samples].to_numpy(dtype=float)
        exp_values = expr_df[exp_samples].to_numpy(dtype=float)

        self._logger.info("开始差异表达统计计算，可能需要一些时间")

        control_count = np.sum(~np.isnan(control_values), axis=1)
        exp_count = np.sum(~np.isnan(exp_values), axis=1)
        min_samples = getattr(self.cfg, "min_samples_per_group", 3)
        valid_mask = (control_count >= min_samples) & (exp_count >= min_samples)

        if not np.any(valid_mask):
            self._logger.error("没有足够的样本值进行差异分析")
            raise ValueError("没有足够的样本值进行差异分析")

        ttest_result = scipy.stats.ttest_ind(
            control_values,
            exp_values,
            axis=1,
            equal_var=False,
            nan_policy='omit'
        )

        if self.analyzer._is_log(expr_df):
            log2fc = np.nanmean(exp_values, axis=1) - np.nanmean(control_values, axis=1)
        else:
            log2fc = np.log2(np.nanmean(exp_values, axis=1) + 1) - np.log2(np.nanmean(control_values, axis=1) + 1)

        gene_labels = self._extract_gene_labels(expr_df)

        result = pd.DataFrame({
            "Gene": gene_labels,
            "log2FC": log2fc,
            "P_value": ttest_result.pvalue
        })
        result = result.loc[valid_mask].reset_index(drop=True)
        _, result['padj'] = fdrcorrection(result['P_value'].fillna(1))
        return result

    def _read_cached_summary(self, path: str) -> Optional[pd.DataFrame]:
        if os.path.exists(path) and not self.cfg.debug:
            df = pd.read_pickle(path)
            if self._summary_has_valid_genes(df):
                self._logger.info(f"发现现存分析结果：{path}，跳过差异分析")
                return df
            self._logger.warning(
                f"现存差异分析结果 {path} 中的 Gene 列似乎为数字索引，将重新计算以获取真实基因标签"
            )
        return None

    def _summary_has_valid_genes(self, df: pd.DataFrame) -> bool:
        if 'Gene' not in df.columns:
            return False
        genes = df['Gene'].dropna().astype(str).str.strip()
        if genes.empty:
            return False
        if genes.str.fullmatch(r'\d+').all():
            return False
        return True

    def _extract_gene_labels(self, expr_df: pd.DataFrame) -> pd.Series:
        """从表达矩阵中提取最优基因标签，优先使用 SYMBOL 等列"""
        label_columns = ['SYMBOL', 'GENE', 'GENENAME', 'ENSEMBL', 'ENTREZID', 'ID_REF', 'TARGETID']
        for col in label_columns:
            if col in expr_df.columns:
                labels = expr_df[col].astype(str).replace({'nan': pd.NA, 'None': pd.NA}).str.strip()
                if not labels.notna().any():
                    continue
                labels = labels.where(labels.notna(), expr_df.index.astype(str))
                if labels.eq(expr_df.index.astype(str)).all():
                    continue
                if labels.str.fullmatch(r'\d+').all():
                    continue
                return labels
        return expr_df.index.astype(str)
