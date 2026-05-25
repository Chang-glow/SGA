import pandas as pd

from modules.calculater import prepare_expr_matrix


class ImmuneStrategy:
    """免疫浸润分析策略，使用 TumorDecon 进行去卷积"""

    def __init__(self, analyzer):
        self.analyzer = analyzer
        self._logger = analyzer._logger
        self.cfg = analyzer.cfg

    def calculate(self) -> pd.DataFrame:
        """执行免疫浸润去卷积分析

        Returns:
            cell frequency DataFrame (samples × cell types, values 0~1)
        """
        pack = self.analyzer.roaming_data()

        # 查找表达矩阵
        expr_df = None
        for key in pack:
            if key in ("meta", "meta_full", "group_info"):
                continue
            item = pack[key]
            if isinstance(item, dict) and "matrix_aligned" in item:
                expr_df = item["matrix_aligned"]
            elif isinstance(item, pd.DataFrame):
                expr_df = item
            if expr_df is not None:
                break

        if expr_df is None:
            raise ValueError("未在数据包中找到表达矩阵，无法进行免疫浸润分析")

        method = getattr(self.cfg, "immune_method", "DeconRNASeq")
        method_lower = method.lower()

        # 将基因名列设为索引（TumorDecon 要求 Hugo_Symbol 索引）
        expr_df = prepare_expr_matrix(expr_df)
        expr_df.index.name = "Hugo_Symbol"

        self._logger.info(
            f"使用 {method} 进行免疫浸润分析"
            f"（{expr_df.shape[0]} 基因 × {expr_df.shape[1]} 样本）"
        )

        # 根据方法加载所需数据
        sig_matrix = None
        up_genes = None
        down_genes = None

        import TumorDecon as td  # 懒加载：仅免疫浸润分析时需要，减少模块导入时的启动耗时

        if method_lower in ("cibersort", "deconrnaseq"):
            self._logger.info("加载 LM22 签名矩阵...")
            sig_matrix = td.read_sig_file()
            self._logger.info(
                f"LM22 签名矩阵: {sig_matrix.shape[0]} 基因 × {sig_matrix.shape[1]} 细胞类型"
            )
            if method_lower == "cibersort":
                self._logger.info(
                    "CIBERSORT 使用 NuSVR，每个样本约需 10-30 秒，请耐心等待..."
                )
        elif method_lower in ("ssgsea", "singscore"):
            self._logger.info("加载 LM22 基因集...")
            up_genes = td.read_geneset()

        # 记录基因重叠情况，方便排查问题
        if sig_matrix is not None:
            overlap = set(expr_df.index) & set(sig_matrix.index)
            self._logger.info(
                f"基因重叠: {len(overlap)}/{len(sig_matrix.index)} (签名矩阵) "
                f"vs {expr_df.shape[0]} (表达矩阵)"
            )

        self._logger.info("正在去卷积...")
        result = td.tumor_deconvolve(
            expr_df,
            method=method,
            sig_matrix=sig_matrix,
            up_genes=up_genes,
            down_genes=down_genes,
        )
        result.index.name = "Sample"
        self._logger.info(f"免疫浸润分析完成，共 {result.shape[1]} 种免疫细胞类型")
        return result
