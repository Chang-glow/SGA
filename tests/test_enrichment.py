"""EnrichStrategy._ensure_gene_symbols() 单元测试。"""
import os
import sys
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.strategies.enrichment import EnrichStrategy


class TestEnsureGeneSymbols:
    """验证 _ensure_gene_symbols() 使用 _ensembl_to_symbol 映射转换 Ensembl ID。"""

    def _make_analyzer(self, pack=None):
        """构建最小 mock analyzer。"""
        analyzer = MagicMock()
        analyzer._meta_matrix_pack = pack
        analyzer._logger = MagicMock()
        from utils.config_manager import Config
        analyzer.cfg = Config(
            tar_gene="POLB", gse_id="GSE12345", analysis_mode="enrich",
            debug=False, data_dir="/tmp",
        )
        return analyzer

    def test_converts_ensembl_using_pack_mapping(self):
        """pack 中有 _ensembl_to_symbol 映射时应转换 Gene 列。"""
        pack = {
            "_ensembl_to_symbol": {
                "ENSMUSG00000028842": "Acta2",
                "ENSMUSG00000037742": "Col1a1",
                "ENSMUSG00000096856": "Tnf",
            },
        }
        analyzer = self._make_analyzer(pack=pack)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame({
            "Gene": [
                "ENSMUSG00000028842",
                "ENSMUSG00000037742",
                "ENSMUSG00000096856",
            ],
            "log2FC": [1.5, -0.8, 2.1],
            "padj": [0.01, 0.02, 0.03],
        })

        result = strategy._ensure_gene_symbols(df)

        assert list(result["Gene"]) == ["Acta2", "Col1a1", "Tnf"]

    def test_falls_back_to_expr_matrix_symbol_column(self):
        """pack 无 _ensembl_to_symbol 但 expr_matrix 有 SYMBOL 列时回退。"""
        import copy
        expr_matrix = pd.DataFrame({
            "SYMBOL": ["Acta2", "Col1a1", "Tnf"],
            "s1": [1.0, 2.0, 3.0],
            "s2": [4.0, 5.0, 6.0],
        }, index=["ENSMUSG00000028842", "ENSMUSG00000037742", "ENSMUSG00000096856"])
        pack = {"expr_matrix": expr_matrix}
        analyzer = self._make_analyzer(pack=pack)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame({
            "Gene": [
                "ENSMUSG00000028842",
                "ENSMUSG00000037742",
                "ENSMUSG00000096856",
            ],
            "log2FC": [1.5, -0.8, 2.1],
            "padj": [0.01, 0.02, 0.03],
        })

        result = strategy._ensure_gene_symbols(df)

        assert list(result["Gene"]) == ["Acta2", "Col1a1", "Tnf"]

    def test_no_pack_no_change(self):
        """pack 为 None 时返回原 DataFrame（不抛异常）。"""
        analyzer = self._make_analyzer(pack=None)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame({
            "Gene": [
                "ENSMUSG00000028842",
                "ENSMUSG00000037742",
            ],
            "log2FC": [1.5, -0.8],
            "padj": [0.01, 0.02],
        })

        result = strategy._ensure_gene_symbols(df)

        assert list(result["Gene"]) == [
            "ENSMUSG00000028842",
            "ENSMUSG00000037742",
        ]

    def test_symbol_genes_not_converted(self):
        """基因符号为主的 Gene 列不触发转换。"""
        pack = {
            "_ensembl_to_symbol": {
                "ENSMUSG00000028842": "Acta2",
            },
        }
        analyzer = self._make_analyzer(pack=pack)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame({
            "Gene": ["Acta2", "Col1a1", "Tnf"],
            "log2FC": [1.5, -0.8, 2.1],
            "padj": [0.01, 0.02, 0.03],
        })

        result = strategy._ensure_gene_symbols(df)

        # 全是基因符号，不应触发转换
        assert list(result["Gene"]) == ["Acta2", "Col1a1", "Tnf"]

    def test_partial_conversion_keeps_unmapped(self):
        """映射中不存在的 Ensembl ID 保持原值。"""
        pack = {
            "_ensembl_to_symbol": {
                "ENSMUSG00000028842": "Acta2",
            },
        }
        analyzer = self._make_analyzer(pack=pack)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame({
            "Gene": [
                "ENSMUSG00000028842",
                "ENSMUSG99999999999",
            ],
            "log2FC": [1.5, -0.8],
            "padj": [0.01, 0.02],
        })

        result = strategy._ensure_gene_symbols(df)

        assert list(result["Gene"]) == ["Acta2", "ENSMUSG99999999999"]

    def test_empty_df_pass_through(self):
        """空 DataFrame 不抛异常。"""
        pack = {
            "_ensembl_to_symbol": {
                "ENSMUSG00000028842": "Acta2",
            },
        }
        analyzer = self._make_analyzer(pack=pack)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame()

        result = strategy._ensure_gene_symbols(df)

        assert result.empty

    def test_no_gene_column_pass_through(self):
        """无 Gene 列的 DataFrame 不抛异常。"""
        pack = {
            "_ensembl_to_symbol": {
                "ENSMUSG00000028842": "Acta2",
            },
        }
        analyzer = self._make_analyzer(pack=pack)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame({"log2FC": [1.5], "padj": [0.01]})

        result = strategy._ensure_gene_symbols(df)

        assert list(result.columns) == ["log2FC", "padj"]

    def test_ensembl_mapping_with_expr_matrix_symbol_as_index(self):
        """复现 bug: expr_matrix 的 SYMBOL 不在 columns 中(被 set_index 移走)
        _ensembl_to_symbol 映射仍能正常工作。"""
        # 模拟 data_packer 将 SYMBOL set_index 后的 expr_matrix
        # 此时 expr_matrix.columns 中没有 "SYMBOL"
        expr_matrix = pd.DataFrame(
            np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]),
            index=["Acta2", "Col1a1", "Tnf"],  # SYMBOL 已是 index
            columns=["s1", "s2"],
        )

        pack = {
            "expr_matrix": expr_matrix,
            "_ensembl_to_symbol": {
                "ENSMUSG00000028842": "Acta2",
                "ENSMUSG00000037742": "Col1a1",
                "ENSMUSG00000096856": "Tnf",
            },
        }
        analyzer = self._make_analyzer(pack=pack)
        strategy = EnrichStrategy(analyzer)

        df = pd.DataFrame({
            "Gene": [
                "ENSMUSG00000028842",
                "ENSMUSG00000037742",
                "ENSMUSG00000096856",
            ],
            "log2FC": [1.5, -0.8, 2.1],
            "padj": [0.01, 0.02, 0.03],
        })

        result = strategy._ensure_gene_symbols(df)

        # 应通过 _ensembl_to_symbol 转换（而不是走旧 expr_matrix SYMBOL 列路径）
        assert list(result["Gene"]) == ["Acta2", "Col1a1", "Tnf"]
