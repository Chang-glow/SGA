"""CorrelationStrategy 单元测试。"""
import os
import sys
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from modules.strategies.correlation import CorrelationStrategy


class TestCorrelationPackIteration:
    """验证 calculate() 在遍历 pack dict 时跳过非 DataFrame 值。"""

    def test_skips_non_dataframe_pack_values(self):
        """_organism, group_info 等非 DataFrame 值不应传递到 fetch_gene_vector。"""
        tar_gene = "POLB"

        # 构建一个包含非 DataFrame 值的 minimal pack
        expr_df = pd.DataFrame(
            np.random.randn(5, 3),
            index=[tar_gene, "GENE2", "GENE3", "GENE4", "GENE5"],
            columns=["GSM1", "GSM2", "GSM3"],
        )
        pack = {
            "meta": pd.DataFrame({"group": ["Control", "Control", "Experiment"]}),
            "_organism": "mouse",
            "group_info": {"group_col": "col", "mapping": {}},
            "_batch_exp_groups": [("MCD", "MCD 1")],
            "_ensembl_to_symbol": {"ENSG1": "GENE1"},
            "expr_data": expr_df,
        }

        # Mock analyzer — _get_hfm_dict 返回空 dict 避免遍历 marker genes
        analyzer = MagicMock()
        analyzer.cfg.tar_gene = tar_gene
        analyzer._meta_matrix_pack = pack
        analyzer._get_hfm_dict.return_value = {}

        strategy = CorrelationStrategy(analyzer)

        with patch("modules.strategies.correlation.fetch_gene_vector") as mock_fgv:
            mock_fgv.return_value = pd.Series(
                [1.0, 2.0, 3.0], index=["GSM1", "GSM2", "GSM3"]
            )
            strategy.calculate()

        # 断言：所有传给 fetch_gene_vector 的第一个参数都是 DataFrame
        non_df_calls = []
        for call in mock_fgv.call_args_list:
            df_arg = call[0][0]
            if not isinstance(df_arg, pd.DataFrame):
                non_df_calls.append((type(df_arg), df_arg))
        assert (
            not non_df_calls
        ), f"fetch_gene_vector called with non-DataFrame args: {non_df_calls}"

        # 确认至少为 expr_data 调用了一次
        assert mock_fgv.call_count >= 1, (
            "fetch_gene_vector 至少应该为 expression DataFrame 调用一次"
        )
