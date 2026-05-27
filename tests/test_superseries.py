"""SuperSeries 合并 / Ensembl 映射 / diff & hilo 适配 单元测试。"""
import os
import sys
import pytest
import pandas as pd
import numpy as np
from unittest.mock import MagicMock, patch, PropertyMock

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from utils.config_manager import Config
from modules.data_packer import DataPacker


# ── helpers ──────────────────────────────────────────────

def _make_cfg(**overrides) -> Config:
    defaults = dict(
        tar_gene="POLB", gse_id="GSE12345", analysis_mode="diff",
        debug=False, group_select_col="source_name_ch1",
        control_label=["control"], exp_label=["CCl4"], exp_type="Fibrosis",
        group_memory_enabled=False, group_memory_use=False, data_dir="/tmp",
    )
    defaults.update(overrides)
    return Config(**defaults)


def _make_expr_df(index_values, columns, data=1.0):
    """构建最小表达矩阵，每行值为 data。"""
    n_rows = len(index_values)
    n_cols = len(columns)
    return pd.DataFrame(
        np.full((n_rows, n_cols), data, dtype=float),
        index=index_values,
        columns=columns,
    )


# ── normalize_gene_index ─────────────────────────────────

class TestNormalizeGeneIndex:
    """验证 Ensembl 版本号剥离和内联去重。"""

    def test_strip_version_suffix(self):
        from modules.calculater import normalize_gene_index
        df = _make_expr_df(
            ["ENSMUSG00000028842.15", "ENSMUSG00000037742", "ENSMUSG00000096856.7"],
            ["s1", "s2"],
        )
        result = normalize_gene_index(df)
        assert list(result.index) == ["ENSMUSG00000028842", "ENSMUSG00000037742", "ENSMUSG00000096856"]

    def test_no_version_no_change(self):
        from modules.calculater import normalize_gene_index
        df = _make_expr_df(["ENSMUSG00000028842", "ENSMUSG00000037742"], ["s1"])
        result = normalize_gene_index(df)
        assert list(result.index) == list(df.index)

    def test_dedup_after_strip_keeps_higher_mean(self):
        from modules.calculater import normalize_gene_index
        df = _make_expr_df(
            ["ENSMUSG00000028842.15", "ENSMUSG00000028842.7"],
            ["s1", "s2"],
        )
        df.iloc[0, :] = 100.0  # higher expr
        df.iloc[1, :] = 10.0
        result = normalize_gene_index(df)
        assert list(result.index) == ["ENSMUSG00000028842"]
        assert result.iloc[0, 0] == 100.0

    def test_no_duplicates_after_strip(self):
        from modules.calculater import normalize_gene_index
        df = _make_expr_df(["ENSMUSG00000028842.15", "ENSMUSG00000037742.3"], ["s1"])
        result = normalize_gene_index(df)
        assert len(result) == 2


# ── map_ensembl_to_symbol ────────────────────────────────

class TestMapEnsemblToSymbol:
    """验证 Ensembl 索引检测和 mygene 映射。"""

    def test_non_ensembl_index_passthrough(self):
        from modules.calculater import map_ensembl_to_symbol
        df = _make_expr_df(["Acta2", "Col1a1", "Tnf"], ["s1"])
        result = map_ensembl_to_symbol(df)
        assert "SYMBOL" not in result.columns
        assert list(result.index) == list(df.index)

    def test_mixed_index_below_threshold_passthrough(self):
        from modules.calculater import map_ensembl_to_symbol
        # 只有 2 个 Ensembl 格式（18 个非 Ensembl），不到 80% 阈值
        idx = ["Acta2"] * 18 + ["ENSMUSG99999999999", "ENSMUSG99999999998"]
        df = _make_expr_df(idx, ["s1"])
        result = map_ensembl_to_symbol(df)
        assert "SYMBOL" not in result.columns

    def test_ensembl_index_adds_symbol_column(self):
        from modules.calculater import map_ensembl_to_symbol
        idx = ["ENSMUSG00000051951"]  # known: Xkr4
        df = _make_expr_df(idx, ["s1"])
        result = map_ensembl_to_symbol(df)
        assert "SYMBOL" in result.columns
        # 应映射为基因符号或回退为 Ensembl ID
        assert result["SYMBOL"].iloc[0] in ("Xkr4", "ENSMUSG00000051951")

    def test_version_suffix_stripped_before_query(self):
        from modules.calculater import map_ensembl_to_symbol
        df = _make_expr_df(["ENSMUSG00000051951.5"], ["s1"])
        result = map_ensembl_to_symbol(df)
        assert "SYMBOL" in result.columns
        # 版本后缀应被剥离，查询应成功
        assert result["SYMBOL"].iloc[0] == "Xkr4"

    def test_unknown_ensembl_fills_with_original_id(self):
        from modules.calculater import map_ensembl_to_symbol
        # 混合已知和未知 Ensembl ID，未知的应 fallback 为原始 ID
        df = _make_expr_df(["ENSMUSG00000051951", "ENSMUSG99999999999"], ["s1"])
        result = map_ensembl_to_symbol(df)
        assert result["SYMBOL"].iloc[0] == "Xkr4"
        assert result["SYMBOL"].iloc[1] == "ENSMUSG99999999999"


# ── DiffStrategy SuperSeries merge ───────────────────────

class TestDiffSuperseriesMerge:
    """验证 diff 策略的 SuperSeries 多矩阵合并。"""

    def _make_pack(self, expr_matrices, meta, meta_full, group_info=True):
        """构建模拟 pack。"""
        pack = {"meta": meta, "meta_full": meta_full}
        pack.update(expr_matrices)
        if group_info:
            pack["group_info"] = {"group_col": "group", "mapping": {
                "Control": ["GSM01", "GSM02"],
                "Fibrosis": ["GSM03", "GSM04", "GSM05", "GSM06"],
            }}
        return pack

    def test_single_matrix_no_superseries_path(self):
        """单矩阵不走 SuperSeries 路径。"""
        from modules.strategies.difference import DiffStrategy
        analyzer = MagicMock()
        analyzer.cfg = _make_cfg(min_samples_per_group=1)
        analyzer._logger = MagicMock()
        analyzer._meta_matrix_pack = None

        meta = pd.DataFrame(
            {"group": ["Control", "Control", "Fibrosis", "Fibrosis"]},
            index=["GSM01", "GSM02", "GSM03", "GSM04"],
        )
        expr = _make_expr_df(["GeneA", "GeneB"], ["GSM01", "GSM02", "GSM03", "GSM04"])
        pack = {"meta": meta, "key1.txt": expr}
        analyzer._meta_matrix_pack = pack
        analyzer.roaming_data = MagicMock(return_value=pack)

        strategy = DiffStrategy(analyzer)
        strategy._align_expression_matrix = MagicMock(return_value=expr)
        strategy._prefilter_expression_matrix = MagicMock(return_value=expr)
        strategy._calculate_diff = MagicMock(return_value=pd.DataFrame({"Gene": ["GeneA"]}))

        result = strategy.calculate()
        assert len(result) == 1

    def test_superseries_merge_called(self):
        """有 / 的 expr key 触发 SuperSeries 合并。"""
        from modules.strategies.difference import DiffStrategy
        analyzer = MagicMock()
        analyzer.cfg = _make_cfg(min_samples_per_group=1)
        analyzer._logger = MagicMock()
        analyzer._meta_matrix_pack = None

        meta = pd.DataFrame(
            {"group": ["Control", "Control", "Fibrosis", "Fibrosis", "Fibrosis", "Fibrosis"]},
            index=["GSM01", "GSM02", "GSM03", "GSM04", "GSM05", "GSM06"],
        )
        meta_full = pd.DataFrame(
            {"series_id": ["GSE1,GSEparent", "GSE1,GSEparent", "GSE2,GSEparent",
                           "GSE1,GSEparent", "GSE1,GSEparent", "GSE2,GSEparent"]},
            index=["GSM01", "GSM02", "GSM04", "GSM03", "GSM05", "GSM06"],
        )
        expr1 = _make_expr_df(["GeneA", "GeneB", "GeneC"], ["GSM01", "GSM02", "GSM03", "GSM04"])
        expr2 = _make_expr_df(["GeneB", "GeneC", "GeneD"], ["GSM05", "GSM06"])

        pack = {
            "meta": meta,
            "meta_full": meta_full,
            "GSE1/GSE1_counts.txt": expr1,
            "GSE2/GSE2_counts.txt": expr2,
        }
        analyzer._meta_matrix_pack = pack
        analyzer.roaming_data = MagicMock(return_value=pack)

        strategy = DiffStrategy(analyzer)
        strategy._align_expression_matrix = MagicMock(side_effect=lambda df, samples, pack, **kw: df)
        strategy._prefilter_expression_matrix = MagicMock(side_effect=lambda df: df)
        strategy._set_gene_index = MagicMock(side_effect=lambda df: df)
        strategy._normalize_gene_index = MagicMock(side_effect=lambda df: df)
        strategy._calculate_diff = MagicMock(return_value=pd.DataFrame({"Gene": ["GeneB"]}))

        result = strategy.calculate()
        # _align_expression_matrix 应被调用两次（每个子系列一次）
        assert strategy._align_expression_matrix.call_count == 2
        assert len(result) == 1

    def test_batch_exp_groups_excluded_from_expr_keys(self):
        """_batch_exp_groups 不会被当作表达矩阵。"""
        from modules.strategies.difference import DiffStrategy
        analyzer = MagicMock()
        analyzer.cfg = _make_cfg()
        analyzer._logger = MagicMock()
        analyzer._meta_matrix_pack = None

        meta = pd.DataFrame(
            {"group": ["Control", "Fibrosis"]},
            index=["GSM01", "GSM02"],
        )
        expr = _make_expr_df(["GeneA"], ["GSM01", "GSM02"])
        pack = {
            "meta": meta,
            "_batch_exp_groups": [("label", "value")],
            "key1.txt": expr,
        }
        analyzer._meta_matrix_pack = pack
        analyzer.roaming_data = MagicMock(return_value=pack)

        strategy = DiffStrategy(analyzer)
        strategy._align_expression_matrix = MagicMock(return_value=expr)
        strategy._prefilter_expression_matrix = MagicMock(return_value=expr)
        strategy._calculate_diff = MagicMock(return_value=pd.DataFrame())

        result = strategy.calculate()
        # _batch_exp_groups 不应被当作表达矩阵处理
        assert strategy._align_expression_matrix.call_count == 1

    def test_organism_and_ensembl_to_symbol_excluded_from_expr_keys(self):
        """_organism 和 _ensembl_to_symbol 不会被当作表达矩阵。"""
        from modules.strategies.difference import DiffStrategy
        analyzer = MagicMock()
        analyzer.cfg = _make_cfg()
        analyzer._logger = MagicMock()
        analyzer._meta_matrix_pack = None

        meta = pd.DataFrame(
            {"group": ["Control", "Fibrosis"]},
            index=["GSM01", "GSM02"],
        )
        expr = _make_expr_df(["GeneA"], ["GSM01", "GSM02"])
        pack = {
            "meta": meta,
            "_organism": "mouse",
            "_ensembl_to_symbol": {"ENSMUSG00000051951": "Xkr4"},
            "key1.txt": expr,
        }
        analyzer._meta_matrix_pack = pack
        analyzer.roaming_data = MagicMock(return_value=pack)

        strategy = DiffStrategy(analyzer)
        strategy._align_expression_matrix = MagicMock(return_value=expr)
        strategy._prefilter_expression_matrix = MagicMock(return_value=expr)
        strategy._calculate_diff = MagicMock(return_value=pd.DataFrame())

        strategy.calculate()
        # _organism 和 _ensembl_to_symbol 不应被当作表达矩阵处理
        # 只有 key1.txt 是真正的表达矩阵，应被传入 _align_expression_matrix
        assert strategy._align_expression_matrix.call_count == 1
        # 验证传入的是真正的表达矩阵（DataFrame），而不是 _organism 字符串
        call_df = strategy._align_expression_matrix.call_args[0][0]
        assert isinstance(call_df, pd.DataFrame), (
            f"_align_expression_matrix 应收到 DataFrame，但收到 {type(call_df).__name__}"
        )

    def test_control_samples_filtered_after_merge(self):
        """合并后 control/exp samples 过滤到矩阵中实际存在的列。"""
        from modules.strategies.difference import DiffStrategy
        analyzer = MagicMock()
        analyzer.cfg = _make_cfg(min_samples_per_group=1)
        analyzer._logger = MagicMock()
        analyzer._meta_matrix_pack = None

        # meta 有 8 个样本，但两个子矩阵合并后只有 6 列
        meta = pd.DataFrame(
            {"group": ["Control", "Control", "Control", "Fibrosis",
                       "Fibrosis", "Fibrosis", "Fibrosis", "Fibrosis"]},
            index=["GSM01", "GSM02", "GSM_MISSING", "GSM03",
                   "GSM04", "GSM05", "GSM06", "GSM_MISSING2"],
        )
        meta_full = pd.DataFrame(
            {"series_id": ["GSE1,GSEp", "GSE1,GSEp", "GSE1,GSEp",
                           "GSE1,GSEp", "GSE1,GSEp", "GSE2,GSEp",
                           "GSE2,GSEp", "GSE2,GSEp"]},
            index=["GSM01", "GSM02", "GSM_MISSING", "GSM03",
                   "GSM04", "GSM05", "GSM06", "GSM_MISSING2"],
        )
        expr1 = _make_expr_df(["GeneA"], ["GSM01", "GSM02", "GSM03", "GSM04"])
        expr2 = _make_expr_df(["GeneA"], ["GSM05", "GSM06"])

        pack = {
            "meta": meta,
            "meta_full": meta_full,
            "GSE1/GSE1.txt": expr1,
            "GSE2/GSE2.txt": expr2,
        }
        analyzer._meta_matrix_pack = pack
        analyzer.roaming_data = MagicMock(return_value=pack)

        strategy = DiffStrategy(analyzer)
        strategy._align_expression_matrix = MagicMock(side_effect=lambda df, samples, pack, **kw: df)
        strategy._prefilter_expression_matrix = MagicMock(side_effect=lambda df: df)
        strategy._set_gene_index = MagicMock(side_effect=lambda df: df)
        strategy._normalize_gene_index = MagicMock(side_effect=lambda df: df)
        strategy._calculate_diff = MagicMock(return_value=pd.DataFrame())

        strategy.calculate()
        call_args = strategy._calculate_diff.call_args
        ctrl = call_args[0][1]
        exp = call_args[0][2]
        # 缺失样本应已被过滤掉
        assert "GSM_MISSING" not in ctrl
        assert "GSM_MISSING2" not in exp


# ── DataPacker._merge_superseries_matrices (hilo) ────────

class TestHiloSuperseriesMerge:
    """验证 hilo 分组的 SuperSeries 合并。"""

    def test_merge_with_version_suffix(self):
        """不同子系列的 Ensembl 版本号应被剥离后合并。"""
        cfg = _make_cfg(analysis_mode="hilo", tar_gene="GeneB")
        # mock _fetch_target_gene_vector 以避免 mygene 网络调用
        meta = pd.DataFrame({
            "series_id": ["GSE1,GSEp", "GSE1,GSEp", "GSE2,GSEp"],
            "title": ["a", "b", "c"],
        }, index=["GSM01", "GSM02", "GSM03"])
        gse = MagicMock()
        gse.phenotype_data = meta

        # 两个子系列，带版本号的 Ensembl ID
        downloaded = {
            "GSE1/GSE1.txt": "/tmp/gse1.txt",
            "GSE2/GSE2.txt": "/tmp/gse2.txt",
        }

        packer = DataPacker(cfg, gse, downloaded)

        # mock 文件加载
        df1 = _make_expr_df(
            ["ENSMUSG00000028842.15", "ENSMUSG00000037742"],
            ["c1", "c2"],
        )
        df2 = _make_expr_df(
            ["ENSMUSG00000028842.3", "ENSMUSG00000037742"],
            ["c3"],
        )

        with patch.object(packer, '_load_expression_file', side_effect=[df1, df2]), \
             patch.object(packer, '_fetch_target_gene_vector', return_value=pd.Series([1.0, 2.0], index=["GSM01", "GSM03"])), \
             patch.object(packer, '_apply_hilo_threshold'):
            packer._prepare_hilo_group(meta)
            # 去版本号后合并应找到共同基因，不抛异常

    def test_non_superseries_uses_single_file(self):
        """无 / 的 key 走原单文件路径。"""
        cfg = _make_cfg(analysis_mode="hilo", tar_gene="POLB")
        meta = pd.DataFrame(
            {"source_name_ch1": ["control", "control", "CCl4", "CCl4"]},
            index=["GSM01", "GSM02", "GSM03", "GSM04"],
        )
        gse = MagicMock()
        gse.phenotype_data = meta

        downloaded = {"single.txt": "/tmp/single.txt"}
        packer = DataPacker(cfg, gse, downloaded)

        expr_df = _make_expr_df(["POLB", "Acta2"], ["GSM01", "GSM02", "GSM03", "GSM04"])

        with patch.object(packer, '_load_expression_file', return_value=expr_df), \
             patch.object(packer, '_select_hilo_expression_file', return_value="/tmp/single.txt"), \
             patch.object(packer, '_apply_hilo_threshold'):
            packer._prepare_hilo_group(meta)
            # 不应抛异常

    def test_superseries_detected_when_keys_contain_slash(self):
        """有 / 的 key + series_id 列 → 走 SuperSeries 路径。"""
        cfg = _make_cfg(analysis_mode="hilo", tar_gene="GeneX")
        meta = pd.DataFrame({
            "series_id": ["GSE1,GSEp", "GSE1,GSEp"],
            "title": ["a", "b"],
        }, index=["GSM01", "GSM02"])
        gse = MagicMock()
        gse.phenotype_data = meta

        downloaded = {"GSE1/sub.txt": "/tmp/sub.txt"}
        packer = DataPacker(cfg, gse, downloaded)

        merged_mock = _make_expr_df(["GeneX"], ["GSM01", "GSM02"])

        with patch.object(packer, '_merge_superseries_matrices', return_value=merged_mock) as mock_merge, \
             patch.object(packer, '_fetch_target_gene_vector', return_value=pd.Series([1.0, 2.0], index=["GSM01", "GSM02"])), \
             patch.object(packer, '_apply_hilo_threshold'):
            packer._prepare_hilo_group(meta)
            mock_merge.assert_called_once()


class TestBoxplotPartialOverlap:
    """箱线图 _prepare_diff_data: expr 列与 meta 索引部分重叠时不抛异常且产出有效 p 值。"""

    def test_partial_overlap_produces_valid_pvalue(self):
        """复现 GSE329774 场景：部分样本不在 meta.index 中，剩余样本正常计算。"""
        from modules.fig_plotter import FigurePlotter
        from utils.loggers import get_logger

        FigurePlotter.__abstractmethods__ = set()  # bypass ABC for testing _prepare_diff_data

        cfg = _make_cfg(analysis_mode="diff", tar_gene="GeneX", multi_gene=None)
        meta = pd.DataFrame({
            "group": ["Control", "Control", "Fibrosis", "Fibrosis"],
        }, index=["GSM01", "GSM02", "GSM03", "GSM04"])

        # 表达矩阵：6 个样本列，其中 GSM05/GSM06 不在 meta 中
        expr_df = pd.DataFrame(
            {"SYMBOL": ["GeneX"],
             "GSM01": [2.0], "GSM02": [2.2], "GSM03": [5.0], "GSM04": [5.5],
             "GSM05": [3.0], "GSM06": [4.0]},
        )

        plotter = FigurePlotter(cfg)
        plotter._meta_matrix_pack = {"meta": meta}
        plotter._logger = get_logger("test")
        # patch _get_expr_matrix to return our controlled data
        import types
        plotter._get_expr_matrix = types.MethodType(lambda self: expr_df, plotter)

        result = plotter._prepare_diff_data(gene="GeneX")
        assert result is not None, "部分重叠时不应返回 None"
        assert not np.isnan(result["p_value"]), f"p_value 不应为 NaN"
        # 应有 2 Control + 2 Fibrosis = 4 样本参与 t-test
        x = result["x"]
        assert (x == "Control").sum() == 2, f"Control 应为 2，实际: {(x == 'Control').sum()}"
        assert (x == "Fibrosis").sum() == 2, f"Fibrosis 应为 2，实际: {(x == 'Fibrosis').sum()}"

    def test_all_samples_in_meta_no_warning(self):
        """全部样本在 meta 中 → 不走交集修剪，行为不变。"""
        from modules.fig_plotter import FigurePlotter
        from utils.loggers import get_logger

        FigurePlotter.__abstractmethods__ = set()

        cfg = _make_cfg(analysis_mode="diff", tar_gene="GeneX", multi_gene=None)
        meta = pd.DataFrame({
            "group": ["Control", "Control", "Fibrosis", "Fibrosis"],
        }, index=["GSM01", "GSM02", "GSM03", "GSM04"])

        expr_df = pd.DataFrame(
            {"SYMBOL": ["GeneX"],
             "GSM01": [1.0], "GSM02": [1.2], "GSM03": [4.0], "GSM04": [4.5]},
        )

        plotter = FigurePlotter(cfg)
        plotter._meta_matrix_pack = {"meta": meta}
        plotter._logger = get_logger("test")
        import types
        plotter._get_expr_matrix = types.MethodType(lambda self: expr_df, plotter)

        result = plotter._prepare_diff_data(gene="GeneX")
        assert result is not None
        assert not np.isnan(result["p_value"])
        assert (result["x"] == "Control").sum() == 2
        assert (result["x"] == "Fibrosis").sum() == 2


class TestModeMeta:
    """验证全局模式映射表覆盖所有分析模式，消除 if/elif 重复。"""

    def test_mode_meta_covers_all_handled_modes(self):
        """_MODE_META 应包含 corr/diff/hilo/enrich/immune 五个模式。"""
        from modules.fig_plotter import FigurePlotter

        assert hasattr(FigurePlotter, "_MODE_META"), (
            "FigurePlotter 应有 _MODE_META 类属性"
        )
        meta = FigurePlotter._MODE_META
        for mode in ["corr", "diff", "hilo", "enrich", "immune"]:
            assert mode in meta, f"_MODE_META 缺少模式: {mode}"
            entry = meta[mode]
            assert "attr" in entry, f"{mode} 缺少 'attr' 字段"
            assert "summary_suffix" in entry, f"{mode} 缺少 'summary_suffix' 字段"

    def test_diff_and_hilo_share_same_table_attr(self):
        """diff 和 hilo 应共享 _gene_diff_table（hilo 内部复用 DiffStrategy）。"""
        from modules.fig_plotter import FigurePlotter

        if not hasattr(FigurePlotter, "_MODE_META"):
            pytest.skip("_MODE_META 尚未定义")
        meta = FigurePlotter._MODE_META
        assert meta["diff"]["attr"] == meta["hilo"]["attr"], (
            f"diff 和 hilo 应共享同一 table attr，"
            f"实际 diff={meta['diff']['attr']}, hilo={meta['hilo']['attr']}"
        )


class TestDataPlotterHilo:
    """DataPlotter._load_data() 不应拒绝 hilo 模式。"""

    def test_load_data_handles_hilo(self):
        """hilo 模式走 DataPlotter 内存路径不抛 ValueError。"""
        from modules.fig_plotter import DataPlotter
        from utils.config_manager import Config, DataHandler

        cfg = Config(
            tar_gene="POLB", gse_id="GSE12345", analysis_mode="hilo",
            debug=False, data_dir="/tmp",
        )
        data = DataHandler()
        data.gene_diff_table = pd.DataFrame({
            "Gene": ["A"], "log2FC": [1.0], "P_value": [0.01], "padj": [0.05],
        })

        plotter = DataPlotter(cfg, data)
        try:
            plotter._load_data()
        except ValueError as e:
            if "未知分析模式" in str(e):
                pytest.fail(f"DataPlotter._load_data() 不应拒绝 hilo 模式: {e}")
            raise
        assert plotter._gene_diff_table is not None, "hilo 数据应加载到 _gene_diff_table"


class TestCreateFallbackToDataPlotter:
    """FigurePlotter.create() 在 PKL 缺失时应回退到 DataPlotter。"""

    def test_missing_table_pkl_falls_back_to_data_plotter(self):
        """pack 存在但 table PKL 不存在 -> 返回 DataPlotter 使用内存数据。"""
        import tempfile
        from modules.fig_plotter import FigurePlotter
        from utils.config_manager import Config, DataHandler

        with tempfile.TemporaryDirectory() as tmpdir:
            gse_dir = os.path.join(tmpdir, "GSE12345")
            pkl_dir = os.path.join(gse_dir, "pkl")
            os.makedirs(pkl_dir, exist_ok=True)
            # 创建 pack 文件（存在）
            pack_path = os.path.join(
                pkl_dir, "GSE12345_diff_group_processed_pack.pkl"
            )
            pd.to_pickle({"meta": pd.DataFrame()}, pack_path)
            # table PKL 不存在

            cfg = Config(
                tar_gene="POLB", gse_id="GSE12345", analysis_mode="diff",
                debug=False, data_dir=tmpdir,
            )
            data = DataHandler()
            data.meta_matrix_pack = {"meta": pd.DataFrame()}
            data.gene_diff_table = pd.DataFrame({"Gene": ["A"]})

            plotter = FigurePlotter.create(cfg, data)
            from modules.fig_plotter import DataPlotter
            assert isinstance(plotter, DataPlotter), (
                f"PKL 缺失时应回退到 DataPlotter，实际: {type(plotter).__name__}"
            )

    def test_both_missing_returns_data_plotter(self):
        """pack 和 table 都不存在 -> 返回 DataPlotter（原有行为不变）。"""
        from modules.fig_plotter import FigurePlotter
        from utils.config_manager import Config, DataHandler

        cfg = Config(
            tar_gene="POLB", gse_id="GSE99999", analysis_mode="diff",
            debug=False, data_dir="/nonexistent/path",
        )
        data = DataHandler()
        plotter = FigurePlotter.create(cfg, data)
        from modules.fig_plotter import DataPlotter
        assert isinstance(plotter, DataPlotter), (
            f"两者都缺失时应返回 DataPlotter，实际: {type(plotter).__name__}"
        )


class TestSkipCheckHasHilo:
    """main.py _SKIP_CHECK 应包含 hilo 模式。"""

    def test_skip_check_contains_hilo(self):
        """验证 _SKIP_CHECK 有 "hilo" 键。"""
        from main import _SKIP_CHECK
        assert "hilo" in _SKIP_CHECK, (
            f"_SKIP_CHECK 缺少 hilo 模式，当前: {list(_SKIP_CHECK.keys())}"
        )


class TestSummaryValidGenes:
    """_summary_has_valid_genes() 应拒绝 Ensembl ID 为主的 Gene 列。"""

    def test_ensembl_genes_are_invalid(self):
        from modules.strategies.difference import DiffStrategy
        from unittest.mock import MagicMock
        strategy = DiffStrategy(MagicMock())
        df = pd.DataFrame({
            "Gene": ["ENSMUSG00000028238", "ENSMUSG00000044534", "Polb"],
            "log2FC": [1.0, -0.5, 2.0],
            "P_value": [0.01, 0.02, 0.03],
            "padj": [0.04, 0.05, 0.06],
        })
        # 66% Ensembl ID > 50% → 应判定无效
        assert not strategy._summary_has_valid_genes(df), (
            "Ensembl ID 占比 > 50% 时应判定为无效"
        )

    def test_symbol_genes_are_valid(self):
        from modules.strategies.difference import DiffStrategy
        from unittest.mock import MagicMock
        strategy = DiffStrategy(MagicMock())
        df = pd.DataFrame({
            "Gene": ["Polb", "Tp53", "Il6"],
            "log2FC": [1.0, -0.5, 2.0],
            "P_value": [0.01, 0.02, 0.03],
            "padj": [0.04, 0.05, 0.06],
        })
        assert strategy._summary_has_valid_genes(df), (
            "基因符号应判定为有效"
        )

    def test_empty_gene_column_is_invalid(self):
        from modules.strategies.difference import DiffStrategy
        from unittest.mock import MagicMock
        strategy = DiffStrategy(MagicMock())
        df = pd.DataFrame({
            "Gene": [],
            "log2FC": [],
            "P_value": [],
            "padj": [],
        })
        assert not strategy._summary_has_valid_genes(df)

    def test_missing_gene_column_is_invalid(self):
        from modules.strategies.difference import DiffStrategy
        from unittest.mock import MagicMock
        strategy = DiffStrategy(MagicMock())
        df = pd.DataFrame({
            "log2FC": [1.0],
            "P_value": [0.01],
            "padj": [0.04],
        })
        assert not strategy._summary_has_valid_genes(df)

    def test_numeric_genes_are_invalid(self):
        from modules.strategies.difference import DiffStrategy
        from unittest.mock import MagicMock
        strategy = DiffStrategy(MagicMock())
        df = pd.DataFrame({
            "Gene": ["12345", "67890"],
            "log2FC": [1.0, -0.5],
            "P_value": [0.01, 0.02],
            "padj": [0.04, 0.05],
        })
        assert not strategy._summary_has_valid_genes(df)
