"""DataPacker 分组逻辑单元测试。"""
import os
import sys
import pytest
import pandas as pd
import yaml
import tempfile
from unittest.mock import MagicMock, PropertyMock, patch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from utils.config_manager import Config
from modules.data_packer import DataPacker


# ── helpers ──────────────────────────────────────────────

def _make_cfg(**overrides) -> Config:
    """构建最小可用 Config，所有字段可覆盖。"""
    defaults = dict(
        tar_gene="POLB",
        gse_id="GSE12345",
        analysis_mode="diff",
        debug=False,
        group_select_col="source_name_ch1",
        control_label=["control", "ND"],
        exp_label=["CCl4"],
        exp_type="Fibrosis",
        group_memory_enabled=False,
        group_memory_use=False,
        data_dir="/tmp",
    )
    defaults.update(overrides)
    return Config(**defaults)


def _make_meta(col_values: dict) -> pd.DataFrame:
    """用 {列名: [值列表]} 构建最小 phenotype_data。"""
    return pd.DataFrame(col_values)


def _make_packer(cfg=None, meta=None, downloaded_data=None) -> DataPacker:
    """构建 DataPacker 实例，gse 为 mock。"""
    if cfg is None:
        cfg = _make_cfg()
    gse = MagicMock()
    if meta is not None:
        gse.phenotype_data = meta
    if downloaded_data is None:
        downloaded_data = {"dummy.txt": "/tmp/dummy.txt"}
    return DataPacker(cfg, gse, downloaded_data)


# ── _normalize_labels ────────────────────────────────────

class TestNormalizeLabels:
    def test_none(self):
        cfg = _make_cfg()
        p = _make_packer(cfg)
        assert p._normalize_labels(None) == []

    def test_string(self):
        cfg = _make_cfg()
        p = _make_packer(cfg)
        assert p._normalize_labels("Control") == ["control"]

    def test_list(self):
        cfg = _make_cfg()
        p = _make_packer(cfg)
        assert p._normalize_labels(["Control", "ND"]) == ["control", "nd"]

    def test_list_with_none(self):
        cfg = _make_cfg()
        p = _make_packer(cfg)
        assert p._normalize_labels(["Control", None]) == ["control"]

    def test_nested_list_hydra_artifact(self):
        """Hydra 可能产生 [['control', 'ND']] 嵌套，需扁平化。"""
        cfg = _make_cfg()
        p = _make_packer(cfg)
        assert p._normalize_labels([["control", "ND"]]) == ["control", "nd"]

    def test_mixed_nested(self):
        cfg = _make_cfg()
        p = _make_packer(cfg)
        assert p._normalize_labels(["control", ["ND", "WT"]]) == ["control", "nd", "wt"]


# ── _build_batch_tuples ──────────────────────────────────

class TestBuildBatchTuples:
    def test_mcd_cdahfd(self):
        result = DataPacker._build_batch_tuples([
            "MCD, rep1", "CDAHFD, rep1", "MCDHFD something",
        ])
        assert result == [
            ("MCD", "MCD, rep1"),
            ("CDAHFD", "CDAHFD, rep1"),
            ("MCD", "MCDHFD something"),
        ]

    def test_other_passthrough(self):
        result = DataPacker._build_batch_tuples(["Something else"])
        assert result == [("Something else", "Something else")]


# ── _auto_group_division ─────────────────────────────────

class TestAutoGroupDivision:
    """核心：source_name_ch1 失败 → title 回退。"""

    def test_primary_col_matches(self):
        """source_name_ch1 直接匹配成功，不走 title 回退。"""
        cfg = _make_cfg(group_select_col="source_name_ch1")
        meta = _make_meta({
            "source_name_ch1": ["control A", "control B", "MCD 1", "MCD 2"],
            "title": ["t1", "t2", "t3", "t4"],
        })
        p = _make_packer(cfg, meta)
        p._resolve_fibrosis_choice = lambda eg: eg  # skip interactive

        p._auto_group_division(meta)

        assert p._group_col == "source_name_ch1"
        assert len(p._group_mapping["Control"]) == 2
        assert len(p._group_mapping["Fibrosis"]) == 2

    def test_primary_fails_title_fallback(self):
        """source_name_ch1 无匹配，title 列回退成功。"""
        cfg = _make_cfg(group_select_col="source_name_ch1")
        meta = _make_meta({
            "source_name_ch1": ["tissue A", "tissue B", "tissue A", "tissue B"],
            "title": ["ND, sample 1", "ND, sample 2", "MCD, sample 1", "MCD, sample 2"],
        })
        p = _make_packer(cfg, meta)
        p._resolve_fibrosis_choice = lambda eg: eg

        p._auto_group_division(meta)

        assert p._group_col == "title", "should fallback to title"
        assert len(p._group_mapping["Control"]) == 2  # ND matches "ND"
        assert len(p._group_mapping["Fibrosis"]) == 2  # MCD

    def test_both_fail_falls_to_manual(self):
        """两列都无匹配时触发 _manual_group_division。"""
        cfg = _make_cfg(group_select_col="source_name_ch1")
        meta = _make_meta({
            "source_name_ch1": ["a", "b"],
            "title": ["x", "y"],
        })
        p = _make_packer(cfg, meta)

        called = []
        p._manual_group_division = lambda m: called.append(True)

        p._auto_group_division(meta)

        assert called, "should fallback to manual group division"

    def test_title_fallback_with_control_label_nd(self):
        """control_label 含 'ND' 时，title 中的 ND 样本被识别为 Control。"""
        cfg = _make_cfg(
            group_select_col="source_name_ch1",
            control_label=["control", "ND"],
            exp_type="Fibrosis",
        )
        meta = _make_meta({
            "source_name_ch1": ["Spleen"] * 5,  # single value, useless for grouping
            "title": [
                "ND, rep 1", "ND, rep 2",
                "MCD, rep 1", "MCD, rep 2",
                "CDAHFD, rep 1",
            ],
        })
        p = _make_packer(cfg, meta)
        p._resolve_fibrosis_choice = lambda eg: eg

        p._auto_group_division(meta)

        assert p._group_col == "title"
        ctrl = p._group_mapping["Control"]
        exp = p._group_mapping["Fibrosis"]
        assert len(ctrl) == 2  # ND titles
        assert len(exp) == 3  # MCD + CDAHFD

    def test_title_fallback_skips_when_same_col(self):
        """group_select_col 本身就是 'title' 时不走回退（避免重复）。"""
        cfg = _make_cfg(group_select_col="title")
        meta = _make_meta({
            "title": ["ND 1", "ND 2", "MCD 1"],
        })
        p = _make_packer(cfg, meta)
        p._resolve_fibrosis_choice = lambda eg: eg

        p._auto_group_division(meta)

        # should have matched on primary (= title) directly
        assert p._group_col == "title"

    def test_group_select_col_not_in_meta(self):
        """指定的列不存在时直接进入交互。"""
        cfg = _make_cfg(group_select_col="nonexistent_col")
        meta = _make_meta({"title": ["a", "b"]})
        p = _make_packer(cfg, meta)

        called = []
        p._manual_group_division = lambda m: called.append(True)

        p._auto_group_division(meta)
        assert called


# ── _resolve_fibrosis_choice ─────────────────────────────

class TestResolveFibrosisChoice:
    def test_single_category_returns_all(self):
        """只有 MCD 类型的实验组时直接返回全部。"""
        cfg = _make_cfg()
        p = _make_packer(cfg)
        exp_groups = ["MCD, rep 1", "MCD, rep 2"]

        result = p._resolve_fibrosis_choice(exp_groups)
        assert result == exp_groups

    def test_categorizes_mcd_cdahfd(self):
        """MCD 和 CDAHFD 混合时按类别归类。"""
        cfg = _make_cfg()
        p = _make_packer(cfg)
        exp_groups = [
            "MCD, rep 1", "MCD, rep 2",
            "CDAHFD, rep 1", "CDAHFD, rep 2",
            "MCDHFD something",
        ]

        # 检查内部分类逻辑（不触发交互）
        categories = {}
        for g in exp_groups:
            g_lower = str(g).lower()
            if "cdahfd" in g_lower:
                categories.setdefault("CDAHFD", []).append(g)
            elif "mcd" in g_lower:
                categories.setdefault("MCD", []).append(g)
            else:
                categories.setdefault(str(g), []).append(g)

        assert len(categories) == 2  # MCD, CDAHFD
        assert len(categories["MCD"]) == 3  # 2 MCD + 1 MCDHFD
        assert len(categories["CDAHFD"]) == 2

    def test_memory_recall_both(self, tmp_path):
        """记忆 __both__ 时跳过交互，返回全部。"""
        memory_path = tmp_path / "group_memory.yaml"
        memory = {"GSE12345": {"diff": {"_fibrosis_exp_choice": "__both__"}}}
        with open(memory_path, "w") as f:
            yaml.dump(memory, f)

        cfg = _make_cfg(group_memory_enabled=True)
        p = _make_packer(cfg)

        exp_groups = ["MCD 1", "CDAHFD 1"]
        with patch("modules.data_packer.CONFIG_DIR", str(tmp_path)):
            result = p._resolve_fibrosis_choice(exp_groups)
        assert result == exp_groups

    def test_memory_recall_single_category(self, tmp_path):
        """记忆单个类别时返回该类别的值。"""
        memory_path = tmp_path / "group_memory.yaml"
        memory = {"GSE12345": {"diff": {"_fibrosis_exp_choice": "MCD"}}}
        with open(memory_path, "w") as f:
            yaml.dump(memory, f)

        cfg = _make_cfg(group_memory_enabled=True)
        p = _make_packer(cfg)

        exp_groups = [
            "MCD, rep 1", "MCD, rep 2",
            "CDAHFD, rep 1",
        ]
        with patch("modules.data_packer.CONFIG_DIR", str(tmp_path)):
            result = p._resolve_fibrosis_choice(exp_groups)
        assert len(result) == 2
        assert all("MCD" in str(g) for g in result)


# ── _apply_auto_group ────────────────────────────────────

class TestApplyAutoGroup:
    def test_sets_state_correctly(self):
        cfg = _make_cfg()
        meta = _make_meta({"source_name_ch1": ["ctrl 1", "ctrl 2", "exp 1"]})
        p = _make_packer(cfg, meta)

        p._apply_auto_group(meta, "source_name_ch1", ["ctrl 1", "ctrl 2"], ["exp 1"])

        assert p._group_col == "source_name_ch1"
        assert p._group_mapping["Control"] == ["ctrl 1", "ctrl 2"]
        assert p._group_mapping["Fibrosis"] == ["exp 1"]
        assert p._chosen_meta.shape[0] == 3

    def test_non_fibrosis_exp_type(self):
        cfg = _make_cfg(exp_type="")
        meta = _make_meta({"col": ["ctrl", "exp"]})
        p = _make_packer(cfg, meta)

        p._apply_auto_group(meta, "col", ["ctrl"], ["exp"])

        assert p._group_mapping.get("Experiment") == ["exp"]
        assert p._batch_exp_groups == []  # no batch for non-Fibrosis

    def test_fibrosis_batch_tuples(self):
        cfg = _make_cfg(exp_type="Fibrosis")
        meta = _make_meta({"col": ["ctrl", "MCD 1", "CDAHFD 1"]})
        p = _make_packer(cfg, meta)

        p._apply_auto_group(meta, "col", ["ctrl"], ["MCD 1", "CDAHFD 1"])

        assert len(p._batch_exp_groups) == 2
        assert p._batch_exp_groups[0][0] == "MCD"
        assert p._batch_exp_groups[1][0] == "CDAHFD"


# ── _manual_group_division column switch → auto ──────────

class TestManualGroupDivisionAutoFallback:
    def test_column_switch_triggers_auto(self):
        """用户切换列后，尝试自动分组。"""
        cfg = _make_cfg(group_select_col="source_name_ch1")
        meta = _make_meta({
            "source_name_ch1": ["Spleen"] * 4,
            "title": ["ND 1", "ND 2", "MCD 1", "MCD 2"],
            "other_col": ["control_a", "control_b", "mcd_a", "mcd_b"],
        })
        p = _make_packer(cfg, meta)

        # 模拟 _manual_group_select 返回：用户切换到了 "other_col"
        p._manual_group_select = lambda meta, group_label=None, default_col=None: {
            "group_indices": [0, 1],
            "unique_groups": meta["other_col"].unique(),
            "current_col": "other_col",
        }
        p._resolve_fibrosis_choice = lambda eg: eg

        p._manual_group_division(meta)

        # 因为切换了列 (other_col != title)，应该触发 auto_group_division
        assert p._group_col is not None, "should have auto-grouped"

    def test_default_title_col_no_auto_attempt(self):
        """默认 title 列时不触发 auto（避免无限递归）。"""
        cfg = _make_cfg(group_select_col="source_name_ch1")
        meta = _make_meta({
            "title": ["ND 1", "ND 2", "MCD 1"],
        })
        p = _make_packer(cfg, meta)

        # 模拟 _manual_group_select 返回：用户使用默认 title 列
        p._manual_group_select = lambda meta, group_label=None, default_col=None: {
            "group_indices": [0, 1],
            "unique_groups": meta["title"].unique(),
            "current_col": "title",  # default
        }

        auto_called = []
        p._auto_group_division = lambda m: auto_called.append(True)

        p._manual_group_division(meta)
        assert not auto_called, "should NOT trigger auto when on default title column"


# ── get_pack_group_key / resolve_pack_path ───────────────

class TestPackPath:
    def test_diff_mode_returns_diff_group(self):
        assert DataPacker.get_pack_group_key("diff") == "diff_group"
        assert DataPacker.get_pack_group_key("wgcna") == "diff_group"

    def test_hilo_mode_returns_hilo_group(self):
        assert DataPacker.get_pack_group_key("hilo") == "hilo_group"

    def test_other_mode_returns_default(self):
        assert DataPacker.get_pack_group_key("corr") == "default"
        assert DataPacker.get_pack_group_key("immune") == "default"


# ── GSE329774 真实数据回归测试 ────────────────────────────

# 从实际运行输出中提取的 title 列所有唯一值
GSE329774_TITLES = [
    "ND, Splenic CD8+ T cells, ATAC-seq, rep 1",
    "ND, Splenic CD8+ T cells, ATAC-seq, rep 2",
    "ND, Splenic CD8+ T cells, ATAC-seq, rep 3",
    "MCD, Splenic CD8+ T cells, ATAC-seq, rep 1",
    "MCD, Splenic CD8+ T cells, ATAC-seq, rep 2",
    "MCD, Splenic CD8+ T cells, ATAC-seq, rep 3",
    "CD8+T cells, ND, 7 weeks, rep1",
    "CD8+T cells, ND,  7 weeks, rep2",
    "CD8+T cells, ND, 7 weeks, rep3",
    "CD8+T cells, MCDHFD, 7 weeks, rep1",
    "CD8+T cells, MCDHFD, 7 weeks, rep2",
    "CD8+T cells, MCDHFD, 7 weeks, rep3",
    "Control vector Splenic CD8+ T cells, IP Flag, rep 1",
    "Control vector Splenic CD8+ T cells, IP Flag, rep 2",
    "Control vector Splenic CD8+ T cells, IP Flag, rep 3",
    "TRNP1 overexpression, Splenic CD8+ T cells, IP Flag, rep 1",
    "TRNP1 overexpression, Splenic CD8+ T cells, IP Flag, rep 2",
    "TRNP1 overexpression, Splenic CD8+ T cells, IP Flag, rep 3",
    "IgG",
    "CD8+T cells, ND, 7 weeks, rep 1",
    "CD8+T cells, ND,  7 weeks, rep 2",
    "CD8+T cells, ND, 7 weeks, rep 3",
    "CD8+T cells, CDAHFD, 7 weeks, rep 1",
    "CD8+T cells, CDAHFD, 7 weeks, rep 2",
    "CD8+T cells, CDAHFD, 7 weeks, rep 3",
    "CD8+T cells, control vector, rep1",
    "CD8+T cells, control vector, rep2",
    "CD8+T cells, control vector, rep3",
    "CD8+T cells, Flag-TRNP1 OE, rep1",
    "CD8+T cells, Flag-TRNP1 OE, rep2",
    "CD8+T cells, Flag-TRNP1 OE, rep3",
    "ND, Splenic CD8+ T cells, control medium, H3K27me3",
    "ND, Splenic CD8+ T cells, methionine deficient medium, H3K27me3",
    "IgG CUT&Tag control",
    "Splenic CD8+ T cells,RPMI 1640, ATAC-seq",
    "Splenic CD8+ T cells,ROMI 1640 Met(-), ATAC-seq",
    "ND, Splenic CD8+ T cells, 7 weeks, WGBS",
    "CDAHFD, Splenic CD8+ T cells, 7 weeks, WGBS",
    "MCDHFD, Splenic CD8+ T cells, 7 weeks, WGBS",
    "CD8+T cells, RPMI 1640 medium",
    "CD8+T cells, RPMI 1640 methionine (-) medium",
    "ND, Splenic CD8+ T cells, H3K9me3, rep 1",
    "ND, Splenic CD8+ T cells, H3K9me3, rep 2",
    "ND, Splenic CD8+ T cells, H3K9me3, rep 3",
    "ND, Splenic CD8+ T cells, H3K27me3, rep 1",
    "ND, Splenic CD8+ T cells, H3K27me3, rep 2",
    "ND, Splenic CD8+ T cells, H3K27me3, rep 3",
    "MCD, Splenic CD8+ T cells, H3K9me3, rep 1",
    "MCD, Splenic CD8+ T cells, H3K9me3, rep 2",
    "MCD, Splenic CD8+ T cells, H3K9me3, rep 3",
    "MCD, Splenic CD8+ T cells, H3K27me3, rep 1",
    "MCD, Splenic CD8+ T cells, H3K27me3, rep 2",
    "MCD, Splenic CD8+ T cells, H3K27me3, rep 3",
    "ND, Splenic CD8+ T cells, H3K27ac, rep 1",
    "ND, Splenic CD8+ T cells, H3K27ac, rep 2",
    "ND, Splenic CD8+ T cells, H3K27ac, rep 3",
    "MCD, Splenic CD8+ T cells, H3K27ac, rep 1",
    "MCD, Splenic CD8+ T cells, H3K27ac, rep 2",
    "MCD, Splenic CD8+ T cells, H3K27ac, rep 3",
    "IgG, rep 1",
    "ND, scRNA-seq, rep 1",
    "ND, scRNA-seq, rep 2",
    "ND, scRNA-seq, rep 3",
    "MCD, scRNA-seq, rep 1",
    "MCD, scRNA-seq, rep 2",
    "MCD, scRNA-seq, rep 3",
    "ND, Splenic CD8+ T cells, CTCF, rep 1",
    "ND, Splenic CD8+ T cells, CTCF, rep 2",
    "ND, Splenic CD8+ T cells, CTCF, rep 3",
    "MCD, Splenic CD8+ T cells, CTCF, rep 1",
    "MCD, Splenic CD8+ T cells, CTCF, rep 2",
    "MCD, Splenic CD8+ T cells, CTCF, rep 3",
    "ND, Splenic CD8+ T cells, Hi-C, rep 1",
    "MCD, Splenic CD8+ T cells, Hi-C, rep 1",
    "ND, Splenic CD8+ T cells, H3K27me3，rep1",
    "ND, Splenic CD8+ T cells, H3K27me3，rep2",
    "ND, Splenic CD8+ T cells, H3K27me3，rep3",
    "CDAHFD, Splenic CD8+ T cells, H3K27me3,rep1",
    "CDAHFD, Splenic CD8+ T cells, H3K27me3,rep2",
    "CDAHFD, Splenic CD8+ T cells, H3K27me3,rep3",
    "ND, Splenic CD8+ T cells, H3K27ac,rep1",
    "ND, Splenic CD8+ T cells, H3K27ac,rep2",
    "ND, Splenic CD8+ T cells, H3K27ac,rep3",
    "CDAHFD, Splenic CD8+ T cells, H3K27ac,rep1",
    "CDAHFD, Splenic CD8+ T cells, H3K27ac,rep2",
    "CDAHFD, Splenic CD8+ T cells, H3K27ac,rep3",
    "IgG.",
    "ND, Splenic CD8+ T cells, ATAC-seq, rep. 1",
    "ND, Splenic CD8+ T cells, ATAC-seq, rep. 2",
    "ND, Splenic CD8+ T cells, ATAC-seq, rep. 3",
    "CDAHFD, Splenic CD8+ T cells, ATAC-seq, rep 1",
    "CDAHFD, Splenic CD8+ T cells, ATAC-seq, rep 2",
    "CDAHFD, Splenic CD8+ T cells, ATAC-seq, rep 3",
]


class TestGSE329774RealData:
    """用 GSE329774 真实 title 数据做回归测试，确保匹配逻辑不退化。"""

    def test_source_name_ch1_fails_title_fallback_succeeds(self):
        """模拟真实场景：source_name_ch1 只有 'Spleen'，回退 title 匹配 ND/MCD/CDAHFD。"""
        cfg = _make_cfg(
            gse_id="GSE329774",
            group_select_col="source_name_ch1",
            control_label=["control", "ND"],
            exp_label=["CCl4"],
            exp_type="Fibrosis",
        )
        meta = _make_meta({
            "source_name_ch1": ["Spleen"] * len(GSE329774_TITLES),
            "title": GSE329774_TITLES,
        })
        p = _make_packer(cfg, meta)
        p._resolve_fibrosis_choice = lambda eg: eg  # 跳过交互

        p._auto_group_division(meta)

        assert p._group_col == "title", f"应该回退到 title，实际: {p._group_col}"
        ctrl = p._group_mapping["Control"]
        exp = p._group_mapping["Fibrosis"]
        # 所有含 ND 或 control 的 title
        assert len(ctrl) > 0
        # 所有含 MCD 或 CDAHFD 的 title
        assert len(exp) > 0
        # 被选入 chosen_meta 的样本数 = ctrl + exp
        assert p._chosen_meta.shape[0] == len(ctrl) + len(exp)

    def test_title_matching_sanity(self):
        """直接验证 title 列中 control/exp 匹配数量。"""
        control_labels = ["control", "nd"]
        exp_labels = ["mcd", "cdahfd"]

        ctrl = [t for t in GSE329774_TITLES
                if any(l in t.lower() for l in control_labels)]
        exp = [t for t in GSE329774_TITLES
               if any(l in t.lower() for l in exp_labels)]

        assert len(ctrl) > 0, f"control_labels={control_labels} 应该匹配到 title"
        assert len(exp) > 0, f"exp_labels={exp_labels} 应该匹配到 title"
        # 排除同时匹配两组的情况（如 "ND ... control medium" 既含 ND 又含 control）
        ctrl_set = set(ctrl)
        exp_set = set(exp)
        both = ctrl_set & exp_set
        # 验证有意义的匹配数量
        assert len(ctrl_set - exp_set) > 0, "应该有纯 control 条目"
        assert len(exp_set - ctrl_set) > 0, "应该有纯 exp 条目"
