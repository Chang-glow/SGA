import os

import hydra
import logging
import time
import pandas as pd

from modules import DataLoader, Analyzer, FigurePlotter, DataPacker
from modules.calculater import normalize_tar_gene_from_pack, detect_gene_case_convention, normalize_gene_symbol
from modules.strategies import NO_GEO_MODES

from utils import loggers, Config, DataHandler, FIGURE_DIR, CONFIG_DIR, parse_user_input, BASE_DIR, relpath


_RESULT_ATTR = {
    "corr": "gene_corr_table",
    "diff": "gene_diff_table",
    "hilo": "gene_hilo_table",
    "immune": "gene_immune_table",
    "wgcna": "gene_wgcna_table",
    "enrich": "gene_enrich_table",
}

_SKIP_CHECK = {
    "corr": ("correlation", "基因相关性"),
    "diff": ("differential", "基因差异"),
    "hilo": ("highlow", "基因高低表达"),
    "immune": ("immune", "免疫浸润"),
    "enrich": ("enrichment", "基因富集"),
    "wgcna": ("wgcna", "WGCNA"),
}


def _resolve_config(cfg: Config) -> None:
    """解析配置冲突与自动推导。"""
    cfg.gse_id = str(cfg.gse_id)
    if cfg.gse_id.isdigit():
        cfg.gse_id = f"GSE{cfg.gse_id}"
    # multi_gene 优先：自动清 tar_gene
    if cfg.multi_gene and cfg.tar_gene:
        logger = logging.getLogger(__name__)
        logger.warning(
            f"tar_gene ({cfg.tar_gene!r}) 和 multi_gene 同时存在，"
            f"multi_gene 优先，已自动清空 tar_gene。"
        )
        cfg.tar_gene = ""

    if cfg.multi_gene and cfg.analysis_mode in NO_GEO_MODES:
        if os.path.isfile(cfg.multi_gene):
            stem = os.path.splitext(os.path.basename(cfg.multi_gene))[0]
        else:
            stem = "multi_gene_analysis"
        cfg.gse_id = stem


def _sync_organism_from_pack_path(cfg: Config, pack_path: str) -> None:
    """从缓存 pack 文件读取 _organism 并同步到 cfg.organism。"""
    try:
        import pandas as pd
        import os as _os
        if not _os.path.exists(pack_path):
            return
        pack = pd.read_pickle(pack_path)
        if isinstance(pack, dict) and '_organism' in pack:
            cfg.organism = pack['_organism']
    except Exception:
        pass


def _init(cfg: Config) -> tuple:
    """初始化日志、数据处理器"""
    logging.getLogger().setLevel(logging.INFO)
    logger = loggers.get_logger()
    logger.info("初始化中...")
    data = DataHandler()
    logger.info("初始化完成")
    return data, logger


def _normalize_tar_gene_from_raw(cfg: Config, downloaded_data: dict) -> None:
    """从原始下载文件中提取表达矩阵，检测基因命名规则并规范化 cfg.tar_gene。

    在 pack 构建前调用，保障 hilo 模式的 _prepare_hilo_group 使用规范化后的基因名。
    """
    import pandas as pd

    if not cfg.tar_gene or cfg.analysis_mode == "enrich" or not downloaded_data:
        return
    for file_path in downloaded_data.values():
        if file_path and os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path, sep="\t", compression="infer", index_col=0)
                convention = detect_gene_case_convention(df.index)
                new_gene = normalize_gene_symbol(cfg.tar_gene, convention)
                if new_gene != cfg.tar_gene:
                    logger = logging.getLogger(__name__)
                    logger.info(
                        f"基因符号已规范化: {cfg.tar_gene!r} → {new_gene!r}"
                        f"（检测到 {convention} 命名规则）"
                    )
                    cfg.tar_gene = new_gene
                return
            except Exception:
                continue


# 主函数
@hydra.main(version_base=None, config_path=CONFIG_DIR, config_name="config")
def main(cfg: Config):
    # data_dir 由 Hydra 注入 YAML 值，需手动转为绝对路径
    if not os.path.isabs(cfg.data_dir):
        cfg.data_dir = os.path.abspath(os.path.join(BASE_DIR, cfg.data_dir))

    _resolve_config(cfg)

    data, logger = _init(cfg)

    while True:
        data_dir = os.path.join(cfg.data_dir, cfg.gse_id)

        # 阶段1: 数据获取与清洗
        if "1" in str(cfg.process):
            _t0 = time.time()
            if not (cfg.multi_gene and cfg.analysis_mode in NO_GEO_MODES) and data.meta_matrix_pack is None:
                data_pack_path = DataPacker.resolve_pack_path(data_dir, cfg.gse_id, cfg.analysis_mode)

                if os.path.exists(data_pack_path) and not cfg.force:
                    logger.info(f"发现数据包：{data_pack_path}，跳过下载与清洗")
                    data.meta_matrix_pack = pd.read_pickle(data_pack_path)
                    if isinstance(data.meta_matrix_pack, dict) and '_organism' in data.meta_matrix_pack:
                        cfg.organism = data.meta_matrix_pack['_organism']
                else:
                    logger.info("开始获取数据并清洗")
                    loader = DataLoader(cfg)
                    downloaded_data = loader.loader()
                    # 基因大小写规范化（pack 构建前，保障 hilo 模式）
                    _normalize_tar_gene_from_raw(cfg, downloaded_data)
                    packer = DataPacker(cfg, loader.gse, downloaded_data)
                    data.meta_matrix_pack = packer.build_pack()
                    # 从 pack 自动检测的 organism 同步到 cfg（影响 enrich/WGCNA/基因规范化）
                    if data.meta_matrix_pack.get('_organism'):
                        cfg.organism = data.meta_matrix_pack['_organism']
                    # 基因大小写规范化（pack 构建后，表达矩阵已正确索引）
                    normalize_tar_gene_from_pack(cfg, data.meta_matrix_pack)
                    logger.info("数据获取完成")
            elif cfg.multi_gene and cfg.analysis_mode in NO_GEO_MODES:
                logger.info("multi_gene 直接输入，无需 GEO 数据，跳过下载与清洗。")
            logger.info(f"阶段1 完成，耗时 {time.time() - _t0:.1f}s")
        else:
            logger.info("跳过阶段1（数据下载与清洗）")

        # 阶段2: 分析
        if "2" in str(cfg.process):
            _t0 = time.time()
            _no_batch = not getattr(cfg, "_batch_suffix", "")
            _skip = False
            if not cfg.force and _no_batch:
                _pattern = _SKIP_CHECK.get(cfg.analysis_mode)
                if _pattern:
                    _pkl = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_{_pattern[0]}_summary.pkl")
                    if os.path.exists(_pkl):
                        logger.info(f"发现{_pattern[1]}分析结果，跳过计算")
                        _skip = True
            if not _skip:
                logger.info("开始分析")
                normalize_tar_gene_from_pack(cfg, data.meta_matrix_pack)
                calculater = Analyzer.create(cfg, data)
                result = calculater.calculate()
                attr = _RESULT_ATTR.get(cfg.analysis_mode)
                if attr is None:
                    logger.error(f"未知分析模式：{cfg.analysis_mode}，结果未保存")
                    raise ValueError(f"未知分析模式：{cfg.analysis_mode}，结果未保存")
                if not (cfg.analysis_mode == "enrich" and result is None):
                    setattr(data, attr, result)
                logger.info("分析完成")
            logger.info(f"阶段2 完成，耗时 {time.time() - _t0:.1f}s")
        else:
            logger.info("跳过阶段2（分析计算）")

        # 阶段3: 画图
        if "3" in str(cfg.process):
            _t0 = time.time()
            logger.info("开始绘图")
            normalize_tar_gene_from_pack(cfg, data.meta_matrix_pack)
            plotter = FigurePlotter.create(cfg, data)
            plotter.plotter()
            logger.info(f"绘图结果保存在 {relpath(FIGURE_DIR)}")
            logger.info(f"阶段3 完成，耗时 {time.time() - _t0:.1f}s")
        else:
            logger.info("跳过阶段3（画图）")

        # 批处理：Fibrosis "两者都分析" — 切到下一个实验组，continue 复用阶段 2/3
        if data.meta_matrix_pack and data.meta_matrix_pack.get("_batch_exp_groups"):
            next_label, next_exp = data.meta_matrix_pack["_batch_exp_groups"][0]
            object.__setattr__(cfg, "_batch_suffix", f"_{next_label}")
            logger.info(f"--- 批处理：对实验组 {next_label} ({next_exp}) 进行分析 ---")

            pack = data.meta_matrix_pack
            meta_full = pack["meta_full"]
            group_col = pack["group_info"]["group_col"]
            control_values = pack["group_info"]["mapping"]["Control"]
            exp_type = cfg.exp_type or "Experiment"
            chosen_meta = meta_full[
                meta_full[group_col].isin(control_values + [next_exp])
            ].copy()
            chosen_meta["group"] = chosen_meta[group_col].apply(
                lambda v: "Control" if v in control_values else exp_type
            )
            chosen_meta = chosen_meta[chosen_meta["group"].notna()]
            pack["meta"] = chosen_meta
            pack["group_info"]["mapping"][exp_type] = [next_exp]
            pack["_batch_exp_groups"] = pack["_batch_exp_groups"][1:]
            continue

        # 清除批处理标记，避免交互菜单切换模式后残留
        object.__setattr__(cfg, "_batch_suffix", "")

        # 交互菜单
        menu = (
            "\n" + "=" * 48 + "\n"
            "分析完成！请选择下一步操作：\n"
            "  1. 相关性分析 (corr)\n"
            "  2. 差异分析 (diff)\n"
            "  3. 高低表达分析 (hilo)\n"
            "  4. 富集分析 (enrich)\n"
            "  5. 免疫浸润分析 (immune)\n"
            "  6. WGCNA数据准备 (wgcna)\n"
            "  7. 切换新数据集 \n(tips.若要修改GSE_ID和tar_gene以外的其他配置请重启本项目)\n"
            "---\n"
            "  0. 退出\n"
            + "=" * 48
        )
        if not cfg.debug:
            print(menu)
        else:
            break

        choice = parse_user_input(
            prompt="请输入选项序号: ",
            max_index=7,
            whitelist="0",
        )

        if choice == "0":
            logger.info("程序退出。")
            break

        if not choice or not isinstance(choice, list):
            if not cfg.debug:
                print("无效输入，请重新选择...")
                continue
            else:
                break

        option = choice[0]
        mode_map = {1: "corr", 2: "diff", 3: "hilo", 4: "enrich", 5: "immune", 6: "wgcna"}

        if option in mode_map:
            cfg.analysis_mode = mode_map[option]
            continue

        if option == 7:
            new_gse = input("请输入新的 GSE ID(如GSE123456): ").strip()
            if not new_gse:
                print("GSE ID 不能为空，返回菜单。")
                continue
            new_gene = input("请输入新的目标基因(如TP53): ").strip()
            if not new_gene:
                print("目标基因不能为空，返回菜单。")
                continue
            cfg.gse_id = new_gse
            cfg.tar_gene = new_gene
            cfg.multi_gene = ""
            data = DataHandler()
            logger.info(f"已切换至 GSE: {new_gse}, 目标基因: {new_gene}")
            continue


if __name__ == '__main__':
    main()
