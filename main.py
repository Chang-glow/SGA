import os

import hydra
import logging
import time

from modules import DataLoader, Analyzer, FigurePlotter
from modules.data_packer import DataPacker

from utils import loggers, Config, DataHandler, FIGURE_DIR, CONFIG_DIR, parse_user_input
from utils.paths import BASE_DIR


def _resolve_config(cfg: Config) -> None:
    """Hydra 会覆盖 __post_init__ 的修改，因此将自动推导逻辑放在此处。"""
    # multi_gene 优先：自动清 tar_gene
    if cfg.multi_gene and cfg.tar_gene:
        logger = logging.getLogger(__name__)
        logger.warning(
            f"tar_gene ({cfg.tar_gene!r}) 和 multi_gene 同时存在，"
            f"multi_gene 优先，已自动清空 tar_gene。"
        )
        cfg.tar_gene = ""

    # multi_gene 自动推导 gse_id
    if cfg.multi_gene:
        if os.path.isfile(cfg.multi_gene):
            stem = os.path.splitext(os.path.basename(cfg.multi_gene))[0]
        else:
            stem = "multi_gene_analysis"
        cfg.gse_id = stem


def _init(cfg: Config) -> tuple:
    """初始化日志、数据处理器"""
    logging.getLogger().setLevel(logging.INFO)
    logger = loggers.get_logger()
    logger.info("---欢迎使用本项目---")
    if not cfg.debug:
        time.sleep(3)
    logger.info("初始化中...")
    data = DataHandler()
    logger.info("初始化完成")
    if not cfg.debug:
        time.sleep(3)
    return data, logger


# 主函数
@hydra.main(version_base=None, config_path=CONFIG_DIR, config_name="config")
def main(cfg: Config):
    # 相对路径基于项目根目录解析（Hydra 用 YAML 值覆盖了 __post_init__ 的解析结果）
    if not os.path.isabs(cfg.data_dir):
        cfg.data_dir = os.path.abspath(os.path.join(BASE_DIR, cfg.data_dir))

    _resolve_config(cfg)

    data, logger = _init(cfg)

    while True:
        data_dir = os.path.join(cfg.data_dir, cfg.gse_id)

        # 阶段1: 数据获取与清洗
        if "1" in str(cfg.process):
            if cfg.analysis_mode != "enrich" and data.meta_matrix_pack is None:
                data_pack_path = DataPacker.resolve_pack_path(data_dir, cfg.gse_id, cfg.analysis_mode)

                if os.path.exists(data_pack_path) and not cfg.debug:
                    logger.info(f"发现数据包：{data_pack_path}，跳过下载与清洗")
                else:
                    logger.info("开始获取数据并清洗")
                    loader = DataLoader(cfg)
                    downloaded_data = loader.loader()
                    packer = DataPacker(cfg, loader.gse, downloaded_data)
                    data.meta_matrix_pack = packer.build_pack()
                    logger.info("数据获取完成")
            elif cfg.analysis_mode == "enrich":
                logger.info("富集分析模式：跳过数据下载与清洗阶段。")
            if not cfg.debug:
                time.sleep(1)
        else:
            logger.info("跳过阶段1（数据下载与清洗）")

        # 阶段2: 分析
        if "2" in str(cfg.process):
            if cfg.analysis_mode == "corr" and os.path.exists(os.path.join(data_dir, "pkl", f"{cfg.gse_id}_correlation_summary.pkl")) and not cfg.debug:
                logger.info("发现基因相关性分析结果，跳过计算")
            elif cfg.analysis_mode == "diff" and os.path.exists(os.path.join(data_dir, "pkl", f"{cfg.gse_id}_differential_summary.pkl")) and not cfg.debug:
                logger.info("发现基因差异分析结果，跳过计算")
            elif cfg.analysis_mode == "immune" and os.path.exists(os.path.join(data_dir, "pkl", f"{cfg.gse_id}_immune_summary.pkl")) and not cfg.debug:
                logger.info("发现免疫浸润分析结果，跳过计算")
            elif cfg.analysis_mode == "enrich" and os.path.exists(os.path.join(data_dir, "pkl", f"{cfg.gse_id}_enrichment_summary.pkl")) and not cfg.debug:
                logger.info("发现基因富集分析结果，跳过计算")
            elif cfg.analysis_mode == "wgcna" and os.path.exists(os.path.join(data_dir, "pkl", f"{cfg.gse_id}_wgcna_summary.pkl")) and not cfg.debug:
                logger.info("发现WGCNA分析结果，跳过计算")
            else:
                logger.info("开始分析")
                calculater = Analyzer.create(cfg, data)
                result = calculater.calculate()
                if cfg.analysis_mode == "corr":
                    data.gene_corr_table = result
                elif cfg.analysis_mode == "diff":
                    data.gene_diff_table = result
                elif cfg.analysis_mode == "hilo":
                    data.gene_hilo_table = result
                elif cfg.analysis_mode == "immune":
                    data.gene_immune_table = result
                elif cfg.analysis_mode == "wgcna":
                    data.gene_wgcna_table = result
                elif cfg.analysis_mode == "enrich":
                    if result is not None:
                        data.gene_enrich_table = result
                else:
                    logger.error(f"未知分析模式：{cfg.analysis_mode}，结果未保存")
                    raise ValueError(f"未知分析模式：{cfg.analysis_mode}，结果未保存")
                logger.info("分析完成")
            if not cfg.debug:
                time.sleep(1)
        else:
            logger.info("跳过阶段2（分析计算）")

        # 阶段3: 画图
        if "3" in str(cfg.process):
            logger.info("开始绘图")
            plotter = FigurePlotter.create(cfg, data)
            plotter.plotter()
            logger.info(f"绘图结果保存在{FIGURE_DIR}")
        else:
            logger.info("跳过阶段3（画图）")

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
