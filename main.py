import os

import hydra
import logging
import time

from modules import DataLoader, Analyzer, FigurePlotter

from utils import loggers, Config, DataHandler, FIGURE_DIR, CONFIG_DIR

# 创建配置文件
if not os.path.exists(os.path.join(CONFIG_DIR, "config.yaml")):
    config_path = os.path.join(CONFIG_DIR, "config.yaml")
    with open(config_path, "w") as f:
        with open(os.path.join(CONFIG_DIR, "config.yaml.template"), "r") as template:
            f.write(template.read())

# 主函数
@hydra.main(version_base=None, config_path=CONFIG_DIR, config_name="config")
def main(cfg: Config):
    # 初始化
    logging.getLogger().setLevel(logging.INFO)
    logger = loggers.get_logger()
    logger.info("---欢迎使用本项目---")
    time.sleep(3)
    logger.info("初始化中...")
    data = DataHandler()
    data_dir = os.path.join(cfg.data_dir, cfg.gse_id)
    logger.info("初始化完成")
    time.sleep(3)

    # 数据获取与清洗
    data_pack_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_processed_pack.pkl")

    # 判断是否存在缓存
    if os.path.exists(data_pack_path) and not cfg.debug:
        logger.info(f"发现数据包：{data_pack_path}，跳过下载与清洗")
    else:
        logger.info("开始获取数据并清洗")
        loader = DataLoader(cfg)
        data.meta_matrix_pack = loader.loader()
        logger.info("数据获取完成")
    time.sleep(3)

    # 分析
    if cfg.analysis_mode == "corr" and os.path.exists(os.path.join(data_dir, "pkl", f"{cfg.gse_id}_correlation_summary.pkl")) and not cfg.debug:
        logger.info("发现基因相关性分析结果，跳过计算")
    elif cfg.analysis_mode == "diff" and os.path.exists(os.path.join(data_dir, "pkl", f"{cfg.gse_id}_differential_summary.pkl")) and not cfg.debug:
        logger.info("发现基因差异分析结果，跳过计算")
    else:
        logger.info("开始分析")
        calculater = Analyzer.create(cfg, data)
        result = calculater.calculate()
        if cfg.analysis_mode == "corr":
            data.gene_corr_table = result
        elif cfg.analysis_mode == "diff":
            data.gene_diff_table = result
        else:
            logger.error(f"未知分析模式：{cfg.analysis_mode}，结果未保存")
            raise ValueError(f"未知分析模式：{cfg.analysis_mode}，结果未保存")
        logger.info("分析完成")
    time.sleep(3)

    # 绘图
    logger.info("开始绘图")
    plotter = FigurePlotter.create(cfg, data)
    plotter.plotter()
    logger.info(f"绘图结果保存在{FIGURE_DIR}")


if __name__ == '__main__':
    main()
