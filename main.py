import os
import sys

import hydra
import logging
import time
import yaml

from modules import DataLoader, Analyzer, FigurePlotter
from modules.data_packer import DataPacker

from utils import loggers, Config, DataHandler, FIGURE_DIR, CONFIG_DIR, parse_user_input


def _read_cfg_version(config_path: str) -> str | None:
    """读取配置文件中的版本数据"""
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if line.startswith("version:"):
                    _, value = line.split(":", 1)
                    value = value.strip()
                    if not value:
                        return None
                    if (value.startswith('"') and value.endswith('"')) or (
                        value.startswith("'") and value.endswith("'")
                    ):
                        return value[1:-1]
                    return value
    except FileNotFoundError:
        return None
    return None


def _print_cfg_version_and_exit() -> None:
    """输出版本信息"""
    config_path = os.path.join(CONFIG_DIR, "config.yaml")
    if not os.path.exists(config_path):
        config_path = os.path.join(CONFIG_DIR, "config.yaml.template")

    version = _read_cfg_version(config_path)
    if version is None:
        print("No version field defined in config.")
        sys.exit(0)

    print(version)
    sys.exit(0)

def _print_custom_help_and_exit() -> None:
    """读取 conf/help.yml 并格式化打印帮助信息"""
    help_path = os.path.join(CONFIG_DIR, "help.yml")
    if not os.path.exists(help_path):
        print("帮助文件未找到: conf/help.yml")
        sys.exit(1)

    with open(help_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    width = 72
    print("=" * width)
    print("SGA (Simple GEO Analyzer) — 配置项参考")
    print("=" * width)
    print()
    print(data.get("usage", "").strip())
    print()

    for section in data.get("sections", []):
        print(f"[{section['title']}]")
        print("-" * 48)
        for field in section.get("fields", []):
            key = field.get("key", "")
            ftype = field.get("type", "")
            default = field.get("default", "")
            desc = field.get("desc", "").strip()
            choices = field.get("choices")
            note = field.get("note")

            print(f"  {key}")
            print(f"    类型: {ftype}")
            print(f"    默认: {default}")
            if choices:
                print(f"    可选: {', '.join(choices)}")
            for line in desc.split("\n"):
                print(f"    {line.strip()}")
            if note:
                print(f"    注意: {note}")
            print()
    sys.exit(0)


# 如果命令为查看版本信息则直接输出
if len(sys.argv) > 1 and sys.argv[1] in {"version", "--version", "-v"}:
    _print_cfg_version_and_exit()

# 如果命令为查看帮助则直接输出
if len(sys.argv) > 1 and sys.argv[1] in {"help", "--help", "-h"}:
    _print_custom_help_and_exit()

# 创建/补全配置文件（缺失字段从 template 追加）
_config_path = os.path.join(CONFIG_DIR, "config.yaml")
_template_path = os.path.join(CONFIG_DIR, "config.yaml.template")

if not os.path.exists(_config_path):
    with open(_config_path, "w") as f:
        with open(_template_path, "r") as t:
            f.write(t.read())
else:
    with open(_template_path, "r") as t:
        template_data = yaml.safe_load(t)
    with open(_config_path, "r") as c:
        config_data = yaml.safe_load(c) or {}
    # 将缺失的顶层键追写到文件末尾
    missing = [k for k in template_data if k not in config_data and k != "hydra"]
    if missing:
        with open(_config_path, "a") as f:
            f.write("\n")
            for key in missing:
                val = template_data[key]
                if isinstance(val, list):
                    f.write(f"{key}:\n")
                    for item in val:
                        f.write(f"  - {item}\n")
                elif isinstance(val, bool):
                    f.write(f"{key}: {str(val).lower()}\n")
                elif isinstance(val, str):
                    f.write(f'{key}: "{val}"\n')
                else:
                    f.write(f"{key}: {val}\n")


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
    data, logger = _init(cfg)

    while True:
        data_dir = os.path.join(cfg.data_dir, cfg.gse_id)

        # 阶段1: 数据获取与清洗
        if "1" in str(cfg.process):
            if cfg.analysis_mode != "enrich" and data.meta_matrix_pack is None:
                data_pack_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_processed_pack.pkl")

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
            "  6. 切换新数据集 \n(tips.若要修改GSE_ID和tar_gene以外的其他配置请重启本项目)\n"
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
            max_index=6,
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
        mode_map = {1: "corr", 2: "diff", 3: "hilo", 4: "enrich", 5: "immune"}

        if option in mode_map:
            cfg.analysis_mode = mode_map[option]
            continue

        if option == 5:
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
            data = DataHandler()
            logger.info(f"已切换至 GSE: {new_gse}, 目标基因: {new_gene}")
            continue


if __name__ == '__main__':
    main()
