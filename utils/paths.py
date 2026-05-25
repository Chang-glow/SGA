import os

# 路径配置
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULT_DIR = os.path.join(BASE_DIR, "res")
FIGURE_DIR = os.path.join(RESULT_DIR, "figures")
CONFIG_DIR = os.path.join(BASE_DIR, "conf")
LOGGER_DIR = os.path.join(BASE_DIR, "error_logs")


def relpath(path: str) -> str:
    """将绝对路径转为相对于项目根目录的路径，用于日志输出"""
    try:
        return os.path.relpath(path, BASE_DIR)
    except ValueError:
        return path


def dirs_init():
    for d in [BASE_DIR, DATA_DIR, RESULT_DIR, FIGURE_DIR, CONFIG_DIR, LOGGER_DIR]:
        os.makedirs(d, exist_ok=True)


dirs_init()
