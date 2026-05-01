import os
import hashlib
from typing import Optional

import pandas as pd


def safe_filepath(filepath: str) -> str:
    """若文件已存在，在文件名后加编号 (1), (2), ... 直到不冲突

    Args:
        filepath: 原始文件路径

    Returns:
        不冲突的文件路径
    """
    if not os.path.exists(filepath):
        return filepath

    dirpath = os.path.dirname(filepath)
    name, ext = os.path.splitext(os.path.basename(filepath))
    counter = 1
    while True:
        new_name = f"{name} ({counter}){ext}"
        new_path = os.path.join(dirpath, new_name)
        if not os.path.exists(new_path):
            return new_path
        counter += 1


def file_md5(filepath: str) -> str:
    """计算文件的 MD5 哈希"""
    h = hashlib.md5()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


def df_content_hash(df: pd.DataFrame) -> str:
    """计算 DataFrame 内容哈希（基于 pandas 内置哈希，与行顺序无关）"""
    return hashlib.md5(
        pd.util.hash_pandas_object(df, index=True).values.tobytes()
    ).hexdigest()


def resolve_save_path(filepath: str, content_hash: Optional[str] = None) -> Optional[str]:
    """解析保存路径，防止覆盖同名文件且内容相同时跳过

    Args:
        filepath: 期望的保存路径
        content_hash: 待保存内容的哈希值。为 None 时不做比对，直接递增编号

    Returns:
        安全路径，或 None 表示内容相同无需保存
    """
    if not os.path.exists(filepath):
        return filepath

    if content_hash is not None:
        existing_hash = file_md5(filepath)
        if existing_hash == content_hash:
            return None

    return safe_filepath(filepath)
