import os
import re
import sys
import shutil
import tempfile
import logging
import requests
import GEOparse
from typing import Optional, Dict

import pandas as pd
import yaml

from utils import loggers, Config, parse_user_input, safe_filepath, CONFIG_DIR

logging.getLogger("GEOparse").setLevel(logging.WARNING)


class DataLoader:
    """从GEO下载数据并储存"""
    _logger = loggers.get_logger()

    def __init__(self, cfg: Config):
        """初始化数据获取对象

        Args:
            cfg: 基本配置项
        """
        self.cfg: Config = cfg
        self._data: Optional[Dict[str, pd.DataFrame]] = None
        self._download_data = None
        self._group_mapping = {}
        self._group_col = None
        self.gse = None

    def loader(self) -> dict:
        """
        用于调用数据的API

        Returns:
            downloaded_data: 下载后的文件路径字典
        """
        if self._download_data is None:
            gse = self._get_gse()
            self.gse = gse
            self._download_data = self._user_selection_flow(gse)
            if not self._download_data:
                raise RuntimeError("未获取到补充矩阵文件，请检查 GEO 数据集是否包含匹配的 matrix/count 文件，并确保按提示选择文件")
            return self._download_data
        return self._download_data

    def download_geo_data(self, url: str) -> Optional[str]:
        """从GEO下载所需数据

        Args:
            url: 下载链接

        Returns:
            local_path: 下载文件所在位置
        """
        # 确定本地位置
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        file_name = os.path.basename(url)
        local_path = os.path.join(data_dir, file_name)

        # 检验是否存在文件
        if os.path.exists(local_path):
            DataLoader._logger.info(f"文件{file_name} 已存在，跳过下载。")
            return local_path

        # 转换请求头
        if url.startswith("ftp://"):
            url = url.replace("ftp://", "https://", 1)

        # 下载逻辑
        DataLoader._logger.info(f"正在从NCBI下载{file_name}...")
        try:
            response = requests.get(url, stream=True, timeout=30)
            response.raise_for_status()

            self._save_response_to_file(response, local_path)
            DataLoader._logger.info(f"{local_path}下载完成")
            return local_path
        except Exception as e:
            DataLoader._logger.error(f"下载失败: {e}")
            if os.path.exists(local_path):
                os.remove(local_path)
            return None

    def _save_response_to_file(self, response, local_path):
        total_length = response.headers.get('content-length')
        downloaded = 0
        if total_length is not None:
            total_length = int(total_length)

        with open(local_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if not chunk:
                    continue
                f.write(chunk)
                downloaded += len(chunk)
                self._print_download_progress(downloaded, total_length)

        sys.stdout.write("\n")

    def _print_download_progress(self, downloaded, total_length):
        if total_length is not None:
            percent = int(downloaded / total_length * 100)
            bar_len = 30
            filled = int(bar_len * percent / 100)
            bar = '#' * filled + '-' * (bar_len - filled)
            sys.stdout.write(
                f"\r下载进度: [{bar}] {percent}% ({downloaded // 1024}KB/{total_length // 1024}KB)"
            )
        else:
            sys.stdout.write(f"\r已下载 {downloaded // 1024}KB")
        sys.stdout.flush()

    def _get_gse(self) -> GEOparse.GEOTypes.GSE:
        """获取GEO数据包

        SOFT 文件保留在 data/{GSE_ID}/，GEOparse 调用期间 CWD 隔离到 /tmp，
        防止平台注释等中间文件污染项目根目录。

        Returns:
            gse: 下载的数据包
        """
        gse_id = self.cfg.gse_id
        dest_dir = os.path.normpath(os.path.join(self.cfg.data_dir, self.cfg.gse_id))

        try:
            is_str = isinstance(self.cfg.gse_id, str) and isinstance(dest_dir, str)
            if not is_str:
                raise TypeError("请输入字符串而非其他类型参数")

            is_gse = re.match(r'^GSE\d+$', gse_id)
            if not is_gse:
                raise ValueError("请输入正确的GSE ID！")

            # CWD 隔离到 /tmp，GEOparse 任何非 destdir 写入都不会污染项目目录
            tmpdir = tempfile.mkdtemp(prefix="sga_tmp_")
            orig_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                DataLoader._logger.info(f"开始下载或调用现存{gse_id}_family.soft.gz")
                DataLoader._logger.info("正在检索远程服务器或本地缓存...")
                gse = GEOparse.get_GEO(geo=gse_id, destdir=dest_dir)
            finally:
                os.chdir(orig_cwd)
                shutil.rmtree(tmpdir, ignore_errors=True)

            if not gse:
                raise Exception("出现未知错误，请检查\n1、GSE编号是否正确\n2、下载地址是否正确/有权限写入")

            return gse

        except ValueError as e:
            DataLoader._logger.error(f"【输入错误】:{e}")

        except TypeError as e:
            DataLoader._logger.error(f"【类型错误】:{e}")

        except Exception as e:
            DataLoader._logger.error(f"【未知错误】:{e}")
            raise
    
    def _user_selection_flow(self, gse) -> None:
        """用户交互下载数据，支持记忆矩阵选择

        Args:
            gse: 下载的GES对象
        """
        try:
            # 从gse中读取补充文件列表
            sp_files = gse.metadata.get('supplementary_file', [])
            if not sp_files:
                DataLoader._logger.warning("未发现补充文件")

            # 简单筛选去除明显不是目标文件内容
            candidates = [f for f in sp_files if (
                '.matrix' in f.lower()
                or '.count' in f.lower()
                or '.txt' in f.lower()
            ) and 'readme' not in f.lower()]

            # 无候选补充文件时尝试从 SOFT 提取内嵌表达矩阵
            if not candidates:
                DataLoader._logger.warning("未在补充文件中发现 matrix/count/txt 文件，尝试从 SOFT 文件提取")
                soft_result = self._extract_expression_from_soft(gse)
                if soft_result is not None:
                    return soft_result
                raise FileNotFoundError("未发现补充文件且 SOFT 文件中也没有表达数据")

            selected_idx = None

            # 尝试从记忆中加载矩阵选择
            if getattr(self.cfg, "group_memory_use", False):
                saved = self._load_matrix_memory(candidates)
                if saved is not None:
                    selected_idx = saved
                    DataLoader._logger.info("已从记忆中恢复矩阵选择，跳过交互。")

            # 无记忆时交互选择
            if selected_idx is None:
                print("\n--- 发现以下疑似矩阵文件 ---")
                for i, url in enumerate(candidates):
                    print(f"[{i}] {os.path.basename(url)}")

                selected_idx = parse_user_input(
                    prompt="请输入需要的矩阵序号(如1:8,11):",
                    max_index=len(candidates) - 1,
                )

                # 保存当前选择
                if getattr(self.cfg, "group_memory_enabled", False):
                    self._save_matrix_memory(candidates, selected_idx)

            selected_urls = [candidates[i] for i in selected_idx]

            # 下载选中的文件
            downloaded_data = {}
            for url in selected_urls:
                file_path = self.download_geo_data(url)

                # 只有当路径不为None时才存入字典，防止后面read_csv报错
                if file_path:
                    downloaded_data[os.path.basename(url)] = file_path
                else:
                    DataLoader._logger.warning(f"文件 {os.path.basename(url)} 下载失败，将不会被加载。")

            if not downloaded_data:
                raise Exception("出现未知错误，请检查\n1、GSE编号是否正确\n2、下载地址是否正确/有权限写入")

            return downloaded_data

        except FileNotFoundError as e:
            DataLoader._logger.error(f"【文件未找到】:{e}")
        except Exception as e:
            DataLoader._logger.error(f"【未知错误】:{e}")

    def _extract_expression_from_soft(self, gse) -> Optional[dict]:
        """从 SOFT 文件内嵌的 GSM 数据表提取基因级表达矩阵

        适用于芯片数据集（SOFT 中每个 GSM 内嵌表达值表），
        不适用于 RNA-seq（Sample_data_row_count = 0，无内嵌数据）。

        Returns:
            {filename: filepath} 或 None
        """
        try:
            first_gsm = next(iter(gse.gsms.values()))
            if first_gsm.table is None or first_gsm.table.empty:
                DataLoader._logger.info("SOFT 文件中无内嵌表达数据")
                return None

            gpl_counts = {}
            for gsm in gse.gsms.values():
                gpl_name = gsm.metadata.get('platform_id', ['unknown'])[0]
                gpl_counts[gpl_name] = gpl_counts.get(gpl_name, 0) + 1
            main_gpl_name = max(gpl_counts, key=gpl_counts.get)
            gpl = gse.gpls.get(main_gpl_name)
            if gpl is None or gpl.table is None:
                DataLoader._logger.warning(f"无法获取平台 {main_gpl_name} 的注释表")
                return None

            gene_col = None
            for candidate in ['GeneSymbol', 'Gene Symbol', 'SYMBOL', 'GENE_SYMBOL']:
                if candidate in gpl.table.columns:
                    gene_col = candidate
                    break
            if gene_col is None:
                DataLoader._logger.warning(
                    f"平台 {main_gpl_name} 注释表中未找到基因名列，可用列: {list(gpl.table.columns)}"
                )
                return None

            # 构建探针→基因符号映射
            gpl_df = gpl.table.set_index('ID')
            probe_to_gene = gpl_df[gene_col].dropna()
            DataLoader._logger.info(
                f"探针→基因映射: {len(probe_to_gene)} 个探针有基因注释"
            )

            # 逐列构建表达矩阵，避免 pivot_and_annotate 的内存峰值
            n_samples = len(gse.gsms)
            DataLoader._logger.info(f"从 {n_samples} 个样本提取表达值...")
            data_series = {}
            for i, (gsm_name, gsm) in enumerate(gse.gsms.items()):
                if gsm.table is None or gsm.table.empty:
                    continue
                if i % 100 == 0:
                    DataLoader._logger.info(f"  提取进度: {i + 1}/{n_samples}")
                tbl = gsm.table.set_index('ID_REF')
                if 'VALUE' in tbl.columns:
                    data_series[gsm_name] = tbl['VALUE']

            if not data_series:
                DataLoader._logger.warning("未能从任何 GSM 提取到表达值")
                return None

            expr_matrix = pd.DataFrame(data_series)
            DataLoader._logger.info(f"探针级矩阵: {expr_matrix.shape[0]} 探针 × {expr_matrix.shape[1]} 样本")

            # 映射探针到基因，丢弃无注释的行
            mapped_genes = expr_matrix.index.map(probe_to_gene)
            valid = mapped_genes.notna()
            expr_matrix = expr_matrix.loc[valid]
            mapped_genes = mapped_genes[valid]
            DataLoader._logger.info(f"有基因注释的探针: {len(expr_matrix)} 行")

            if expr_matrix.empty:
                DataLoader._logger.warning("去除无基因名注释的探针后矩阵为空")
                return None

            # 聚合成基因级
            expr_matrix = expr_matrix.astype(float)
            gene_names = pd.Index(mapped_genes, name='_gene')
            gene_matrix = expr_matrix.groupby(gene_names).mean()
            DataLoader._logger.info(f"基因级矩阵: {gene_matrix.shape[0]} 基因 × {gene_matrix.shape[1]} 样本")

            data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
            os.makedirs(data_dir, exist_ok=True)
            file_name = f"{self.cfg.gse_id}_soft_extracted.txt"
            file_path = safe_filepath(os.path.join(data_dir, file_name))
            gene_matrix.to_csv(file_path, sep="\t")
            DataLoader._logger.info(f"已保存至 {file_path}")

            return {file_name: file_path}

        except Exception as e:
            DataLoader._logger.error(f"从 SOFT 提取表达矩阵失败: {e}")
            return None

    def _save_matrix_memory(self, candidates: list, selected_idx: list) -> None:
        """保存矩阵选择到记忆文件（与分组记忆共用同一文件）"""
        memory_path = os.path.join(CONFIG_DIR, "group_memory.yaml")
        memory = {}
        if os.path.exists(memory_path):
            with open(memory_path, "r", encoding="utf-8") as f:
                memory = yaml.safe_load(f) or {}
        gse = self.cfg.gse_id
        memory.setdefault(gse, {})["matrices"] = [
            os.path.basename(candidates[i]) for i in selected_idx
        ]
        with open(memory_path, "w", encoding="utf-8") as f:
            yaml.dump(memory, f, allow_unicode=True, default_flow_style=False)
        DataLoader._logger.info("矩阵选择已保存至记忆文件")

    def _load_matrix_memory(self, candidates: list) -> Optional[list]:
        """从记忆文件加载矩阵选择，若记忆中的文件名仍在 candidates 中则返回对应索引"""
        memory_path = os.path.join(CONFIG_DIR, "group_memory.yaml")
        if not os.path.exists(memory_path):
            return None
        with open(memory_path, "r", encoding="utf-8") as f:
            memory = yaml.safe_load(f) or {}
        gse = self.cfg.gse_id
        saved = memory.get(gse, {}).get("matrices")
        if not saved:
            return None
        candidate_basenames = [os.path.basename(c) for c in candidates]
        indices = []
        for name in saved:
            if name in candidate_basenames:
                indices.append(candidate_basenames.index(name))
            else:
                DataLoader._logger.warning(
                    f"记忆中的矩阵文件 {name} 在当前候选列表中不存在，将回退到手动选择。"
                )
                return None
        return indices if indices else None


if __name__ == "__main__":
    test_gse_id = "GSE300437"
    test_tar_gene = "Polb"
    test_cfg = Config(tar_gene=test_tar_gene, gse_id=test_gse_id)
    loader = DataLoader(test_cfg)
    data_pack = loader.loader()
