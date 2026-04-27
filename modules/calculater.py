import os
import scipy

import numpy as np
import pandas as pd

from abc import ABC, abstractmethod
from statsmodels.stats.multitest import fdrcorrection
from typing import Optional

from utils import loggers, Config, DataHandler


def fetch_gene_vector(df, tar_gene) -> pd.Series:
    """提取目标基因向量

    Args:
        df: 被提取的数据库DataFrame
        tar_gene: 提取的目标基因

    Returns:
        vector: 目标基因所属列转换的向量
    """
    logger = loggers.get_logger()

    if not isinstance(df, pd.DataFrame):
        logger.error(f"传入的df不是DataFrame，实际类型: {type(df)}，值: {df}")
        return pd.Series(dtype=float)

    target_gene_upper = str(tar_gene).upper()
    vector = None

    # 1. 尝试从索引查找
    index_upper = df.index.astype(str).str.upper()
    if target_gene_upper in index_upper.values:
        # 找到所有匹配的行（处理重复索引）
        vector = df.loc[df.index[index_upper == target_gene_upper]]

    # 2. 如果索引没找到，尝试从列查找
    else:
        potential_columns = [col for col in df.columns if 'SYMBOL' in col.upper() or 'GENE' in col.upper()]
        for col in potential_columns:
            col_values_upper = df[col].astype(str).str.upper()
            if target_gene_upper in col_values_upper.values:
                # 提取匹配的所有行
                vector = df[col_values_upper == target_gene_upper]
                break  # 找到了就跳出列循环

    # 3. 后处理：清洗聚合
    if vector is not None:
        # 只保留数值列，剔除掉注释列
        vector = vector.select_dtypes(include=[np.number])

        # 处理多行情况
        if isinstance(vector, pd.DataFrame):
            if vector.empty:
                logger.warning(f"基因 {tar_gene} 匹配行为空或无数值数据")
                return pd.Series(dtype=float)
            vector = vector.mean(axis=0)

        return vector

    logger.warning(f"未能在矩阵中找到基因: {tar_gene}")
    return pd.Series(dtype=float)


class Analyzer(ABC):
    """分析从GEO下载的测序数据

    Attributes:
        hfm_dict: 肝纤维化常见标志物字典,全称hepatic_fibrosis_marker_dict
        cfg: 基础配置
    """
    hfm_dict = {
        "Classic": ["Acta2", "Vim", "Col1a1", "Col3a1"],
        "Inflammation": ["Il6", "Tnfa", "Il4", "Il1b"],
        "Signaling_Advanced": ["Tem1", "Arrb1", "Gas6", "Axl", "Pdgfb"],
        "Apoptosis": ["Fas", "Fasl", "Bcl2", "Trp53"],
        "Hedgehog": ["Ptch1", "Smo"]
    }
    _logger = loggers.get_logger()

    def __init__(self, cfg: Config):
        """初始化分析对象

        Args:
            cfg: 基础配置
        """
        self.cfg = cfg
        self._meta_matrix_pack: Optional[dict] = None
        self._corr_result: Optional[pd.DataFrame] = None
        self._diff_result: Optional[pd.DataFrame] = None

    @classmethod
    def create(cls, cfg: Config, data: DataHandler):
        """根据cfg检查使用哪个子类"""
        data_dir = os.path.join(cfg.data_dir, cfg.gse_id)

        pack_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_processed_pack.pkl")
        if os.path.exists(pack_path) and not cfg.debug:
            cls._logger.info(f"发现数据包：{pack_path}，将从数据包中分析")
            return FileAnalyzer(cfg)
        elif data.meta_matrix_pack:
            cls._logger.info("开始从数据包中分析")
            return DataAnalyzer(cfg, data)


    def calculate(self) -> pd.DataFrame:
        """用于调用数据的API"""        
        if self.cfg.analysis_mode == "corr":
            self._logger.info("分析模式设定为相关性分析，将执行相关性分析")
            if self._corr_result is not None and not self._corr_result.empty:
                return self._corr_result
            try:
                return self.corr_analyzer()
            except Exception:
                self._logger.exception("分析失败")
                raise
        elif self.cfg.analysis_mode == "diff":
            self._logger.info("分析模式设定为差异分析，将执行差异分析")
            if self._diff_result is not None and not self._diff_result.empty:
                return self._diff_result            
            try:
                return self.diff_analyzer()
            except Exception:
                self._logger.exception("分析失败")
                raise
        else:
            self._logger.error(f"不支持的分析模式: {self.cfg.analysis_mode}，请检查配置文件")
            raise ValueError(f"不支持的分析模式: {self.cfg.analysis_mode}")

    @staticmethod
    def pearson_analyze(vec1: pd.Series, vec2: pd.Series) -> tuple[Optional[float], Optional[float]]:
        """分析r、p值

        Args:
            vec1: 用于计算相关性的向量1,为目标基因数据向量
            vec2: 用于计算相关性的向量2,为标识物数据向量

        Returns:
            r: 两基因的相关性
            p: 两基因的相关性p-value
        """
        # 保证顺序性
        v1, v2 = vec1.align(vec2, join='inner')

        # 剔除缺失值
        mask = v1.notna() & v2.notna()
        v1, v2 = v1[mask].astype(float), v2[mask].astype(float)

        # 保证长度相同且大于3
        if len(v1) < 3 or v1.std() == 0 or v2.std() == 0:
            return None, None

        r, p = scipy.stats.pearsonr(v1, v2)
        return r, p

    def _is_log(self, df: pd.DataFrame) -> bool:
        """判断数据是否已经 log 转化"""
        numeric_df = df.select_dtypes(include=[np.number])
        if numeric_df.empty:
            return True  # 无数值列，视为已转化（无法判断）
        max_val = numeric_df.max().max()
        # 阈值可根据实际数据调整，默认为50
        return max_val <= self.cfg.log_threshold

    def _clean_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """清洗单个矩阵,去除NaN,Inf和自动监测log"""
        df = df.copy()

        # Inf -> NaN
        df.replace([np.inf, -np.inf], np.nan, inplace=True)

        # drop掉NaN行
        origin_rows = df.shape[0]
        df.dropna(axis=0, how='any', inplace=True)
        if df.shape[0] < origin_rows:
            self._logger.info(f"删除了{origin_rows - df.shape[0]}行NaN/Inf数据")

        # 判断log并转化没有log的矩阵
        if not self._is_log(df):
            self._logger.warning("原始数据未log转换,将执行log2(x+1)转换")
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            df[numeric_cols] = np.log2(df[numeric_cols] + 1)
        else:
            self._logger.debug("数据已log转换,跳过该步骤")

        return df

    def _get_sample_columns(self, df: pd.DataFrame) -> list:
        """识别表达矩阵中的样本列，排除注释列和检测 p-value 列"""
        sample_columns = []
        for col in df.columns:
            if isinstance(col, str):
                lowered = col.lower()
                if any(keyword in lowered for keyword in [
                    'ensembl', 'entrezid', 'symbol', 'genename', 'probeid',
                    'id_ref', 'targetid', 'gene', 'description'
                ]):
                    continue
                if 'detection' in lowered and 'pval' in lowered:
                    continue
                sample_columns.append(col)
            else:
                sample_columns.append(col)
        return sample_columns

    def _rename_expr_columns_by_meta_order(self, expr_df: pd.DataFrame, full_meta: pd.DataFrame) -> pd.DataFrame:
        """尝试按元数据顺序或元数据值将表达矩阵列名映射为 GSM 样本名"""
        sample_columns = self._get_sample_columns(expr_df)
        if len(sample_columns) == len(full_meta.index):
            rename_map = {old: new for old, new in zip(sample_columns, full_meta.index.astype(str))}
            expr_df = expr_df.rename(columns=rename_map)
            return expr_df

        # 如果列名与某个元数据字段值直接匹配，则基于值映射
        for meta_col in ['geo_accession', 'title', 'source_name_ch1', 'label_ch1']:
            if meta_col not in full_meta.columns:
                continue
            col_values = full_meta[meta_col].astype(str).tolist()
            if set(sample_columns).issubset(set(col_values)):
                rename_map = {}
                for sample_col in sample_columns:
                    idx = col_values.index(sample_col)
                    rename_map[sample_col] = str(full_meta.index[idx])
                expr_df = expr_df.rename(columns=rename_map)
                return expr_df

        return expr_df

    def _clean_pack(self, raw_pack: dict) -> dict:
        if self._meta_matrix_pack is None:
            cleaned_pack = {}
            for name, item in raw_pack.items():
                if name in {"meta", "meta_full"}:
                    cleaned_pack[name] = item.copy()
                elif isinstance(item, dict):
                    if 'matrix_aligned' in item and isinstance(item['matrix_aligned'], pd.DataFrame):
                        cleaned_pack[name] = self._clean_dataframe(item['matrix_aligned'])
                    else:
                        # 保留非矩阵形式的字典元数据，如 group_info 或其他配置
                        cleaned_pack[name] = item.copy()
                elif isinstance(item, pd.DataFrame):
                    cleaned_pack[name] = self._clean_dataframe(item)
                else:
                    cleaned_pack[name] = item
            return cleaned_pack

    def roaming_data(self) -> dict:
        """缓存读取数据

        Returns:
            self._meta_matrix_pack: 读取到的{文件名：矩阵}字典，存为属性
        """
        if self._meta_matrix_pack is None:
            self._meta_matrix_pack = self._load_data()
        return self._meta_matrix_pack

    @property
    def significant(self):
        """最显著值"""
        if self.cfg.analysis_mode == "corr":
            if self._corr_result is None:
                df = self.calculate()
            else:
                df = self._corr_result
            significant_findings = df[df['P_value'] < 0.05].sort_values('R')
        elif self.cfg.analysis_mode == "diff":
            if self._diff_result is None:
                df = self.calculate()
            else:
                df = self._diff_result
            significant_findings = df[df['padj'] < 0.05].sort_values('log2FC')
        else:
            self._logger.error(f"不支持的分析模式: {self.cfg.analysis_mode}，请检查配置文件")
            raise ValueError(f"不支持的分析模式: {self.cfg.analysis_mode}")
        return significant_findings

    def _corr_calculater(self, meta_matrix_pack: dict) -> Optional[pd.DataFrame]:
        """计算目标基因与常见标识物的相关性

        Args:
            meta_matrix_pack: 打包的测序矩阵数据，格式为{文件名: 对应矩阵数据}

        Returns:
            gene_corr_table: DataFrame格式的分析结果
        """
        # 读取配置
        tar_gene = self.cfg.tar_gene

        # 提取pack里的数据计算
        results_list = []
        df: pd.DataFrame
        for name, df in meta_matrix_pack.items():
            # 避开元数据矩阵
            if name == "meta":
                continue

            self._logger.info(f"--- 当前处理数据{name} ---")
            self._logger.info(f"提取目标基因{tar_gene}的数据中...")
            target_vec = fetch_gene_vector(df, tar_gene)

            self._logger.info("将以常见标识物分类进行计算并储存")
            for category, gene_list in self.hfm_dict.items():
                for gene in gene_list:
                    self._logger.debug(f"提取标识物基因{gene}的数据中...")
                    marker_vec = fetch_gene_vector(df, gene)

                    r, p = None, None
                    if not marker_vec.empty:
                        self._logger.debug("数据提取完成，计算相关性中...")
                        r, p = self.pearson_analyze(target_vec, marker_vec)

                    if r is not None:
                        self._logger.debug("相关性计算完成！")
                        results_list.append({
                            "Matrix": name,
                            "Category": category,
                            "Gene": gene,
                            "R": r,
                            "P_value": p
                        })

        # 转化为DataFrame，以便画图
        if results_list:
            self._logger.info("相关性计算完成！")
            gene_corr_table = pd.DataFrame(results_list)
        else:
            return None

        return gene_corr_table

    def corr_analyzer(self):
        """调用相关性分析主pipeline,串联数据读取和分析得出相关性

        Returns:
            gene_corr_table: 目标基因和常见标识基因相关性分析结果及数据
        """
        # 查找缓存
        data_dir = self.cfg.data_dir
        gse_id = self.cfg.gse_id

        gene_corr_path = os.path.join(data_dir, "pkl", f"{gse_id}_correlation_summary.pkl")
        if os.path.exists(gene_corr_path) and self.cfg.analysis_mode == "corr" and not self.cfg.debug:
            self._logger.info(f"发现现存分析结果：{gene_corr_path}，跳过相关性分析")
            self._corr_result = pd.read_pickle(gene_corr_path)
            return self._corr_result
        
        # 读取数据
        if not self._meta_matrix_pack:
            self._logger.info("正在读取数据...")
            self.roaming_data()
        meta_matrix_pack = self._meta_matrix_pack
        if meta_matrix_pack:
            self._logger.info("读取成功，将继续分析")

        # 计算相关性
        result_df = self._corr_calculater(meta_matrix_pack)
        if result_df is not None and not result_df.empty:
            self._corr_result = result_df
            if self.cfg.storage:
                self._data_storage(result_df, "pkl")
                self._data_storage(result_df, "csv")
            return result_df
        else:
            self._logger.error("结果矩阵为空")
            raise

    def diff_analyzer(self) -> Optional[pd.DataFrame]:
        """差异分析主pipeline"""
        # 查找缓存
        data_dir = self.cfg.data_dir
        gse_id = self.cfg.gse_id

        gene_diff_path = os.path.join(data_dir, "pkl", f"{gse_id}_differential_summary.pkl")
        if os.path.exists(gene_diff_path) and self.cfg.analysis_mode == "diff" and not self.cfg.debug:
            self._logger.info(f"发现现存分析结果：{gene_diff_path}，跳过差异分析")
            self._diff_result = pd.read_pickle(gene_diff_path)
            return self._diff_result
        
        # 读取数据
        if not self._meta_matrix_pack:
            self._logger.info("正在读取数据...")
            self.roaming_data()
        pack = self._meta_matrix_pack
        meta = pack.get("meta")

        if 'group' in meta.columns and not meta['group'].isnull().all():
            self._logger.info("读取成功，将继续分析")
        else:
            self._logger.error("未找到分组信息，无法进行差异分析")
            raise ValueError("未找到分组信息，无法进行差异分析")
        
        # 根据分组进行差异分析
        # 1. 选择表达矩阵
        expr_keys = [k for k, v in pack.items() if k not in {"meta", "meta_full", "group_info"}]
        if not expr_keys:
            raise KeyError("No expression matrix found in data pack")
        expr_df = pack[expr_keys[0]]
        if isinstance(expr_df, dict) and 'matrix_aligned' in expr_df:
            expr_df = expr_df['matrix_aligned']

        # 2. 根据group拆分
        control_samples = meta[meta['group'] == 'Control'].index.tolist()
        fib_samples = meta[meta['group'] == 'Fibrosis'].index.tolist()
        if not control_samples or not fib_samples:
            self._logger.error("分组信息不完整,无法找到Control或Fibrosis组的样本")
            raise ValueError("分组信息不完整,无法找到Control或Fibrosis组的样本")
        
        # 3.对齐矩阵样本
        common_samples = list(dict.fromkeys(control_samples + fib_samples))

        # 如果样本名在行索引中而不是列名中，则转置矩阵
        if set(common_samples).issubset(expr_df.index):
            expr_df = expr_df.T

        missing_samples = [s for s in common_samples if s not in expr_df.columns]
        if missing_samples:
            self._logger.warning(f"样本名未在表达矩阵列中找到: {missing_samples[:10]}{'...' if len(missing_samples)>10 else ''}")
            if 'meta_full' in pack:
                expr_df = self._rename_expr_columns_by_meta_order(expr_df, pack['meta_full'])
                missing_samples = [s for s in common_samples if s not in expr_df.columns]
                if not missing_samples:
                    expr_df = expr_df.loc[:, common_samples]
                    self._logger.info("已按原始元数据顺序将表达矩阵列映射为 GSM 样本名")
                else:
                    self._logger.error("样本名与原始元数据顺序无法完全对应，无法重新映射表达矩阵")
                    raise KeyError("样本名与表达矩阵列不匹配，请检查数据和元数据")
            else:
                self._logger.error("无法从表达矩阵列名中推断Control/Fibrosis分组，请检查数据")
                raise KeyError("样本名与表达矩阵列不匹配，请检查数据和元数据")
        else:
            expr_df = expr_df.loc[:, common_samples]

        control_samples = expr_df[control_samples]
        fib_samples = expr_df[fib_samples]

        # 4.计算统计量
        results_list = []
        for gene in expr_df.index:
            c_vals = control_samples.loc[gene].dropna()
            f_vals = fib_samples.loc[gene].dropna()
            if len(c_vals) < 3 or len(f_vals) < 3:
                continue
            # t 检验
            _, p_val = scipy.stats.ttest_ind(c_vals, f_vals, equal_var=False)
            # log2FC 均值差
            if self._is_log(expr_df):
                log2fc = f_vals.mean() - c_vals.mean()
            else:
                log2fc = np.log2(f_vals.mean() + 1) - np.log2(c_vals.mean() + 1)
            results_list.append({
                "Gene": gene,
                "log2FC": log2fc,
                "P_value": p_val
            })

        # 5. FDR校正
        _diff_result = pd.DataFrame(results_list)
        _, _diff_result['padj'] = fdrcorrection(_diff_result['P_value'].fillna(1))

        self._diff_result = _diff_result
        if self.cfg.storage:
            self._data_storage(_diff_result, "pkl")
            self._data_storage(_diff_result, "csv")
        return _diff_result

    def _data_storage(self, result_df: pd.DataFrame, save_format: str):
        """存储DataFrame数据至pkl和csv"""
        # 读取配置
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        gse_id = self.cfg.gse_id

        # 映射文件格式向格式名与后缀名
        if save_format == "pkl":
            fmt = "pickle"
            ext = "pkl"
        elif save_format == "csv":
            fmt = ext = "csv"
        else:
            raise ValueError(f"不支持的格式:{save_format}")
        
        if self.cfg.analysis_mode == "corr":
            file_name = f"{gse_id}_correlation_summary.{ext}"
        else:
            file_name = f"{gse_id}_differential_summary.{ext}"

        storage_path = os.path.join(data_dir, ext, file_name)
        os.makedirs(os.path.dirname(storage_path), exist_ok=True)

        save_mode = getattr(result_df, f"to_{fmt}")

        try:
            if fmt == "csv":
                save_mode(storage_path, index=False)
            else:
                save_mode(storage_path)
            self._logger.info(f"{fmt}文件已保存至{storage_path}")
        except PermissionError as e:
            self._logger.error(f"无法写入{fmt}，错误：{e}")
        except OSError as e:
            self._logger.error(f"系统错误：{e}")
        except Exception as e:
            self._logger.error(f"存储{fmt}时发生未知错误：{e}")

    @abstractmethod
    def _load_data(self) -> dict:
        pass


class DataAnalyzer(Analyzer):
    """基于直接读取DataFrame数据的分析流程"""
    def __init__(self, cfg: Config, data: DataHandler):
        super().__init__(cfg)
        self.meta_matrix_pack = data.meta_matrix_pack
        self._logger.info("将从打包的pack中读取数据")

    def _load_data(self) -> dict:
        self.meta_matrix_pack = self._clean_pack(self.meta_matrix_pack)
        return self.meta_matrix_pack


class FileAnalyzer(Analyzer):
    """基于读取pkl文件的数据分析流程"""
    def __init__(self, cfg: Config):
        super().__init__(cfg)
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        self.data_path = os.path.join(data_dir, "pkl", f"{self.cfg.gse_id}_processed_pack.pkl")
        self._logger.info("将从打包的pkl中读取数据")

    def read_pkl(self) -> dict:
        """读取pkl文件

        Returns:
            data_dict: 从pkl文件中提取的数据包
        """
        try:
            if not isinstance(self.data_path, str):
                raise TypeError("请输入正确的路径")
            whole_path = self.data_path

            if not os.path.exists(whole_path):
                raise FileNotFoundError(f"当前查找路径：{whole_path},该路径未找到文件")

            # 读取pkl中字典
            data_dict = pd.read_pickle(whole_path)

            if not isinstance(data_dict, dict):
                raise TypeError("读取到的结果并非字典")
            if not data_dict:
                raise ValueError("文件为空，未查询到有效数据")
            return data_dict

        except TypeError as e:
            self._logger.error(f"【类型错误】：{e}")
            raise
        except FileNotFoundError as e:
            self._logger.error(f"【路径错误】：{e}")
        except ValueError as e:
            self._logger.warning(f"【数据错误】：{e}")
        except Exception as e:
            self._logger.error(f"【未知错误】：{e}")

    def _load_data(self) -> dict:
        if not self._meta_matrix_pack:
            raw_pack = self.read_pkl()
            self._meta_matrix_pack = self._clean_pack(raw_pack)
        return self._meta_matrix_pack


if __name__ == "__main__":
    test_gse_id = "GSE300437"
    test_tar_gene = "Polb"
    test_cfg = Config(tar_gene=test_tar_gene, gse_id=test_gse_id)
    test_analyzer = FileAnalyzer(test_cfg)
    test_gene_corr = test_analyzer.calculate()
    print(test_analyzer.significant)
