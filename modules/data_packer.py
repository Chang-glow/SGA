import os
import yaml
from typing import Optional, Dict

import pandas as pd

from utils import loggers, Config, parse_user_input, safe_filepath, CONFIG_DIR, relpath


class DataPacker:
    """负责将 GEO 下载结果分组并打包成分析器可用的数据包"""
    _logger = loggers.get_logger()

    def __init__(self, cfg: Config, gse, downloaded_data: dict):
        self.cfg = cfg
        self.gse = gse
        self.downloaded_data = downloaded_data
        self._chosen_meta = None
        self._group_mapping = {}
        self._group_col = None

    def build_pack(self) -> dict:
        """构建数据包"""
        if self.gse is None:
            raise ValueError("GSE对象为空，无法构建数据包")

        if not isinstance(self.downloaded_data, dict) or not self.downloaded_data:
            raise ValueError("下载数据为空，无法构建数据包")

        raw_meta = self.gse.phenotype_data
        if raw_meta is None or raw_meta.empty:
            raise ValueError("未能读取到 GEO 元数据")

        if self.cfg.analysis_mode in ("diff", "wgcna"):
            self._prepare_diff_group(raw_meta)
        elif self.cfg.analysis_mode == "hilo":
            self._prepare_hilo_group(raw_meta)
        else:
            self._prepare_default_meta(raw_meta)

        chosen_meta = self._chosen_meta
        if chosen_meta is None or chosen_meta.empty:
            raise RuntimeError("数据分组失败，无法构建数据包")

        chosen_meta = chosen_meta.copy()
        if self.cfg.analysis_mode in ("diff", "wgcna") and self._group_col is not None and self._group_mapping:
            chosen_meta["group"] = chosen_meta[self._group_col].apply(self._map_group)
            chosen_meta = chosen_meta[chosen_meta["group"].notna()]

        meta_matrix_pack = {
            "meta": chosen_meta,
            "meta_full": raw_meta.copy()
        }

        for datafile_name, file_path in self.downloaded_data.items():
            self._logger.info(f"正在将 {datafile_name} 加载至 DataFrame...")
            if not file_path:
                self._logger.warning(f"文件 {datafile_name} 的路径为空，跳过处理。")
                continue
            if not os.path.exists(file_path):
                self._logger.warning(f"找不到本地文件: {file_path}")
                continue

            df_temp = pd.read_csv(file_path, sep="\t", compression="infer", index_col=0)
            meta_matrix_pack[datafile_name] = df_temp

            if self.cfg.strict_mode:
                common_samples = chosen_meta.index.intersection(df_temp.columns)
                matrix_aligned = df_temp.loc[:, common_samples]
                meta_aligned = chosen_meta.loc[common_samples]
                meta_matrix_pack[datafile_name] = {
                    "matrix_aligned": matrix_aligned,
                    "meta_aligned": meta_aligned
                }

        if self._group_col is not None and self._group_mapping:
            meta_matrix_pack["group_info"] = {
                "group_col": self._group_col,
                "mapping": self._group_mapping
            }

        if self.cfg.storage:
            self._save_pack(meta_matrix_pack)

        return meta_matrix_pack

    def _prepare_default_meta(self, meta: pd.DataFrame) -> None:
        """拷贝元数据以免影响数据源"""
        self._chosen_meta = meta.copy()
        self._group_col = None
        self._group_mapping = {}

    def _prepare_diff_group(self, meta: pd.DataFrame) -> None:
        """差异分析分组（支持分组记忆）"""
        if getattr(self.cfg, "group_memory_use", False):
            saved = self._load_group_memory()
            if saved is not None:
                group_col = saved["group_select_col"]
                if group_col not in meta.columns:
                    self._logger.warning(
                        f"记忆的分组列 '{group_col}' 在当前元数据中不存在，回退到交互式选择"
                    )
                else:
                    control_values = saved["control_values"]
                    exp_values = saved["exp_values"]
                    available = set(meta[group_col].unique())
                    missing_ctrl = [v for v in control_values if v not in available]
                    missing_exp = [v for v in exp_values if v not in available]
                    if missing_ctrl or missing_exp:
                        self._logger.warning(
                            f"记忆的分组值在当前元数据中不存在"
                            + (f"（Control 缺失: {missing_ctrl}）" if missing_ctrl else "")
                            + (f"（Exp 缺失: {missing_exp}）" if missing_exp else "")
                            + "，回退到交互式选择"
                        )
                    else:
                        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"
                        self._group_col = group_col
                        self._chosen_meta = meta[meta[group_col].isin(control_values + exp_values)]
                        self._group_mapping = {"Control": control_values, exp_type: exp_values}
                        self._logger.info(
                            f"已应用记忆分组（GSE={self.cfg.gse_id}, "
                            f"列={group_col}, Control={len(control_values)}样本, "
                            f"Exp={len(exp_values)}样本）"
                        )
                        return

        if not self.cfg.group_select_col:
            selected_res = self._manual_group_select(meta)
            selected_group_indices = selected_res["group_indices"]
            unique_groups = selected_res["unique_groups"]
            current_col = selected_res["current_col"]
            if isinstance(selected_group_indices, list):
                target_groups = [unique_groups[i] for i in selected_group_indices]
                condition = meta[current_col].isin(target_groups)
                self._chosen_meta = meta[condition]
                return
        self._auto_group_division(meta)

    def _prepare_hilo_group(self, meta: pd.DataFrame) -> None:
        """高低表达分析分组"""
        expr_file = self._select_hilo_expression_file()
        expr_df = self._load_expression_file(expr_file)
        expr_df = self._rename_expr_columns_by_meta_order(expr_df, meta)

        common_samples = [col for col in expr_df.columns if col in meta.index]
        if not common_samples:
            raise ValueError("表达矩阵列名与元数据索引无法匹配，无法执行 hilo 分组")

        target_vec = self._fetch_target_gene_vector(expr_df, self.cfg.tar_gene)
        if target_vec.empty:
            raise ValueError(f"未能在表达矩阵中找到目标基因：{self.cfg.tar_gene}")

        expr_df = expr_df.loc[:, common_samples]

        numeric_values = pd.to_numeric(target_vec, errors="coerce")
        if numeric_values.dropna().empty:
            raise ValueError(f"目标基因 {self.cfg.tar_gene} 的表达值无有效数值")

        threshold = numeric_values.median()
        self._logger.info(f"hilo 分组阈值设定为目标基因 {self.cfg.tar_gene} 的中值 {threshold}")

        group_labels = numeric_values.apply(lambda x: "High" if x >= threshold else "Low")
        group_meta = meta.loc[common_samples].copy()
        group_meta["group"] = [group_labels.get(sample) for sample in group_meta.index]
        self._chosen_meta = group_meta
        self._group_col = "group"
        self._group_mapping = {
            "Low": group_labels[group_labels == "Low"].index.tolist(),
            "High": group_labels[group_labels == "High"].index.tolist()
        }

    def _map_group(self, value):
        """映射组名和数据"""
        if value in self._group_mapping.get("Control", []):
            return "Control"
        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"
        if value in self._group_mapping.get(exp_type, []):
            return exp_type
        return None

    def _normalize_labels(self, label):
        """标准化标签"""
        if label is None:
            return []
        if isinstance(label, (list, tuple)):
            return [str(item).lower() for item in label if item is not None]
        return [str(label).lower()]

    def _auto_group_division(self, meta: pd.DataFrame) -> None:
        """按配置自动化分组"""
        group_select_col = self.cfg.group_select_col
        control_label = self.cfg.control_label
        exp_label = self.cfg.exp_label

        control_labels = self._normalize_labels(control_label)
        exp_labels = self._normalize_labels(exp_label)

        if group_select_col not in meta.columns:
            self._logger.warning(f"未在元数据中找到指定的分组列 '{group_select_col}'，进入交互式分组流程")
            self._manual_group_division(meta)
            return

        unique_groups = meta[group_select_col].unique()
        control_groups = [
            g for g in unique_groups
            if any(label in str(g).lower() for label in control_labels)
        ]
        exp_groups = [
            g for g in unique_groups
            if any(label in str(g).lower() for label in exp_labels)
        ]

        if control_groups and exp_groups:
            self._chosen_meta = meta[meta[group_select_col].isin(control_groups + exp_groups)]
            exp_type = self.cfg.exp_type if self.cfg.exp_type else 'Experiment'
            self._logger.info(f"自动分组成功, 控制组标签: {control_groups}, 实验组标签: {exp_groups}")
            self._group_mapping = {
                'Control': control_groups,
                exp_type: exp_groups
            }
            self._group_col = group_select_col
            return

        self._logger.warning("未能自动识别到组别标签，进入交互式分组流程")
        self._manual_group_division(meta)

    def _manual_group_division(self, meta: pd.DataFrame) -> None:
        """交互式分组"""
        res1 = self._manual_group_select(meta, group_label="Control")
        control_values = [res1["unique_groups"][i] for i in res1["group_indices"]]
        group_col = res1["current_col"]

        res2 = self._manual_group_select(meta, group_label="Experiment", default_col=group_col)
        exp_values = [res2["unique_groups"][i] for i in res2["group_indices"]]

        if res2["current_col"] != group_col:
            self._logger.warning(
                f"分组列在两次选择中不一致，第一次选择了'{group_col}'，第二次选择了'{res2['current_col']}'，请确保选择的列包含所需的分组信息"
            )

        if not control_values or not exp_values:
            self._logger.error("未选择任何分组标签，无法进行分组")
            raise ValueError("至少需要选择一个控制组标签和一个实验组标签")

        exp_type = self.cfg.exp_type if self.cfg.exp_type else 'Experiment'
        self._group_col = group_col
        self._chosen_meta = meta[meta[group_col].isin(control_values + exp_values)]
        self._group_mapping = {
            'Control': control_values,
            exp_type: exp_values
        }

        if getattr(self.cfg, "group_memory_enabled", False):
            self._save_group_memory(group_col, control_values, exp_values)

    def _manual_group_select(self, meta: pd.DataFrame, group_label: str = None, default_col: str = None) -> dict:
        """交互式分组选择"""
        if default_col and default_col in meta.columns:
            current_col = default_col
        else:
            # 默认使用 "title" 列，如果不存在则使用第一列
            current_col = "title" if "title" in meta.columns else meta.columns[0]

        # 双层循环，外层循环用于重新选择列，内层循环用于选择分组标签
        while True:
            unique_groups = meta[current_col].unique()
            if group_label:
                print(f"\n --- 正在为[{group_label}]选择分组 ---")
            print(f"\n当前查看列:[{current_col}]发现以下样本分组描述")
            for i, group_name in enumerate(unique_groups):
                print(f"[{i}] {group_name}")

            selected_group_indices = parse_user_input(
                prompt=f"请输入 {group_label} 组需要的内容序号(如1:8,11,输入'm'重新选择列):",
                max_index=len(unique_groups) - 1,
                whitelist="m"
            )

            if selected_group_indices == "m":
                print("\n--- 所有元数据列 ---")
                for i, col in enumerate(meta.columns):
                    print(f"[{i}] {col}")

                selected_col = parse_user_input(
                    prompt="请选择(可能)包含分组信息的列序号:",
                    max_index=len(meta.columns) - 1
                )
                current_col = meta.columns[selected_col[0]]
                continue

            return {
                "group_indices": selected_group_indices,
                "unique_groups": unique_groups,
                "current_col": current_col
            }

    def _save_group_memory(self, group_col: str, control_values: list, exp_values: list) -> None:
        """保存分组选择到记忆文件"""
        memory_path = os.path.join(CONFIG_DIR, "group_memory.yaml")
        memory = {}
        if os.path.exists(memory_path):
            with open(memory_path, "r", encoding="utf-8") as f:
                memory = yaml.safe_load(f) or {}

        gse = self.cfg.gse_id
        mode = self.cfg.analysis_mode
        memory.setdefault(gse, {})[mode] = {
            "group_select_col": group_col,
            "control_values": list(control_values),
            "exp_values": list(exp_values),
        }

        with open(memory_path, "w", encoding="utf-8") as f:
            yaml.dump(memory, f, allow_unicode=True, default_flow_style=False)
        self._logger.info(f"分组记忆已保存至 {relpath(memory_path)}")

    def _load_group_memory(self) -> Optional[dict]:
        """从记忆文件加载分组选择，若无记忆则返回 None"""
        memory_path = os.path.join(CONFIG_DIR, "group_memory.yaml")
        if not os.path.exists(memory_path):
            return None
        with open(memory_path, "r", encoding="utf-8") as f:
            memory = yaml.safe_load(f) or {}
        gse = self.cfg.gse_id
        mode = self.cfg.analysis_mode
        return memory.get(gse, {}).get(mode)

    def _select_hilo_expression_file(self) -> str:
        """选择用于 hilo 分组的表达矩阵文件"""
        if len(self.downloaded_data) == 1:
            return next(iter(self.downloaded_data.values()))
        first_file = next(iter(self.downloaded_data.values()))
        self._logger.warning(
            "检测到多个下载文件，hilo 分组将使用第一个文件进行目标基因中值计算。"
        )
        return first_file

    def _load_expression_file(self, file_path: str) -> pd.DataFrame:
        """加载表达矩阵文件并返回 DataFrame"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"表达矩阵文件不存在: {file_path}")
        return pd.read_csv(file_path, sep="\t", compression="infer", index_col=0)

    def _fetch_target_gene_vector(self, df: pd.DataFrame, tar_gene: str) -> pd.Series:
        """从表达矩阵中提取目标基因的表达向量，支持基因名在行索引或特定列中"""
        if not isinstance(df, pd.DataFrame):
            return pd.Series(dtype=float)

        target_gene_upper = str(tar_gene).upper()
        vector = None

        index_upper = df.index.astype(str).str.upper()
        # 首先尝试在行索引中匹配目标基因
        if target_gene_upper in index_upper.values:
            vector = df.loc[df.index[index_upper == target_gene_upper]]
        else:
            potential_columns = [col for col in df.columns if 'SYMBOL' in str(col).upper() or 'GENE' in str(col).upper()]
            for col in potential_columns:
                col_values_upper = df[col].astype(str).str.upper()
                if target_gene_upper in col_values_upper.values:
                    vector = df[col_values_upper == target_gene_upper]
                    break

        # 如果找到的 vector 是 DataFrame，尝试将其转换为 Series
        if vector is not None:
            numeric_df = vector.select_dtypes(include=[float, int]).copy()
            annotation_keywords = [
                'ENSEMBL', 'ENTREZID', 'SYMBOL', 'GENE', 'PROBEID',
                'ID_REF', 'TARGETID', 'DESCRIPTION', 'GENENAME'
            ]
            sample_columns = [
                col for col in numeric_df.columns
                if not any(keyword in str(col).upper() for keyword in annotation_keywords)
            ]
            vector = numeric_df[sample_columns]
            if isinstance(vector, pd.DataFrame):
                if vector.empty:
                    return pd.Series(dtype=float)
                return vector.mean(axis=0)
            return vector

        return pd.Series(dtype=float)

    def _get_sample_columns(self, df: pd.DataFrame) -> list:
        """尝试从表达矩阵中识别样本列，排除常见的注释列和非数值列"""
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
        """尝试通过元数据的索引顺序重新映射表达矩阵的列名，以匹配样本名"""
        sample_columns = self._get_sample_columns(expr_df)
        if len(sample_columns) == len(full_meta.index):
            rename_map = {old: new for old, new in zip(sample_columns, full_meta.index.astype(str))}
            return expr_df.rename(columns=rename_map)

        for meta_col in ['geo_accession', 'title', 'source_name_ch1', 'label_ch1']:
            if meta_col not in full_meta.columns:
                continue
            col_values = full_meta[meta_col].astype(str).tolist()
            if set(sample_columns).issubset(set(col_values)):
                rename_map = {
                    sample_col: str(full_meta.index[col_values.index(sample_col)])
                    for sample_col in sample_columns
                }
                return expr_df.rename(columns=rename_map)
        return expr_df

    def _save_pack(self, pack: dict) -> str:
        """保存处理后的数据包"""
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        gse_id = self.cfg.gse_id
        group_key = self.get_pack_group_key(self.cfg.analysis_mode)
        save_path = safe_filepath(os.path.join(data_dir, "pkl", f"{gse_id}_{group_key}_processed_pack.pkl"))
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        pd.to_pickle(pack, save_path)
        self._logger.info(f"{gse_id}_{group_key}_processed_pack.pkl已存储完成！")
        return save_path

    @staticmethod
    def get_pack_group_key(analysis_mode: str) -> str:
        """相同分组逻辑的模式共用 pack，返回分组键名"""
        if analysis_mode in ("diff", "wgcna"):
            return "diff_group"
        elif analysis_mode == "hilo":
            return "hilo_group"
        return "default"

    @staticmethod
    def resolve_pack_path(data_dir: str, gse_id: str, analysis_mode: str) -> str:
        """查找 pack 文件路径：新命名优先，回退旧命名（向后兼容）"""
        group_key = DataPacker.get_pack_group_key(analysis_mode)
        new_path = os.path.join(data_dir, "pkl", f"{gse_id}_{group_key}_processed_pack.pkl")
        if os.path.exists(new_path):
            return new_path
        old_path = os.path.join(data_dir, "pkl", f"{gse_id}_{analysis_mode}_processed_pack.pkl")
        if os.path.exists(old_path):
            return old_path
        # 最早命名（无模式前缀），如 GSE143318_processed_pack.pkl
        legacy_path = os.path.join(data_dir, "pkl", f"{gse_id}_processed_pack.pkl")
        if os.path.exists(legacy_path):
            return legacy_path
        return new_path
