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
        self._batch_exp_groups = []  # 批处理队列: [(短标签, 原始值), ...]
        self._merged_expr = None  # SuperSeries 合并后的表达矩阵（含 SYMBOL 列）
        self._ensembl_to_symbol = None  # Ensembl ID → 基因符号映射

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

        # SuperSeries: diff/wgcna 模式下在此合并子系列矩阵（hilo 已在 _prepare_hilo_group 中完成）
        if self._merged_expr is None:
            is_superseries = (
                any('/' in k for k in self.downloaded_data)
                and 'series_id' in raw_meta.columns
            )
            if is_superseries:
                merged = self._merge_superseries_matrices(raw_meta)
                from modules.calculater import map_ensembl_to_symbol
                self._merged_expr = map_ensembl_to_symbol(merged)
                if 'SYMBOL' in self._merged_expr.columns:
                    self._ensembl_to_symbol = self._merged_expr['SYMBOL'].dropna().to_dict()
                    self._merged_expr = self._merged_expr.set_index('SYMBOL')
                self._logger.info("SuperSeries 合并矩阵已生成，基因符号映射完成")

        chosen_meta = self._chosen_meta
        if chosen_meta is None or chosen_meta.empty:
            raise RuntimeError("数据分组失败，无法构建数据包")

        chosen_meta = chosen_meta.copy()
        if self.cfg.analysis_mode in ("diff", "wgcna") and self._group_col is not None and self._group_mapping:
            chosen_meta["group"] = chosen_meta[self._group_col].apply(self._map_group)
            chosen_meta = chosen_meta[chosen_meta["group"].notna()]
        elif self._group_col is not None and self._group_mapping:
            chosen_meta["group"] = chosen_meta[self._group_col].apply(self._map_group)

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

        if self._batch_exp_groups:
            meta_matrix_pack["_batch_exp_groups"] = self._batch_exp_groups

        if self._merged_expr is not None:
            meta_matrix_pack["expr_matrix"] = self._merged_expr

        meta_matrix_pack["_organism"] = self._detect_organism(raw_meta)
        if self._ensembl_to_symbol is not None:
            meta_matrix_pack["_ensembl_to_symbol"] = self._ensembl_to_symbol

        if self.cfg.storage:
            self._save_pack(meta_matrix_pack)

        return meta_matrix_pack

    def _prepare_default_meta(self, meta: pd.DataFrame) -> None:
        """拷贝元数据，并尝试从已有的 diff/hilo pack 中继承分组信息"""
        self._chosen_meta = meta.copy()
        self._group_col = None
        self._group_mapping = {}
        self._inherit_group_from_sibling()

    def _inherit_group_from_sibling(self) -> None:
        """从同目录下已有的 diff_group / hilo_group pack 中继承分组信息"""
        import pickle
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        for sibling_key in ("diff_group", "hilo_group"):
            sibling_path = os.path.join(
                data_dir, "pkl", f"{self.cfg.gse_id}_{sibling_key}_processed_pack.pkl"
            )
            if not os.path.exists(sibling_path):
                continue
            try:
                with open(sibling_path, 'rb') as f:
                    sibling_pack = pickle.load(f)
            except Exception:
                continue
            group_info = sibling_pack.get("group_info") if isinstance(sibling_pack, dict) else None
            if not group_info:
                continue
            group_col = group_info["group_col"]
            if group_col not in self._chosen_meta.columns:
                continue
            self._group_col = group_col
            self._group_mapping = group_info["mapping"]
            self._logger.info(
                f"从已有的 {sibling_key} pack 中继承了分组信息"
                f"（{self._group_col}: {list(self._group_mapping.keys())}）"
            )
            return

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
        is_superseries = (
            any('/' in k for k in self.downloaded_data)
            and 'series_id' in meta.columns
        )
        if is_superseries:
            merged = self._merge_superseries_matrices(meta)
            from modules.calculater import map_ensembl_to_symbol
            merged = map_ensembl_to_symbol(merged)
            if 'SYMBOL' in merged.columns:
                self._ensembl_to_symbol = merged['SYMBOL'].dropna().to_dict()
                merged = merged.set_index('SYMBOL')
            self._merged_expr = merged  # 存入 pack 供后续分析和绘图使用
            target_vec = self._fetch_target_gene_vector(merged, self.cfg.tar_gene)
            if target_vec.empty:
                raise ValueError(f"未能在合并矩阵中找到目标基因：{self.cfg.tar_gene}")
            common_samples = [c for c in merged.columns if c in meta.index]
            if not common_samples:
                raise ValueError("合并矩阵中无有效样本，无法执行 hilo 分组")
            self._apply_hilo_threshold(meta, target_vec, common_samples)
            return

        expr_file = self._select_hilo_expression_file()
        expr_df = self._load_expression_file(expr_file)
        expr_df = self._rename_expr_columns_by_meta_order(expr_df, meta)

        common_samples = [col for col in expr_df.columns if col in meta.index]
        if not common_samples:
            raise ValueError("表达矩阵列名与元数据索引无法匹配，无法执行 hilo 分组")

        from modules.calculater import map_ensembl_to_symbol
        expr_df = map_ensembl_to_symbol(expr_df)
        if 'SYMBOL' in expr_df.columns:
            self._ensembl_to_symbol = expr_df['SYMBOL'].dropna().to_dict()
            expr_df = expr_df.set_index('SYMBOL')
        self._merged_expr = expr_df  # 存入 pack 供后续分析和绘图使用
        target_vec = self._fetch_target_gene_vector(expr_df, self.cfg.tar_gene)
        if target_vec.empty:
            raise ValueError(f"未能在表达矩阵中找到目标基因：{self.cfg.tar_gene}")

        self._apply_hilo_threshold(meta, target_vec, common_samples)

    def _merge_superseries_matrices(self, meta: pd.DataFrame) -> pd.DataFrame:
        """SuperSeries: 合并所有子系列表达矩阵并归一化基因索引"""
        merged = None
        for datafile_name, file_path in self.downloaded_data.items():
            sub_series_id = datafile_name.split('/')[0] if '/' in datafile_name else None
            if sub_series_id is None:
                continue

            sub_meta = meta[meta['series_id'].str.contains(sub_series_id, na=False)]
            if sub_meta.empty:
                self._logger.warning(f"未找到子系列 {sub_series_id} 的元数据样本，跳过")
                continue

            expr_df = self._load_expression_file(file_path)
            expr_df = self._rename_expr_columns_by_meta_order(expr_df, sub_meta)

            from modules.calculater import normalize_gene_index
            expr_df = normalize_gene_index(expr_df)

            merged = expr_df if merged is None else pd.concat([merged, expr_df], axis=1, join='inner')

        if merged is None or merged.empty:
            raise ValueError("SuperSeries 子系列矩阵合并失败，无法执行 hilo 分组")
        self._logger.info(f"hilo SuperSeries 合并: {merged.shape[0]} 基因 × {merged.shape[1]} 样本")
        return merged

    def _apply_hilo_threshold(self, meta: pd.DataFrame, target_vec: pd.Series, common_samples: list) -> None:
        """按目标基因中值将样本分为 High/Low 组"""
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
        """标准化标签为小写字符串列表，自动扁平化嵌套结构。"""
        if label is None:
            return []
        if isinstance(label, str):
            return [label.lower()]
        # 扁平化：Hydra 可能将 YAML 列表包裹成嵌套列表
        flat = []
        for item in label:
            if item is None:
                continue
            if isinstance(item, (list, tuple)):
                flat.extend(str(x).lower() for x in item if x is not None)
            else:
                flat.append(str(item).lower())
        return flat

    def _auto_group_division(self, meta: pd.DataFrame) -> None:
        """按配置自动化分组"""
        group_select_col = self.cfg.group_select_col
        control_label = self.cfg.control_label
        exp_label = self.cfg.exp_label

        control_labels = self._normalize_labels(control_label)
        exp_labels = self._normalize_labels(exp_label)

        # 肝纤维化模式：实验组优先匹配 MCD / CDAHFD
        if self.cfg.exp_type == "Fibrosis":
            exp_labels = self._normalize_labels(["MCD", "CDAHFD"])

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

        # 同时匹配到 MCD 和 CDAHFD 时的处理（记忆 + 交互选择）
        if self.cfg.exp_type == "Fibrosis" and len(exp_groups) > 1:
            exp_groups = self._resolve_fibrosis_choice(exp_groups)

        if control_groups and exp_groups:
            self._apply_auto_group(meta, group_select_col, control_groups, exp_groups)
            return

        # 首选列匹配失败，回退尝试 'title' 列（与交互模式默认列一致）
        fallback_col = "title"
        if group_select_col != fallback_col and fallback_col in meta.columns:
            fb_unique = meta[fallback_col].unique()
            fb_controls = [g for g in fb_unique if any(label in str(g).lower() for label in control_labels)]
            fb_exps = [g for g in fb_unique if any(label in str(g).lower() for label in exp_labels)]
            if fb_controls and fb_exps:
                self._logger.info(
                    f"在 '{group_select_col}' 中未匹配到分组，"
                    f"从 '{fallback_col}' 列自动识别成功"
                )
                if self.cfg.exp_type == "Fibrosis" and len(fb_exps) > 1:
                    fb_exps = self._resolve_fibrosis_choice(fb_exps)
                if fb_exps:
                    self._apply_auto_group(meta, fallback_col, fb_controls, fb_exps)
                    return

        self._logger.warning("未能自动识别到组别标签，进入交互式分组流程")
        self._manual_group_division(meta)

    def _apply_auto_group(self, meta: pd.DataFrame, group_select_col: str,
                          control_groups: list, exp_groups: list) -> None:
        """应用自动识别到的分组设置。"""
        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"
        self._group_col = group_select_col
        self._chosen_meta = meta[meta[group_select_col].isin(control_groups + exp_groups)]
        self._group_mapping = {"Control": control_groups, exp_type: exp_groups}

        if self.cfg.exp_type == "Fibrosis":
            self._batch_exp_groups = self._build_batch_tuples(exp_groups)
            if len(exp_groups) == 2:
                self._logger.info(f"检测到两个 Fibrosis 实验组: {exp_groups}，将分别进行对比分析")

        if getattr(self.cfg, "group_memory_enabled", False):
            self._save_group_memory(group_select_col, control_groups, exp_groups)

        self._logger.info(
            f"自动分组完成: 列='{group_select_col}', "
            f"Control={control_groups}, Exp={exp_groups}"
        )

    def _resolve_fibrosis_choice(self, exp_groups: list) -> list:
        """Fibrosis 模式匹配到多个实验组时，按类别（MCD/CDAHFD）通过记忆或交互选择确认。"""
        # 将原始分组值归类为 MCD / CDAHFD / 其他
        categories: dict[str, list] = {}
        for g in exp_groups:
            g_lower = str(g).lower()
            if "cdahfd" in g_lower:
                categories.setdefault("CDAHFD", []).append(g)
            elif "mcd" in g_lower:
                categories.setdefault("MCD", []).append(g)
            else:
                categories.setdefault(str(g), []).append(g)

        cat_names = list(categories.keys())
        if len(cat_names) <= 1:
            return exp_groups

        # 尝试从记忆中恢复
        memory_path = os.path.join(CONFIG_DIR, "group_memory.yaml")
        saved_choice = None
        if os.path.exists(memory_path):
            with open(memory_path, "r", encoding="utf-8") as f:
                memory = yaml.safe_load(f) or {}
            gse = self.cfg.gse_id
            mode = self.cfg.analysis_mode
            saved_choice = memory.get(gse, {}).get(mode, {}).get("_fibrosis_exp_choice")

        if saved_choice:
            if saved_choice == "__both__":
                return exp_groups
            if saved_choice in categories:
                self._logger.info(
                    f"已从记忆中恢复 Fibrosis 实验组选择: {saved_choice}"
                    f"（{len(categories[saved_choice])} 个样本）"
                )
                return categories[saved_choice]

        # 交互选择 — 按类别而非逐个标题展示
        print(f"\n检测到多个 Fibrosis 实验组类别：")
        for i, name in enumerate(cat_names):
            print(f"  [{i}] 仅分析 {name}（{len(categories[name])} 个样本）")
        both_idx = len(cat_names)
        print(f"  [{both_idx}] 两者都分析（先后进行对比分析）")

        choice = parse_user_input(
            prompt=f"请输入选项序号 (0-{both_idx}):",
            max_index=both_idx,
        )

        if not choice:
            return []

        if choice[0] == both_idx:
            if getattr(self.cfg, "group_memory_enabled", False):
                self._save_fibrosis_exp_choice("__both__")
            return exp_groups

        chosen_cat = cat_names[choice[0]]
        if getattr(self.cfg, "group_memory_enabled", False):
            self._save_fibrosis_exp_choice(chosen_cat)
        return categories[chosen_cat]

    def _manual_group_division(self, meta: pd.DataFrame) -> None:
        """交互式分组，若用户切换列则优先尝试自动分组。"""
        res1 = self._manual_group_select(meta, group_label="Control")
        group_col = res1["current_col"]

        # 用户在交互中切换了列 → 对新列尝试自动分组，失败再回退手动
        if group_col != "title" and not getattr(self, "_manual_auto_attempted", False):
            self._manual_auto_attempted = True
            try:
                self.cfg.group_select_col = group_col
                self._auto_group_division(meta)
                return
            finally:
                self._manual_auto_attempted = False

        control_values = [res1["unique_groups"][i] for i in res1["group_indices"]]

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

    @staticmethod
    def _build_batch_tuples(exp_groups: list) -> list:
        """将原始分组值映射为 (短标签, 原始值) 二元组，用于批处理文件命名和过滤。"""
        result = []
        for g in exp_groups:
            g_lower = str(g).lower()
            if "cdahfd" in g_lower:
                result.append(("CDAHFD", g))
            elif "mcd" in g_lower:
                result.append(("MCD", g))
            else:
                result.append((g, g))
        return result

    def _save_fibrosis_exp_choice(self, choice_value: str) -> None:
        """保存 Fibrosis 实验组选择到记忆文件（单选值或 '__both__'）。"""
        memory_path = os.path.join(CONFIG_DIR, "group_memory.yaml")
        memory = {}
        if os.path.exists(memory_path):
            with open(memory_path, "r", encoding="utf-8") as f:
                memory = yaml.safe_load(f) or {}
        gse = self.cfg.gse_id
        mode = self.cfg.analysis_mode
        memory.setdefault(gse, {}).setdefault(mode, {})["_fibrosis_exp_choice"] = choice_value
        with open(memory_path, "w", encoding="utf-8") as f:
            yaml.dump(memory, f, allow_unicode=True, default_flow_style=False)
        self._logger.info(f"Fibrosis 实验组选择已记忆: {choice_value}")

    def rebuild_group_for_batch(self, prev_pack: dict, next_label: str, next_value) -> dict:
        """用批处理队列中的下一个实验组重建 pack 的分组信息。

        Args:
            prev_pack: 上一个 pack（需包含 meta_full 和表达矩阵）
            next_label: 短标签，如 "MCD"
            next_value: 原始分组值，用于过滤 meta

        Returns:
            新 pack，_batch_exp_groups 已弹出当前项
        """
        meta_full = prev_pack["meta_full"].copy()
        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"

        control_values = self._group_mapping.get("Control", [])
        self._group_mapping = {"Control": control_values, exp_type: [next_value]}

        chosen_meta = meta_full[
            meta_full[self._group_col].isin(control_values + [next_value])
        ].copy()
        chosen_meta["group"] = chosen_meta[self._group_col].apply(self._map_group)
        chosen_meta = chosen_meta[chosen_meta["group"].notna()]

        new_pack = {"meta": chosen_meta, "meta_full": meta_full}
        for key, val in prev_pack.items():
            if key in ("meta", "meta_full", "group_info", "_batch_exp_groups"):
                continue
            new_pack[key] = val

        new_pack["group_info"] = {
            "group_col": self._group_col,
            "mapping": self._group_mapping,
        }
        remaining = prev_pack.get("_batch_exp_groups", [])
        if len(remaining) > 1:
            new_pack["_batch_exp_groups"] = remaining[1:]

        if self.cfg.storage:
            self._save_pack(new_pack)
        return new_pack

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
        from modules.calculater import get_sample_columns
        return get_sample_columns(df)

    def _rename_expr_columns_by_meta_order(self, expr_df: pd.DataFrame, full_meta: pd.DataFrame):
        from modules.calculater import rename_expr_columns_by_meta_order
        df, _success = rename_expr_columns_by_meta_order(expr_df, full_meta)
        return df

    def _save_pack(self, pack: dict) -> str:
        """保存处理后的数据包"""
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        gse_id = self.cfg.gse_id
        group_key = self.get_pack_group_key(self.cfg.analysis_mode)
        canonical = os.path.join(data_dir, "pkl", f"{gse_id}_{group_key}_processed_pack.pkl")
        save_path = canonical if self.cfg.force else safe_filepath(canonical)
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
    def _detect_organism(meta: pd.DataFrame) -> str:
        """从 GEO 元数据的 organism_ch1 / taxid_ch1 列自动检测物种。

        Returns:
            'human', 'mouse', 或 'rat'（默认 'human'）
        """
        if 'organism_ch1' in meta.columns:
            values = meta['organism_ch1'].dropna().astype(str).str.lower()
            if not values.empty:
                sample = values.iloc[0]
                if 'mus musculus' in sample:
                    return 'mouse'
                if 'homo sapiens' in sample:
                    return 'human'
                if 'rattus norvegicus' in sample:
                    return 'rat'
        if 'taxid_ch1' in meta.columns:
            taxids = meta['taxid_ch1'].dropna().astype(str)
            if not taxids.empty:
                taxid = taxids.iloc[0]
                taxid_map = {'10090': 'mouse', '9606': 'human', '10116': 'rat'}
                if taxid in taxid_map:
                    return taxid_map[taxid]
        return 'human'

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
