import os
from abc import ABC, abstractmethod
import numpy as np
from typing import Optional
from scipy.stats import ttest_ind, mannwhitneyu, spearmanr

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from utils import Config, DataHandler, loggers, safe_filepath, FIGURE_DIR, RESULT_DIR, parse_tar_genes
from modules.calculater import fetch_gene_vector
from modules.data_packer import DataPacker


class FigurePlotter(ABC):
    """绘图类，用于将基因相关性分析结果汇成含拟合线与误差线的散点图

    Attributes:
        cfg: 基础配置
    """
    _logger = loggers.get_logger()

    def __init__(self, cfg: Config):
        self.cfg = cfg
        self._gene_corr_table: Optional[pd.DataFrame] = None
        self._gene_diff_table: Optional[pd.DataFrame] = None
        self._gene_enrich_table: Optional[pd.DataFrame] = None
        self._gene_immune_table: Optional[pd.DataFrame] = None
        self._meta_matrix_pack: Optional[dict] = None
        # 内部判断是否画图
        self._plotter = False

    @classmethod
    def create(cls, cfg: Config, data: DataHandler):
        """根据cfg检查使用哪个子类"""
        data_dir = os.path.join(cfg.data_dir, cfg.gse_id)
        data_pack_path = DataPacker.resolve_pack_path(data_dir, cfg.gse_id, cfg.analysis_mode)
        gene_corr_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_correlation_summary.pkl")
        gene_diff_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_differential_summary.pkl")
        gene_hilo_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_highlow_summary.pkl")
        gene_enrich_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_enrichment_summary.pkl")
        gene_immune_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_immune_summary.pkl")

        is_pack = os.path.exists(data_pack_path)
        if cfg.analysis_mode == "corr":
            is_table = os.path.exists(gene_corr_path)
        elif cfg.analysis_mode == "diff":
            is_table = os.path.exists(gene_diff_path)
        elif cfg.analysis_mode == "hilo":
            is_table = os.path.exists(gene_hilo_path)
        elif cfg.analysis_mode == "enrich":
            is_table = os.path.exists(gene_enrich_path)
        elif cfg.analysis_mode == "immune":
            is_table = os.path.exists(gene_immune_path)
        elif cfg.analysis_mode == "wgcna":
            is_table = True  # WGCNA 不生成图表，跳过
        else:
            cls._logger.warning(f"未知分析模式：{cfg.analysis_mode}，将默认使用数据包画图")
            is_table = False

        if is_pack and is_table:
            return FilePlotter(cfg)
        elif (is_pack and not is_table) or (is_table and not is_pack):
            lack = "pack" if is_table else "table"
            cls._logger.warning(f"{lack}数据缺失，分析可能出错")
            return FilePlotter(cfg)
        else:
            return DataPlotter(cfg, data)

    @abstractmethod
    def fig_plotter(self) -> None:
        pass

    @abstractmethod
    def _load_data(self) -> None:
        pass

    def plotter(self):
        """画图API"""
        if self._plotter:
            pass
        else:
            return self.fig_plotter()

    def _get_vecs(self, df: pd.DataFrame, gene: str) -> tuple[pd.Series, pd.Series]:
        """取得目标基因和对比基因的向量"""
        # 读取配置
        tar_gene = self.cfg.tar_gene

        # 调用外部函数取向量
        x_vec = fetch_gene_vector(df, tar_gene=tar_gene)
        y_vec = fetch_gene_vector(df, tar_gene=gene)
        return x_vec, y_vec
    
    def figplotter(self) -> None:
        """绘图主流程，包含条件筛选和循环画图"""
        if self.cfg.analysis_mode == "corr":
            self._logger.info("分析模式设定为相关性分析,将绘制相关性散点图")
            self.bar_plotter()
            self.scatter_plotter()
        elif self.cfg.analysis_mode == "diff":
            self._logger.info("分析模式设定为差异分析,将绘制差异箱线图和热图")
            self.box_plotter()
            self.heatmap_plotter()
        elif self.cfg.analysis_mode == "hilo":
            self._logger.info("分析模式设定为高低表达分析,将绘制箱线图、火山图和热图")
            self.box_plotter()
            self.volcano_plotter()
            self.heatmap_plotter()
        elif self.cfg.analysis_mode == "enrich":
            self._logger.info("分析模式设定为富集分析,将绘制富集条形图、气泡图和 GO 合集条形图")
            self.enrich_bar_plotter()
            self.enrich_dot_plotter()
            self.enrich_go_combined_bar()
            self._save_enrich_plotting_csv()
        elif self.cfg.analysis_mode == "immune":
            self._logger.info("分析模式设定为免疫浸润分析，将绘制堆叠柱状图、组间箱线图、相关性热图")
            self.immune_stacked_bar()
            self.immune_box_plots()
            self.immune_corr_heatmap()
        elif self.cfg.analysis_mode == "wgcna":
            self._logger.info(
                "WGCNA 模式：图表（模块-性状热图、eigengene 图、GO 富集等）"
                "已由阶段2 PyWGCNA 管道生成至 res/figures/wgcna/，阶段3跳过。"
            )
            return
        else:
            self._logger.error(f"未知分析模式：{self.cfg.analysis_mode}，无法绘图")
            raise ValueError(f"未知分析模式：{self.cfg.analysis_mode}")

    def bar_plotter(self) -> None:
        """柱状图pipeline"""
        # 加载数据
        self._load_data()
        self._logger.info("构建柱状图绘图数据中...")
        targets, colors = self._mapping_corr_to_bar_data()
        self._logger.info("绘制相关性分析柱状图中...")
        self._bar_plot(targets["Gene"], targets["R"], colors)
        bar_name = self._build_fig_name("correlation_barplot.png")
        if os.path.exists(self._build_fig_path(bar_name)):
            self._logger.info("相关性分析柱状图绘制完成！")
        else:
            self._logger.error("相关性分析柱状图绘制失败！")
            raise FileNotFoundError("相关性分析柱状图绘制失败！")

    def _mapping_corr_to_bar_data(self) -> tuple[pd.DataFrame, list]:
        """将相关性分析结果映射成柱状图数据"""
        # 读取配置
        p_thr = self.cfg.p_threshold
        
        # 构建P值条件
        p_condition = self._gene_corr_table["P_value"] < p_thr
        self._logger.info(f"将以p值阈值为{p_thr}为条件筛选因子")
        targets = self._gene_corr_table[p_condition]

        # 排序映射颜色
        targets = targets.sort_values(by="R", ascending=False)
        colors = ["red" if x > 0 else "blue" for x in targets["R"]]
        
        return targets, colors

    def _bar_plot(self, x: pd.Series, y: pd.Series, colors: list) -> None:
        """画柱状图并存储"""
        plt.figure(figsize=(10, 6))
        sns.barplot(x=x, y=y, hue=x, palette=colors, legend=False)
        plt.title(f"Correlation of {self.cfg.tar_gene} with other genes", fontsize=12)
        plt.xlabel("Gene")
        plt.ylabel("Correlation Coefficient (R)")
        plt.xticks(rotation=90)
        fig_name = self._build_fig_name('correlation_barplot.png')
        self._save_plot(fig_name)

    def scatter_plotter(self):
        """散点图pipeline"""
        self._logger.info("正在构建相关性筛选条件...")
        p_condition, sign_condition = self._mapping_corr_condition()
        self._filter_corr_data_to_plot(p_condition, sign_condition)

    def _mapping_corr_condition(self) -> tuple[pd.Series, pd.Series, str]:
        # 读取配置
        p_thr, signs = self.cfg.p_threshold, self.cfg.signs

        # 构建P值条件
        p_condition = self._gene_corr_table["P_value"] < p_thr
        # 构建相关性条件(可多选)
        sign_condition = None
        for sign in signs:
            if sign == "negative":
                cond = self._gene_corr_table["R"] < 0
            elif sign == "positive":
                cond = self._gene_corr_table["R"] > 0
            else:
                self._logger.warning(f"将忽略未知符号：{sign}")
                continue
            sign_condition = cond if sign_condition is None else (sign_condition | cond)

        # 构建日志描述映射
        sign_map = {'negative': '负相关 (R < 0)', 'positive': '正相关 (R > 0)'}
        if len(signs) == 1:
            sign_desc = sign_map[signs[0]]
        else:
            sign_desc = "或".join([sign_map[s] for s in signs])

        if p_condition is None or sign_condition is None:
            self._logger.error("配置项缺失有效的p值阈值或相关性取向")
            return
        
        self._logger.info(f"将以\n1,p值阈值为{p_thr}\n2,{sign_desc}相关为条件筛选因子")
        return p_condition, sign_condition

    def _filter_corr_data_to_plot(self, p_condition: Optional[pd.Series] = None, sign_condition: Optional[pd.Series] = None):
        # 加载数据
        self._load_data()

        # 根据符号构建筛选条件
        targets = self._gene_corr_table[p_condition & sign_condition]

        self._logger.info("绘图中...")
        for _, row in targets.iterrows():
            matrix_name = row['Matrix']
            gene_name = row['Gene']
            self._logger.debug(f"当前绘图基因 {gene_name}")

            df = self._meta_matrix_pack[matrix_name]
            x_vec, y_vec = self._get_vecs(df, gene_name)

            self._corr_plot(x_vec, y_vec, row)

        self._logger.info("绘图完成！")

    def _corr_plot(self, x: pd.Series, y: pd.Series, info: pd.DataFrame, plot_type: str = "scatter") -> None:
        """画散点图并存储"""
        # 读取配置
        tar_gene = self.cfg.tar_gene

        plt.figure(figsize=(6, 6))
        
        # 画散点图并自动添加回归线
        sns.regplot(x=x, y=y, ci=95,
                    scatter_kws={'alpha': 0.6, 's': 80, 'color': '#34495e'},
                    line_kws={'color': '#c0392b', 'lw': 2})

        # 标注相关系数和p-value
        matrix_info = os.path.splitext(os.path.splitext(info['Matrix'])[0])[0]
        gene_info = info['Gene']
        title_str = f"{matrix_info}\n{tar_gene} vs {gene_info}\nR={info['R']:.3f}, P={info['P_value']:.4e}"
        plt.title(title_str, fontsize=10)
        plt.xlabel(f"{tar_gene} Expression")
        plt.ylabel(f"{gene_info} Expression")
        fig_name = self._build_fig_name(f"{matrix_info}_{gene_info}_corr_scatter.png")
        self._save_plot(fig_name)

    def _gene_token(self) -> str:
        """返回用于标题/文件名的基因标识。
        多基因时返回 "N genes"，单基因时返回基因名。"""
        genes = parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene)
        return genes[0] if len(genes) == 1 else f"{len(genes)} genes"

    def _build_fig_name(self, suffix: str, gene: str = None) -> str:
        """构建包含分析模式的结果文件名。
        gene: 指定时用该基因名；不指定且多基因时用 multi_{N}genes"""
        if gene:
            token = gene
        elif len(parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene)) > 1:
            token = f"multi_{len(parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene))}genes"
        else:
            token = parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene)[0]
        return f"{self.cfg.gse_id}_{token}_{self.cfg.analysis_mode}_{suffix}"

    def _build_fig_path(self, fig_name: str) -> str:
        mode_dir = os.path.join(FIGURE_DIR, self.cfg.analysis_mode)
        os.makedirs(mode_dir, exist_ok=True)
        return os.path.join(mode_dir, fig_name)

    def _resolve_save_path(self, path: str) -> str:
        """根据 overwrite_figures 决定是否迭代编号"""
        if self.cfg.overwrite_figures:
            return path
        return safe_filepath(path)

    def _save_plot(self, fig_name: str) -> None:
        # 保存图片
        fig_path = self._resolve_save_path(self._build_fig_path(fig_name))
        dpi = getattr(self.cfg, "figure_dpi", 300)
        plt.savefig(fig_path, dpi=dpi, bbox_inches='tight')
        plt.close()

    def box_plotter(self) -> None:
        """箱线图pipeline — 支持多基因"""
        if not self._meta_matrix_pack:
            self._load_data()

        for gene in parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene):
            self._logger.info(f"正在准备 {gene} 的差异分析数据...")
            data_dict = self._prepare_diff_data(gene=gene)
            self._logger.info(f"{gene} 差异分析数据准备完成，正在绘图...")
            fig_name = self._build_fig_name("boxplot.png", gene=gene)
            self._save_box_plot(
                data_dict['x'],
                data_dict['y'],
                data_dict['p_value'],
                data_dict['title'],
                fig_name=fig_name,
            )
        n = len(parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene))
        self._logger.info(f"差异分析箱线图绘制完成！(共 {n} 个基因)")

    def _get_expr_matrix(self) -> pd.DataFrame:
        """从 _meta_matrix_pack 中提取表达矩阵 DataFrame"""
        for key, val in self._meta_matrix_pack.items():
            if key in {'meta', 'meta_full'}:
                continue
            # 兼容 strict_mode 嵌套结构
            if isinstance(val, dict) and 'matrix_aligned' in val:
                return val['matrix_aligned']
            elif isinstance(val, pd.DataFrame):
                return val
        raise KeyError("No expression matrix found in _meta_matrix_pack")

    def volcano_plotter(self) -> None:
        """火山图 — 竖线 ±1，下调蓝色 / 上调红色，标注显著基因名"""
        if not self._meta_matrix_pack:
            self._load_data()

        if self._gene_diff_table is None or self._gene_diff_table.empty:
            self._logger.warning("差异分析结果为空，无法绘制火山图")
            return

        self._logger.info("正在准备火山图数据...")
        data = self._gene_diff_table.copy()

        if "padj" not in data.columns or "log2FC" not in data.columns:
            self._logger.error("差异结果缺少 padj 或 log2FC 列")
            return

        data = data.dropna(subset=["padj", "log2FC"])
        p_thr = self.cfg.p_threshold
        log2fc_thr = getattr(self.cfg, "log2fc_threshold", 1.0)
        if log2fc_thr <= 0:
            log2fc_thr = 1.0

        x = data["log2FC"].values
        y = -np.log10(data["padj"].clip(lower=1e-300).values)

        sig_up = (data["padj"] < p_thr) & (data["log2FC"] >= log2fc_thr)
        sig_down = (data["padj"] < p_thr) & (data["log2FC"] <= -log2fc_thr)
        not_sig = ~(sig_up | sig_down)

        fig, ax = plt.subplots(figsize=(9, 7))

        ax.scatter(x[not_sig], y[not_sig], c="#bdbdbd", s=6, alpha=0.4,
                   edgecolors="none", label="NS", zorder=1)
        ax.scatter(x[sig_down], y[sig_down], c="#3498db", s=10, alpha=0.7,
                   edgecolors="none", label=f"Down (log2FC ≤ -{log2fc_thr})", zorder=2)
        ax.scatter(x[sig_up], y[sig_up], c="#e74c3c", s=10, alpha=0.7,
                   edgecolors="none", label=f"Up (log2FC ≥ {log2fc_thr})", zorder=2)

        # 标注 top N 显著基因
        label_n = getattr(self.cfg, "volcano_top_n_labels", 15)
        top_genes = data[sig_up | sig_down].nsmallest(label_n, "padj")
        for _, row in top_genes.iterrows():
            ax.annotate(
                row["Gene"], (row["log2FC"], -np.log10(max(row["padj"], 1e-300))),
                fontsize=6, alpha=0.85,
                xytext=(3, 3), textcoords="offset points",
            )

        y_ref = -np.log10(p_thr)
        ax.axhline(y_ref, linestyle="--", color="#333333", alpha=0.5, linewidth=0.8)
        ax.axvline(-log2fc_thr, linestyle="--", color="#333333", alpha=0.3, linewidth=0.6)
        ax.axvline(log2fc_thr, linestyle="--", color="#333333", alpha=0.3, linewidth=0.6)

        ax.set_xlabel("log2FC", fontsize=11)
        ax.set_ylabel("-log10(padj)", fontsize=11)
        ax.set_title(f"Volcano plot — {self.cfg.tar_gene}", fontsize=12, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right", framealpha=0.85)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()

        fig_name = self._build_fig_name("volcano.png")
        self._save_plot(fig_name)
        self._logger.info("火山图绘制完成！")

    def heatmap_plotter(self) -> None:
        """热图 pipeline — 聚类 + 分组颜色条 + Z-score + RdBu_r 配色"""
        if not self._meta_matrix_pack:
            self._load_data()

        if self._gene_diff_table is None or self._gene_diff_table.empty:
            self._logger.warning("差异分析结果为空，无法绘制热图")
            return

        self._logger.info("正在准备热图数据...")
        data = self._gene_diff_table.copy()
        expr_df = self._get_expr_matrix()
        meta = self._meta_matrix_pack.get('meta')
        if meta is None or 'group' not in meta.columns:
            self._logger.warning("未找到分组信息，无法绘制热图")
            return

        full_meta = self._meta_matrix_pack.get('meta_full')
        if full_meta is not None:
            expr_df = self._rename_expr_columns_by_meta_order(expr_df, full_meta)

        top_n = getattr(self.cfg, "heatmap_top_n_genes", 20)
        top_genes = data.sort_values('padj').head(top_n)['Gene'].tolist()
        target_genes = parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene)
        missing = [g for g in target_genes if g not in top_genes]
        if missing:
            top_genes = top_genes[:top_n - len(missing)] + missing
        heatmap_df = self._match_top_genes_rows(expr_df, top_genes)
        if heatmap_df.empty:
            self._logger.warning('未找到用于热图的基因表达数据')
            return

        sample_columns = self._get_sample_columns(heatmap_df)
        if not sample_columns:
            self._logger.warning('未能识别热图中的样本列')
            return

        # Z-score 标准化（按行）
        heatmap_values = heatmap_df[sample_columns].astype(float)
        heatmap_values = heatmap_values.sub(heatmap_values.mean(axis=1), axis=0)
        heatmap_values = heatmap_values.div(heatmap_values.std(axis=1).replace(0, 1), axis=0)

        # 筛选同时在 meta 和表达矩阵中存在的样本
        available_samples = [s for s in meta.index if s in heatmap_values.columns]
        if not available_samples:
            self._logger.warning('热图样本与表达数据不匹配，跳过热图绘制')
            return
        heatmap_values = heatmap_values.loc[:, available_samples]

        # 构建分组颜色条 — Low 蓝色 High 红色，显式映射避免 sorted 顺序反转
        group_palette = {"Low": "#3498db", "High": "#e74c3c"}
        all_groups = set(meta.loc[available_samples, 'group'].unique())
        if all_groups - set(group_palette):
            extra_colors = ["#2ecc71", "#f39c12"]
            for i, g in enumerate(sorted(all_groups - set(group_palette))):
                group_palette[g] = extra_colors[i % len(extra_colors)]
        col_colors = pd.Series(
            meta.loc[available_samples, 'group'].map(group_palette).values,
            index=available_samples,
            name="",
        )

        n_genes = heatmap_values.shape[0]
        fig_height = max(6, n_genes * 0.35)

        g = sns.clustermap(
            heatmap_values,
            cmap="RdBu_r",
            center=0,
            row_cluster=True,
            col_cluster=True,
            col_colors=col_colors,
            yticklabels=True,
            xticklabels=True,
            figsize=(12, fig_height),
            dendrogram_ratio=(0.10, 0.08),
        )

        # clustermap 内部 tight_layout 会覆盖 cbar_pos → 先 subplots_adjust
        # 留出空间，再 set_position 手动设 cbar 位置
        g.figure.subplots_adjust(top=0.90)
        g.ax_cbar.set_position([0.03, 0.80, 0.015, 0.15])

        # 在 col_colors 条上标注组别名 (Low / High)
        col_order = g.dendrogram_col.reordered_ind
        ordered_series = meta.loc[available_samples, 'group'].iloc[col_order]
        ordered_series = ordered_series.reset_index(drop=True)
        group_changes = ordered_series.ne(ordered_series.shift())
        for i in range(len(ordered_series)):
            if group_changes.iloc[i]:
                group_name = str(ordered_series.iloc[i])
                g.ax_col_colors.text(
                    i + 0.5, 0.5, group_name,
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold",
                )

        g.ax_col_colors.set_visible(True)
        g.ax_heatmap.set_xlabel("Sample")
        g.ax_heatmap.set_ylabel("")

        # 标题置顶，避免与列树状图重叠
        g.figure.suptitle(
            f"{self.cfg.tar_gene} — Top {min(20, len(top_genes))} DEGs",
            fontsize=10, fontweight="bold", y=0.97,
        )

        fig_name = self._build_fig_name("heatmap.png")
        fig_path = self._resolve_save_path(self._build_fig_path(fig_name))
        g.savefig(fig_path, dpi=getattr(self.cfg, "figure_dpi", 300), bbox_inches="tight")
        plt.close(g.figure)
        self._logger.info("热图绘制完成（聚类 + 分组注释）")

    def _rename_expr_columns_by_meta_order(self, expr_df: pd.DataFrame, full_meta: pd.DataFrame) -> pd.DataFrame:
        """尝试通过元数据顺序将表达矩阵的样本列映射为 GSM 样本名"""
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

    def _set_expression_index_by_gene_label(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """如果表达矩阵行索引不是基因名，则使用注释列设置基因索引"""
        label_columns = ['SYMBOL', 'GENE', 'GENENAME', 'ENSEMBL', 'ENTREZID', 'ID_REF', 'TARGETID']
        for col in label_columns:
            if col in expr_df.columns:
                labels = expr_df[col].astype(str).replace({'nan': pd.NA, 'None': pd.NA})
                if labels.notna().any():
                    expr_copy = expr_df.copy()
                    expr_copy.index = labels.where(labels.notna(), expr_df.index.astype(str))
                    return expr_copy
        return expr_df

    def _match_top_genes_rows(self, expr_df: pd.DataFrame, top_genes: list) -> pd.DataFrame:
        """根据候选基因名称或数字索引匹配表达矩阵行"""
        if expr_df.index.isin(top_genes).any():
            return expr_df.loc[expr_df.index.isin(top_genes)]

        if expr_df.index.dtype.kind in 'iu' or expr_df.index.dtype.kind == 'f':
            numeric_genes = [int(g) for g in top_genes if isinstance(g, (int, float)) or (isinstance(g, str) and g.isdigit())]
            if numeric_genes:
                matched = expr_df.loc[expr_df.index.isin(numeric_genes)]
                if not matched.empty:
                    matched = self._set_expression_index_by_gene_label(matched)
                    return matched

        labeled_expr = self._set_expression_index_by_gene_label(expr_df)
        if labeled_expr.index.isin(top_genes).any():
            return labeled_expr.loc[labeled_expr.index.isin(top_genes)]

        if 'ENTREZID' in expr_df.columns:
            matched = expr_df.loc[expr_df['ENTREZID'].astype(str).isin([str(g) for g in top_genes])]
            if not matched.empty:
                matched = self._set_expression_index_by_gene_label(matched)
                return matched

        return expr_df.loc[expr_df.index.isin(top_genes)]

    def _infer_group_labels_from_sample_names(self, sample_names: pd.Index, exp_type: str = 'Fibrosis') -> pd.Series:
        """从表达矩阵列名推断组标签
        
        Args:
            sample_names: 表达矩阵的列索引，通常是样本名
            exp_type: 实验组的标签前缀，默认为'Fibrosis'
        
        Returns:
            pd.Series: 包含组标签的Series，索引与输入的sample_names相同
        
        Assumes:
            - Control组样本名以'N'开头（如N1, N2, ...）
            - 实验组样本名以'D'开头（如D1, D2, ...）
        """
        if self.cfg.exp_type:
            exp_type = self.cfg.exp_type
        labels = []
        for name in sample_names:
            if not isinstance(name, str) or len(name) < 2:
                labels.append(None)
                continue
            prefix = name[1].upper()
            if prefix == 'N':
                labels.append('Control')
            elif prefix == 'D':
                labels.append(exp_type)
            else:
                labels.append(None)
        return pd.Series(labels, index=sample_names)

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

    def _prepare_diff_data(self, exp_type: str = 'Fibrosis', gene: str = None) -> dict:
        """准备差异分析的箱线图数据"""
        if self.cfg.exp_type:
            exp_type = self.cfg.exp_type

        if gene is None:
            gene = parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene)[0]

        meta = self._meta_matrix_pack['meta']
        expr_df = self._get_expr_matrix()

        # 提取基因表达向量
        y = fetch_gene_vector(expr_df, tar_gene=gene)

        # 分组标签
        group_col = 'group' if 'group' in meta.columns else 'group_label'
        x = None
        try:
            x = meta.loc[y.index, group_col]
        except KeyError:
            self._logger.warning("样本名与元数据索引不匹配，尝试按原始元数据顺序映射组标签")
            full_meta = self._meta_matrix_pack.get('meta_full')
            sample_columns = self._get_sample_columns(expr_df)

            if full_meta is not None and len(sample_columns) == len(full_meta.index):
                rename_map = {old: new for old, new in zip(sample_columns, full_meta.index.astype(str))}
                y = y.rename(index=rename_map)
                x = meta.reindex(y.index)[group_col]
            else:
                x = self._infer_group_labels_from_sample_names(y.index)

        # 过滤缺失值
        valid = x.notna() & y.notna()
        x, y = x[valid], y[valid]

        if self.cfg.analysis_mode == "hilo":
            group_a, group_b = "Low", "High"
            title = f'{gene} expression in Low vs High'
        else:
            group_a, group_b = "Control", exp_type
            title = f'{gene} expression in Control vs {exp_type}'

        group_vals = [y[x == group_a], y[x == group_b]]
        _, p_value = ttest_ind(*group_vals, equal_var=False)
        self._logger.info(f"准备差异分析数据完成,组间t检验p值:{p_value:.4e}")

        return {
            'x': x, 'y': y, 'p_value': p_value,
            'title': title
        }

    @staticmethod
    def _format_p_value(p: float) -> str:
        """格式化 p 值：>= 0.05 → ns，否则保留 3 位小数"""
        if p >= 0.05:
            return "ns"
        return f"{p:.3f}"

    @staticmethod
    def _significance_stars(p: float) -> str:
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return ""

    def _save_box_plot(self, x, y, p_value, title, exp_type='Fibrosis',
                       fig_name=None, ylabel="Expression", ax=None) -> None:
        """画箱线图并存储，支持传入外部 Axes 用于子图网格"""
        is_immune = self.cfg.analysis_mode == "immune"
        palette = {}
        if self.cfg.analysis_mode == "hilo":
            palette = {"Low": "#3498db", "High": "#e74c3c"}
        else:
            if self.cfg.exp_type:
                exp_type = self.cfg.exp_type
            if is_immune:
                palette = {"Control": "#7FB3D8", f"{exp_type}": "#E07B7B"}
            else:
                palette = {"Control": "#3498db", f"{exp_type}": "#e74c3c"}

        standalone = ax is None
        if standalone:
            plt.figure(figsize=(6, 6))
            ax = plt.gca()

        sns.boxplot(x=x, y=y, hue=x, palette=palette, legend=False, ax=ax,
                    showfliers=False)
        sns.stripplot(x=x, y=y, color='black', alpha=0.4, jitter=0.3, ax=ax)

        # 标题：immune 模式显示 p_adj 格式化值，其他模式保持原样
        if is_immune:
            ax.set_title(title, fontsize=10)
        else:
            ax.set_title(f"{title}\nP-value: {p_value:.4e}", fontsize=10)

        ax.set_xlabel("Group")
        ax.set_ylabel(ylabel)

        # immune 模式画显著性横线
        if is_immune:
            stars = self._significance_stars(p_value)
            if stars:
                ylim = ax.get_ylim()
                y_top = ylim[1] * 1.12
                bar_y = ylim[1] * 1.04
                ax.plot([0, 0, 1, 1], [bar_y * 0.98, bar_y, bar_y, bar_y * 0.98],
                        'k-', linewidth=0.8, clip_on=False)
                ax.text(0.5, bar_y * 1.015, stars, ha="center", va="bottom",
                        fontsize=9, fontweight="bold")
                ax.set_ylim(top=y_top)

        if standalone:
            if fig_name is None:
                fig_name = self._build_fig_name("boxplot.png")
            self._save_plot(fig_name)


    def _stacked_bar_plot(self, df, groups, title, fig_name,
                           xlabel="Sample", ylabel="Cell Fraction (%)"):
        """通用堆叠柱状图：样本 × 类别 (fraction)，按分组排列"""
        # 按分组 + 样本名排序
        group_order = ["Control", self.cfg.exp_type if self.cfg.exp_type else "Experiment"]
        group_rank = {g: i for i, g in enumerate(group_order)}
        sort_df = pd.DataFrame({
            "sample": groups.index,
            "group_rank": groups.map(lambda g: group_rank.get(g, 99)),
            "name": groups.index.astype(str),
        })
        sort_df = sort_df.sort_values(["group_rank", "name"])
        sorted_samples = [s for s in sort_df["sample"] if s in df.index]
        df_sorted = df.loc[sorted_samples]

        n_types = df_sorted.shape[1]
        n_samples = df_sorted.shape[0]
        colors = plt.cm.tab20(range(min(n_types, 20)))
        if n_types > 20:
            colors = list(colors) + list(plt.cm.tab20b(range(n_types - 20)))

        fig, ax = plt.subplots(figsize=(max(14, n_samples * 0.42), 7))
        bottom = np.zeros(n_samples)
        bars = []
        for i, col in enumerate(df_sorted.columns):
            values = df_sorted[col].values * 100
            bar = ax.bar(range(n_samples), values, bottom=bottom, color=colors[i],
                        edgecolor="white", linewidth=0.3, label=col, width=0.85)
            bars.append(bar)
            bottom += values

        # 分组分隔线
        prev_group = None
        for i, s in enumerate(sorted_samples):
            g = groups.get(s)
            if prev_group is not None and g != prev_group:
                ax.axvline(i - 0.5, color="#333333", linestyle="--", linewidth=0.8, alpha=0.5)
            prev_group = g

        # 分组标注（放在 x 轴下方，避免与柱体重叠）
        group_positions = {}
        for i, s in enumerate(sorted_samples):
            g = groups.get(s, "Unknown")
            group_positions.setdefault(g, []).append(i)
        for g, positions in group_positions.items():
            mid = (positions[0] + positions[-1]) / 2
            ax.text(mid, -0.10, g, ha="center", va="top", fontsize=9,
                    fontweight="bold", transform=ax.get_xaxis_transform(),
                    clip_on=False)

        ax.set_xticks(range(n_samples))
        ax.set_xticklabels([s.replace("GSM", "") for s in sorted_samples],
                          rotation=45 if n_samples > 15 else 0, ha="right", fontsize=7)
        ax.set_xlabel("")
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_ylim(0, 105)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.margins(x=0.01)

        # 图例放图外右侧
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7,
                  ncol=2 if n_types > 15 else 1, frameon=False)

        fig.tight_layout()
        fig_path = self._resolve_save_path(self._build_fig_path(fig_name))
        fig.savefig(fig_path, dpi=getattr(self.cfg, "figure_dpi", 300), bbox_inches="tight")
        plt.close(fig)
        self._logger.info(f"堆叠柱状图已保存至 {fig_path}")

    def _ensure_immune_groups(self, immune_df):
        """确保有分组信息，immune 模式回退到配置匹配或样本名推断"""
        meta = self._meta_matrix_pack.get("meta")
        sample_names = immune_df.index
        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"

        # 1) meta 中有 group 列且索引可匹配
        if meta is not None and "group" in meta.columns:
            groups = meta["group"].reindex(sample_names)
            if groups.notna().any():
                return groups

        # 2) 用 control_label / exp_label 匹配 group_select_col
        if meta is not None and self.cfg.group_select_col:
            col = self.cfg.group_select_col
            if col in meta.columns:
                meta_vals = meta[col].astype(str).str.lower()
                ctrl_lower = [l.lower() for l in self.cfg.control_label]
                exp_lower = [l.lower() for l in self.cfg.exp_label]
                label_map = {}
                for v in meta_vals.unique():
                    if any(c in v for c in ctrl_lower):
                        label_map[v] = "Control"
                    elif any(e in v for e in exp_lower):
                        label_map[v] = exp_type
                if label_map:
                    groups = pd.Series(
                        [label_map.get(str(meta.loc[s, col]).lower(), None) for s in sample_names],
                        index=sample_names
                    )
                    valid = groups.notna()
                    if valid.any():
                        return groups

        # 3) 尝试按顺序映射：表达矩阵列顺序 = meta 行顺序
        if meta is not None and len(sample_names) == len(meta):
            ordered_meta = meta.copy()
            ordered_meta.index = sample_names
            if "group" in ordered_meta.columns:
                return ordered_meta["group"]

        # 4) 样本名首字符推断 (N→Control, D/AH→Fibrosis)
        self._logger.info("分组信息缺失，使用样本名首字符推断分组")
        labels = []
        for name in sample_names:
            s = str(name).upper()
            if s.startswith("N"):
                labels.append("Control")
            elif s.startswith(("D", "AH")):
                labels.append(exp_type)
            else:
                labels.append(None)
        result = pd.Series(labels, index=sample_names)
        return result

    def immune_stacked_bar(self):
        """免疫浸润堆叠柱状图"""
        immune_df = self._gene_immune_table
        if immune_df is None or immune_df.empty:
            self._logger.warning("免疫浸润结果为空，跳过堆叠柱状图")
            return

        self._logger.info("构建免疫浸润堆叠柱状图数据...")
        groups = self._ensure_immune_groups(immune_df)
        valid = groups.notna()
        if not valid.any():
            self._logger.warning("未能确定任何样本的分组，跳过堆叠柱状图")
            return
        immune_df = immune_df.loc[valid[valid].index]
        groups = groups[valid]

        # 检测数据是否接近比例（行和 ≈ 1），偏离时警告用户
        if getattr(self.cfg, "plot_data_warning", True):
            row_sums = immune_df.sum(axis=1)
            median_sum = row_sums.median()
            if median_sum < 0.9 or median_sum > 1.1:
                self._logger.warning(
                    f"数据未经过比例归一化（各样本行和中位数={median_sum:.3f}，"
                    f"偏离 1），堆叠柱状图的 Y 轴 \"Cell Fraction (%)\" 可能"
                    f"无意义。若您的反卷积方法输出富集分数而非比例，这是正常现象。"
                    f"设置 plot_data_warning: false 可关闭此提示。"
                )

        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"
        title = f"Immune Cell Composition — {self.cfg.tar_gene} ({self.cfg.gse_id})"
        fig_name = self._build_fig_name("stacked_bar.png")
        self._stacked_bar_plot(immune_df, groups, title, fig_name)

    def _draw_immune_box_grid(self, immune_df, groups, cell_types, p_values, fig_name):
        """绘制免疫箱线图网格（通用：全量版 / 过滤版复用）"""
        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"
        n = len(cell_types)
        n_cols = getattr(self.cfg, "immune_boxplot_n_cols", 5)
        n_rows = (n + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(n_cols * 3.5, n_rows * 3.2))
        axes_flat = axes.flatten() if isinstance(axes, np.ndarray) else [axes]

        for i, ct in enumerate(cell_types):
            ax = axes_flat[i]
            y = immune_df[ct].dropna()
            x = groups.reindex(y.index)
            valid_xy = x.notna() & y.notna()
            x, y = x[valid_xy], y[valid_xy]

            short_name = ct[:40] + "..." if len(ct) > 43 else ct
            self._save_box_plot(x, y, p_values[ct], short_name,
                                exp_type=exp_type, ylabel="Fraction", ax=ax)

            # 轴标签精简：仅最底行显示 x 标签，仅最左列显示 y 标签
            row, col = divmod(i, n_cols)
            if row < n_rows - 1:
                ax.set_xlabel("")
            if col > 0:
                ax.set_ylabel("")

        # 隐藏多余子图
        for j in range(i + 1, len(axes_flat)):
            axes_flat[j].set_visible(False)

        fig.suptitle(f"Immune Cell Fractions — {self.cfg.tar_gene}",
                     fontsize=13, fontweight="bold", y=1.01)
        fig.tight_layout()
        fig_path = self._resolve_save_path(self._build_fig_path(fig_name))
        fig.savefig(fig_path, dpi=getattr(self.cfg, "figure_dpi", 300), bbox_inches="tight")
        plt.close(fig)
        self._logger.info(f"免疫浸润箱线图已保存至 {fig_path}")

    def immune_box_plots(self):
        """免疫细胞组间差异箱线图（子图网格）"""
        immune_df = self._gene_immune_table
        if immune_df is None or immune_df.empty:
            self._logger.warning("免疫浸润结果为空，跳过箱线图")
            return

        self._logger.info("构建免疫浸润组间差异箱线图数据...")
        groups = self._ensure_immune_groups(immune_df)
        valid = groups.notna()
        if not valid.any():
            self._logger.warning("未能确定任何样本的分组，跳过箱线图")
            return
        immune_df = immune_df.loc[valid[valid].index]
        groups = groups[valid]

        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Experiment"
        cell_types = list(immune_df.columns)
        n = len(cell_types)

        # 1) 计算所有细胞类型的 p 值 + Bonferroni 校正
        self._logger.info(f"对 {n} 种免疫细胞进行 Mann-Whitney U 检验 (Bonferroni 校正)...")
        p_raw = {}
        for ct in cell_types:
            y = immune_df[ct].dropna()
            x = groups.reindex(y.index)
            valid_xy = x.notna() & y.notna()
            x_v, y_v = x[valid_xy], y[valid_xy]
            ctrl_vals = y_v[x_v == "Control"]
            exp_vals = y_v[x_v == exp_type]
            if len(ctrl_vals) >= 3 and len(exp_vals) >= 3:
                _, p_val = mannwhitneyu(ctrl_vals, exp_vals, alternative="two-sided")
            else:
                p_val = 1.0
            p_raw[ct] = p_val

        # Bonferroni 校正
        p_adj = {ct: min(p * n, 1.0) for ct, p in p_raw.items()}

        # 2) 按 p_adj 升序排列
        cell_types_sorted = sorted(cell_types, key=lambda ct: p_adj[ct])

        # 3) 全量版
        self._draw_immune_box_grid(immune_df, groups, cell_types_sorted,
                                   p_adj, self._build_fig_name("boxplot.png"))

        # 4) 过滤版：剔除两边组中位数均低于阈值的细胞类型
        low_thr = getattr(self.cfg, "immune_low_abundance_threshold", 0.001)
        cell_types_filtered = []
        for ct in cell_types_sorted:
            ctrl_med = immune_df.loc[groups == "Control", ct].median()
            exp_med = immune_df.loc[groups == exp_type, ct].median()
            if (pd.notna(ctrl_med) and ctrl_med >= low_thr) or \
               (pd.notna(exp_med) and exp_med >= low_thr):
                cell_types_filtered.append(ct)

        if len(cell_types_filtered) < len(cell_types_sorted):
            self._logger.info(
                f"过滤低丰度细胞: {len(cell_types_sorted)} → {len(cell_types_filtered)}"
            )
            self._draw_immune_box_grid(immune_df, groups, cell_types_filtered,
                                       p_adj, self._build_fig_name("boxplot_filtered.png"))
        else:
            self._logger.info("所有细胞类型均通过过滤阈值，不重复输出过滤版")

    def immune_corr_heatmap(self):
        """多基因 × 免疫细胞丰度的 Spearman 相关性热图"""
        immune_df = self._gene_immune_table
        if immune_df is None or immune_df.empty:
            self._logger.warning("免疫浸润结果为空，跳过相关性热图")
            return

        self._logger.info("构建免疫浸润相关性热图数据...")
        try:
            expr_df = self._get_expr_matrix()
        except KeyError:
            self._logger.warning("未找到表达矩阵，跳过相关性热图")
            return

        genes = parse_tar_genes(self.cfg.tar_gene, self.cfg.multi_gene)
        cell_types = list(immune_df.columns)
        skipped = []
        rows_r = {}  # gene → {cell_type: R}
        rows_p = {}  # gene → {cell_type: P_value}

        for gene in genes:
            gene_vec = fetch_gene_vector(expr_df, tar_gene=gene)
            if gene_vec.empty:
                self._logger.warning(f"未找到 {gene} 的表达数据，跳过")
                skipped.append(gene)
                continue

            common = immune_df.index.intersection(gene_vec.index)
            if len(common) < 5:
                self._logger.warning(f"{gene} 交集样本不足 ({len(common)}), 跳过")
                skipped.append(gene)
                continue

            im_aligned = immune_df.loc[common]
            ge_aligned = gene_vec.loc[common].astype(float)

            r_dict, p_dict = {}, {}
            for ct in cell_types:
                y = im_aligned[ct].astype(float).dropna()
                idx = y.index.intersection(ge_aligned.dropna().index)
                if len(idx) < 5:
                    r_dict[ct] = np.nan
                    p_dict[ct] = 1.0
                    continue
                r, p = spearmanr(ge_aligned.loc[idx], y.loc[idx])
                r_dict[ct] = r
                p_dict[ct] = p
            rows_r[gene] = r_dict
            rows_p[gene] = p_dict

        if not rows_r:
            self._logger.warning("无有效基因，跳过热图")
            return

        # 组装矩阵: genes × cell_types, 按均值 |R| 降序排列列
        corr_df = pd.DataFrame(rows_r).T  # genes × cell_types
        col_order = corr_df.abs().mean().sort_values(ascending=False).index
        corr_df = corr_df[col_order]
        p_df = pd.DataFrame(rows_p).T[col_order]

        n_genes = len(corr_df)
        n_cells = len(corr_df.columns)

        # 标注文本
        annot = corr_df.round(3).astype(str)
        for gene in corr_df.index:
            for ct in corr_df.columns:
                p = p_df.loc[gene, ct]
                if pd.isna(p) or p >= 0.05:
                    continue
                if p < 0.001:
                    annot.loc[gene, ct] += "\n***"
                elif p < 0.01:
                    annot.loc[gene, ct] += "\n**"
                else:
                    annot.loc[gene, ct] += "\n*"

        title_gene = genes[0] if n_genes == 1 else f"{n_genes} genes"
        fig, ax = plt.subplots(figsize=(max(10, n_cells * 0.5),
                                        max(2.8, n_genes * 0.5)))
        sns.heatmap(corr_df, annot=annot.values if n_genes > 0 else annot.values,
                    fmt="", annot_kws={"fontsize": 6.5, "fontweight": "bold"},
                    cmap="RdBu_r", center=0, vmin=-1, vmax=1,
                    linewidths=0.8, ax=ax,
                    cbar_kws={"label": "Spearman R", "shrink": 0.8})
        ax.set_title(f"{title_gene} vs Immune Cell Fractions — Spearman Correlation",
                     fontsize=11, fontweight="bold")
        ax.set_ylabel("")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=45, labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
        fig.tight_layout()
        fig_name = self._build_fig_name("heatmap.png")
        fig_path = self._resolve_save_path(self._build_fig_path(fig_name))
        fig.savefig(fig_path, dpi=getattr(self.cfg, "figure_dpi", 300), bbox_inches="tight")
        plt.close(fig)
        self._logger.info(f"免疫相关性热图已保存至 {fig_path}"
                          + (f" (跳过基因: {skipped})" if skipped else ""))

        # 另存 CSV
        csv_dir = os.path.join(RESULT_DIR, "csv")
        os.makedirs(csv_dir, exist_ok=True)
        csv_token = genes[0] if n_genes == 1 else f"multi_{n_genes}genes"
        csv_name = f"{self.cfg.gse_id}_{csv_token}_immune_correlation.csv"
        csv_path = self._resolve_save_path(os.path.join(csv_dir, csv_name))
        out_df = corr_df.copy()
        out_df.index.name = "Gene"
        out_df.to_csv(csv_path)
        self._logger.info(f"免疫相关性数据已保存至 {csv_path}")

    def _truncate_label(self, label: str, max_len: int = 55) -> str:
        if len(label) <= max_len:
            return label
        return label[:max_len - 3].rstrip() + "..."

    def _extract_gene_count(self, overlap: str) -> int:
        try:
            return int(str(overlap).split("/")[0])
        except (ValueError, IndexError, AttributeError):
            return 1

    def _prepare_enrich_data(self, gs_df: pd.DataFrame, top_term: int = 15, max_label_len: int = None) -> pd.DataFrame:
        df = gs_df.copy()
        if max_label_len is not None:
            df["Term"] = df["Term"].astype(str).apply(lambda s: self._truncate_label(s, max_len=max_label_len))
        df["Gene_Count"] = df["Overlap"].apply(self._extract_gene_count)

        # Combined Score 更具区分度，-log10(adjusted P) 作为备选
        if "Combined Score" in df.columns and df["Combined Score"].nunique() > 1:
            df["Score"] = df["Combined Score"]
            xlabel = "Combined Score"
        else:
            df["Score"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
            xlabel = "-log10(Adjusted P-value)"

        df = df.sort_values("Score", ascending=True)
        df = df.tail(top_term)
        return df

    def _enrich_save(self, fig, fig_path: str):
        fig_path = self._resolve_save_path(fig_path)
        fig.savefig(fig_path, dpi=getattr(self.cfg, "figure_dpi", 300), bbox_inches="tight", facecolor="white")
        plt.close(fig)

    def enrich_bar_plotter(self) -> None:
        """富集分析柱状图 — 自定义 matplotlib 实现"""
        if self._gene_enrich_table is None:
            self._load_data()

        if self._gene_enrich_table is None or self._gene_enrich_table.empty:
            self._logger.error("富集分析结果为空，无法绘制柱状图。")
            return

        gene_sets = self._gene_enrich_table["Gene_set"].unique()
        for gs in gene_sets:
            if gs.startswith("GO"):
                continue  # GO 由 enrich_go_combined_bar 统一处理
            gs_df = self._gene_enrich_table[self._gene_enrich_table["Gene_set"] == gs]
            if gs_df.empty:
                continue
            self._logger.info(f"正在为基因集 {gs} 绘制柱状图...")
            try:
                top_n = getattr(self.cfg, "enrich_plot_top_terms", 15)
                df = self._prepare_enrich_data(gs_df, top_term=top_n, max_label_len=55)
                xlabel = "Combined Score" if "Combined Score" in gs_df.columns and gs_df["Combined Score"].nunique() > 1 else "-log10(Adjusted P-value)"

                values = df["Score"].values
                norm = plt.Normalize(values.min(), values.max())
                colors = plt.cm.viridis_r(norm(values))

                fig, ax = plt.subplots(figsize=(10, 0.45 * len(df) + 1.5))
                ax.barh(range(len(df)), values, color=colors, edgecolor="white", linewidth=0.5)
                ax.set_yticks(range(len(df)))
                ax.set_yticklabels(df["Term"].values, fontsize=9)
                ax.set_xlabel(xlabel, fontsize=11)
                ax.set_title(f"{self._gene_token()} — {gs}", fontsize=12, fontweight="bold")
                ax.invert_yaxis()

                if self.cfg.p_threshold < 1:
                    ref = -np.log10(self.cfg.p_threshold) if "Adjusted" in xlabel else None
                    if ref is not None:
                        ax.axvline(ref, linestyle="--", color="#333333", alpha=0.5, linewidth=0.8)

                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                fig.tight_layout()

                fig_name = self._build_fig_name(f"{gs}_barplot.png")
                fig_path = self._build_fig_path(fig_name)
                self._enrich_save(fig, fig_path)
            except Exception as e:
                self._logger.error(f"{gs} 柱状图绘制失败: {e}")

    def enrich_dot_plotter(self) -> None:
        """富集分析气泡图 — 自定义 matplotlib 实现"""
        if self._gene_enrich_table is None:
            self._load_data()

        if self._gene_enrich_table is None or self._gene_enrich_table.empty:
            self._logger.error("富集分析结果为空，无法绘制气泡图。")
            return

        gene_sets = self._gene_enrich_table["Gene_set"].unique()
        for gs in gene_sets:
            if gs.startswith("GO"):
                continue  # GO 由 enrich_go_combined_bar 统一处理
            gs_df = self._gene_enrich_table[self._gene_enrich_table["Gene_set"] == gs]
            if gs_df.empty:
                continue
            self._logger.info(f"正在为基因集 {gs} 绘制气泡图...")
            try:
                top_n = getattr(self.cfg, "enrich_plot_top_terms", 15)
                df = self._prepare_enrich_data(gs_df, top_term=top_n, max_label_len=55)

                gene_counts = df["Gene_Count"].values
                min_size, max_size = 40, 350
                if gene_counts.max() > gene_counts.min():
                    sizes = min_size + (max_size - min_size) * (
                        (gene_counts - gene_counts.min()) / (gene_counts.max() - gene_counts.min())
                    )
                else:
                    sizes = np.full_like(gene_counts, (min_size + max_size) // 2, dtype=float)

                values = df["Score"].values
                norm = plt.Normalize(values.min(), values.max())

                fig, ax = plt.subplots(figsize=(9, 0.45 * len(df) + 1.5))
                scatter = ax.scatter(
                    values,
                    range(len(df)),
                    s=sizes,
                    c=values,
                    cmap="viridis_r",
                    norm=norm,
                    alpha=0.85,
                    edgecolors="white",
                    linewidth=0.5,
                    zorder=3,
                )
                ax.set_yticks(range(len(df)))
                ax.set_yticklabels(df["Term"].values, fontsize=9)
                ax.invert_yaxis()

                xlabel = "Combined Score" if "Combined Score" in gs_df.columns and gs_df["Combined Score"].nunique() > 1 else "-log10(Adjusted P-value)"
                ax.set_xlabel(xlabel, fontsize=11)
                ax.set_title(f"{self._gene_token()} — {gs}", fontsize=12, fontweight="bold")

                if self.cfg.p_threshold < 1 and "Adjusted" in xlabel:
                    ref = -np.log10(self.cfg.p_threshold)
                    ax.axvline(ref, linestyle="--", color="#333333", alpha=0.5, linewidth=0.8)

                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)

                cbar = fig.colorbar(scatter, ax=ax, fraction=0.03, pad=0.04)
                cbar.ax.set_title('Combined Score' if "Combined Score" in xlabel else 'Score', fontsize=9)

                # 基因数量图例
                legend_counts = [int(np.percentile(gene_counts, p)) for p in [0, 50, 100]]
                legend_counts = sorted(set(max(1, c) for c in legend_counts))
                legend_elements = []
                for cnt in legend_counts:
                    s = min_size + (max_size - min_size) * ((cnt - gene_counts.min()) / (gene_counts.max() - gene_counts.min())) if gene_counts.max() > gene_counts.min() else (min_size + max_size) // 2
                    legend_elements.append(
                        plt.scatter([], [], s=s, c="#555555", alpha=0.7, edgecolors="white", linewidth=0.5, label=str(cnt))
                    )
                leg = ax.legend(
                    handles=legend_elements, 
                    title="Gene Count", 
                    loc="upper right",
                    bbox_to_anchor=(1.2, 1),
                    labelspacing=1.5,
                    frameon=False, 
                    fancybox=True, 
                    fontsize=8, 
                    title_fontsize=9
                    )
                leg.get_frame().set_alpha(0.85)

                fig.tight_layout()

                fig_name = self._build_fig_name(f"{gs}_dotplot.png")
                fig_path = self._build_fig_path(fig_name)
                self._enrich_save(fig, fig_path)
            except Exception as e:
                self._logger.error(f"{gs} 气泡图绘制失败: {e}")

    def enrich_go_combined_bar(self) -> None:
        top_term = getattr(self.cfg, "enrich_plot_top_terms", 15)
        """GO 三合一子图 — BP/MF/CC 各画条形图，垂直拼接，右侧图例"""
        if self._gene_enrich_table is None:
            self._load_data()

        enrich = self._gene_enrich_table
        if enrich is None or enrich.empty:
            self._logger.error("富集分析结果为空，无法绘制 GO 合集图。")
            return

        go_sets = [gs for gs in enrich["Gene_set"].unique() if gs.startswith("GO")]
        if not go_sets:
            self._logger.info("未找到 GO 富集结果，跳过 GO 合集图。")
            return

        go_order = {"Biological": 0, "Cellular": 1, "Molecular": 2}
        go_sets.sort(key=lambda gs: next(
            (go_order[k] for k in go_order if k in gs), 99
        ))
        go_abbr = {"Biological": "BP", "Cellular": "CC", "Molecular": "MF"}

        n = len(go_sets)
        fig, axes = plt.subplots(n, 1, figsize=(10, 2.8 * n), sharex=False)
        if n == 1:
            axes = [axes]

        score_is_combined = (
            "Combined Score" in enrich.columns
            and enrich["Combined Score"].nunique() > 1
        )
        xlabel = "Combined Score" if score_is_combined else "-log10(Adjusted P-value)"

        # 全局归一化颜色映射
        all_scores = []
        for gs in go_sets:
            gs_df = enrich[enrich["Gene_set"] == gs]
            df = self._prepare_enrich_data(gs_df, top_term=top_term, max_label_len=65) 
            all_scores.extend(df["Score"].values)
        global_norm = plt.Normalize(min(all_scores), max(all_scores))

        for ax, gs in zip(axes, go_sets):
            gs_df = enrich[enrich["Gene_set"] == gs]
            df = self._prepare_enrich_data(gs_df, top_term=top_term, max_label_len=65) 

            # 清洗 GO ID 后缀 (如 "(GO:0001234)")
            labels = df["Term"].astype(str).str.replace(r"\s*\(GO[^)]*\)?", "", regex=True)
            labels = labels.str.strip()

            values = df["Score"].values
            colors = plt.cm.viridis_r(global_norm(values))

            ax.barh(range(len(df)), values, color=colors, edgecolor="white", linewidth=0.5)
            ax.set_yticks(range(len(df)))
            ax.set_yticklabels(labels.values, fontsize=9)
            ax.invert_yaxis()

            abbr = next((go_abbr[k] for k in go_abbr if k in gs), gs[:4])
            ax.text(1.02, 0.5, abbr, transform=ax.transAxes, fontsize=12,
                    fontweight="bold", va="center", ha="left")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

        # 底部子图显示统一 xlabel
        axes[-1].set_xlabel(xlabel, fontsize=10)

        # 右侧图例
        cbar_ax = fig.add_axes([1.0, 0.3, 0.015, 0.4])
        sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=global_norm)
        cbar = fig.colorbar(sm, cax=cbar_ax, orientation="vertical", fraction=0.03, pad=0.04)
        cbar.ax.set_title('Score', fontsize=9)

        fig.suptitle(
            f"{self._gene_token()} — GO Enrichment Combined",
            fontsize=14, fontweight="bold", y=1.01,
        )
        fig.subplots_adjust(left=0.60, right=0.90, top=0.93, bottom=0.05)

        fig_name = self._build_fig_name("go_combined_bar.png")
        fig_path = self._build_fig_path(fig_name)
        self._enrich_save(fig, fig_path)
        self._logger.info("GO 合集条形图绘制完成！")

    def _save_enrich_plotting_csv(self) -> None:
        """将画图用的 prepared enrichment data 合并保存为 CSV。"""
        enrich = self._gene_enrich_table
        if enrich is None or enrich.empty:
            return

        gene_sets = enrich["Gene_set"].unique()
        frames = []
        for gs in gene_sets:
            gs_df = enrich[enrich["Gene_set"] == gs]
            if gs_df.empty:
                continue
            prepared = self._prepare_enrich_data(gs_df)
            prepared["Gene_set"] = gs
            frames.append(prepared)

        if not frames:
            return

        plot_df = pd.concat(frames, ignore_index=True)
        csv_dir = os.path.join(RESULT_DIR, "csv")
        os.makedirs(csv_dir, exist_ok=True)
        fname = f"{self.cfg.gse_id}_{self._gene_token()}_enrich_plotting_data.csv"
        out_path = self._resolve_save_path(os.path.join(csv_dir, fname))
        plot_df.to_csv(out_path, index=False)
        self._logger.info(f"画图数据已保存至 {out_path}")


class DataPlotter(FigurePlotter):
    """直接从内存中调用数据画图"""
    def __init__(self, cfg: Config, data: DataHandler):
        """初始化

        Args:
            cfg: 基础配置
            data: 数据传递类，包括相关性分析/差异分析DataFrame和筛选后的原始基因DataFrame数据
        """
        super().__init__(cfg)
        self.gene_corr_table: Optional[pd.DataFrame] = data.gene_corr_table
        self.gene_diff_table: Optional[pd.DataFrame] = data.gene_diff_table
        self.gene_enrich_table: Optional[pd.DataFrame] = data.gene_enrich_table
        self.gene_immune_table: Optional[pd.DataFrame] = data.gene_immune_table
        self.meta_matrix_pack: dict = data.meta_matrix_pack

    def fig_plotter(self):
        """筛选所需目标并画图"""
        # 读取索引和数据仓库
        self._logger.info("读取索引和数据中...")
        if not self._gene_immune_table and not self._gene_corr_table and not self._meta_matrix_pack:
            self._load_data()
        self._logger.info("索引和数据读取成功！")
        self.figplotter()

    def _load_data(self):
        if self.cfg.analysis_mode == "corr":
            if self.gene_corr_table is not None and not self.gene_corr_table.empty:
                self._gene_corr_table = self.gene_corr_table
        elif self.cfg.analysis_mode == "diff":
            if self.gene_diff_table is not None and not self.gene_diff_table.empty:
                self._gene_diff_table = self.gene_diff_table
        elif self.cfg.analysis_mode == "enrich":
            if self.gene_enrich_table is not None and not self.gene_enrich_table.empty:
                self._gene_enrich_table = self.gene_enrich_table
        elif self.cfg.analysis_mode == "immune":
            if self.gene_immune_table is not None and not self.gene_immune_table.empty:
                self._gene_immune_table = self.gene_immune_table
        else:
            self._logger.error(f"未知分析模式：{self.cfg.analysis_mode}，无法加载数据")
            raise ValueError(f"未知分析模式：{self.cfg.analysis_mode}")

        if self.meta_matrix_pack:
            self._meta_matrix_pack = self.meta_matrix_pack


class FilePlotter(FigurePlotter):
    """通过读取文件数据画图"""
    def __init__(self, cfg: Config):
        super().__init__(cfg)

    def fig_plotter(self):
        """筛选所需目标并画图"""
        self._logger.info("从pkl中读取索引和数据中...")
        need_load = (
            self.cfg.analysis_mode in ("corr", "diff", "hilo")
            and (not self._gene_corr_table or not self._meta_matrix_pack)
        )
        if need_load or (
            self.cfg.analysis_mode in ("enrich", "immune")
            and (not self._gene_enrich_table if self.cfg.analysis_mode == "enrich"
                 else not self._gene_immune_table)
        ):
            self._load_data()
        self._logger.info("索引和数据读取成功！")
        self.figplotter()

    def _load_data(self):
        """从文件中加载数据"""
        # 读取配置
        data_dir, gse_id = os.path.join(self.cfg.data_dir, self.cfg.gse_id), self.cfg.gse_id

        # 执行加载
        summary_name = {
            "corr": f"{gse_id}_correlation_summary.pkl",
            "diff": f"{gse_id}_differential_summary.pkl",
            "hilo": f"{gse_id}_highlow_summary.pkl",
            "enrich": f"{gse_id}_enrichment_summary.pkl",
            "immune": f"{gse_id}_immune_summary.pkl",
        }.get(self.cfg.analysis_mode)

        if summary_name is None:
            self._logger.error(f"未知分析模式：{self.cfg.analysis_mode}，无法加载数据")
            raise ValueError(f"未知分析模式：{self.cfg.analysis_mode}")

        summary_path = os.path.join(data_dir, "pkl", summary_name)
        data_pack_path = DataPacker.resolve_pack_path(data_dir, gse_id, self.cfg.analysis_mode)

        if self.cfg.analysis_mode == "corr":
            self._gene_corr_table = pd.read_pickle(summary_path)
        elif self.cfg.analysis_mode == "enrich":
            self._gene_enrich_table = pd.read_pickle(summary_path)
            return
        elif self.cfg.analysis_mode == "immune":
            self._gene_immune_table = pd.read_pickle(summary_path)
        else:
            self._gene_diff_table = pd.read_pickle(summary_path)

        self._meta_matrix_pack = pd.read_pickle(data_pack_path)


if __name__ == "__main__":
    test_gse_id = "GSE300437"
    test_tar_gene = "Polb"
    test_cfg = Config(tar_gene=test_tar_gene, gse_id=test_gse_id)
    test_plotter = FilePlotter(test_cfg)
    if test_plotter.plotter:
        print("Done!")
