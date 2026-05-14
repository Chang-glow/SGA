import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import seaborn as sns
from scipy import stats as sp_stats

from modules.calculater import fetch_gene_vector, prepare_expr_matrix
from utils import RESULT_DIR, FIGURE_DIR


class WgcnaStrategy:
    """WGCNA 策略：数据准备 + PyWGCNA 完整分析管道"""

    def __init__(self, analyzer):
        self.analyzer = analyzer
        self._logger = analyzer._logger
        self.cfg = analyzer.cfg

    def calculate(self) -> pd.DataFrame:
        pack = self.analyzer.roaming_data()

        expr_df = None
        meta = None
        meta_full = None
        for key in pack:
            if key == "meta":
                meta = pack[key]
                continue
            if key == "meta_full":
                meta_full = pack[key]
                continue
            if key == "group_info":
                continue
            item = pack[key]
            if isinstance(item, dict) and "matrix_aligned" in item:
                expr_df = item["matrix_aligned"]
            elif isinstance(item, pd.DataFrame):
                expr_df = item
            if expr_df is not None:
                break

        if expr_df is None:
            raise ValueError("未在数据包中找到表达矩阵，无法进行WGCNA数据准备")

        if meta_full is not None:
            expr_df = self.analyzer._rename_expr_columns_by_meta_order(expr_df, meta_full)

        # 筛选前提取目标基因表达向量
        tar_expr = fetch_gene_vector(expr_df, self.cfg.tar_gene)
        if tar_expr.empty:
            self._logger.warning(f"未找到目标基因 {self.cfg.tar_gene} 的表达数据")
            tar_expr = None

        expr_df = prepare_expr_matrix(expr_df)
        self._logger.info(f"表达矩阵: {expr_df.shape[0]} 基因 × {expr_df.shape[1]} 样本")

        n_genes = getattr(self.cfg, "wgcna_top_n_genes", 3000)
        filtered_df = self._filter_top_genes(expr_df, tar_expr, n_genes)

        # 转置: 样本(行) × 基因(列)
        datExpr = filtered_df.T
        datExpr.index.name = "Sample"
        datExpr.columns = [f"X{c}" if c[0].isdigit() else c for c in datExpr.columns]

        datTraits = self._build_traits(datExpr, meta, tar_expr)

        # 按生物分组重排序样本：Control 在前，Treatment 在后
        if "Group" in datTraits.columns and "Control" in datTraits["Group"].values:
            control_mask = datTraits["Group"] == "Control"
            new_order = (list(datTraits.index[control_mask])
                         + list(datTraits.index[~control_mask]))
            datTraits = datTraits.loc[new_order]
            datExpr = datExpr.loc[new_order]
            self._logger.info(
                f"样本已按分组重排序: {control_mask.sum()} Control + "
                f"{(~control_mask).sum()} Treatment"
            )

        self._save_csv(datExpr, "datExpr")
        self._save_csv(datTraits, "datTraits")

        # ── PyWGCNA 完整分析管道 ──
        from PyWGCNA import WGCNA

        expr_csv_path = os.path.join(RESULT_DIR, "csv",
                                     f"{self.cfg.gse_id}_{self.cfg.tar_gene}_datExpr.csv")
        wgcna_output = os.path.join(FIGURE_DIR, "wgcna") + os.sep

        self._logger.info("启动 PyWGCNA 分析管道...")
        wgcna = WGCNA(
            name=f"{self.cfg.gse_id}_{self.cfg.tar_gene}",
            species=self.cfg.organism,
            geneExpPath=expr_csv_path,
            sampleInfo=datTraits,
            save=True,
            outputPath=wgcna_output,
        )
        # 必须设置 gene_name 列，否则 analyseWGCNA 会跳过所有 GO 富集
        wgcna.datExpr.var["gene_name"] = wgcna.datExpr.var.index

        self._logger.info("运行 WGCNA 模块检测（预处理 → 软阈值 → TOM → 聚类）...")
        wgcna.runWGCNA()

        modules = wgcna.datExpr.var["moduleColors"].unique().tolist()
        self._logger.info(f"检测到 {len(modules)} 个模块: {modules}")

        # 设置 metadata 颜色映射（plotModuleEigenGene 依赖此字典）
        for col in datTraits.columns:
            unique_vals = datTraits[col].dropna().unique()
            if len(unique_vals) <= 10:
                palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                           "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                           "#bcbd22", "#17becf"]
                cmap = {val: palette[i % len(palette)] for i, val in enumerate(sorted(unique_vals))}
            else:
                cmap = ScalarMappable(cmap="viridis", norm=mpl.colors.Normalize())
            wgcna.setMetadataColor(col, cmap)

        self._logger.info("分析模块-性状关系、计算 kME、GO 富集...")
        metadata_cols = datTraits.columns.tolist()

        # 1) 模块-性状关系热图（带显著性星号标注）
        self._module_trait_heatmap_with_stars(
            wgcna, metadata_cols,
            f"{self.cfg.gse_id}_{self.cfg.tar_gene}_module-traitRelationships",
        )

        # 2) 计算 signed kME（模块 membership）
        wgcna.CalculateSignedKME()

        # 3) 各模块 eigengene 图 & GO 富集
        # 猴子补丁：gseapy DotPlot.scatter 将尺寸图例移到图下方，避免与颜色条重叠
        import gseapy.plot as _gp_plot_module
        _original_dotplot_scatter = _gp_plot_module.DotPlot.scatter

        def _patched_scatter(self, outer_ring: bool = False):
            from gseapy.plot import UnitData

            df = self.data.assign(
                area=self.data["Hits_ratio"]
                * (self.scale * plt.rcParams["lines.markersize"]) ** 2
            )
            colmap = df[self.colname].astype(int)
            vmin = np.percentile(colmap.min(), 2)
            vmax = np.percentile(colmap.max(), 98)
            ax = self.get_ax()
            x, xlabel = self.set_x()
            y = self.y
            if all(df[x].map(self.isfloat)):
                df = df.sort_values(by=x)
            xunits = UnitData(self.get_x_order()) if self.x_order else None
            yunits = UnitData(self.get_y_order()) if self.y_order else None

            if outer_ring:
                smax = df["area"].max()
                ax.scatter(x=x, y=y, s=smax * 1.6, edgecolors="none",
                           c="black", data=df, marker=self.marker,
                           xunits=xunits, yunits=yunits, zorder=0)
                ax.scatter(x=x, y=y, s=smax * 1.3, edgecolors="none",
                           c="white", data=df, marker=self.marker,
                           xunits=xunits, yunits=yunits, zorder=1)

            sc = ax.scatter(x=x, y=y, data=df, s="area", edgecolors="none",
                            c=self.colname, cmap=self.cmap, vmin=vmin, vmax=vmax,
                            marker=self.marker, xunits=xunits, yunits=yunits, zorder=2)
            ax.set_xlabel(xlabel, fontsize=14, fontweight="bold")
            ax.xaxis.set_tick_params(labelsize=14)
            ax.yaxis.set_tick_params(labelsize=16)
            ax.set_axisbelow(True)
            ax.grid(axis="y", zorder=-1)
            ax.margins(x=0.25)

            # 尺寸图例放在图下方，避免与右侧颜色条重叠
            handles, labels = sc.legend_elements(
                prop="sizes", num=3, color="gray",
                func=lambda s: 100 * s
                / (plt.rcParams["lines.markersize"] * self.scale) ** 2,
            )
            ax.legend(
                handles, labels,
                title="% Genes\nin set",
                bbox_to_anchor=(0.5, -0.18),
                loc="upper center",
                frameon=False,
                labelspacing=2,
                ncol=len(labels),
            )
            ax.set_title(self.title, fontsize=20, fontweight="bold")

            cbar = self.fig.colorbar(
                sc, shrink=0.25, aspect=10,
                anchor=(0.0, 0.2), location="right",
            )
            cbar.ax.yaxis.set_tick_params(color="white", direction="in",
                                          left=True, right=True)
            cbar.ax.set_title(self.cbar_title, loc="left", fontweight="bold")
            cbar.ax.title.set_position((0, 1.05))
            for key, spine in cbar.ax.spines.items():
                spine.set_visible(False)
            return ax

        _gp_plot_module.DotPlot.scatter = _patched_scatter

        for module in modules:
            try:
                wgcna.plotModuleEigenGene(module, metadata_cols, show=False)
            except Exception as e:
                self._logger.warning(f"模块 {module} eigengene 绘图失败: {e}")
            try:
                wgcna.functional_enrichment_analysis(
                    type="GO", moduleName=module,
                    file_name=f"{self.cfg.gse_id}_{self.cfg.tar_gene}_{module}_GO",
                )
            except Exception as e:
                self._logger.warning(f"模块 {module} GO 富集失败: {e}")

        # 恢复原始 DotPlot.scatter 方法
        _gp_plot_module.DotPlot.scatter = _original_dotplot_scatter

        # 提取基因模块归属 + kME 值
        var_df = wgcna.datExpr.var[["moduleColors", "moduleLabels"]].copy()
        var_df.index.name = "Gene"
        kme_df = wgcna.signedKME.copy() if wgcna.signedKME is not None else None
        if kme_df is not None:
            var_df = var_df.join(kme_df)
        summary_df = var_df.reset_index()
        self._save_summary_csv(summary_df, "wgcna_module_genes")

        # 提取各模块 top 10 hub genes
        hub_parts = []
        for module in modules:
            hub = wgcna.top_n_hub_genes(module, n=10)
            if hub is not None and not hub.empty:
                hub.insert(0, "Module", module)
                hub_parts.append(hub)
        if hub_parts:
            hub_df = pd.concat(hub_parts)
            self._save_summary_csv(hub_df, "wgcna_hub_genes")

        wgcna.saveWGCNA()
        self._logger.info("PyWGCNA 分析管道完成")
        self._logger.info(
            "eigengene 热图中样本排列顺序: Control 组在前, Treatment 组在后。"
            "详细样本 ID 对应关系请查看 res/csv/ 目录下的 datExpr.csv"
        )

        return summary_df

    def _filter_top_genes(self, df, tar_expr, n_genes=3000):
        """筛选 top N 基因：优先用 diff/hilo 缓存的 padj，回退按方差"""
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        gse_id = self.cfg.gse_id

        for mode_key in ("differential", "highlow"):
            cache_path = os.path.join(data_dir, "pkl", f"{gse_id}_{mode_key}_summary.pkl")
            if os.path.exists(cache_path) and not self.cfg.debug:
                try:
                    cached = pd.read_pickle(cache_path)
                    if "padj" in cached.columns and "Gene" in cached.columns:
                        top_genes = (
                            cached.dropna(subset=["padj"])
                            .sort_values("padj")["Gene"]
                            .dropna().unique()
                        )
                        top_genes = [g for g in top_genes if g in df.index]
                        selected = top_genes[:n_genes]
                        self._logger.info(
                            f"使用缓存 {mode_key} 结果的 padj 筛选了 {len(selected)} 个基因"
                        )
                        result = df.loc[selected]
                        return self._ensure_tar_gene(result, df, tar_expr)
                except Exception as e:
                    self._logger.debug(f"读取缓存 {cache_path} 失败: {e}")

        gene_var = df.var(axis=1, numeric_only=True)
        top_genes = gene_var.sort_values(ascending=False).head(n_genes).index
        self._logger.info(f"使用方差筛选了 {len(top_genes)} 个基因")
        result = df.loc[top_genes]
        return self._ensure_tar_gene(result, df, tar_expr)

    def _ensure_tar_gene(self, filtered_df, full_df, tar_expr):
        """确保目标基因在筛选结果中"""
        if tar_expr is None:
            return filtered_df
        tar_gene = self.cfg.tar_gene.upper()
        index_upper = [str(i).upper() for i in filtered_df.index]
        if tar_gene not in index_upper:
            for idx in full_df.index:
                if str(idx).upper() == tar_gene:
                    filtered_df = pd.concat([filtered_df, full_df.loc[[idx]]])
                    self._logger.info(f"目标基因 {self.cfg.tar_gene} 被筛选掉，已强制加回")
                    break
        return filtered_df

    def _build_traits(self, datExpr, meta, tar_expr):
        """构建性状 DataFrame：Group(字符串标签) + 目标基因表达量"""
        samples = datExpr.index
        datTraits = pd.DataFrame(index=samples)

        exp_type = self.cfg.exp_type if self.cfg.exp_type else "Fibrosis"

        if meta is not None and "group" in meta.columns:
            group_series = meta["group"].reindex(samples)
            datTraits["Group"] = group_series.apply(
                lambda x: "Control" if x == "Control" else exp_type
            )
            n_control = (datTraits["Group"] == "Control").sum()
            n_exp = (datTraits["Group"] == exp_type).sum()
            self._logger.info(f"Group列: Control={n_control}, {exp_type}={n_exp}")
        else:
            self._logger.warning("未找到分组信息，Group列将全部设为'Unknown'")
            datTraits["Group"] = "Unknown"

        if tar_expr is not None and not tar_expr.empty:
            datTraits[f"{self.cfg.tar_gene}_expression"] = tar_expr.reindex(samples).values
        else:
            self._logger.warning(f"未找到目标基因 {self.cfg.tar_gene} 的表达数据")
            datTraits[f"{self.cfg.tar_gene}_expression"] = float("nan")

        return datTraits

    def _module_trait_heatmap_with_stars(self, wgcna, metaData, file_name):
        """模块-性状关系热图，带显著性星号标注（p<0.05=* p<0.01=** p<0.001=***）"""
        datTraits_inner = wgcna.getDatTraits(metaData)

        moduleTraitCor = pd.DataFrame(index=wgcna.MEs.columns,
                                      columns=datTraits_inner.columns, dtype="float")
        moduleTraitPvalue = pd.DataFrame(index=wgcna.MEs.columns,
                                         columns=datTraits_inner.columns, dtype="float")
        for i in wgcna.MEs.columns:
            for j in datTraits_inner.columns:
                tmp = sp_stats.pearsonr(wgcna.MEs[i], datTraits_inner[j])
                moduleTraitCor.loc[i, j] = tmp[0]
                moduleTraitPvalue.loc[i, j] = tmp[1]

        if set(metaData) == set(wgcna.datExpr.obs.columns.tolist()):
            wgcna.moduleTraitCor = moduleTraitCor
            wgcna.moduleTraitPvalue = moduleTraitPvalue

        figsize = (max(20, int(moduleTraitPvalue.shape[0] * 1.8)),
                   moduleTraitPvalue.shape[1] * 1.8)
        fig, ax = plt.subplots(figsize=figsize, facecolor="white")

        xlabels = []
        for label in wgcna.MEs.columns:
            xlabels.append(
                label[2:].capitalize() + "("
                + str(sum(wgcna.datExpr.var["moduleColors"] == label[2:])) + ")"
            )
        ylabels = datTraits_inner.columns

        tmp_cor = moduleTraitCor.T.round(decimals=2)
        tmp_pvalue = moduleTraitPvalue.T.round(decimals=3)

        def _pv_stars(p):
            if p < 0.001:
                return "***"
            elif p < 0.01:
                return "**"
            elif p < 0.05:
                return "*"
            return ""

        labels = np.asarray([
            "{0:.2f}\n({1:.3f}){2}".format(cor, pvalue, _pv_stars(pvalue))
            for cor, pvalue in zip(tmp_cor.values.flatten(),
                                   tmp_pvalue.values.flatten())
        ]).reshape(moduleTraitCor.T.shape)

        sns.set(font_scale=1.5)
        res = sns.heatmap(moduleTraitCor.T, annot=labels, fmt="",
                          cmap="RdBu_r", vmin=-1, vmax=1, ax=ax,
                          annot_kws={"size": 18, "weight": "bold"},
                          xticklabels=xlabels, yticklabels=ylabels)
        res.set_xticklabels(res.get_xmajorticklabels(),
                            fontsize=20, fontweight="bold", rotation=90)
        res.set_yticklabels(res.get_ymajorticklabels(),
                            fontsize=20, fontweight="bold")
        plt.yticks(rotation=0)
        ax.set_title(f"Module-trait Relationships heatmap for {wgcna.name}",
                     fontsize=30, fontweight="bold")
        ax.set_facecolor("white")
        fig.tight_layout()
        plt.close(fig)
        if wgcna.save:
            fig.savefig(
                f"{wgcna.outputPath}figures/{file_name}.{wgcna.figureType}",
                bbox_inches="tight",
            )

    def _save_summary_csv(self, df, suffix):
        """保存总结 CSV 到 res/csv/（index=False）"""
        csv_dir = os.path.join(RESULT_DIR, "csv")
        os.makedirs(csv_dir, exist_ok=True)
        filename = f"{self.cfg.gse_id}_{self.cfg.tar_gene}_{suffix}.csv"
        filepath = os.path.join(csv_dir, filename)
        df.to_csv(filepath, index=False)
        self._logger.info(f"已保存: {filepath} ({df.shape[0]} 行 × {df.shape[1]} 列)")

    def _save_csv(self, df, suffix):
        """保存 DataFrame 到 res/csv/"""
        csv_dir = os.path.join(RESULT_DIR, "csv")
        os.makedirs(csv_dir, exist_ok=True)
        filename = f"{self.cfg.gse_id}_{self.cfg.tar_gene}_{suffix}.csv"
        filepath = os.path.join(csv_dir, filename)
        df.to_csv(filepath, index=True, index_label="Sample")
        self._logger.info(f"已保存: {filepath} ({df.shape[0]} 行 × {df.shape[1]} 列)")
