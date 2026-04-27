import os
from abc import ABC, abstractmethod
from typing import Optional
from scipy.stats import ttest_ind

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from utils.paths import FIGURE_DIR
from modules.calculater import fetch_gene_vector
from utils import Config, DataHandler, loggers


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
        self._meta_matrix_pack: Optional[dict] = None
        # 内部判断是否画图
        self._plotter = False

    @classmethod
    def create(cls, cfg: Config, data: DataHandler):
        """根据cfg检查使用哪个子类"""
        data_dir = os.path.join(cfg.data_dir, cfg.gse_id)
        data_pack_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_processed_pack.pkl")
        gene_corr_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_correlation_summary.pkl")
        gene_diff_path = os.path.join(data_dir, "pkl", f"{cfg.gse_id}_differential_summary.pkl")

        is_pack = os.path.exists(data_pack_path)
        if cfg.analysis_mode == "corr":
            is_table = os.path.exists(gene_corr_path)
        elif cfg.analysis_mode == "diff":
            is_table = os.path.exists(gene_diff_path)
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
    def fig_plotter(self):
        pass

    @abstractmethod
    def _load_data(self):
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
    
    def figplotter(self):
        """绘图主流程，包含条件筛选和循环画图"""
        if self.cfg.analysis_mode == "corr":
            self._logger.info("分析模式设定为相关性分析,将绘制相关性散点图")
            self.scatter_plotter()
        elif self.cfg.analysis_mode == "diff":
            self._logger.info("分析模式设定为差异分析,将绘制差异箱线图")
            self.box_plotter()
        else:
            self._logger.error(f"未知分析模式：{self.cfg.analysis_mode}，无法绘图")
            raise ValueError(f"未知分析模式：{self.cfg.analysis_mode}")

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
        tar_gene, _ = self.cfg.tar_gene, os.path.join(self.cfg.data_dir, self.cfg.gse_id)

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
        fig_name = f"{matrix_info}_{gene_info}.png"
        self._save_plot(fig_name)

    def _save_plot(self, fig_name: str) -> None:
        # 保存图片
        fig_path = os.path.join(FIGURE_DIR, fig_name)
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close()

    def box_plotter(self) -> None:
        """箱线图pipeline"""
        # 加载数据
        if not self._meta_matrix_pack:
            self._load_data()

        # 准备数据并画图
        self._logger.info("正在准备差异分析数据...")
        data_dict = self._prepare_diff_data()
        self._logger.info("差异分析数据准备完成，正在绘图...")
        self._save_box_plot(
            data_dict['x'], 
            data_dict['y'], 
            data_dict['p_value'], 
            data_dict['title']
            )
        self._logger.info("差异分析箱线图绘制完成！")

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

    def _infer_group_labels_from_sample_names(self, sample_names):
        """从表达矩阵列名推断组标签"""
        labels = []
        for name in sample_names:
            if not isinstance(name, str) or len(name) < 2:
                labels.append(None)
                continue
            prefix = name[1].upper()
            if prefix == 'N':
                labels.append('Control')
            elif prefix == 'D':
                labels.append('Fibrosis')
            else:
                labels.append(None)
        return pd.Series(labels, index=sample_names)

    def _prepare_diff_data(self) -> dict:
        """准备差异分析的箱线图数据"""
        # 获取数据矩阵
        gene = self.cfg.tar_gene
        meta = self._meta_matrix_pack['meta']
        expr_df = self._get_expr_matrix()

        # 提取基因表达向量
        y = fetch_gene_vector(expr_df, tar_gene=gene)

        # 分组标签
        group_col = 'group' if 'group' in meta.columns else 'group_label'
        try:
            x = meta.loc[y.index, group_col]
        except KeyError:
            self._logger.warning("样本名与元数据索引不匹配，尝试按原始元数据顺序映射组标签")
            full_meta = self._meta_matrix_pack.get('meta_full')
            if full_meta is not None and len(y.index) == len(full_meta.index):
                x = pd.Series(full_meta[group_col].values, index=y.index)
            else:
                x = self._infer_group_labels_from_sample_names(y.index)

        # 过滤缺失值
        valid = x.notna() & y.notna()
        x, y = x[valid], y[valid]

        # check
        group_vals = [y[x == "Control"], y[x == "Fibrosis"]]
        _, p_value = ttest_ind(*group_vals, equal_var=False)
        self._logger.info(f"准备差异分析数据完成,组间t检验p值:{p_value:.4e}")

        return {
            'x': x, 'y': y, 'p_value': p_value,
            'title': f'{gene} expression in Control vs Fibrosis'
        }

    def _save_box_plot(self, x, y, p_value, title) -> None:
        """画箱线图并存储"""
        plt.figure(figsize=(6, 6))
        sns.boxplot(x=x, y=y, palette={'Control':'#3498db', 'Fibrosis':'#e74c3c'})
        sns.stripplot(x=x, y=y, color='black', alpha=0.6, jitter=True)
        plt.title(f"{title}\nP-value: {p_value:.4e}", fontsize=10)
        plt.xlabel("Group")
        plt.ylabel("Expression")
        fig_name = f"{self.cfg.gse_id}_{self.cfg.tar_gene}_boxplot.png"
        self._save_plot(fig_name)


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
        self.meta_matrix_pack: dict = data.meta_matrix_pack

    def fig_plotter(self):
        """筛选所需目标并画图"""
        # 读取索引和数据仓库
        self._logger.info("读取索引和数据中...")
        if not self._gene_corr_table and not self._meta_matrix_pack:
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
        if not self._gene_corr_table or not self._meta_matrix_pack:
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
        }.get(self.cfg.analysis_mode)

        if summary_name is None:
            self._logger.error(f"未知分析模式：{self.cfg.analysis_mode}，无法加载数据")
            raise ValueError(f"未知分析模式：{self.cfg.analysis_mode}")

        summary_path = os.path.join(data_dir, "pkl", summary_name)
        data_pack_path = os.path.join(data_dir, "pkl", f"{gse_id}_processed_pack.pkl")

        if self.cfg.analysis_mode == "corr":
            self._gene_corr_table = pd.read_pickle(summary_path)
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
