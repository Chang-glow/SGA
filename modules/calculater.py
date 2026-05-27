import os
import re
import scipy

import numpy as np
import pandas as pd

from abc import ABC, abstractmethod
from typing import Optional, List

from utils import loggers, Config, DataHandler, RESULT_DIR, resolve_save_path, df_content_hash, relpath
from modules.data_packer import DataPacker


_PROBE_RE = re.compile(
    r"^(ILMN_|AFFY-|A_\d|CUST_|GI_|NM_|NR_|XM_|XR_"
    r"|TC\d|Hs\.|Mm\.|Rn\.|DDB_|ENS[A-Z]*[EPT]|FBgn_|agc_|GE_)"
)

_SAMPLE_COL_EXCLUDE_KEYWORDS = [
    'ensembl', 'entrezid', 'symbol', 'genename', 'probeid',
    'id_ref', 'targetid', 'gene', 'description',
]


def get_sample_columns(df: pd.DataFrame) -> list:
    """识别表达矩阵中的样本列，排除注释列和检测 p-value 列"""
    sample_columns = []
    for col in df.columns:
        if isinstance(col, str):
            lowered = col.lower()
            if any(kw in lowered for kw in _SAMPLE_COL_EXCLUDE_KEYWORDS):
                continue
            if 'detection' in lowered and 'pval' in lowered:
                continue
            sample_columns.append(col)
        else:
            sample_columns.append(col)
    return sample_columns


def rename_expr_columns_by_meta_order(expr_df: pd.DataFrame, full_meta: pd.DataFrame,
                                       warn_on_fail: bool = True) -> tuple:
    """将表达矩阵列名映射为 GSM 样本名。优先按位置匹配，其次按元数据字段值匹配。

    Returns:
        (expr_df, True)  — 重命名成功
        (expr_df, False) — 重命名失败（列名未变），已 emit warning（若 warn_on_fail=True）
    """
    sample_columns = get_sample_columns(expr_df)
    if len(sample_columns) == len(full_meta.index):
        rename_map = {old: new for old, new in zip(sample_columns, full_meta.index.astype(str))}
        return expr_df.rename(columns=rename_map), True

    for meta_col in ['geo_accession', 'title', 'source_name_ch1', 'label_ch1']:
        if meta_col not in full_meta.columns:
            continue
        col_values = full_meta[meta_col].astype(str).tolist()
        if set(sample_columns).issubset(set(col_values)):
            rename_map = {
                sample_col: str(full_meta.index[col_values.index(sample_col)])
                for sample_col in sample_columns
            }
            return expr_df.rename(columns=rename_map), True
    if warn_on_fail:
        import logging
        logging.warning(
            f"表达矩阵列名未能对齐为 GSM ID：样本列数({len(sample_columns)}) "
            f"!= meta 行数({len(full_meta.index)})，"
            f"且样本列名与 meta 各字段值均不匹配，列名保持不变"
        )
    return expr_df, False
_NON_CODING_RES = [
    re.compile(r"^MIR", re.IGNORECASE),
    re.compile(r"^SNOR[A-Z]\d", re.IGNORECASE),
    re.compile(r"^SCARNA\d", re.IGNORECASE),
    re.compile(r"^RNU\d", re.IGNORECASE),
    re.compile(r"^RNVU\d", re.IGNORECASE),
    re.compile(r"^LINC\d", re.IGNORECASE),
    re.compile(r"^CASC\d", re.IGNORECASE),
    re.compile(r"^FLJ\d", re.IGNORECASE),
    re.compile(r"^RNA5", re.IGNORECASE),
    re.compile(r"^RN7SKP\d", re.IGNORECASE),
    re.compile(r"-AS\d", re.IGNORECASE),
    re.compile(r"-IT\d+$", re.IGNORECASE),
    re.compile(r"-DT\d*$", re.IGNORECASE),
    re.compile(r"^LOC\d", re.IGNORECASE),
    re.compile(r"^RP[LS]\d*[A-Z]*P\d", re.IGNORECASE),
    re.compile(r"[A-Z][A-Z0-9]*[A-Z]P\d+$", re.IGNORECASE),
    re.compile(r"\d+P\d"),
    re.compile(r"\d+P$"),
]

_TUPLE_PATTERNS = {
    "mirna": [re.compile(r"^MIR", re.IGNORECASE)],
}


def prepare_expr_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """基因名列设为索引，仅保留数值列，按基因去重取均值"""
    df = df.copy()

    gene_col = None
    for col in df.columns:
        if col in ("Gene", "Symbol", "SYMBOL", "GeneSymbol", "gene_symbol",
                   "Hugo_Symbol", "hugo", "Gene.Symbol", "GENE"):
            gene_col = col
            break
    if gene_col is None:
        for col in df.columns:
            if not pd.api.types.is_numeric_dtype(df[col]):
                gene_col = col
                break

    if gene_col is not None:
        df = df.set_index(gene_col)

    # 剔除残留的注释列（ENTREZID 等），仅保留样本表达列
    _ANNO_KEYWORDS = [
        'ensembl', 'entrezid', 'symbol', 'genename', 'probeid',
        'id_ref', 'targetid', 'gene', 'description',
    ]
    sample_cols = [
        c for c in df.columns
        if not (isinstance(c, str) and any(k in c.lower() for k in _ANNO_KEYWORDS))
    ]
    df = df[sample_cols]

    # 过滤基因索引中的注释类名称（如 ENTREZID 行）
    df = df[~df.index.astype(str).str.lower().isin(
        ["entrezid", "ensembl", "symbol", "probeid", "genename"]
    )]

    numeric_cols = df.select_dtypes(include=["number"]).columns
    df = df[numeric_cols]

    if df.index.duplicated().any():
        df = df.groupby(df.index).mean()

    return df


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
        numeric_df = vector.select_dtypes(include=[np.number])
        annotation_keywords = [
            'ENSEMBL', 'ENTREZID', 'SYMBOL', 'GENE', 'PROBEID',
            'ID_REF', 'TARGETID', 'DESCRIPTION', 'GENENAME'
        ]
        sample_columns = [
            col for col in numeric_df.columns
            if not any(keyword in str(col).upper() for keyword in annotation_keywords)
        ]
        vector = numeric_df[sample_columns]

        # 处理多行情况
        if isinstance(vector, pd.DataFrame):
            if vector.empty:
                logger.warning(f"基因 {tar_gene} 匹配行为空或无数值数据")
                return pd.Series(dtype=float)
            vector = vector.mean(axis=0)

        return vector

    logger.warning(f"未能在矩阵中找到基因: {tar_gene}")
    return pd.Series(dtype=float)


def detect_gene_case_convention(index: pd.Index, sample_size: int = 100) -> str:
    """从表达矩阵行索引采样，检测基因命名规则。

    Returns:
        'human' — 全大写（如 ACTA2, TNF）
        'mouse' — 首字母大写其余小写（如 Acta2, Tnf）
    """
    candidates: List[str] = []
    for g in index:
        s = str(g).strip()
        if not s or s.lower() in ('nan', 'none', ''):
            continue
        if s[0].isdigit():
            continue
        if len(s) < 2:
            continue
        candidates.append(s)
        if len(candidates) >= sample_size:
            break

    if not candidates:
        return 'human'

    upper_count = 0
    title_count = 0
    for g in candidates:
        alpha_chars = [c for c in g if c.isalpha()]
        if not alpha_chars:
            continue
        if all(c.isupper() for c in alpha_chars):
            upper_count += 1
        first_alpha_idx = next((i for i, c in enumerate(g) if c.isalpha()), None)
        if first_alpha_idx is not None:
            tail_alpha = [c for c in g[first_alpha_idx + 1:] if c.isalpha()]
            if (tail_alpha and g[first_alpha_idx].isupper()
                    and all(c.islower() for c in tail_alpha)):
                title_count += 1

    return 'human' if upper_count >= title_count else 'mouse'


def normalize_gene_symbol(gene: str, convention: str) -> str:
    """将单个基因符号转为指定命名规则格式。"""
    if not gene or not gene.strip():
        return gene
    if convention == 'human':
        return gene.upper()
    elif convention == 'mouse':
        return gene[0].upper() + gene[1:].lower()
    return gene


def resolve_gene_in_matrix(gene: str, matrix_index: pd.Index, convention: str,
                           homolog_map: dict | None = None) -> str:
    """在表达矩阵中查找基因的最佳匹配，返回矩阵中实际存在的基因名。

    依次尝试：原始名 → 规范化名 → 同源映射（来自 cfg.homolog_map）。均不匹配返回空字符串。
    """
    idx_set = {str(i) for i in matrix_index}
    # 1. 原始名
    if gene in idx_set:
        return gene
    # 2. 规范化名
    normalized = normalize_gene_symbol(gene, convention)
    if normalized != gene and normalized in idx_set:
        return normalized
    # 3. 同源映射（cfg.homolog_map: {human: mouse} or {mouse: human}）
    if homolog_map:
        # 查找: gene → mapped
        mapped = homolog_map.get(gene) or homolog_map.get(gene.upper()) or homolog_map.get(
            gene[0].upper() + gene[1:].lower() if gene else "")
        if mapped and mapped in idx_set:
            return mapped
        # 反向查找: mapped → gene
        reverse = {v: k for k, v in homolog_map.items()}
        rev_mapped = reverse.get(gene) or reverse.get(gene.upper()) or reverse.get(
            gene[0].upper() + gene[1:].lower() if gene else "")
        if rev_mapped and rev_mapped in idx_set:
            return rev_mapped
    return ""


def normalize_tar_gene_from_pack(cfg, meta_matrix_pack: dict) -> None:
    """从 pack 中提取表达矩阵，检测基因命名规则，原地修改 cfg.tar_gene。

    跳过 enrich 模式（无表达矩阵）或无 tar_gene 的情况。
    """
    import logging

    if not cfg.tar_gene:
        return
    if cfg.analysis_mode == "enrich":
        return
    if not meta_matrix_pack:
        return

    expr_df = None
    for key, val in meta_matrix_pack.items():
        if key in ('meta', 'meta_full', 'group_info', '_batch_exp_groups'):
            continue
        if isinstance(val, pd.DataFrame):
            expr_df = val
            break
        if isinstance(val, dict) and 'matrix_aligned' in val:
            expr_df = val['matrix_aligned']
            break

    if expr_df is None or expr_df.empty:
        return

    # 优先使用 SYMBOL 列检测命名规则（Ensembl 索引无法反映基因符号规范）
    index_for_detection = expr_df.index
    for col in ('SYMBOL', 'GENE', 'GENENAME'):
        if col in expr_df.columns:
            symbols = expr_df[col].dropna().astype(str).str.strip()
            symbols = symbols[symbols.str.fullmatch(r'[A-Za-z][A-Za-z0-9]+')]
            if len(symbols) > 10:
                index_for_detection = pd.Index(symbols)
                break

    # 若矩阵索引无法反映命名规则（如 Ensembl ID），回退到 pack 自动检测或 organism 配置
    idx_sample = expr_df.index.astype(str)
    if idx_sample.str.match(r'^ENS[A-Z]*G\d+').mean() >= 0.5:
        organism = meta_matrix_pack.get('_organism')
        if not organism:
            meta_full = meta_matrix_pack.get('meta_full')
            if meta_full is not None:
                from modules.data_packer import DataPacker
                organism = DataPacker._detect_organism(meta_full)
        if not organism:
            organism = getattr(cfg, 'organism', 'human')
        cfg.organism = organism
        convention = 'mouse' if organism.lower() == 'mouse' else 'human'
    else:
        convention = detect_gene_case_convention(index_for_detection)

    new_gene = normalize_gene_symbol(cfg.tar_gene, convention)
    if new_gene != cfg.tar_gene:
        logger = logging.getLogger(__name__)
        logger.info(
            f"基因符号已规范化: {cfg.tar_gene!r} → {new_gene!r}"
            f"（检测到 {convention} 命名规则）"
        )
        cfg.tar_gene = new_gene


def normalize_gene_index(expr_df: pd.DataFrame) -> pd.DataFrame:
    """去除基因索引的 Ensembl 版本号后缀，按表达均值去重"""
    clean = expr_df.index.astype(str).str.replace(r'\.\d+$', '', regex=True)
    if clean.equals(expr_df.index.astype(str)):
        return expr_df
    expr_df = expr_df.copy()
    expr_df.index = clean
    if expr_df.index.duplicated().any():
        numeric_cols = expr_df.select_dtypes(include=["number"]).columns
        expr_df["_mean_expr_"] = expr_df[numeric_cols].mean(axis=1)
        expr_df = expr_df.sort_values("_mean_expr_", ascending=False)
        expr_df = expr_df[~expr_df.index.duplicated(keep='first')]
        expr_df = expr_df.drop(columns=["_mean_expr_"])
    return expr_df


def map_ensembl_to_symbol(expr_df: pd.DataFrame) -> pd.DataFrame:
    """若行索引以 Ensembl ID 为主，通过 mygene 查询添加 SYMBOL 列"""
    idx_str = expr_df.index.astype(str)
    sample = idx_str[:min(20, len(idx_str))]
    ensembl_re = re.compile(r'^ENS[A-Z]*G\d+')
    if sum(1 for x in sample if ensembl_re.match(str(x))) < len(sample) * 0.8:
        return expr_df

    import mygene
    logger = loggers.get_logger()
    mg = mygene.MyGeneInfo()
    # biothings_client 默认 timeout=None（永不超时），强制设 30s 防止网络卡死
    import httpx as _httpx
    mg._set_http_client()
    mg.http_client.timeout = _httpx.Timeout(30.0)

    # 抑制第三方库日志噪音
    import logging as _logging
    for name in ('biothings.client', 'httpx', 'urllib3'):
        _logging.getLogger(name).setLevel(_logging.ERROR)

    clean_ids = idx_str.str.replace(r'\.\d+$', '', regex=True)
    unique_ids = clean_ids.unique().tolist()
    total_batches = (len(unique_ids) + 999) // 1000

    symbol_map = {}
    batch_size = 1000
    for i in range(0, len(unique_ids), batch_size):
        batch = unique_ids[i:i + batch_size]
        batch_no = i // batch_size + 1
        logger.debug(f"查询基因符号 {batch_no}/{total_batches} ...")
        try:
            for r in mg.querymany(batch, scopes='ensembl.gene', fields='symbol', species='all'):
                if 'symbol' in r:
                    symbol_map[r['query']] = r['symbol']
        except Exception:
            logger.warning(f"mygene 查询超时或失败 (batch {batch_no}/{total_batches})，跳过")
            continue

    # 恢复日志级别
    for name in ('biothings.client', 'httpx', 'urllib3'):
        _logging.getLogger(name).setLevel(_logging.NOTSET)

    if not symbol_map:
        return expr_df

    # 添加 SYMBOL 列，未命中则保留原始 Ensembl ID
    expr_df = expr_df.copy()
    expr_df['SYMBOL'] = clean_ids.map(symbol_map)
    mapped = expr_df['SYMBOL'].notna()
    logger.debug(
        f"Ensembl→Symbol 映射完成: {mapped.sum()}/{len(expr_df)} 有符号, "
        f"{len(expr_df) - mapped.sum()} 未命中"
    )
    # 未命中用原始 Ensembl ID 填充
    missing = expr_df['SYMBOL'].isna()
    expr_df.loc[missing, 'SYMBOL'] = idx_str[missing].astype(str)
    return expr_df


def map_mouse_genes_to_human_orthologs(mouse_genes, cache_dir=None, logger=None):
    """通过 mygene HomoloGene 数据库将小鼠基因符号映射为人类同源基因符号。

    两步查询：1) 小鼠符号 → HomoloGene 组 → 人类 Entrez Gene ID
              2) 人类 Entrez Gene ID → 人类基因符号

    Returns:
        dict: {mouse_symbol: [human_symbol, ...]}，无同源基因的 key 不出现在结果中
    """
    import json
    import mygene
    import httpx as _httpx
    import logging as _logging

    if logger is None:
        logger = loggers.get_logger()

    if cache_dir is None:
        from utils.paths import BASE_DIR
        cache_dir = os.path.join(BASE_DIR, "data", "cache", "orthologs")
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, "mouse_to_human.json")

    # 加载缓存
    cache = {}
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'r') as f:
                cache = json.load(f)
            cache.pop("_version", None)
        except (json.JSONDecodeError, OSError):
            logger.warning("同源基因缓存损坏，将重新构建")

    # 过滤出缓存中不存在的基因
    genes_to_query = [g for g in mouse_genes if str(g) not in cache]
    if not genes_to_query:
        logger.debug(f"所有 {len(mouse_genes)} 个小鼠基因已在同源缓存中")
        return {g: cache[str(g)] for g in mouse_genes if str(g) in cache}

    logger.info(f"正在查询 mygene 人鼠同源基因映射（{len(genes_to_query)} 个新基因）...")

    mg = mygene.MyGeneInfo()
    mg._set_http_client()
    mg.http_client.timeout = _httpx.Timeout(30.0)
    for name in ('biothings.client', 'httpx', 'urllib3'):
        _logging.getLogger(name).setLevel(_logging.ERROR)

    # Step 1: 小鼠基因符号 → HomoloGene 组 → 人类 Entrez Gene ID
    mouse_to_human_ids = {}
    batch_size = 1000
    total_batches = (len(genes_to_query) + batch_size - 1) // batch_size
    for i in range(0, len(genes_to_query), batch_size):
        batch = genes_to_query[i:i + batch_size]
        batch_no = i // batch_size + 1
        try:
            results = mg.querymany(batch, scopes='symbol', species='mouse', fields='homologene')
            for r in results:
                if 'homologene' not in r:
                    continue
                human_ids = [g[1] for g in r['homologene']['genes'] if g[0] == 9606]
                if human_ids:
                    mouse_to_human_ids[r['query']] = human_ids
        except Exception:
            logger.warning(f"mygene 同源查询失败 (batch {batch_no}/{total_batches})，跳过")

    for name in ('biothings.client', 'httpx', 'urllib3'):
        _logging.getLogger(name).setLevel(_logging.NOTSET)

    if not mouse_to_human_ids:
        logger.warning("mygene 同源查询无结果，将仅使用缓存中的同源映射")
        return {g: cache[str(g)] for g in mouse_genes if str(g) in cache}

    # Step 2: 人类 Entrez Gene ID → 基因符号
    all_human_ids = list(set(hid for ids in mouse_to_human_ids.values() for hid in ids))
    id_to_symbol = {}
    for i in range(0, len(all_human_ids), batch_size):
        batch = all_human_ids[i:i + batch_size]
        try:
            results = mg.querymany(batch, scopes='entrezgene', species='human', fields='symbol')
            for r in results:
                if 'symbol' in r:
                    id_to_symbol[str(r['query'])] = r['symbol']
        except Exception:
            logger.warning(f"mygene 符号查询失败 (batch {i // batch_size + 1})，跳过")

    # 构建映射并更新缓存
    new_entries = 0
    for mouse_gene, human_ids in mouse_to_human_ids.items():
        symbols = []
        for hid in human_ids:
            sym = id_to_symbol.get(str(hid))
            if sym:
                symbols.append(sym)
        if symbols:
            # 去重（多个 Entrez ID 可能对应同一符号，如不同转录本）
            symbols = list(dict.fromkeys(symbols))
            cache[str(mouse_gene)] = symbols
            new_entries += 1

    # 持久化缓存
    cache_out = {"_version": 1, **cache}
    try:
        with open(cache_file, 'w') as f:
            json.dump(cache_out, f, ensure_ascii=False)
    except OSError:
        logger.debug("同源基因缓存写入失败（不影响分析）")

    # 返回本次查询的所有基因的映射（含缓存命中的）
    result = {}
    for g in mouse_genes:
        if str(g) in cache:
            result[g] = cache[str(g)]

    logger.debug(
        f"同源基因映射完成: 新增 {new_entries} 条, 总计 {len(cache)} 条缓存, "
        f"{len(result)}/{len(mouse_genes)} 个基因匹配到人类同源基因"
    )
    return result


def apply_ortholog_conversion(expr_df, ortholog_map, logger=None):
    """将表达矩阵 index 从小鼠基因符号替换为人类同源基因符号。

    - 一对多同源（1 mouse → N human）：复制表达行
    - 多对一同源（N mouse → 1 human）：groupby.mean() 合并
    - 无同源基因（ortholog_map 中值为 [] 或不存在）：丢弃该行
    """
    if not ortholog_map:
        return expr_df

    if logger is None:
        logger = loggers.get_logger()

    expr_df = expr_df.copy()
    original_index = expr_df.index.astype(str)

    new_rows = []         # [(human_symbol, row_values_array)]
    multi_hit_count = 0
    dropped = []

    for gene in original_index:
        human_symbols = ortholog_map.get(gene)
        if human_symbols is None or len(human_symbols) == 0:
            dropped.append(gene)
            continue
        row_vals = expr_df.loc[gene].values
        for hs in human_symbols:
            new_rows.append((hs, row_vals))
        if len(human_symbols) > 1:
            multi_hit_count += 1

    if dropped:
        logger.debug(
            f"无同源基因被移除 ({len(dropped)} 个): "
            f"{dropped[:10]}{'...' if len(dropped) > 10 else ''}"
        )

    new_index = [r[0] for r in new_rows]
    new_data = [r[1] for r in new_rows]
    result = pd.DataFrame(new_data, index=new_index, columns=expr_df.columns)

    # 多对一：按 index 合并（groupby mean）
    if result.index.duplicated().any():
        result = result.groupby(result.index).mean()

    logger.info(
        f"同源转换: {len(original_index)} 基因 → {len(result)} 基因 "
        f"(一对多 {multi_hit_count}, 无同源移除 {len(dropped)})"
    )
    return result


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

    def _get_hfm_dict(self) -> dict:
        """返回自定义或默认标志物基因集"""
        custom = getattr(self.cfg, "custom_marker_dict", None)
        if custom and any(custom.values()):
            return custom
        return self.hfm_dict
    _logger = loggers.get_logger()

    def __init__(self, cfg: Config):
        """初始化分析对象

        Args:
            cfg: 基础配置
        """
        self.cfg = cfg
        self._meta_matrix_pack: Optional[dict] = None
        self._analysis_result: Optional[pd.DataFrame] = None
        self._corr_result: Optional[pd.DataFrame] = None
        self._diff_result: Optional[pd.DataFrame] = None
        self._hilo_result: Optional[pd.DataFrame] = None
        self._immune_result: Optional[pd.DataFrame] = None
        self._wgcna_result: Optional[pd.DataFrame] = None
        self._strategy = self._get_strategy()

    @classmethod
    def create(cls, cfg: Config, data: DataHandler):
        """根据cfg检查数据传入方式,创建分析对象"""
        data_dir = os.path.join(cfg.data_dir, cfg.gse_id)

        pack_path = DataPacker.resolve_pack_path(data_dir, cfg.gse_id, cfg.analysis_mode)

        if cfg.analysis_mode == "enrich":
            if data.meta_matrix_pack:
                return DataAnalyzer(cfg, data)
            return FileAnalyzer(cfg)

        if os.path.exists(pack_path) and not cfg.force:
            cls._logger.info(f"发现数据包：{pack_path}，将从数据包中分析")
            return FileAnalyzer(cfg)
        elif data.meta_matrix_pack:
            cls._logger.info("开始从数据包中分析")
            return DataAnalyzer(cfg, data)
        else:
            if cfg.force and "1" not in str(cfg.process):
                hint = (
                    "force=true 会跳过数据包缓存读取，且 process 不含 '1' 跳过了数据下载，"
                    "两者同时生效导致无数据可用。请设置 force=false 以读取缓存，或在 process 中加入 '1' 以重新下载数据。"
                )
            elif cfg.force:
                hint = "force=true 已跳过数据包缓存读取，且未加载到内存数据。请设置 force=false 后重试。"
            else:
                hint = "未找到数据包缓存，请确保先执行阶段1（数据下载与清洗）。"
            raise RuntimeError(hint)


    def calculate(self) -> pd.DataFrame:
        """用于调用数据的API"""
        mode = self.cfg.analysis_mode
        self._logger.info(f"分析模式设定为{mode}，将执行分析")

        # 如果已经计算过且结果非空，直接返回缓存结果
        if self._analysis_result is not None and not self._analysis_result.empty:
            return self._analysis_result

        try:
            # 执行分析策略
            result = self._strategy.calculate()
            if result is None or result.empty:
                self._logger.error("结果矩阵为空")
                if mode == "enrich":
                    return None
                raise ValueError("结果矩阵为空")

            # 在 padj/P-value 过滤前保存 raw CSV
            if mode in ("diff", "hilo") and self.cfg.storage and result is not None and not result.empty:
                self._data_storage(result, "csv", suffix="raw")

            # 对 diff / hilo 结果清洗探针、padj 过滤、限制基因数
            result = self._clean_diff_results(result, mode)

            # tar_tuple 基因类别过滤
            if self.cfg.tar_tuple and "Gene" in result.columns:
                result = self._apply_tuple_filter(result)

            # 根据分析模式将结果存储到对应属性，并执行存储
            self._analysis_result = result
            if mode == "corr":
                self._corr_result = result
            elif mode == "diff":
                self._diff_result = result
            elif mode == "hilo":
                self._hilo_result = result
            elif mode == "immune":
                self._immune_result = result
            elif mode == "wgcna":
                self._wgcna_result = result
            elif mode == "enrich":
                pass
            if self.cfg.storage:
                self._data_storage(result, "pkl")
                self._data_storage(result, "csv")

            return result
        
        except Exception:
            self._logger.exception("分析失败")
            raise

    def _get_strategy(self):
        """根据分析模式选择策略类"""
        mode = self.cfg.analysis_mode
        if mode == "corr":
            from .strategies.correlation import CorrelationStrategy
            return CorrelationStrategy(self)
        elif mode == "diff":
            from .strategies.difference import DiffStrategy
            return DiffStrategy(self)
        elif mode == "hilo":
            from .strategies.highlow import HighLowStrategy
            return HighLowStrategy(self)
        elif mode == "enrich":
            from .strategies.enrichment import EnrichStrategy
            return EnrichStrategy(self)
        elif mode == "immune":
            from .strategies.immune import ImmuneStrategy
            return ImmuneStrategy(self)
        elif mode == "wgcna":
            from .strategies.wgcna import WgcnaStrategy
            return WgcnaStrategy(self)
        else:
            self._logger.error(f"不支持的分析模式: {mode}，请检查配置文件")
            raise ValueError(f"不支持的分析模式: {mode}")

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
            return True  # 无数值列，视为已转化
        max_val = numeric_df.max().max()
        # 阈值可根据实际数据调整，默认为50
        return max_val <= self.cfg.log_threshold

    def _clean_dataframe(self, df: pd.DataFrame, skip_log: bool = False) -> pd.DataFrame:
        """清洗单个矩阵,去除NaN,Inf和自动监测log"""
        df = df.copy()

        # Inf -> NaN
        df.replace([np.inf, -np.inf], np.nan, inplace=True)

        # drop掉NaN行
        origin_rows = df.shape[0]
        df.dropna(axis=0, how='any', inplace=True)
        if df.shape[0] < origin_rows:
            self._logger.info(f"删除了{origin_rows - df.shape[0]}行NaN/Inf数据")

        if skip_log:
            return df

        # 判断log并转化没有log的矩阵
        if not self._is_log(df):
            self._logger.warning("原始数据未log转换,将执行log2(x+1)转换")
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            df[numeric_cols] = np.log2(df[numeric_cols] + 1)
        else:
            self._logger.debug("数据已log转换,跳过该步骤")

        return df

    @staticmethod
    def _tpm_convert(df: pd.DataFrame) -> pd.DataFrame:
        """将表达矩阵按样本标准化至 TPM 尺度（每列总和 × 1M，用于免疫浸润分析）"""
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if numeric_cols.empty:
            return df
        col_sums = df[numeric_cols].sum(axis=0).replace(0, np.nan)
        df = df.copy()
        df[numeric_cols] = df[numeric_cols].div(col_sums, axis=1) * 1e6
        return df

    def _apply_tuple_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """按 tar_tuple 配置对结果做基因类别正向筛选"""
        patterns = _TUPLE_PATTERNS.get(self.cfg.tar_tuple)
        if patterns is None:
            self._logger.warning(f"未知的 tar_tuple: {self.cfg.tar_tuple}，跳过过滤")
            return df

        gene_col = df["Gene"].astype(str).str.strip()
        keep = pd.Series(False, index=df.index)
        for pat in patterns:
            keep = keep | gene_col.str.match(pat, na=False)
        df = df.loc[keep]
        self._logger.info(f"tar_tuple='{self.cfg.tar_tuple}' 过滤后保留 {len(df)} 个基因")
        return df

    def _clean_diff_results(self, df: pd.DataFrame, mode: str) -> pd.DataFrame:
        """清洗 diff / hilo 结果：剔除探针 ID、按显著性过滤、限制基因数

        strict_filter=True:  padj < p_threshold（严格 FDR 校正）
        strict_filter=False: P_value < p_threshold 且 |log2FC| > log2fc_threshold（宽松）
        """
        if mode not in ("diff", "hilo") or df is None or df.empty:
            return df

        before = len(df)

        if "Gene" in df.columns:
            gene_col = df["Gene"].astype(str).str.strip()
            keep = (
                gene_col.notna()
                & (gene_col != "")
                & (gene_col != "nan")
                & (gene_col != "None")
                & ~gene_col.str.fullmatch(r"\d+")
            )
            # tar_tuple 模式跳过探针/非编码剔除，由 _apply_tuple_filter 正向筛选
            if not self.cfg.tar_tuple:
                keep = keep & ~gene_col.str.match(_PROBE_RE)
                for ncrna_re in _NON_CODING_RES:
                    keep = keep & ~gene_col.str.contains(ncrna_re, regex=True, na=False)

            blacklist = getattr(self.cfg, "gene_blacklist", [])
            if blacklist:
                black_upper = [g.upper() for g in blacklist]
                keep = keep & ~gene_col.str.upper().isin(black_upper)

            df = df.loc[keep]

        strict = getattr(self.cfg, "strict_filter", True)
        p_thr = self.cfg.p_threshold

        if strict:
            if "padj" in df.columns:
                p_before = len(df)
                df = df[df["padj"].notna() & (df["padj"] < p_thr)]
                self._logger.info(f"严格模式 padj < {p_thr}: {p_before} → {len(df)}")
        else:
            if "P_value" in df.columns:
                p_before = len(df)
                if mode == "hilo":
                    df = df[df["P_value"].notna() & (df["P_value"] < p_thr)]
                    self._logger.info(f"hilo宽松模式 P_value < {p_thr}: {p_before} → {len(df)}")
                elif "log2FC" in df.columns:
                    fc_thr = getattr(self.cfg, "log2fc_threshold", 0.5)
                    df = df[
                        df["P_value"].notna() & (df["P_value"] < p_thr)
                        & df["log2FC"].notna() & (df["log2FC"].abs() > fc_thr)
                    ]
                    self._logger.info(f"宽松模式 P_value < {p_thr} & |log2FC| > {fc_thr}: {p_before} → {len(df)}")

        max_genes = getattr(self.cfg, "max_output_genes", 0)
        if max_genes > 0 and len(df) > max_genes:
            sort_col = "padj" if strict else "P_value"
            df = df.sort_values(sort_col).head(max_genes)

        removed = before - len(df)
        if removed:
            self._logger.info(f"diff/hilo 结果清洗完成，剔除 {removed} 行，保留 {len(df)} 行。")
        return df

    def _get_sample_columns(self, df: pd.DataFrame) -> list:
        return get_sample_columns(df)

    def _rename_expr_columns_by_meta_order(self, expr_df: pd.DataFrame, full_meta: pd.DataFrame):
        return rename_expr_columns_by_meta_order(expr_df, full_meta)

    def _clean_pack(self, raw_pack: dict) -> dict:
        """清洗整个数据包,处理其中的DataFrame矩阵并保留非矩阵元数据"""
        if self._meta_matrix_pack is None:
            is_immune = self.cfg.analysis_mode == "immune"
            skip_log = self.cfg.analysis_mode in ("immune", "wgcna")
            cleaned_pack = {}
            for name, item in raw_pack.items():
                if name in {"meta", "meta_full"}:
                    cleaned_pack[name] = item.copy()
                elif isinstance(item, dict):
                    if 'matrix_aligned' in item and isinstance(item['matrix_aligned'], pd.DataFrame):
                        df = self._clean_dataframe(item['matrix_aligned'], skip_log=skip_log)
                        if is_immune:
                            self._logger.info("正在进行 TPM 转换（免疫浸润分析预处理）...")
                            df = self._tpm_convert(df)
                        cleaned_pack[name] = df
                    else:
                        # 保留非矩阵形式的字典元数据，如 group_info 或其他配置
                        cleaned_pack[name] = item.copy()
                elif isinstance(item, pd.DataFrame):
                    df = self._clean_dataframe(item, skip_log=skip_log)
                    if is_immune:
                        self._logger.info("正在进行 TPM 转换（免疫浸润分析预处理）...")
                        df = self._tpm_convert(df)
                    cleaned_pack[name] = df
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
        """最显著值,corr模式为按相关性R排序,diff模式以及hilo模式为按log2FC排序"""
        # 根据分析模式返回显著结果
        if self._analysis_result is None or self._analysis_result.empty:
             self._logger.warning("分析结果为空，无法提取显著值")
             return None
        
        if self.cfg.analysis_mode == "corr":
            if self._corr_result is None:
                df = self.calculate()
            else:
                df = self._corr_result
            significant_findings = df[df['P_value'] < self.cfg.p_threshold].sort_values('R')
        
        elif self.cfg.analysis_mode in ["diff", "hilo"]:
            if self._diff_result is None and self._hilo_result is None:
                df = self.calculate()
            else:
                df = self._diff_result if self._diff_result is not None else self._hilo_result
            significant_findings = df[df['padj'] < self.cfg.p_threshold].sort_values('log2FC')

        elif self.cfg.analysis_mode == "enrich":
            if self._analysis_result is None or self._analysis_result.empty:
                df = self.calculate()
            else:
                df = self._analysis_result
            significant_findings = df[df.get('Adjusted P-value', pd.Series(dtype=float)) < self.cfg.p_threshold].sort_values('Adjusted P-value')
        
        else:
            self._logger.error(f"不支持的分析模式: {self.cfg.analysis_mode}，请检查配置文件")
            raise ValueError(f"不支持的分析模式: {self.cfg.analysis_mode}")
        
        return significant_findings

    def _data_storage(self, result_df: pd.DataFrame, save_format: str, suffix: str = None):
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
        
        # 映射分析模式向文件名关键词
        mode_mapping = {
            "corr": "correlation",
            "diff": "differential",
            "hilo": "highlow",
            "enrich": "enrichment",
            "immune": "immune",
            "wgcna": "wgcna"
        }

        # 构建文件名和路径
        keyword = mode_mapping.get(self.cfg.analysis_mode, 'unknown')
        if suffix and save_format == "csv":
            file_name = f"{gse_id}_{keyword}_summary_{suffix}.{ext}"
        else:
            file_name = f"{gse_id}_{keyword}_summary.{ext}"

        if save_format == "csv":
            storage_path = os.path.join(RESULT_DIR, save_format, file_name)
        else:
            storage_path = os.path.join(data_dir, ext, file_name)
        os.makedirs(os.path.dirname(storage_path), exist_ok=True)

        # 防覆盖 + 内容去重
        storage_path = resolve_save_path(storage_path, df_content_hash(result_df))
        if storage_path is None:
            self._logger.info(f"{fmt} 文件内容与已有文件相同，跳过保存")
            return

        save_mode = getattr(result_df, f"to_{fmt}")

        # 尝试保存文件，并捕获可能的异常
        try:
            if fmt == "csv":
                save_mode(storage_path, index=False)
            else:
                save_mode(storage_path)
            self._logger.info(f"{fmt}文件已保存至 {relpath(storage_path)}")
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
        self._meta_matrix_pack = data.meta_matrix_pack
        self._logger.info("将从打包的pack中读取数据")

    def _load_data(self) -> dict:
        self._meta_matrix_pack = self._clean_pack(self._meta_matrix_pack)
        return self._meta_matrix_pack


class FileAnalyzer(Analyzer):
    """基于读取pkl文件的数据分析流程"""
    def __init__(self, cfg: Config):
        super().__init__(cfg)
        data_dir = os.path.join(self.cfg.data_dir, self.cfg.gse_id)
        self.data_path = DataPacker.resolve_pack_path(data_dir, self.cfg.gse_id, self.cfg.analysis_mode)
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
