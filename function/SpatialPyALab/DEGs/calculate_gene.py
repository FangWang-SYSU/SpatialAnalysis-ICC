import pandas as pd
import numpy as np
from scipy.sparse import issparse

def compute_gene_expression_ratio(adata, group_key="cell_type", min_cells=1):
    """
    计算每个基因在每个细胞类型中的表达比例。

    Parameters:
        adata: AnnData 对象，要求 adata.obs[group_key] 包含细胞类型信息
        group_key: str，用于分组的列名
        min_cells: int，最少包含细胞数的组才会被考虑

    Returns:
        pd.DataFrame: 行是基因，列是细胞类型，值是表达比例（0~1）
    """
    gene_names = adata.var_names
    group_labels = adata.obs[group_key].unique()

    result = pd.DataFrame(index=gene_names)

    for label in group_labels:
        cells_in_group = adata.obs[group_key] == label
        sub_adata = adata[cells_in_group]

        if sub_adata.n_obs < min_cells:
            continue

        if issparse(sub_adata.X):
            expr = sub_adata.X.toarray()
        else:
            expr = sub_adata.X

        # 非零表达细胞数除以总数
        expr_ratio = (expr > 0).sum(axis=0) / sub_adata.n_obs
        result[label] = expr_ratio

    return result

def filter_genes_by_expression_ratio(expr_ratio_df, min_ratio=0.01, fold_change=2):
    """
    筛选满足以下条件的基因：
    - 在至少一个细胞类型中表达比例 > min_ratio
    - 且该比例为所有其他细胞类型表达比例的 fold_change 倍以上

    Parameters:
        expr_ratio_df: DataFrame，行是基因，列是细胞类型，值是表达比例
        min_ratio: float，最低表达比例
        fold_change: float，倍数阈值

    Returns:
        筛选后基因名的列表
    """
    ratio_dict = {}

    for group in expr_ratio_df.columns:
        selected_genes = set()
        for gene, row in expr_ratio_df.iterrows():
            expr_value = row[group]
            if expr_value <= min_ratio:
                continue
            other_ratios = row.drop(index=group)
            if all(expr_value >= fold_change * other for other in other_ratios):
                selected_genes.add(gene)
        ratio_dict[group] = selected_genes

    return ratio_dict


import pandas as pd

def extract_significant_degs(adata, logfc_thresh=1.0, pval_thresh=0.01):
    """
    从 adata.uns['rank_genes_groups'] 中提取显著差异基因。
    
    Parameters:
        adata: AnnData 对象
        logfc_thresh: float，logFC 阈值
        pval_thresh: float，p 值阈值

    Returns:
        pd.DataFrame：包含 group, gene, logfc, pval 的显著基因表
    """
    res = adata.uns['rank_genes_groups']
    groups = res['names'].dtype.names
    all_degs = []

    for group in groups:
        genes = res['names'][group]
        logfc = res['logfoldchanges'][group]
        pvals = res['pvals'][group]

        df = pd.DataFrame({
            'group': group,
            'gene': genes,
            'logfc': logfc,
            'pval': pvals
        })

        sig_df = df[(df['logfc'] > logfc_thresh) & (df['pval'] < pval_thresh)]
        all_degs.append(sig_df)

    return pd.concat(all_degs, ignore_index=True)



def intersect_gene_sets(deg_df, ratio_dict):
    """
    将差异表达基因表与表达比例筛选结果取交集。

    参数：
        deg_df: pd.DataFrame，包含 group, gene, logfc, pval
        ratio_dict: dict，{cell_type: set of genes}

    返回：
        dict: {cell_type: sorted list of intersected genes}
    """
    deg_dict = (
        deg_df.groupby("group")["gene"]
        .apply(set)
        .to_dict()
    )

    intersect_dict = {
        group: sorted(deg_dict.get(group, set()) & ratio_dict.get(group, set()))
        for group in deg_dict.keys()
    }

    return intersect_dict






import numpy as np

def compute_cell_scores(
    adata, 
    gene_dict, 
    score_names=None, 
    use_raw=True, 
    method="mean"  # 可选："mean" 或 "nonzero_mean"
):
    if score_names is None:
        score_names = gene_dict.keys()

    for label in score_names:
        genes = gene_dict[label]
        genes_present = [g for g in genes if g in adata.var_names]

        if not genes_present:
            print(f"Warning: No genes for {label} found in adata.var_names")
            continue

        if use_raw and adata.raw is not None:
            expr = adata.raw[:, genes_present].X
        else:
            expr = adata[:, genes_present].X

        if method == "mean":
            if isinstance(expr, np.ndarray):
                mean_score = expr.mean(axis=1)
            else:
                mean_score = expr.mean(axis=1).A1  # sparse matrix to flat array
            adata.obs[f"{label}_mean"] = mean_score

        elif method == "nonzero_mean":
            if not isinstance(expr, np.ndarray):
                expr = expr.toarray()  # 稀疏转密集

            nonzero_counts = (expr > 0).sum(axis=1)
            expr_sum = expr.sum(axis=1)
            mean_expr = np.divide(expr_sum, nonzero_counts, out=np.zeros_like(expr_sum), where=nonzero_counts > 0)
            adata.obs[f"{label}_nonzero_mean"] = mean_expr

        else:
            raise ValueError(f"Unsupported method '{method}'. Choose 'mean' or 'nonzero_mean'.")

    return adata
