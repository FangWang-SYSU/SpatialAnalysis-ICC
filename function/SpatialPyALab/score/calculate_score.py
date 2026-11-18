import numpy as np
import pandas as pd
import scipy.sparse as sp

def add_expression_score(adata, gene_list, score_name, method='mean'):
    """
    给定一组基因，计算每个spot/cell的表达得分（总和或均值），添加到obs中。

    Parameters:
    - adata: AnnData 对象
    - gene_list: 要计算的基因列表
    - score_name: 得分列名
    - method: 'sum' 或 'mean'，决定是否对基因数量归一化
    """
    # 过滤掉不在 var_names 中的基因
    genes_present = [g for g in gene_list if g in adata.var_names]
    if not genes_present:
        raise ValueError(f"None of the genes in {gene_list} are found in adata.")

    # 提取表达矩阵（稀疏 or 密集）
    expr = adata[:, genes_present].X
    if sp.issparse(expr):
        score = expr.sum(axis=1).A1  # 转为一维 array
    else:
        score = expr.sum(axis=1)

    if method == 'mean':
        score = score / len(genes_present)

    adata.obs[score_name] = score


def compute_region_scores(adata, score_keys, region_key='region_label'):
    
    region_scores_dict = {}
    for key in score_keys:
        region_scores_df = (
            adata.obs
            .groupby(region_key)[key]
            .mean()
            .sort_index()
            .reset_index()
            .rename(columns={key: 'score', region_key: 'region'})
        )
        region_scores_dict[key] = region_scores_df
    return region_scores_dict


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def plot_region_scores_bar(score_dict, fig, axs):
    axs = axs.flatten()

    for ax, (score_name, scores) in zip(axs, score_dict.items()):
        region_df = pd.DataFrame({
            'Region': [f"region_{i+1}" for i in range(len(scores))],
            'Score': scores
        })

        # 按 Score 升序排序
        region_df = region_df.sort_values(by='Score', ascending=True).reset_index(drop=True)

        sns.barplot(
            data=region_df,
            y='Region',
            x='Score',
            palette=sns.color_palette("YlOrRd", len(region_df)),
            ax=ax
        )

        for i, row in region_df.iterrows():
            ax.text(row['Score'] + 0.01, i, f"{row['Score']:.2f}",
                    va='center', ha='left', fontsize=8)

        ax.set_title(f"{score_name} per Region")
        ax.set_xlabel("Average Score")
        ax.set_ylabel("Region")
        ax.set_xlim(0, max(region_df['Score']) + 0.1)
