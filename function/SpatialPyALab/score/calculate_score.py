import numpy as np
import pandas as pd
import scipy.sparse as sp

def add_expression_score(adata, gene_list, score_name, method='mean'):
    genes_present = [g for g in gene_list if g in adata.var_names]
    if not genes_present:
        raise ValueError(f"None of the genes in {gene_list} are found in adata.")
    expr = adata[:, genes_present].X
    if sp.issparse(expr):
        score = expr.sum(axis=1).A1 
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
