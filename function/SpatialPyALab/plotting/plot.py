import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
from skimage import measure

from rich.console import Console
from rich.table import Table
from IPython.display import display

import matplotlib.pyplot as plt

def find_optimal_df(D_max, alpha, max_iter=1000):
    from scipy.stats import chi2

    # 处理自由度为0或v=1的情况
    if D_max < chi2.isf(alpha, 1):
        return 0

    # 倍增法确定初始上界high
    high = 1
    while chi2.isf(alpha, high) <= D_max:
        high *= 2
        # 防止无限循环，设置最大迭代次数
        max_iter -= 1
        if max_iter <= 0:
            raise ValueError("Exceeded maximum iterations while finding upper bound")

    low = high // 2

    # 二分查找
    max_iter = 1000  # 重置最大迭代次数
    while low < high and max_iter > 0:
        mid = (low + high + 1) // 2  # 向上取整
        q = chi2.isf(alpha, mid)
        if q <= D_max:
            low = mid
        else:
            high = mid - 1
        max_iter -= 1

    if max_iter <= 0:
        raise ValueError("Exceeded maximum iterations during binary search")

    return low


def plot_deviance(adata, n_seq = 5000, alpha = 0.001, figsize = (10, 6)):
    from scipy.stats import chi2

    assert (n_seq > 0) and (isinstance(n_seq, int)), "n_seq must be a positive integer"
    assert (alpha > 0) and (alpha < 1), "alpha must be a value between 0 and 1"

    deviance = adata.var['deviance']
    deviance_info = adata.uns['deviance']
    sorted_idx = np.argsort(deviance)[::-1]
    sorted_deviance = deviance[sorted_idx]

    critical_deviance = chi2.isf(alpha, len(deviance)-1)
    max_genes = (sorted_deviance > critical_deviance).sum()

    plt.figure(figsize = figsize)
    plt.scatter(x = np.arange(len(sorted_deviance)), y = sorted_deviance, label = f"{deviance_info['fam'].capitalize()} Deviance")

    plt.axvline(x = max_genes, linestyle="--", color="red")
    plt.text(max_genes+10, max(sorted_deviance) * 0.95, 
             f"DF: {len(deviance)-1}", color="gray", ha="center", fontsize=9)
    for i in range(1, len(sorted_deviance) // n_seq + 1):
        pos = i * 5000
        opt_df = find_optimal_df(sorted_deviance[pos-1], alpha)
        if pos <= len(sorted_deviance):
            plt.axvline(x=pos - 1, linestyle="--", color="red")
            plt.text(pos, max(sorted_deviance) * 0.95, 
                     f"DF: {opt_df}", color="gray", ha="center", fontsize=9)

    plt.xlabel("Ranked Genes")
    plt.ylabel(f"{deviance_info['fam'].capitalize()} Deviance")
    plt.title("Feature Selection with Deviance")
    plt.legend(loc='upper right')
    plt.show()

# To be Modified
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math

def stat_data(ax, data):
    data = data[~np.isnan(data)]
    if len(data) == 0:
        return
    median = np.median(data)
    mean = np.mean(data)
    mode = data[data > 0].mode().iloc[0] if np.any(data > 0) else 0
    q1 = np.percentile(data, 1)
    q99 = np.percentile(data, 99)

    _ = ax.axvline(median, color='r', linestyle='dashed', linewidth=1)
    _ = ax.text(0.6 * (np.max(data) - np.min(data)) + np.min(data),
                0.8 * ax.get_ylim()[1], f'Median: {median:.2f}',
                fontsize=10, verticalalignment='bottom')
    _ = ax.text(0.6 * (np.max(data) - np.min(data)) + np.min(data),
                0.7 * ax.get_ylim()[1], f'Mode: {mode:.2f}',
                fontsize=10, verticalalignment='bottom')
    _ = ax.text(0.6 * (np.max(data) - np.min(data)) + np.min(data),
                0.6 * ax.get_ylim()[1], f'Mean: {mean:.2f}',
                fontsize=10, verticalalignment='bottom')
    _ = ax.text(0.6 * (np.max(data) - np.min(data)) + np.min(data),
                0.5 * ax.get_ylim()[1], f'1%: {q1:.2f}',
                fontsize=10, verticalalignment='bottom')
    _ = ax.text(0.6 * (np.max(data) - np.min(data)) + np.min(data),
                0.4 * ax.get_ylim()[1], f'99%: {q99:.2f}',
                fontsize=10, verticalalignment='bottom')

def plot_obs_histograms(adata, keys, bins=60, n_cols=4, show=False):
    n_keys = len(keys)
    n_rows = math.ceil(n_keys / n_cols)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axs = np.array(axs).reshape(-1)  

    any_plotted = False

    for i, key in enumerate(keys):
        ax = axs[i]
        if key not in adata.obs.columns:
            ax.set_visible(False)
            continue

        data = adata.obs[key]
        
        sns.histplot(data, kde=False, bins=bins, ax=ax)
        ax.set_title(f"{key} Distribution")
        ax.set_xlabel(key)
        ax.set_ylabel("Cell Count")
        stat_data(ax, data)

    for j in range(n_keys, len(axs)):
        axs[j].set_visible(False)

    plt.tight_layout() 

    if show:
        plt.show()

    return fig

    
# To be Modified
def print_rich_log(log_path, console = Console()):
    # 创建表格
    table = Table(title="STAR Genome Log", show_header=True, header_style="bold magenta")
    table.add_column("Category", style="dim", width=20)
    table.add_column("Details", style="cyan")
    table.add_column("Results", style="green")

    sub_dictionary = {}
    currrent_dictionary = {}
    sub_iloc = 0

    def print_section(description, dict):
        i = 0

        for key, value in dict.items():
            if i == 0:
                table.add_row(description, key, value)
            else:
                table.add_row(None, key, value)
            i += 1
        table.add_section()

    with open(log_path, "r") as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines if line.strip()]

    for line in lines:
        if line.endswith(":"):
            sub_dictionary[sub_iloc] = currrent_dictionary
            currrent_dictionary = {}
            sub_iloc += 1
        else:
            key, value = line.split(" |\t")
            if key.endswith(" on"):
                key = key.split(" on")[0]
            if key.endswith(" %"):
                key = key.split(" %")[0]
            if key.endswith(", %"):
                key = key.split(", %")[0]
            currrent_dictionary[key] = value

    print_section("SUMMARY", sub_dictionary[0])

    print_section("UNIQUE READS", sub_dictionary[1])

    print_section("MULTI-MAPPING READS", sub_dictionary[2])

    print_section("UNMAPPED READS", sub_dictionary[3])

    # 打印表格
    console.print(table)

# To be Modified
def print_rich_adata(adata, console = Console()):
    n_points, n_genes = adata.shape

    # 创建表格
    table = Table(title=f"AnnData Structure: {n_points} points X {n_genes} genes", show_header=True, header_style="bold magenta")
    table.add_column("Category", style="dim", width=20)
    table.add_column("Description", style="cyan")
    table.add_column("Details", style="green")

    # 添加obs信息
    obs_keys = ", ".join(adata.obs.keys())
    table.add_row("obs", "Stats information by cells", obs_keys)
    table.add_section()

    # 添加var信息
    var_keys = ", ".join(adata.var.keys())
    table.add_row("var", "Stats information by genes", var_keys)
    table.add_section()

    # 添加uns信息
    uns_keys = ", ".join(adata.uns.keys())
    table.add_row("uns", "Unstructured annotations", uns_keys)
    table.add_section()

    # 添加obsm信息
    obsm_keys = ", ".join(adata.obsm.keys())
    table.add_row("obsm", "Spatial coordinates", obsm_keys)
    table.add_section()

    # 打印表格
    console.print(table)

# To be Modified
def print_rich_table(puck_path, console = Console()):
    df = pd.read_csv(puck_path, sep = ',', header = 0)
    df = df.drop(columns=["puck_barcode_file", "parsed_barcode_file"])
    df.px_by_um = df.px_by_um.round(5)

    n_rows = df.shape[0]

    # 创建表格
    table = Table(title=f"Puck Barcode Files: {n_rows} tiles", show_header=True, header_style="bold magenta")

    # 添加表头（排除不要的列）
    for col in df.columns:
        if col == 'puck_barcode_file_id':
            table.add_column(col, style = 'cyan', min_width = 25, justify='center')
        elif col == 'px_by_um':
            table.add_column(col, style = 'cyan', min_width = 10, justify='center')
        else:
            table.add_column(col, style = 'cyan', min_width = 14, justify='center')

    # 添加数据到表格
    for row in df.itertuples(index=False):
        table.add_row(*map(str, row))
        table.add_section()

    # 打印表格
    display(table)

def draw_single_cluster(adata, cluster_id, key = 'leiden', fig_size = (6, 10), spot_size = 500, basis = 'spatial'):
    plot_adata = adata.copy()
    orig_palette = adata.uns[f'{key}_colors']
    new_palette = {str(i): orig_palette[i] if i == int(cluster_id) else 'lightgrey' for i in range(len(orig_palette))}
    fig, axs = plt.subplots(figsize = fig_size)
    sc.pl.spatial(plot_adata, color = key, img_key = None, spot_size = spot_size, basis = basis, palette = new_palette, ax = axs)


def compare_two_cluster(adata, target, reference, groupby = None, key = None, x_interval = (-10, 10), y_interval = (0, 75000)):
    if groupby is None:
        groupby = 'leiden'
    if key is None:
        key = 'rank_genes_groups'
    plot = adata[adata.obs[groupby].isin([target, reference])].copy()
    sc.tl.rank_genes_groups(plot, groupby = groupby, groups = [target], reference = reference, method = 'wilcoxon')
    score = pd.DataFrame(plot.uns[key]['scores'])
    names = pd.DataFrame(plot.uns[key]['names'])

    scores_values = score.iloc[:, 0].values
    genes = names.iloc[:, 0].values
    y_values = plot.var.loc[genes, 'deviance'].values.flatten()

    colors = []
    red_indices = []
    blue_indices = []
    purple_indices = []
    black_indices = []
    
    for i in range(len(scores_values)):
        x = scores_values[i]
        y = y_values[i]
        
        if y < y_interval[0]:
            colors.append('black')
            black_indices.append(i)
        elif y > y_interval[1]:
            colors.append('purple')
            purple_indices.append(i)
        elif x > x_interval[1]:
            colors.append('red')
            red_indices.append(i)
        elif x < x_interval[0]:
            colors.append('blue')
            blue_indices.append(i)
        else:
            colors.append('grey')

    plt.figure(figsize=(12, 12))
    scatter = plt.scatter(
        x=scores_values,
        y=y_values,
        c=colors,
        alpha=0.7,
        edgecolors='w',
        linewidths=0.5
    )
    
    def add_labels(indices, color):
        for i in indices:
            plt.text(
                x=scores_values[i] + 0.5,  # 微调x位置防止重叠
                y=y_values[i],
                s=genes[i],
                color=color,
                fontsize=8,
                alpha=0.8,
                verticalalignment='center'
            )
    
    add_labels(red_indices, 'red')
    add_labels(blue_indices, 'blue')
    add_labels(black_indices, 'black')
    add_labels(purple_indices, 'purple')

    # 添加图例
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label=f'x > {x_interval[1]} ({len(red_indices)})',
                    markerfacecolor='red', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label=f'x < {x_interval[0]} ({len(blue_indices)})',
                    markerfacecolor='blue', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label=f'y < {y_interval[0]} ({len(black_indices)})',
                    markerfacecolor='black', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label=f'y > {y_interval[1]} ({len(purple_indices)})',
                    markerfacecolor='purple', markersize=8)
    ]
    
    plt.legend(handles=legend_elements, loc='upper right')
    
    plt.ylabel("Genes Poisson Deviance")
    plt.xlabel("Wilcoxon Rank Test Score")
    plt.title(f"Cluster {target} vs {reference} Differential Expression")
    plt.grid(alpha=0.2)
    plt.show()


def subset_adata_to_mask(adata, mask, spatial_key: str = 'spatial'):
    # Subset adata to the valid coordinates from the mask
    adata = adata[(adata.obsm[spatial_key][:, 0] >= 0) & 
                  (adata.obsm[spatial_key][:, 1] >= 0) &
                  (adata.obsm[spatial_key][:, 0] < mask.shape[0]) & 
                  (adata.obsm[spatial_key][:, 1] < mask.shape[1])].copy()

    # Subset the labels to those in the mask
    labels = mask[adata.obsm[spatial_key][:, 0].astype(int), adata.obsm[spatial_key][:, 1].astype(int)]

    # Set first label as zero for background
    labels[0] = 0

    # Assign label as cell_ID
    adata.obs["cell_ID"] = labels

    # Get centroid and label ID from mask
    props = measure.regionprops_table(mask, properties=["label", "centroid"])
    props = pd.DataFrame(props)

    props_filter = props[props.label.isin(np.unique(adata.obs["cell_ID"]))]
    return adata, props_filter