from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from anndata import AnnData
    from numpy import ndarray
    from scipy.sparse import csr_matrix
    from frozendict import frozendict


from rich.console import Console
from rich.table import Table
import numpy as np




def print_stats_table(stat: frozendict):
    console = Console()
    table = Table(title="Stat Info", show_header=True, header_style="bold magenta")

    def print_section(description, dict):
        i = 0

        for key, value in dict.items():
            if i == 0:
                table.add_row(description, key, value)
            else:
                table.add_row(None, key, value)
            i += 1
        table.add_section()

    table.add_column("Category", style="green", no_wrap=True)
    table.add_column("Description", style="green", no_wrap=True)
    table.add_column("Value", style="bold purple")

    for category, metrics in stat.items():
        print_section(category, metrics)

    # 打印表格
    console.print(table)


import scanpy as sc
import numpy as np


def stat_data(ax, data):
    data = data[~np.isnan(data) & ~np.isinf(data)]
    median = np.median(data)
    mean = np.mean(data)
    percentile_1 = np.percentile(data, 1)
    percentile_99 = np.percentile(data, 99)
    mode = data[data > 0].mode().iloc[0]
    _ = ax.axvline(median, color='r', linestyle='dashed', linewidth=1)
    _ = ax.text(0.6*(np.max(data) - np.min(data)) + np.min(data),
                0.9 * ax.get_ylim()[1], f'Median: {median:.2f}', 
                fontsize=12, verticalalignment='bottom')
    _ = ax.text(0.6*(np.max(data) - np.min(data)) + np.min(data), 
                0.8 * ax.get_ylim()[1], f'Mode: {mode:.2f}', 
                fontsize=12, verticalalignment='bottom')
    _ = ax.text(0.6*(np.max(data) - np.min(data)) + np.min(data), 
                0.7 * ax.get_ylim()[1], f'Mean: {mean:.2f}', 
                fontsize=12, verticalalignment='bottom')
    _ = ax.text(0.6*(np.max(data) - np.min(data)) + np.min(data), 
                0.6 * ax.get_ylim()[1], f'1%: {percentile_1:.2f}', 
                fontsize=12, verticalalignment='bottom')
    _ = ax.text(0.6*(np.max(data) - np.min(data)) + np.min(data), 
                0.5 * ax.get_ylim()[1], f'99%: {percentile_99:.2f}', 
                fontsize=12, verticalalignment='bottom')