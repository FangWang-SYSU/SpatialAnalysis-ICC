import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_areas(areas, bins = 100):
    # 创建画布和坐标轴
    fig, ax = plt.subplots(figsize=(7, 4))

    # 自动计算合适的对数分箱
    min_area = max(1, np.min(areas))  # 确保最小值为1（0会导至对数计算错误）
    max_area = np.max(areas)
    bins = np.logspace(np.log10(min_area), np.log10(max_area), bins)

    # 绘制直方图
    counts, bins, patches = ax.hist(areas, 
                                bins=bins,
                                edgecolor='black',
                                alpha=0.7)

    # 设置对数坐标轴
    ax.set_xscale('log')
    ax.set_yscale('log')

    # 设置坐标轴标签和格式
    ax.set_xlabel('Area (log scale)', fontsize=10)
    ax.set_ylabel('Count (log scale)', fontsize=10)
    ax.set_title('Connected Components Area Distribution', fontsize=10)

    # 优化刻度显示
    ax.tick_params(axis='both', which='both', labelsize=10)
    # ax.grid(True, which='both', linestyle='--', alpha=0.5)

    # 自动调整坐标轴范围
    ax.set_xlim(left=min_area*0.9, right=max_area*1.1)
    ax.set_ylim(bottom=0.9, top=np.max(counts)*10)

    # 添加色标增强可读性
    cmap = plt.cm.viridis
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    for c, p in zip(bin_centers, patches):
        plt.setp(p, 'facecolor', cmap(np.log10(c)/np.log10(max_area)))

    # 添加颜色条加强数据解读
    sm = plt.cm.ScalarMappable(cmap=cmap, 
                            norm=plt.Normalize(vmin=np.log10(min_area), 
                                            vmax=np.log10(max_area)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Log10(Area)', rotation=270, labelpad=15)

    plt.tight_layout()
    plt.show()