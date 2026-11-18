
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.io import mmread, mmwrite
import scipy.sparse as sp
from scipy.spatial import ConvexHull


from scipy.interpolate import splprep, splev
from shapely.ops import unary_union
from shapely.geometry import Polygon, MultiPolygon, MultiPoint, LinearRing
from shapely import STRtree

def split_geoms_to_polygons(alpha_shape):
    """
    将 alpha_shape 拆分为 Polygon 的列表[^1^]。
    
    参数：
    alpha_shape: 由 alphashape.alphashape 函数生成的几何对象。
    
    返回值：
    polygon_list: 包含所有 Polygon 对象的列表。
    """
    polygon_list = []
    
    # 如果 alpha_shape 是 Polygon 类型，直接添加到列表中
    if isinstance(alpha_shape, Polygon):
        polygon_list.append(alpha_shape)
    # 如果 alpha_shape 是 MultiPolygon 类型，提取其中的每个 Polygon
    elif isinstance(alpha_shape, MultiPolygon):
        polygon_list.extend(alpha_shape.geoms)
    # 如果 alpha_shape 是其他类型（如 GeometryCollection），需要进一步处理
    else:
        # 遍历 alpha_shape 的所有几何元素
        for geom in alpha_shape.geoms:
            # 如果几何元素是 Polygon 类型，添加到列表中
            if isinstance(geom, Polygon):
                polygon_list.append(geom)
            # 如果几何元素是 MultiPolygon 类型，提取其中的每个 Polygon
            elif isinstance(geom, MultiPolygon):
                polygon_list.extend(geom.geoms) 
    return polygon_list



def smooth_polygon(polygon, smooth_factor=0.1):
    x = np.array(polygon.exterior.xy[0])
    y = np.array(polygon.exterior.xy[1])
    tck, u = splprep([x, y], s=smooth_factor)
    new_points = np.array(splev(u, tck)).T
    return Polygon(new_points)


def process_polygons(alpha_shape, smooth_factor=0.1, buffer_distance=1):
    """
    处理 alpha_shape, 包括平滑边界和合并重叠多边形。
    """
    # 拆分 alpha_shape 为 Polygon 列表
    polygons = split_geoms_to_polygons(alpha_shape)

    # 平滑每个多边形边界
    smoothed_polygons = [smooth_polygon(poly, smooth_factor) for poly in polygons]

    # 对每个平滑后的多边形进行缓冲区操作
    buffered_polygons = [poly.buffer(buffer_distance) for poly in smoothed_polygons]

    # 合并重叠的多边形
    merged_polygon = unary_union(buffered_polygons)

    # 修复可能的拓扑错误
    if not merged_polygon.is_valid:
        merged_polygon = merged_polygon.buffer(0)

    return merged_polygon

from shapely.geometry import Point, Polygon, MultiPolygon

def filter_polygons_by_cell_count(cell_coords, threshold, polygons):
    """
    根据每个多边形内的细胞数进行筛选，返回包含细胞数超过阈值的多边形。
    
    参数：
        cell_coords: np.ndarray 或 list of [x, y]，细胞坐标。
        threshold: int,保留多边形所需的最小细胞数量。
        polygons: shapely MultiPolygon 或 list of Polygon。
    
    返回：
        list of Polygon, 满足条件的多边形。
    """
    # 转换为 shapely 的 Point 列表
    cell_points = [Point(xy) for xy in cell_coords]

    if isinstance(polygons, MultiPolygon):
        poly_list = list(polygons.geoms)
    else:
        poly_list = list(polygons)

    filtered_polys = []
    for poly in poly_list:
        count = sum(poly.contains(p) for p in cell_points)
        if count >= threshold:
            filtered_polys.append(poly)

    return filtered_polys



def plot_polygons_boundaries(buffered_layers, 
                             ax = None, 
                             edgecolor = 'black', 
                             linewidth = 0.5, 
                             alpha = 1,
                             linestyle ='-',
                             boundary = None):
    if ax is None:
        fig, ax = plt.subplots()
    for layer in buffered_layers:
        if boundary:
            layer = layer.intersection(boundary)
            
        if isinstance(layer, MultiPolygon):
            for geom in layer.geoms:
                xs, ys = geom.exterior.xy
                ax.plot(xs, ys, color = edgecolor, linewidth = linewidth, alpha = alpha,
                       linestyle=linestyle)
                for interior in geom.interiors:
                    xs, ys = interior.xy
                    ax.plot(xs, ys, color = edgecolor, linewidth = linewidth, alpha = alpha,
                           linestyle=linestyle )
        else:
            xs, ys = layer.exterior.xy
            ax.plot(xs, ys, color = edgecolor, linewidth = linewidth, alpha = alpha,
                   linestyle=linestyle)
            for interior in layer.interiors:
                xs, ys = interior.xy
                ax.plot(xs, ys, color = edgecolor, linewidth = linewidth, alpha = alpha,
                       linestyle=linestyle)
    return ax


import numpy as np
from scipy.stats import binned_statistic_2d

def filter_sparse_points(points, x_bins=50, y_bins=50, density_threshold=2):
    """
    根据网格密度过滤稀疏点，只保留聚集区域。

    参数:
        points: ndarray (N, 2)，二维坐标数组。
        x_bins: int，X方向网格数量。
        y_bins: int，Y方向网格数量。
        density_threshold: int，密度阈值，小于等于该值的视为稀疏点。

    返回:
        core_points: ndarray，密度高的点。
        sparse_points: ndarray，密度低的点。
    """
    # 统计每个网格中的点数
    stat, x_edges, y_edges, _ = binned_statistic_2d(
        points[:, 0], points[:, 1], None, statistic='count', bins=[x_bins, y_bins]
    )

    # 找出每个点所在网格的密度
    x_bin_idx = np.digitize(points[:, 0], x_edges) - 1
    y_bin_idx = np.digitize(points[:, 1], y_edges) - 1

    # 避免越界
    valid = (x_bin_idx >= 0) & (x_bin_idx < x_bins) & (y_bin_idx >= 0) & (y_bin_idx < y_bins)
    x_bin_idx = x_bin_idx[valid]
    y_bin_idx = y_bin_idx[valid]
    points_valid = points[valid]

    grid_density = stat[x_bin_idx, y_bin_idx]

    # 按密度阈值筛选
    core_points = points_valid[grid_density > density_threshold]
    sparse_points = points_valid[grid_density <= density_threshold]

    return core_points, sparse_points


from shapely.geometry import Polygon
from shapely.ops import unary_union

def remove_intersections(polygons):
    cleaned_polygons = []
    for i, poly in enumerate(polygons):
        # 组合其他所有 polygon
        others = polygons[:i] + polygons[i+1:]
        others_union = unary_union(others)

        # 去掉交集部分
        diff = poly.difference(others_union)

        # 有可能 diff 仍是多个部分
        if diff.is_empty:
            continue
        elif diff.geom_type == 'Polygon':
            cleaned_polygons.append(diff)
        elif diff.geom_type == 'MultiPolygon':
            cleaned_polygons.extend(diff.geoms)
    return cleaned_polygons



 
def plot_four_panels(
    sparse_points, 
    target_coords, 
    alpha_shape, 
    polys_sorted, 
    final_polys, 
    figsize=(16, 5)
):
    """
    绘制四联图：稀疏点过滤 / Alpha Shape / 平滑多边形 / 筛选多边形

    参数
    ----
    sparse_points : np.ndarray
        稀疏点坐标 (N, 2)
    target_coords : np.ndarray
        目标簇坐标 (M, 2)
    alpha_shape : shapely Polygon or MultiPolygon
        Alpha Shape 结果
    polys_sorted : list of shapely.Polygon
        平滑后的多边形
    final_polys : list of shapely.Polygon
        最终筛选的多边形
    figsize : tuple
        图像大小 (默认 (16,5))
    """

    fig, axs = plt.subplots(1, 4, figsize=figsize)

    # Panel 0: Sparse point filtering
    axs[0].scatter(sparse_points[:, 0], sparse_points[:, 1], s=0.1, c='lightgray', label='Sparse')
    axs[0].scatter(target_coords[:, 0], target_coords[:, 1], s=0.1, c='red', label='Core Cluster')
    axs[0].set_title("Filtered Dense Cluster")
    axs[0].invert_yaxis()
    axs[0].legend()

    # Panel 1: Alpha Shape
    axs[1].scatter(target_coords[:, 0], target_coords[:, 1], color='#41b6c4', s=2, label="Cells")
    if alpha_shape.geom_type == 'Polygon':
        bx, by = alpha_shape.exterior.xy
        axs[1].plot(bx, by, 'r-', label="Alpha Shape")
    elif alpha_shape.geom_type == 'MultiPolygon':
        for geom in alpha_shape.geoms:
            bx, by = geom.exterior.xy
            axs[1].plot(bx, by, 'r-', label="Alpha Shape")
    axs[1].set_title("Alpha Shape")
    axs[1].invert_yaxis()

    # Panel 2: Smoothed Polygons
    for i, poly in enumerate(polys_sorted):
        x, y = poly.exterior.xy
        axs[2].fill(x, y, alpha=0.5)
        axs[2].plot(x, y, color='black', linewidth=1)
        axs[2].text(poly.centroid.x, poly.centroid.y, f"{i+1}", ha='center', va='center', fontsize=8)
    axs[2].set_title(f"Smoothed Polygons\nTotal: {len(polys_sorted)}")
    axs[2].invert_yaxis()

    # Panel 3: Filtered Polygons
    for i, poly in enumerate(final_polys):
        x, y = poly.exterior.xy
        axs[3].fill(x, y, alpha=0.5)
        axs[3].plot(x, y, color='black', linewidth=1)
        axs[3].text(poly.centroid.x, poly.centroid.y, f"{i+1}", ha='center', va='center', fontsize=8)
    axs[3].set_title(f"Filtered Polygons\nTotal: {len(final_polys)}")
    axs[3].invert_yaxis()

    # Common layout
    for ax in axs:
        ax.axis("equal")
        ax.axis("off")

    plt.tight_layout()
    plt.show()
