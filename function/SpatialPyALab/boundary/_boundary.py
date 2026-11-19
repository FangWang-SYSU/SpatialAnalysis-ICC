
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

    polygon_list = []
    if isinstance(alpha_shape, Polygon):
        polygon_list.append(alpha_shape)
    elif isinstance(alpha_shape, MultiPolygon):
        polygon_list.extend(alpha_shape.geoms)
    else:
        for geom in alpha_shape.geoms:
            if isinstance(geom, Polygon):
                polygon_list.append(geom)
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

    polygons = split_geoms_to_polygons(alpha_shape)
    smoothed_polygons = [smooth_polygon(poly, smooth_factor) for poly in polygons]

    buffered_polygons = [poly.buffer(buffer_distance) for poly in smoothed_polygons]
    merged_polygon = unary_union(buffered_polygons)
    if not merged_polygon.is_valid:
        merged_polygon = merged_polygon.buffer(0)

    return merged_polygon

from shapely.geometry import Point, Polygon, MultiPolygon

def filter_polygons_by_cell_count(cell_coords, threshold, polygons):
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
    stat, x_edges, y_edges, _ = binned_statistic_2d(
        points[:, 0], points[:, 1], None, statistic='count', bins=[x_bins, y_bins]
    )

    x_bin_idx = np.digitize(points[:, 0], x_edges) - 1
    y_bin_idx = np.digitize(points[:, 1], y_edges) - 1

    valid = (x_bin_idx >= 0) & (x_bin_idx < x_bins) & (y_bin_idx >= 0) & (y_bin_idx < y_bins)
    x_bin_idx = x_bin_idx[valid]
    y_bin_idx = y_bin_idx[valid]
    points_valid = points[valid]

    grid_density = stat[x_bin_idx, y_bin_idx]

    core_points = points_valid[grid_density > density_threshold]
    sparse_points = points_valid[grid_density <= density_threshold]

    return core_points, sparse_points


from shapely.geometry import Polygon
from shapely.ops import unary_union

def remove_intersections(polygons):
    cleaned_polygons = []
    for i, poly in enumerate(polygons):
        others = polygons[:i] + polygons[i+1:]
        others_union = unary_union(others)
        diff = poly.difference(others_union)
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
