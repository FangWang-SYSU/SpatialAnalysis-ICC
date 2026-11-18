from shapely.geometry import Polygon
from shapely.ops import unary_union

# 每个 poly 向外扩展
def create_outer_layers_per_poly(polys, buffer_step, num_steps, consider_others_expansion=False, safe_margin=0.0):
    results = {}

    for i, poly in enumerate(polys):
        previous = poly
        layers = []

        for step in range(1, num_steps + 1):
            distance = step * buffer_step
            this_buffer = poly.buffer(distance)

            layer_raw = this_buffer.difference(previous)

            others = polys[:i] + polys[i + 1:]
            
            if not consider_others_expansion:
                others_union = unary_union(others).buffer(safe_margin)
                layer_clean = layer_raw.difference(others_union)

            else:
                others_union = unary_union([p.buffer(distance) for p in others])
                layer_clean = layer_raw.difference(others_union)
                
            if not layer_clean.is_empty:
                layers.append((layer_clean, step))

            previous = this_buffer

        results[i] = layers

    return results

# 全局向外扩展：将所有 poly 合并处理，仍返回 dict 格式，key 用 -1 表示整体
def create_outer_layers_global(polys, buffer_step, num_steps):
    union_poly = unary_union(polys)
    buffered_layers = [union_poly]

    for step in range(1, num_steps + 1):
        distance = step * buffer_step
        buffered = union_poly.buffer(distance)
        buffered_layers.append(buffered)

    layers = []
    for i in range(len(buffered_layers) - 1):
        shell = buffered_layers[i + 1].difference(buffered_layers[i])
        if not shell.is_empty:
            layers.append(shell)

    return {-1: layers}


# 每个 poly 向内收缩
def create_inward_layers_per_poly(polys, buffer_step=10, min_area=0.1):
    results = {}

    for i, poly in enumerate(polys):
        inward_layers = []
        previous = poly
        step = 1

        while True:
            distance = -step * buffer_step
            shrunken = poly.buffer(distance)

            if shrunken.is_empty or shrunken.area < min_area:
                if step > 1:  # 第一层
                    if inward_layers:
                            merged = previous
                            inward_layers[-1] = (inward_layers[-1][0].union(merged), inward_layers[-1][1])
                else:
                    inward_layers.append((poly, -1))
                break

            layer = previous.difference(shrunken)
            if not layer.is_empty:
                inward_layers.append((layer,-step))

            previous = shrunken
            step += 1

        results[i] = inward_layers

    return results



# 全局向内收缩：合并后整体收缩
def create_inward_layers_global(polys, buffer_step, min_area):
    union_poly = unary_union(polys)
    buffer_step_inward = -abs(buffer_step)
    buffered_layers = [union_poly]
    step = 1

    while True:
        distance = step * buffer_step_inward
        try:
            shrunken = union_poly.buffer(distance)
            if shrunken.area > min_area:
                buffered_layers.append(shrunken)
            else:
                break
        except:
            break
        step += 1

    layers = []
    for i in range(len(buffered_layers) - 1):
        diff = buffered_layers[i].difference(buffered_layers[i + 1])
        if not diff.is_empty:
            layers.append(diff)

    return {-1: layers[::-1]}

"""
# 主函数：支持 per_poly / global 两种模式，统一输出格式
def create_layers(
    polys,
    operation_mode="per_poly",     # "per_poly" or "global"
    expansion_mode="both",         # "outward", "inward", "both"
    buffer_step=10,
    num_steps=5,
    min_area=0.1,
    consider_others_expansion=False,
    safe_margin=0.0,
    include_center = True
):
    outward_layers = {}
    inward_layers = {}


    if operation_mode == "per_poly":
        if expansion_mode in ("outward", "both"):
            outward_layers = create_outer_layers_per_poly(polys, buffer_step, num_steps, consider_others_expansion, safe_margin)
        
        if expansion_mode in ("inward", "both"):
            inward_layers = create_inward_layers_per_poly(polys, buffer_step, min_area)

    if expansion_mode == "both":
        all_keys = set(outward_layers.keys()) | set(inward_layers.keys())
        combined = {}
        for k in all_keys:
            combined[k] = []
            inward = inward_layers.get(k, [])
            outward = outward_layers.get(k, [])

            inward = sorted(inward, key=lambda x: x[1])
            outward = sorted(outward, key=lambda x: x[1])
 
            combined[k].extend(inward)
            if include_center:
                combined[k].append((polys[k], 0))

            combined[k].extend(outward)
        return combined
    elif expansion_mode == "outward":
        return outward_layers
    elif expansion_mode == "inward":
        return inward_layers
    
"""


from shapely.ops import unary_union

# === 你已有的函数保持不动（outer/global/inward非滑动版） ===
# ...（此处省略你前面的 create_outer_layers_per_poly / create_outer_layers_global /
#      create_inward_layers_per_poly / create_inward_layers_global 定义）...

# ====== 新增：滑动外扩（per_poly） ======
def create_sliding_outer_layers_per_poly(
    polys,
    window_size=50,
    slide_step=10,
    num_steps=None,
    max_dist=None,
    consider_others_expansion=False,
    safe_margin=0.0
):
    if max_dist is None:
        if num_steps is None:
            raise ValueError("sliding outward: 需要 num_steps 或 max_dist 之一")
        max_dist = num_steps * slide_step  # 最后一窗起点

    results = {}
    for i, poly in enumerate(polys):
        layers = []
        others = polys[:i] + polys[i+1:]
        others_union = unary_union(others) if others else None

        for dist in range(0, max_dist + 1, slide_step):
            outer = poly.buffer(dist + window_size)
            inner = poly.buffer(dist)
            layer_raw = outer.difference(inner)

            if others_union:
                if not consider_others_expansion:
                    layer_clean = layer_raw.difference(others_union.buffer(safe_margin))
                else:
                    # buffer(union) 近似等价 union(buffer)（正向buffer时成立），更快
                    layer_clean = layer_raw.difference(others_union.buffer(dist + window_size))
            else:
                layer_clean = layer_raw

            if not layer_clean.is_empty:
                layers.append((layer_clean, float(dist)))  # 向外为正，单位=坐标单位（通常像素或µm）

        results[i] = layers
    return results


# ====== 新增：滑动内收（per_poly） ======
def create_sliding_inward_layers_per_poly(
    polys,
    window_size=50,
    slide_step=10,
    min_area=0.1
):
    results = {}
    for i, poly in enumerate(polys):
        layers = []
        dist = 0  # 从边界向内推进的起点
        while True:
            outer = poly.buffer(-dist)
            if outer.is_empty or outer.area < min_area:
                break

            inner = poly.buffer(-(dist + window_size))
            if inner.is_empty or inner.area < min_area:
                # 最后一层：剩余核心整体作为最终层（厚度可能 < window_size）
                if not outer.is_empty:
                    layers.append((outer, -float(dist)))
                break

            ring = outer.difference(inner)
            if not ring.is_empty:
                layers.append((ring, -float(dist)))  # 向内为负

            dist += slide_step

        results[i] = layers
    return results


# ====== 新增：滑动外扩（global） ======
def create_sliding_outer_layers_global(
    polys,
    window_size=50,
    slide_step=10,
    num_steps=None,
    max_dist=None
):
    union_poly = unary_union(polys)
    if max_dist is None:
        if num_steps is None:
            raise ValueError("sliding outward global: 需要 num_steps 或 max_dist 之一")
        max_dist = num_steps * slide_step

    layers = []
    for dist in range(0, max_dist + 1, slide_step):
        outer = union_poly.buffer(dist + window_size)
        inner = union_poly.buffer(dist)
        ring = outer.difference(inner)
        if not ring.is_empty:
            layers.append((ring, float(dist)))
    return {-1: layers}


# ====== 新增：滑动内收（global） ======
def create_sliding_inward_layers_global(
    polys,
    window_size=50,
    slide_step=10,
    min_area=0.1
):
    union_poly = unary_union(polys)
    layers = []
    dist = 0
    while True:
        outer = union_poly.buffer(-dist)
        if outer.is_empty or outer.area < min_area:
            break

        inner = union_poly.buffer(-(dist + window_size))
        if inner.is_empty or inner.area < min_area:
            if not outer.is_empty:
                layers.append((outer, -float(dist)))
            break

        ring = outer.difference(inner)
        if not ring.is_empty:
            layers.append((ring, -float(dist)))

        dist += slide_step
    return {-1: layers}


# ====== 主入口：支持滑动模式 ======
def create_layers(
    polys,
    operation_mode="per_poly",     # "per_poly" or "global"
    expansion_mode="both",         # "outward", "inward", "both"
    buffer_step=10,
    num_steps=5,
    min_area=0.1,
    consider_others_expansion=False,
    safe_margin=0.0,
    include_center=True,
    # 新增：滑动模式参数
    sliding_mode=False,
    window_size=None,              # 必填（滑动模式）
    slide_step=None,               # 必填（滑动模式）
    max_dist=None                  # 可选（滑动外扩）
):
    outward_layers = {}
    inward_layers  = {}

    if sliding_mode:
        if window_size is None or slide_step is None:
            raise ValueError("sliding_mode=True 时须提供 window_size 与 slide_step")

        if operation_mode == "per_poly":
            if expansion_mode in ("outward", "both"):
                outward_layers = create_sliding_outer_layers_per_poly(
                    polys,
                    window_size=window_size,
                    slide_step=slide_step,
                    num_steps=num_steps,
                    max_dist=max_dist,
                    consider_others_expansion=consider_others_expansion,
                    safe_margin=safe_margin
                )
            if expansion_mode in ("inward", "both"):
                inward_layers = create_sliding_inward_layers_per_poly(
                    polys,
                    window_size=window_size,
                    slide_step=slide_step,
                    min_area=min_area
                )

        elif operation_mode == "global":
            if expansion_mode in ("outward", "both"):
                outward_layers = create_sliding_outer_layers_global(
                    polys,
                    window_size=window_size,
                    slide_step=slide_step,
                    num_steps=num_steps,
                    max_dist=max_dist
                )
            if expansion_mode in ("inward", "both"):
                inward_layers = create_sliding_inward_layers_global(
                    polys,
                    window_size=window_size,
                    slide_step=slide_step,
                    min_area=min_area
                )

    else:
        # === 保持你原来的非滑动逻辑 ===
        if operation_mode == "per_poly":
            if expansion_mode in ("outward", "both"):
                outward_layers = create_outer_layers_per_poly(
                    polys, buffer_step, num_steps,
                    consider_others_expansion, safe_margin
                )
            if expansion_mode in ("inward", "both"):
                inward_layers = create_inward_layers_per_poly(
                    polys, buffer_step=buffer_step, min_area=min_area
                )
        elif operation_mode == "global":
            if expansion_mode in ("outward", "both"):
                outward_layers = create_outer_layers_global(polys, buffer_step, num_steps)
            if expansion_mode in ("inward", "both"):
                inward_layers = create_inward_layers_global(polys, buffer_step, min_area)

    # === 汇总 ===
    if expansion_mode == "both":
        all_keys = set(outward_layers.keys()) | set(inward_layers.keys())
        combined = {}
        for k in all_keys:
            combined[k] = []
            inward  = inward_layers.get(k, [])
            outward = outward_layers.get(k, [])

            # 统一排序：负值(向内) -> 0 -> 正值(向外)
            inward  = sorted(inward,  key=lambda x: x[1])
            outward = sorted(outward, key=lambda x: x[1])

            combined[k].extend(inward)
            if include_center and operation_mode == "per_poly" and isinstance(k, int) and k >= 0:
                combined[k].append((polys[k], 0))
            combined[k].extend(outward)
        return combined
    elif expansion_mode == "outward":
        return outward_layers
    elif expansion_mode == "inward":
        return inward_layers

 
 



def plot_polys_with_layers_grid(polys, 
                                layers_dict, 
                                cols=5, 
                                figsize=(20, 15), 
                                suptitle=None, 
                                edge_width=1,
                                boundary=None,
                                layer_colors=None):
    
    import matplotlib.pyplot as plt
    import math
    from shapely.geometry import Polygon, MultiPolygon
    
    if layer_colors is None:
        layer_colors = [
            '#8dd3c7', '#ffffb3', '#bebada', '#fb8072',
            '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9'
        ]

    n = len(polys)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    axes = axes.flatten()

    for idx in range(n):
        ax = axes[idx]

        # 画其他多边形（灰色）
        for j, poly in enumerate(polys):
            if j != idx:
                poly_to_draw = poly.intersection(boundary) if boundary else poly
                if isinstance(poly_to_draw, (Polygon, MultiPolygon)):
                    _draw_polygon_edges(poly_to_draw, ax, edgecolor='lightgrey', linewidth=edge_width)

        # 画当前 polygon 的各层（带颜色）
        if idx in layers_dict:
            for i, layer in enumerate(layers_dict[idx]):
                color = layer_colors[i % len(layer_colors)]
                layer_to_draw = layer.intersection(boundary) if boundary else layer
                _fill_polygon(layer_to_draw, ax, facecolor=color, edgecolor='black',
                              linewidth=edge_width, alpha=0.6)

        ax.set_title(f"Poly {idx}")
        ax.set_aspect('equal')
        ax.axis('off')
        ax.invert_yaxis()

    # 删除空白子图
    for i in range(n, len(axes)):
        fig.delaxes(axes[i])

    if suptitle:
        fig.suptitle(suptitle, fontsize=20)

    plt.tight_layout()
    plt.show()

def _draw_polygon_edges(poly, ax, edgecolor='black', linewidth=1, alpha=1):
    """画 Polygon 或 MultiPolygon 的边界线"""
    if isinstance(poly, MultiPolygon):
        for p in poly.geoms:
            _draw_single_polygon(p, ax, edgecolor, linewidth, alpha)
    else:
        _draw_single_polygon(poly, ax, edgecolor, linewidth, alpha)

def _draw_single_polygon(p, ax, edgecolor, linewidth, alpha):
    xs, ys = p.exterior.xy
    ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha)
    for interior in p.interiors:
        xs, ys = interior.xy
        ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha)

def _fill_polygon(poly, ax, facecolor, edgecolor, linewidth, alpha):
    """填充 Polygon 或 MultiPolygon"""
    if isinstance(poly, MultiPolygon):
        for p in poly.geoms:
            ax.fill(*p.exterior.xy, facecolor=facecolor, edgecolor=edgecolor,
                    linewidth=linewidth, alpha=alpha)
    else:
        ax.fill(*poly.exterior.xy, facecolor=facecolor, edgecolor=edgecolor,
                linewidth=linewidth, alpha=alpha)



import math
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon

DEFAULT_LAYER_COLORS = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999']

def plot_polys_with_layers_grid_lineonly(polys, layers_dict, cols=5, figsize=(20, 15),
                                          layer_colors=None, edge_width=1, boundary=None,
                                          suptitle=None):
    if layer_colors is None:
        layer_colors = DEFAULT_LAYER_COLORS

    n = len(polys)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    axes = axes.flatten()

    for idx in range(n):
        ax = axes[idx]

        # 画其他 polygon 边界（灰色）
        for j, poly in enumerate(polys):
            if j != idx:
                poly_to_draw = poly.intersection(boundary) if boundary else poly
                _draw_polygon_edges(poly_to_draw, ax, edgecolor='lightgrey', linewidth=edge_width)

        # 画当前 polygon 的每一层
        if idx in layers_dict:
            layers = layers_dict[idx]
            for i, layer in enumerate(layers):
                color = layer_colors[i % len(layer_colors)]
                layer_to_draw = layer.intersection(boundary) if boundary else layer
                _draw_polygon_edges(layer_to_draw, ax, edgecolor=color, linewidth=edge_width)

        ax.set_title(f"Poly {idx}")
        ax.set_aspect('equal')
        ax.axis('off')
        ax.invert_yaxis()

    # 删除多余 subplot
    for i in range(n, len(axes)):
        fig.delaxes(axes[i])

    if suptitle:
        fig.suptitle(suptitle, fontsize=20)

    plt.tight_layout()
    plt.show()

def _draw_polygon_edges(poly, ax, edgecolor='black', linewidth=1, alpha=1):
    if poly is None or poly.is_empty:
        return
    if isinstance(poly, MultiPolygon):
        for p in poly.geoms:
            _draw_single_polygon(p, ax, edgecolor, linewidth, alpha)
    else:
        _draw_single_polygon(poly, ax, edgecolor, linewidth, alpha)

def _draw_single_polygon(p, ax, edgecolor, linewidth, alpha):
    if not p.is_valid or p.is_empty:
        return
    xs, ys = p.exterior.xy
    ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha)
    for interior in p.interiors:
        xs, ys = interior.xy
        ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha)
 
# from shapely.geometry import Point
# import pandas as pd
# from rtree import index

# def assign_points_to_layers_optimized(coords, combined_layers_with_index):
#     """
#     使用空间索引优化，将每个点分配到对应 polygon 的层（内外层）。
    
#     Parameters:
#         coords: numpy.ndarray, shape=(n_points, 2)
#             每个点的空间坐标。
#         combined_layers_with_index: dict
#             每个多边形的层列表，每层是 (Polygon, 层编号)。
#             内层编号为负，外层为正，建议从内到外排序。

#     Returns:
#         pd.DataFrame
#             index: 点的索引
#             columns: region_0, region_1, ...
#             values: 该点在该区域的层编号，若无则为 None
#     """
#     n_points = coords.shape[0]
#     region_ids = [f"region_{k}" for k in combined_layers_with_index.keys()]
#     layer_matrix = pd.DataFrame(index=range(n_points), columns=region_ids)

#     # 为每个 region 构建独立的空间索引
#     spatial_indexes = {}
#     for region_id, layers in combined_layers_with_index.items():
#         idx = index.Index()
#         for i, (poly, layer_num) in enumerate(layers):
#             idx.insert(i, poly.bounds)
#         spatial_indexes[region_id] = idx

#     # 分配点到每个 region 的层
#     for point_idx, (x, y) in enumerate(coords):
#         pt = Point(x, y)

#         for region_key in combined_layers_with_index:
#             col_name = f"region_{region_key}"
#             layers = combined_layers_with_index[region_key]
#             idx = spatial_indexes[region_key]

#             # 先按层编号排序，负数优先，从内到外（绝对值小的在前）
#             sorted_layers = sorted(layers, key=lambda x: (x[1] >= 0, abs(x[1])))

#             for i, (poly, layer_num) in enumerate(sorted_layers):
#                 # 仅考虑可能的候选 polygon
#                 if any(idx.intersection(poly.bounds)):
#                     if poly.contains(pt):
#                         layer_matrix.at[point_idx, col_name] = layer_num
#                         break 

#     return layer_matrix


import numpy as np
import pandas as pd
from shapely.prepared import prep
from shapely.geometry import Point
from rtree import index

def assign_points_nonoverlap_fast_from_obsm(
    adata,
    layers_dict,                     # {region_id: [(geom, layer_num), ...]} ——各层在同一 region 内不重叠
    obsm_key="spatial_um",           # 坐标在 obsm；若是像素请自行做单位换算
    cell_id_col=None,                # None=用 obs_names；或传你自定义的 'cell_id'
    region_name_prefix="region",
    swap_xy=False,                   # 若 obsm 是 (row,col) = (y,x)，设 True
):
    # 1) 取坐标 & 建点索引（R-tree）
    if obsm_key not in adata.obsm:
        raise KeyError(f"{obsm_key} 不在 adata.obsm 中")
    arr = np.asarray(adata.obsm[obsm_key], dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"{obsm_key} 形状应为 (N,2)，当前 {arr.shape}")
    if swap_xy:
        arr = arr[:, [1, 0]]

    # 过滤非法坐标
    valid = np.isfinite(arr).all(axis=1)
    if not valid.any():
        raise ValueError(f"{obsm_key} 没有可用坐标（全是 NaN/Inf？）")
    arr = arr[valid]
    obs_idx_kept = np.flatnonzero(valid)

    # R-tree（点）
    pt_idx = index.Index()
    for j, (x, y) in enumerate(arr):
        pt_idx.insert(j, (float(x), float(y), float(x), float(y)))
    points = [Point(float(x), float(y)) for x, y in arr]

    # 2) 准备输出容器
    cell_ids = (
        adata.obs_names.to_numpy()[obs_idx_kept]
        if cell_id_col is None else
        adata.obs[cell_id_col].to_numpy()[obs_idx_kept]
    )
    region_keys = [f"{region_name_prefix}_{rid}" for rid in layers_dict.keys()]
    layer_mat = pd.DataFrame(index=cell_ids, columns=region_keys, dtype='float32')

    long_records = []

    # 3) 逐 region / 逐层：查候选点 → covers 精筛 → 赋值
    for rid, items in layers_dict.items():
        rkey = f"{region_name_prefix}_{rid}"
        # 每个 region 内唯一归属：先开一个数组存层号，默认 NaN
        assigned = np.full(len(points), np.nan, dtype='float32')

        for geom, layer_num in items:
            if geom.is_empty:
                continue
            pg = prep(geom)  # prepared geometry 加速重复 covers
            minx, miny, maxx, maxy = geom.bounds
            cand_local = list(pt_idx.intersection((minx, miny, maxx, maxy)))
            if not cand_local:
                continue

            for lid in cand_local:
                if not np.isnan(assigned[lid]):     # 已经被本 region 其它层命中（虽然不重叠，一般不会发生）
                    continue
                if pg.covers(points[lid]):
                    assigned[lid] = layer_num

        # 写入 layer_matrix（仅对命中者赋值，其余保持 NaN）
        hit_mask = ~np.isnan(assigned)
        if hit_mask.any():
            hit_cell_ids = cell_ids[hit_mask]
            layer_mat.loc[hit_cell_ids, rkey] = assigned[hit_mask]

            # 生成 long_df 记录
            for cid, ln in zip(hit_cell_ids, assigned[hit_mask]):
                long_records.append({"cell_id": cid, "region_key": rkey, "layer": int(ln)})

    long_df = pd.DataFrame(long_records)
    return layer_mat, long_df
