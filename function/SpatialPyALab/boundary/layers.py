from shapely.geometry import Polygon
from shapely.ops import unary_union

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
                if step > 1: 
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



from shapely.ops import unary_union

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
            raise ValueError("sliding outward: need one of num_steps or max_dist")
        max_dist = num_steps * slide_step  

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
                    layer_clean = layer_raw.difference(others_union.buffer(dist + window_size))
            else:
                layer_clean = layer_raw

            if not layer_clean.is_empty:
                layers.append((layer_clean, float(dist))) 

        results[i] = layers
    return results


def create_sliding_inward_layers_per_poly(
    polys,
    window_size=50,
    slide_step=10,
    min_area=0.1
):
    results = {}
    for i, poly in enumerate(polys):
        layers = []
        dist = 0 
        while True:
            outer = poly.buffer(-dist)
            if outer.is_empty or outer.area < min_area:
                break

            inner = poly.buffer(-(dist + window_size))
            if inner.is_empty or inner.area < min_area:
                if not outer.is_empty:
                    layers.append((outer, -float(dist)))
                break

            ring = outer.difference(inner)
            if not ring.is_empty:
                layers.append((ring, -float(dist))) 

            dist += slide_step

        results[i] = layers
    return results


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
            raise ValueError("sliding outward global: need num_steps or max_dist")
        max_dist = num_steps * slide_step

    layers = []
    for dist in range(0, max_dist + 1, slide_step):
        outer = union_poly.buffer(dist + window_size)
        inner = union_poly.buffer(dist)
        ring = outer.difference(inner)
        if not ring.is_empty:
            layers.append((ring, float(dist)))
    return {-1: layers}


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
    sliding_mode=False,
    window_size=None,
    slide_step=None, 
    max_dist=None  
):
    outward_layers = {}
    inward_layers  = {}

    if sliding_mode:
        if window_size is None or slide_step is None:
            raise ValueError("sliding_mode=True need window_size and slide_step")

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

    if expansion_mode == "both":
        all_keys = set(outward_layers.keys()) | set(inward_layers.keys())
        combined = {}
        for k in all_keys:
            combined[k] = []
            inward  = inward_layers.get(k, [])
            outward = outward_layers.get(k, [])
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
        for j, poly in enumerate(polys):
            if j != idx:
                poly_to_draw = poly.intersection(boundary) if boundary else poly
                if isinstance(poly_to_draw, (Polygon, MultiPolygon)):
                    _draw_polygon_edges(poly_to_draw, ax, edgecolor='lightgrey', linewidth=edge_width)

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

    for i in range(n, len(axes)):
        fig.delaxes(axes[i])

    if suptitle:
        fig.suptitle(suptitle, fontsize=20)

    plt.tight_layout()
    plt.show()

def _draw_polygon_edges(poly, ax, edgecolor='black', linewidth=1, alpha=1):

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

        for j, poly in enumerate(polys):
            if j != idx:
                poly_to_draw = poly.intersection(boundary) if boundary else poly
                _draw_polygon_edges(poly_to_draw, ax, edgecolor='lightgrey', linewidth=edge_width)

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
 
 
import numpy as np
import pandas as pd
from shapely.prepared import prep
from shapely.geometry import Point
from rtree import index

def assign_points_nonoverlap_fast_from_obsm(
    adata,
    layers_dict,                    
    obsm_key="spatial_um",         
    cell_id_col=None,               
    region_name_prefix="region",
    swap_xy=False,              
):
    
    if obsm_key not in adata.obsm:
        raise KeyError(f"{obsm_key} is not in adata.obsm")
    arr = np.asarray(adata.obsm[obsm_key], dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"{obsm_key} shape must be (N,2)，but now {arr.shape}")
    if swap_xy:
        arr = arr[:, [1, 0]]

    valid = np.isfinite(arr).all(axis=1)
    if not valid.any():
        raise ValueError(f"{obsm_key} haven't coor to used")
    arr = arr[valid]
    obs_idx_kept = np.flatnonzero(valid)
    pt_idx = index.Index()
    for j, (x, y) in enumerate(arr):
        pt_idx.insert(j, (float(x), float(y), float(x), float(y)))
    points = [Point(float(x), float(y)) for x, y in arr]

    cell_ids = (
        adata.obs_names.to_numpy()[obs_idx_kept]
        if cell_id_col is None else
        adata.obs[cell_id_col].to_numpy()[obs_idx_kept]
    )
    region_keys = [f"{region_name_prefix}_{rid}" for rid in layers_dict.keys()]
    layer_mat = pd.DataFrame(index=cell_ids, columns=region_keys, dtype='float32')

    long_records = []

    for rid, items in layers_dict.items():
        rkey = f"{region_name_prefix}_{rid}"
        assigned = np.full(len(points), np.nan, dtype='float32')

        for geom, layer_num in items:
            if geom.is_empty:
                continue
            pg = prep(geom)  
            minx, miny, maxx, maxy = geom.bounds
            cand_local = list(pt_idx.intersection((minx, miny, maxx, maxy)))
            if not cand_local:
                continue

            for lid in cand_local:
                if not np.isnan(assigned[lid]): 
                    continue
                if pg.covers(points[lid]):
                    assigned[lid] = layer_num

        hit_mask = ~np.isnan(assigned)
        if hit_mask.any():
            hit_cell_ids = cell_ids[hit_mask]
            layer_mat.loc[hit_cell_ids, rkey] = assigned[hit_mask]
            for cid, ln in zip(hit_cell_ids, assigned[hit_mask]):
                long_records.append({"cell_id": cid, "region_key": rkey, "layer": int(ln)})

    long_df = pd.DataFrame(long_records)
    return layer_mat, long_df
