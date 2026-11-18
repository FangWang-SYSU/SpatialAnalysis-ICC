import numpy as np
import pandas as pd
from shapely.geometry import Point
from rtree import index

# 1) 展平窗口（层）
def flatten_layers(layers_dict):
    rows = []
    for rid, items in layers_dict.items():
        for w_idx, (geom, offset) in enumerate(items):
            side = "inward" if offset < 0 else ("center" if offset == 0 else "outward")
            rows.append({"region_id": rid, "window_id": w_idx, "offset": float(offset),
                         "side": side, "geom": geom})
    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values(["region_id", "offset", "window_id"], inplace=True)
    return df

# 2) 从 obsm 取点并建索引（鲁棒版，返回保留下来的原始行号）
def make_points_from_obsm(adata, obsm_key="spatial_um", unit_scale=1.0, swap_xy=False):
    if obsm_key not in adata.obsm:
        raise KeyError(f"{obsm_key} 不在 adata.obsm 中")

    arr = np.asarray(adata.obsm[obsm_key], dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"{obsm_key} 形状应为 (N,2)，现在是 {arr.shape}")

    xy = arr.copy()
    if swap_xy:
        xy = xy[:, [1, 0]]

    # 过滤非法坐标
    valid_mask = np.isfinite(xy).all(axis=1)
    kept_idx = np.flatnonzero(valid_mask)
    xy = xy[valid_mask]

    if unit_scale != 1.0:
        xy = xy * unit_scale

    # 0 点处理：不要用流式构建，空索引也能创建成功
    pt_idx = index.Index()
    points = []
    for j, (x, y) in enumerate(xy):
        # bounds 顺序: (minx, miny, maxx, maxy)；点用 (x, y, x, y)
        pt_idx.insert(j, (float(x), float(y), float(x), float(y)))
        points.append(Point(float(x), float(y)))

    return points, pt_idx, kept_idx  # kept_idx 是 obs 的原始行号

# 3) 核心：按滑动窗口统计（点可多重归属），带原索引映射
def proportions_by_sliding_windows(
    adata,
    layers_dict,
    obsm_key="spatial_um",
    type_key="CellType_ficture",
    unit_scale=1.0,
    swap_xy=False,
    min_cells=1
):
    layer_df = flatten_layers(layers_dict)
    points, pt_idx, kept_idx = make_points_from_obsm(adata, obsm_key, unit_scale, swap_xy)

    if len(points) == 0:
        raise ValueError(
            f"{obsm_key} 中没有可用坐标（可能被 NaN/Inf 过滤掉了或对象为空）。"
            "请检查：过滤前后细胞数量、坐标是否为数值、是否需要 swap_xy。"
        )

    # 只在保留下来的行上取 cell_type，保证索引对齐
    cell_types_all = adata.obs[type_key].astype("category")
    
    rows = []
    for _, row in layer_df.iterrows():
        geom = row["geom"]
        minx, miny, maxx, maxy = geom.bounds

        # 候选点（在窗口 bbox 内）
        cand_local_ids = list(pt_idx.intersection((minx, miny, maxx, maxy)))
        if not cand_local_ids:
            continue

        # 精筛（covers 包含边界）
        hit_local_ids = [lid for lid in cand_local_ids if geom.covers(points[lid])]
        if len(hit_local_ids) < min_cells:
            continue

        # 映射回 obs 原始行号
        hit_obs_ids = kept_idx[hit_local_ids]
        counts = cell_types_all.iloc[hit_obs_ids].value_counts()
        total  = int(counts.sum())
        for ct, n in counts.items():
            rows.append({
                "region_id": row["region_id"],
                "window_id": row["window_id"],
                "offset": row["offset"],
                "side": row["side"],
                "cell_type": ct,
                "n": int(n),
                "total": total,
                "prop": float(n / total)
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values(["region_id", "offset", "cell_type"], inplace=True)
    return df

 
import numpy as np
import pandas as pd


def aggregate_by_clone(df_prop, region_to_clone, region_name_prefix="region", unknown="Unknown"):
    df = df_prop.copy()

    # ---- 1) 补齐 region_key（字符串），兼容多种输入 ----
    if "region_key" not in df.columns:
        if "region" in df.columns:
            # 你旧版 long_df 风格
            df["region_key"] = df["region"].astype(str)
        elif "region_id" in df.columns:
            df["region_key"] = df["region_id"].map(lambda x: f"{region_name_prefix}_{int(x)}")
        else:
            raise KeyError("df_prop 中既无 'region_key' 也无 'region' / 'region_id'，无法映射到 clone。")

    # 基本字段检查
    needed = {"offset", "side", "cell_type", "n", "total"}
    missing = needed - set(df.columns)
    if missing:
        raise KeyError(f"df_prop 缺少必要列: {missing}")

    # ---- 3) region -> clone 映射 ----
    df["clone"] = df["region_key"].map(region_to_clone)
    # 兼容映射字典里用的是整数 region_id 的情况
    if df["clone"].isna().any() and "region_id" in df.columns:
        df.loc[df["clone"].isna(), "clone"] = df.loc[df["clone"].isna(), "region_id"].map(region_to_clone)
    df["clone"] = df["clone"].fillna(unknown)

    # ---- 4) 加权汇总（Σn / Σtotal），避免简单平均带来的偏差 ----
    grp = (
        df.groupby(["clone", "side", "offset", "cell_type"], as_index=False)[["n", "total"]]
          .sum()
    )
    grp["prop"] = grp["n"] / grp["total"].where(grp["total"] > 0, np.nan)
    grp.sort_values(["clone", "cell_type"], inplace=True)

    return grp

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import ceil

def plot_facets_offset(
    df: pd.DataFrame,
    facet_by: str,  
    color_by:str,                  # "clone" 或 "cell_type"
    clones: list[str] | None = None,  # 可选：只画这些克隆
    cell_types: list[str] | None = None,  # 可选：只画这些细胞类型
    ncols: int = 3,
    sharey: bool = False,
    palette: dict | None = None,       # 针对“颜色维度”的 {name: color}
    figsize = (14,8),
    bbox_to_anchor = (1.1,0.5)
):

    dd = df.copy()
    if clones is not None:
        dd = dd[dd["clone"].isin(clones)]
    if cell_types is not None:
        dd = dd[dd["cell_type"].isin(cell_types)]
    if dd.empty:
        raise ValueError("筛选后数据为空，请检查 clones / cell_types。")

    # ---- 分面维度 & 颜色维度 ----
    facet_vals = dd[facet_by].dropna().unique().tolist()
    color_vals = dd[color_by].dropna().unique().tolist()

    # ---- 调色板：对“颜色维度”固定颜色 ----
    if palette is None:
        base = plt.get_cmap("tab20").colors
        palette = {str(k): base[i % len(base)] for i, k in enumerate(sorted(map(str, color_vals)))}

    # ---- 排序稳定 ----
    dd = dd.sort_values([facet_by, color_by, "offset"]).reset_index(drop=True)

    # ---- 画布布局 ----
    n_panels = len(facet_vals)
    ncols = max(1, min(ncols, n_panels))
    nrows = ceil(n_panels / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,figsize=figsize,sharex=True, sharey=sharey)

    if n_panels == 1:
        axes = np.array([axes])
    axes = axes.reshape(nrows, ncols)

    # ---- 作图 ----
    handles_ref, labels_ref = [], []
    for i, facet_val in enumerate(facet_vals):
        ax = axes[i // ncols, i % ncols]
        sub = dd[dd[facet_by] == facet_val]

        for name in color_vals:
            ss = sub[sub[color_by] == name]
            if ss.empty:
                continue
            x = ss["offset"].to_numpy()
            y = ss["prop"].to_numpy()

            line, = ax.plot(
                x, y, marker="o", linewidth=1.2, markersize=2,
                color=palette[str(name)], label=str(name)
            )
            if str(name) not in labels_ref:
                handles_ref.append(line)
                labels_ref.append(str(name))

        ax.axvline(0, linestyle="--", linewidth=1)
        ax.set_title(f"{facet_val}")
        ax.set_xlabel("Distance from Clone boundary (µm)")
        ax.grid(True, alpha=0.3)


    fig.supylabel("Cell proportion")

    # 删空轴
    for j in range(n_panels, nrows * ncols):
        fig.delaxes(axes[j // ncols, j % ncols])

    # 全局图例
    if handles_ref and labels_ref:
        legend_title = "Cell type" if color_by == "cell_type" else "Clone"
        fig.legend(handles_ref, labels_ref, title=legend_title, 
                   loc="right",bbox_to_anchor=bbox_to_anchor)

    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.show()
    return fig

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from math import ceil

# def _ensure_center_offset(df: pd.DataFrame, window_size: float | None) -> pd.DataFrame:
#     df = df.copy()
#     if "center_offset" not in df.columns:
#         if window_size is None:
#             raise ValueError("df_clone 缺少 center_offset，请传入 window_size 用于自动计算。")
#         df["center_offset"] = df["offset"] + (window_size / 2.0)
#         inward = df["side"].astype(str).str.lower().eq("inward")
#         df.loc[inward, "center_offset"] = df.loc[inward, "offset"] - (window_size / 2.0)
#     return df

# def plot_celltype_facets_by_clone(
#     df_clone: pd.DataFrame,
#     window_size: float | None = None,
#     cell_types: list[str] | None = None,
#     top_k: int = 8,
#     ncols: int = 3,
#     smooth_win: int = 1,
#     sharey: bool = False,
#     dpi: int = 150,
# ):
#     # 基础检查
#     req = {"clone", "cell_type", "offset", "side"}
#     if not req.issubset(df_clone.columns):
#         missing = req - set(df_clone.columns)
#         raise KeyError(f"df_clone 缺少必要列: {missing}")

#     df = df_clone.copy()
#     if "prop" not in df.columns:
#         if {"n", "total"}.issubset(df.columns):
#             df["prop"] = df["n"] / df["total"].replace({0: np.nan})
#         else:
#             raise KeyError("df_clone 既无 'prop' 也无 'n'/'total'，无法作图。")

#     df = _ensure_center_offset(df, window_size)
#     df = df.sort_values(["cell_type", "clone", "center_offset"]).reset_index(drop=True)

#     # 选择 cell types
#     if cell_types is None:
#         cell_types = df["cell_type"].drop_duplicates().tolist()
#     else:
#         cell_types = [ct for ct in cell_types if ct in df["cell_type"].unique()]

#     n_panels = len(cell_types)
#     ncols = max(1, min(ncols, n_panels))
#     nrows = ceil(n_panels / ncols)


#     fig, axes = plt.subplots(nrows=nrows, ncols=ncols, 
#                              figsize=(4.8 * ncols, 3.6 * nrows),sharex=True, sharey=sharey)
#     if n_panels == 1:
#         axes = np.array([axes])
#     axes = axes.reshape(nrows, ncols)


#     # 作图
#     handles_ref, labels_ref = [], []
#     for i, ct in enumerate(cell_types):
#         ax = axes[i // ncols, i % ncols]
#         sub = df[df["cell_type"] == ct]
#         for cl in sub["clone"].dropna().unique():
#             ss = sub[sub["clone"] == cl].sort_values("center_offset")
#             x = ss["center_offset"].to_numpy()
#             y = ss["prop"].to_numpy()

#             # 简单平滑
#             if isinstance(smooth_win, int) and smooth_win >= 3 and smooth_win % 2 == 1 and len(y) >= smooth_win:
#                 k = smooth_win
#                 y = np.convolve(y, np.ones(k)/k, mode="valid")
#                 pad = k // 2
#                 x = x[pad:-pad]

#             line, = ax.plot(x, y, marker="o", linewidth=1.2, markersize=2, label=str(cl))
#             if str(cl) not in labels_ref:
#                 handles_ref.append(line)
#                 labels_ref.append(str(cl))

#         ax.axvline(0, linestyle="--", linewidth=1)
#         ax.set_title(str(ct))
#         ax.set_xlabel("Distance from boundary (µm)")
#         ax.grid(True, alpha=0.3)


#     for ax in axes.ravel():
#         if ax is not None:
#             ax.tick_params(axis='x', labelbottom=True)

#     fig.supylabel("Cell proportion")

#     # 清理多余子图
#     for j in range(n_panels, nrows * ncols):
#         fig.delaxes(axes[j // ncols, j % ncols])

#     # 合并图例 → 固定在右边
#     if handles_ref and labels_ref:
#         fig.legend(handles_ref, labels_ref, title="Clone",
#                    loc="upper right", bbox_to_anchor=(1.06, 0.88),
#                    frameon=False)

#     plt.tight_layout(rect=[0,0,1,1])
#     plt.show()
#     return fig
