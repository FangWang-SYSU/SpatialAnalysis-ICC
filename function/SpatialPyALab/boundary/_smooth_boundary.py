import numpy as np
import pandas as pd
from shapely.geometry import Point
from rtree import index

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

def make_points_from_obsm(adata, obsm_key="spatial_um", unit_scale=1.0, swap_xy=False):

    arr = np.asarray(adata.obsm[obsm_key], dtype=float)

    xy = arr.copy()
    if swap_xy:
        xy = xy[:, [1, 0]]

    valid_mask = np.isfinite(xy).all(axis=1)
    kept_idx = np.flatnonzero(valid_mask)
    xy = xy[valid_mask]

    if unit_scale != 1.0:
        xy = xy * unit_scale

    pt_idx = index.Index()
    points = []
    for j, (x, y) in enumerate(xy):
        pt_idx.insert(j, (float(x), float(y), float(x), float(y)))
        points.append(Point(float(x), float(y)))

    return points, pt_idx, kept_idx  

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
    cell_types_all = adata.obs[type_key].astype("category")
    
    rows = []
    for _, row in layer_df.iterrows():
        geom = row["geom"]
        minx, miny, maxx, maxy = geom.bounds

        cand_local_ids = list(pt_idx.intersection((minx, miny, maxx, maxy)))
        if not cand_local_ids:
            continue

        hit_local_ids = [lid for lid in cand_local_ids if geom.covers(points[lid])]
        if len(hit_local_ids) < min_cells:
            continue

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

    if "region_key" not in df.columns:
        if "region" in df.columns:
            df["region_key"] = df["region"].astype(str)
        elif "region_id" in df.columns:
            df["region_key"] = df["region_id"].map(lambda x: f"{region_name_prefix}_{int(x)}")
        else:
            raise KeyError("df_prop don't have 'region_key' and no 'region' / 'region_id',can not mapping clone。")

    needed = {"offset", "side", "cell_type", "n", "total"}
    missing = needed - set(df.columns)
    if missing:
        raise KeyError(f"df_prop missing : {missing}")

    df["clone"] = df["region_key"].map(region_to_clone)
    if df["clone"].isna().any() and "region_id" in df.columns:
        df.loc[df["clone"].isna(), "clone"] = df.loc[df["clone"].isna(), "region_id"].map(region_to_clone)
    df["clone"] = df["clone"].fillna(unknown)

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
    color_by:str,                
    clones: list[str] | None = None,  
    cell_types: list[str] | None = None, 
    ncols: int = 3,
    sharey: bool = False,
    palette: dict | None = None,
    figsize = (14,8),
    bbox_to_anchor = (1.1,0.5)
):

    dd = df.copy()
    if clones is not None:
        dd = dd[dd["clone"].isin(clones)]
    if cell_types is not None:
        dd = dd[dd["cell_type"].isin(cell_types)]
    if dd.empty:
        raise ValueError("empty")

    facet_vals = dd[facet_by].dropna().unique().tolist()
    color_vals = dd[color_by].dropna().unique().tolist()

    if palette is None:
        base = plt.get_cmap("tab20").colors
        palette = {str(k): base[i % len(base)] for i, k in enumerate(sorted(map(str, color_vals)))}

    dd = dd.sort_values([facet_by, color_by, "offset"]).reset_index(drop=True)

    n_panels = len(facet_vals)
    ncols = max(1, min(ncols, n_panels))
    nrows = ceil(n_panels / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,figsize=figsize,sharex=True, sharey=sharey)

    if n_panels == 1:
        axes = np.array([axes])
    axes = axes.reshape(nrows, ncols)

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

    for j in range(n_panels, nrows * ncols):
        fig.delaxes(axes[j // ncols, j % ncols])

    if handles_ref and labels_ref:
        legend_title = "Cell type" if color_by == "cell_type" else "Clone"
        fig.legend(handles_ref, labels_ref, title=legend_title, 
                   loc="right",bbox_to_anchor=bbox_to_anchor)

    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.show()
    return fig
