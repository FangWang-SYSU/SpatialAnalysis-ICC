import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
import pandas as pd

def merge_ficture_result(adata, class_file, anno_file, anno_name):

    class_df = pd.read_csv(class_file, sep='\t', compression='gzip')
    
    if anno_file is not None:
        anno_df = pd.read_csv(anno_file, sep='\t', header=None, names=['K1', anno_name])
        merged_df = pd.merge(class_df, anno_df, on='K1', how='left')
    else:
        merged_df = class_df.copy()
        merged_df[anno_name] = "Cluster_" + merged_df["K1"].astype(str)

    df = adata.obs.copy()
    
    if 'barcode' in df.columns:
        df.index.name = None
    else: 
        df['barcode'] = df.index.str.split(':').str[0]
    
    df = pd.merge(df, merged_df[['barcode', anno_name]], on='barcode', how='left')
    df.index = df["barcode"]
    adata.obs = df
    return adata



def spot_to_hex(mesh_adata):

    import pandas as pd
    import numpy as np

    indices_joined_spatial_units = mesh_adata.uns['indices_joined_spatial_units']
    rows, cols = indices_joined_spatial_units.nonzero()
    
    aggregated_spots = np.empty((len(rows), 2), dtype=object)
    aggregated_spots[:, 0] = mesh_adata.obs.index[rows]
    aggregated_spots[:, 1] = [mesh_adata.uns['spatial_units_obs_names'][col] for col in cols]
    aggregated_spots_df = pd.DataFrame(aggregated_spots, columns=["Hex_ID", "Spot_ID"])
    
    return aggregated_spots_df


import numpy as np
import pandas as pd
from anndata import AnnData


def spot_to_cellid(spot_adata: AnnData,
                   spatial_key: str = 'spatial_manual_fine',
                   mask_key: str = 'staining_image_mask') -> pd.DataFrame:
    cell_ids = []
    valid_barcodes = []

    coords = spot_adata.obsm[spatial_key]
    pixel_coords = np.round(coords).astype(int)
    barcodes = spot_adata.obs_names.to_numpy()
    mask = spot_adata.uns[mask_key]

    for bc, (x, y) in zip(barcodes, pixel_coords):
        try:
            cell_id = mask[x, y] 
            if cell_id > 0:
                valid_barcodes.append(bc)
                cell_ids.append(cell_id)
        except IndexError:
            continue

    df = pd.DataFrame({
        'Spot_ID': valid_barcodes,
        'cell_id': cell_ids
    })
    return df


def add_ficture_result(
    spot_to_other,
    ficture_anno,
    threshold: float = 0.5,
    label_cols: list = ['CellType','CellType_ficture', 'GLMPCA_leiden'],
    spot_id: str = 'Spot_ID',
    segment_id: str = 'Hex_ID',
    barcode_col: str = 'barcode',
    min_label_spots: int = 3,
    conflict_label: str = 'conflict',
):

    import pandas as pd
    import numpy as np
    ficture_anno = ficture_anno.copy()
    if ficture_anno.index.name == barcode_col and barcode_col in ficture_anno.columns:
        ficture_anno = ficture_anno.reset_index(drop=True)

    if barcode_col not in ficture_anno.columns:
        if ficture_anno.index.name == barcode_col:
            ficture_anno = ficture_anno.reset_index()
        else:
            raise ValueError(
                f"error ficture_anno must include `{barcode_col}`）"
            )
    else:
        if ficture_anno.index.name is not None:
            ficture_anno.index.name = None

    label_cols = [c for c in label_cols if c in ficture_anno.columns]
    if len(label_cols) == 0:
        raise ValueError("error")
    merged_df = pd.merge(
        ficture_anno,
        spot_to_other,
        left_on=barcode_col,
        right_on=spot_id,
        how='left',
        validate='m:1'
    )
    merged_df[segment_id] = merged_df[segment_id].astype("string")
    merged_df_in = merged_df[merged_df[segment_id].notna()].copy()
    if merged_df_in.empty:
        raise ValueError(f"error")

    segs = merged_df_in[segment_id].dropna().unique()
    result = pd.DataFrame({segment_id: segs})

    for label_col in label_cols:
        n_labeled = (
            merged_df_in[label_col].notna()
            .groupby(merged_df_in[segment_id], dropna=True)
            .sum()
            .rename("n_labeled")
            .reset_index()
        )

        label_counts = (
            merged_df_in.groupby(segment_id, dropna=True)[label_col]
            .value_counts(dropna=True)
            .rename("n_label")
            .reset_index()
        )
        label_counts = label_counts.merge(
            n_labeled,
            on=segment_id,
            how="left")
        label_counts["frac"] = label_counts["n_label"] / label_counts["n_labeled"]
        label_counts = label_counts.sort_values(
            by=[segment_id, "n_label", "frac", label_col],
            ascending=[True, False, False, True],
            kind="mergesort"
        )
        best = label_counts.groupby(segment_id, as_index=False).head(1).copy()
        best[label_col] = np.where(
            (best["frac"] >= threshold) & (best["n_label"] >= min_label_spots),
            best[label_col].astype(str),
            conflict_label
        )

        best = best[[segment_id, label_col]]
        result = result.merge(best, on=segment_id, how="left")
        result[label_col] = result[label_col].fillna(conflict_label)

    return result


import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns

def build_cell_adata(
    spot_adata,
    cell_adata,
    spatial_key="spatial_manual_fine",
    mask_key="staining_image_mask",
    min_prop=0.8,
    min_spots=3,
    include_nan=False,
    plot=True
):

    cell_type = spot_adata.obs[['barcode', 'CellType_ficture']].reset_index(drop=True)
    cell_type_clean = cell_type.dropna(how="any") 

    spot_to_cell = spot_to_cellid(spot_adata, spatial_key, mask_key)
    spot_to_cell['cell_id'] = spot_to_cell['cell_id'].astype(str)

    merged_df = pd.merge(cell_type_clean, spot_to_cell, on="barcode", how="inner")

    counts = (
        merged_df.groupby(["cell_id", "CellType_ficture"])
        .size()
        .reset_index(name="n_spots")
    )
    totals = counts.groupby("cell_id")["n_spots"].sum().reset_index(name="total_spots")
    counts = counts.merge(totals, on="cell_id")
    counts["prop"] = counts["n_spots"] / counts["total_spots"]

    max_prop = counts.groupby("cell_id")["prop"].max().rename("max_prop")
    counts = counts.merge(max_prop, on="cell_id")
    major_types = counts.loc[
        counts["prop"] == counts["max_prop"],
        ["cell_id", "CellType_ficture", "n_spots", "max_prop"]
    ].drop_duplicates("cell_id")

    if plot:
        plt.figure(figsize=(5, 4))
        sns.histplot(major_types["max_prop"], bins=30, edgecolor="black", kde=True)
        plt.axvline(min_prop, color="red", linestyle="--", label=f"threshold={min_prop}")
        plt.xlabel("Max cell-type proportion per cell (NaN ignored)")
        plt.ylabel("Cell count")
        plt.title("Distribution of major cell-type proportion")
        plt.legend()
        plt.tight_layout()
        plt.show()

    major_types = major_types[
        (major_types["max_prop"] >= min_prop) & (major_types["n_spots"] >= min_spots)
    ]

    spot_keep = merged_df.merge(
        major_types[["cell_id", "CellType_ficture"]],
        on=["cell_id", "CellType_ficture"],
        how="inner"
    )

    if include_nan:
        spot_nan = pd.merge(cell_type[cell_type["CellType_ficture"].isna()],
                            spot_to_cell, on="barcode", how="inner")
        spot_nan = spot_nan[spot_nan["cell_id"].isin(major_types["cell_id"])]
        spot_keep = pd.concat([spot_keep, spot_nan], ignore_index=True)

    spot_major = spot_adata[spot_keep["barcode"].values, :].copy()
    spot_to_cell_map = spot_keep.set_index("barcode")["cell_id"].reindex(spot_major.obs_names)
    cell_index = pd.Categorical(spot_to_cell_map)
    cell_ids = cell_index.categories.astype(str)

    rows, cols = spot_major.X.nonzero()
    data = spot_major.X.data
    new_rows = cell_index.codes[rows]

    cell_matrix = csr_matrix(
        (data, (new_rows, cols)),
        shape=(len(cell_ids), spot_major.n_vars)
    )
    cell_adata.obs['cell_ID_mask'] = cell_adata.obs['cell_ID_mask'].astype(str)
    mask_indexer = cell_adata.obs["cell_ID_mask"].isin(cell_ids)
    cell_adata_new = cell_adata[mask_indexer].copy()

    order = pd.Categorical(
        cell_adata_new.obs["cell_ID_mask"], categories=cell_ids, ordered=True
    ).codes
    cell_adata_new.X = cell_matrix[order, :]
    cell_type_map = major_types.set_index("cell_id")["CellType_ficture"].to_dict()
    cell_adata_new.obs["major_celltype"] = cell_adata_new.obs["cell_ID_mask"].map(cell_type_map)
    return cell_adata_new

if __name__ == "__main__":
    ...
