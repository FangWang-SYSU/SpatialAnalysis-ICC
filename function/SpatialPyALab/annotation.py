import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix
import pandas as pd



def merge_ficture_result(adata, class_file, anno_file, anno_name):
    """
    合并注释数据并更新 `spot_adata.obs` 中的相关信息。

    Parameters:
    - adata: AnnData 对象，包含 `.obs` 属性。
    - project_path: 项目文件路径，用于定位文件。
    - class_file: 类别文件路径，包含 barcodes 和 K1。
    - anno_file: 注释文件路径，包含 K1 和 Label。

    Returns:
    - updated adata: 更新后的 AnnData 对象。
    """

    class_df = pd.read_csv(class_file, sep='\t', compression='gzip')
    
    if anno_file is not None:
        # 读取注释文件
        anno_df = pd.read_csv(anno_file, sep='\t', header=None, names=['K1', anno_name])
        merged_df = pd.merge(class_df, anno_df, on='K1', how='left')
    else:
        # 没有 anno_file → 生成 NonRef 列
        merged_df = class_df.copy()
        merged_df[anno_name] = "Cluster_" + merged_df["K1"].astype(str)

    # 合并 barcode 和 Label 信息
    df = adata.obs.copy()
    
    if 'barcode' in df.columns:
        df.index.name = None
    else: 
        df['barcode'] = df.index.str.split(':').str[0]
    
    df = pd.merge(df, merged_df[['barcode', anno_name]], on='barcode', how='left')
    df.index = df["barcode"]

    # 更新 adata.obs
    adata.obs = df

    return adata




def spot_to_hex(mesh_adata):
    """
    Extracts and formats aggregated spatial units from mesh_adata.
    
    Parameters:
        mesh_adata: AnnData object with .uns['indices_joined_spatial_units'] and 
                    .uns['spatial_units_obs_names'] populated.
    
    Returns:
        pd.DataFrame: A DataFrame with columns ["Hex_ID", "Spot_ID"] representing
                      the mapping between hexagons and spatial spots.
    """

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
    """
    从 spot_adata 中的空间坐标和 adata_segment 中的 mask 提取 barcode 与 cell_id 的对应关系。

    参数：
        spot_adata: 包含坐标和条码的 AnnData（通常是 spot 级别）
        adata_segment: 包含掩膜图像的 AnnData
        spatial_key: 坐标在 obsm 中的键（默认是 'spatial_manual_fine'）
        mask_key: 掩膜在 uns 中的键（默认是 'staining_image_mask'）

    返回：
        pd.DataFrame，包含两列：barcode 和 cell_id
    """
    cell_ids = []
    valid_barcodes = []

    coords = spot_adata.obsm[spatial_key]
    pixel_coords = np.round(coords).astype(int)
    barcodes = spot_adata.obs_names.to_numpy()
    mask = spot_adata.uns[mask_key]

    for bc, (x, y) in zip(barcodes, pixel_coords):
        try:
            cell_id = mask[x, y]  # 行是 y，列是 x
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
    """
    给 segment（如 Hex / cell）分配标签：由 segment 内 spot 的标签投票决定

    ✅ 投票原则（已按你的要求修改）：
    - 分母 n_labeled = 该 segment 内“该 label_col 非 NA”的 spot 数（只在已标注 spot 中投票）
    - 某 label 需同时满足：
        1) n_label / n_labeled >= threshold
        2) n_label >= min_label_spots
      才能分配为该 label，否则为 conflict_label

    其余逻辑保持不变（merge/检查/输出结构一致）
    """
    import pandas as pd
    import numpy as np

    # ---------- 基本检查 ----------
    if spot_id not in spot_to_other.columns:
        raise ValueError(f"❌ spot_to_other 缺少 spot_id 列: `{spot_id}`")
    if segment_id not in spot_to_other.columns:
        raise ValueError(f"❌ spot_to_other 缺少 segment_id 列: `{segment_id}`")

    # ---------- 规范化 ficture_anno：强制 barcode_col 成为“唯一的列”，避免 pandas merge 歧义 ----------
    ficture_anno = ficture_anno.copy()

    # 情况A：barcode 同时在 index.name 和 columns 里（merge 会报歧义）
    if ficture_anno.index.name == barcode_col and barcode_col in ficture_anno.columns:
        ficture_anno = ficture_anno.reset_index(drop=True)

    # 情况B：barcode 不在列里，但在 index.name 里
    if barcode_col not in ficture_anno.columns:
        if ficture_anno.index.name == barcode_col:
            ficture_anno = ficture_anno.reset_index()
        else:
            raise ValueError(
                f"❌ ficture_anno 必须包含 `{barcode_col}` 列（或其 index.name 为 `{barcode_col}`）"
            )
    else:
        # 情况C：列里有 barcode，但 index 可能也带名字；清掉避免未来坑
        if ficture_anno.index.name is not None:
            ficture_anno.index.name = None

    # ---------- 只保留存在的 label_cols ----------
    label_cols = [c for c in label_cols if c in ficture_anno.columns]
    if len(label_cols) == 0:
        raise ValueError("❌ ficture_anno 中未找到任何 label_cols 列，请检查 label_cols 参数")

    # ---------- 合并：spot 关联到 segment ----------
    merged_df = pd.merge(
        ficture_anno,
        spot_to_other,
        left_on=barcode_col,
        right_on=spot_id,
        how='left',
        validate='m:1'
    )

    # segment_id 统一转 string（避免 1 和 "1" 分裂；<NA> 仍是 NA）
    merged_df[segment_id] = merged_df[segment_id].astype("string")

    # 只统计真正落入 segment 的 spot
    merged_df_in = merged_df[merged_df[segment_id].notna()].copy()
    if merged_df_in.empty:
        raise ValueError(f"❌ merge 后所有 spot 的 `{segment_id}` 都是 NA，检查 spot_id/segment_id 是否对齐")

    # result 框架：每个 segment 一行（这里沿用你原先“来自数据里出现过的 segment”）
    segs = merged_df_in[segment_id].dropna().unique()
    result = pd.DataFrame({segment_id: segs})

    # ---------- 对每个 label_col 做一次分配 ----------
    for label_col in label_cols:
        # 分母：每个 segment 内“已标注（非NA）” spot 数
        n_labeled = (
            merged_df_in[label_col].notna()
            .groupby(merged_df_in[segment_id], dropna=True)
            .sum()
            .rename("n_labeled")
            .reset_index()
        )

        # 分子：每个 segment 内各 label 的 spot 数（自动忽略 NA 标签）
        label_counts = (
            merged_df_in.groupby(segment_id, dropna=True)[label_col]
            .value_counts(dropna=True)
            .rename("n_label")
            .reset_index()
        )

        # 合并分母
        label_counts = label_counts.merge(
            n_labeled,
            on=segment_id,
            how="left"
        )

        # 计算占比（⚠️ 分母改为 n_labeled）
        label_counts["frac"] = label_counts["n_label"] / label_counts["n_labeled"]

        # 找每个 segment 的最优 label（先 n_label，再 frac，再 label 字典序稳定）
        label_counts = label_counts.sort_values(
            by=[segment_id, "n_label", "frac", label_col],
            ascending=[True, False, False, True],
            kind="mergesort"
        )
        best = label_counts.groupby(segment_id, as_index=False).head(1).copy()

        # 应用双阈值：frac>=threshold 且 n_label>=min_label_spots
        best[label_col] = np.where(
            (best["frac"] >= threshold) & (best["n_label"] >= min_label_spots),
            best[label_col].astype(str),
            conflict_label
        )

        best = best[[segment_id, label_col]]
        result = result.merge(best, on=segment_id, how="left")

        # segment 内该 label_col 全是 NA（即 n_labeled=0）时：best 会缺失 -> 填 conflict
        result[label_col] = result[label_col].fillna(conflict_label)

    return result




# def add_ficture_result(spot_to_other: pd.DataFrame,
#                        ficture_anno: pd.DataFrame,
#                        threshold: float = 0.5,
#                        label_cols: list = ['CellType_ficture', 'GLMPCA_leiden'],
#                        spot_id: str ='Spot_ID',
#                        segment_id: str = 'Hex_ID') -> pd.DataFrame:
 
#     """
#     Parameters:
#     - spot_to_other: pd.DataFrame
#     - ficture_anno: pd.DataFrame
#     - threshold: float
#     - id_col: str,Hex ID 所在列的列名,默认是 'Hex_ID'


#     Returns:
#     - result: pd.DataFrame,包含每个 hex 区域的主标签结果
#     """

#     import pandas as pd
    
#     if ficture_anno.index.name == 'barcode':
#         ficture_anno.index.name = None
#     ficture_anno = ficture_anno.reset_index()  # 解决 index/column 冲突


#     #ficture_anno = ficture_anno.rename(columns={'CellType_ficture': 'Label'})
    
#     # 合并得到 Hex_ID
#     merged_df = pd.merge(ficture_anno, spot_to_other, left_on='barcode', right_on=spot_id, how='left')
#     merged_df[segment_id] = merged_df[segment_id].astype(str)

#     result = pd.DataFrame({segment_id: merged_df[segment_id].unique()})

#     for label_col in label_cols:
#         label_counts = merged_df.groupby(merged_df[segment_id])[label_col].value_counts(normalize=True).unstack(fill_value=0)


#         label_assignments = [
#             (hex_id, row.idxmax() if row.max() >= threshold else 'conflict')
#             for hex_id, row in label_counts.iterrows()
#         ]

#         label_df = pd.DataFrame(label_assignments, columns=[segment_id, label_col])
#         result = pd.merge(result, label_df, on=segment_id, how='left')

#     return result


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
    """
    基于 spot_adata 的细胞类型注释，聚合为 cell_adata。
    
    步骤：
    1. 比例计算只考虑非 NaN spots。
    2. 保留的条件：major celltype 占比 >= min_prop 且该 celltype 的 spot 数 >= min_spots。
    3. 聚合时 include_nan 控制是否把 NaN spots 加回。
    4. 最终在 cell_adata.obs 中增加 major_celltype。
    """

    # Step 1: 提取 spot 的 celltype 注释
    cell_type = spot_adata.obs[['barcode', 'CellType_ficture']].reset_index(drop=True)
    cell_type_clean = cell_type.dropna(how="any")  # 非 NaN 部分用于比例计算

    # Step 2: spot -> cell 映射
    spot_to_cell = spot_to_cellid(spot_adata, spatial_key, mask_key)
    spot_to_cell['cell_id'] = spot_to_cell['cell_id'].astype(str)

    # Step 3: 合并（只用非 NaN 部分）
    merged_df = pd.merge(cell_type_clean, spot_to_cell, on="barcode", how="inner")

    # Step 4: 每个 cell 内各类型 spot 数
    counts = (
        merged_df.groupby(["cell_id", "CellType_ficture"])
        .size()
        .reset_index(name="n_spots")
    )
    totals = counts.groupby("cell_id")["n_spots"].sum().reset_index(name="total_spots")
    counts = counts.merge(totals, on="cell_id")
    counts["prop"] = counts["n_spots"] / counts["total_spots"]



    # Step 5: 找到 major celltype
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


    # Step 6: 过滤条件
    major_types = major_types[
        (major_types["max_prop"] >= min_prop) & (major_types["n_spots"] >= min_spots)
    ]

    # Step 7: 保留 major celltype 的 spots
    spot_keep = merged_df.merge(
        major_types[["cell_id", "CellType_ficture"]],
        on=["cell_id", "CellType_ficture"],
        how="inner"
    )

    # Step 8: 如果 include_nan=True，把对应 cell 的 NaN spots 加回
    if include_nan:
        spot_nan = pd.merge(cell_type[cell_type["CellType_ficture"].isna()],
                            spot_to_cell, on="barcode", how="inner")
        spot_nan = spot_nan[spot_nan["cell_id"].isin(major_types["cell_id"])]
        spot_keep = pd.concat([spot_keep, spot_nan], ignore_index=True)

    # Step 9: 子集化 spot_adata
    spot_major = spot_adata[spot_keep["barcode"].values, :].copy()

    # Step 10: 构建稀疏矩阵 cell_matrix
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

    # Step 11: 匹配 cell_adata 并替换 .X
    cell_adata.obs['cell_ID_mask'] = cell_adata.obs['cell_ID_mask'].astype(str)
    mask_indexer = cell_adata.obs["cell_ID_mask"].isin(cell_ids)
    cell_adata_new = cell_adata[mask_indexer].copy()

    order = pd.Categorical(
        cell_adata_new.obs["cell_ID_mask"], categories=cell_ids, ordered=True
    ).codes
    cell_adata_new.X = cell_matrix[order, :]

    # Step 12: 增加 major_celltype
    cell_type_map = major_types.set_index("cell_id")["CellType_ficture"].to_dict()
    cell_adata_new.obs["major_celltype"] = cell_adata_new.obs["cell_ID_mask"].map(cell_type_map)


    return cell_adata_new



# def add_ficture_result(class_file_path, 
#                            annotation_file_path, 
#                            spot_to_hex, 
#                            threshold=0.5):
#     """
#     给 Hex_ID 分配 celltype 标签，并返回处理后的 mesh_adata.obs 格式的数据（DataFrame）。
    
#     Parameters:
#     - class_file_path: str, 含有 barcode-K1 映射的文件路径
#     - annotation_file_path: str, 含有 K1-Label 映射的注释文件路径
#     - aggregated_spots_df: pd.DataFrame, 至少包含 'Hex_ID' 和 'Spot_ID' 列
#     - threshold: float, 分配标签所需的最小比例阈值，默认 0.5
    
#     Returns:
#     - final_label_df: pd.DataFrame, 包含 'Hex_ID' 和 'celltype.ficture' 两列
#     """
#     df_class = pd.read_csv(class_file_path, sep='\t', compression='gzip')
#     anno = pd.read_csv(annotation_file_path, sep='\t', header=None, names=['K1', 'Label'])
#     merged_df = pd.merge(df_class, anno, on='K1', how='left')
#     merged_df = pd.merge(merged_df, spot_to_hex, left_on='barcode', right_on='Spot_ID', how='left')

#     label_counts = merged_df.groupby('Hex_ID')['Label'].value_counts(normalize=True).unstack(fill_value=0)

#     label_assignments = []
#     for hex_id, counts in label_counts.iterrows():
#         max_label = counts.idxmax()
#         max_count = counts.max()
#         if max_count >= threshold:
#             label_assignments.append((hex_id, max_label))
#         else:
#             label_assignments.append((hex_id, 'conflict'))

#     label_df = pd.DataFrame(label_assignments, columns=['Hex_ID', 'CellType_ficture'])

#     return label_df


if __name__ == "__main__":
    ...