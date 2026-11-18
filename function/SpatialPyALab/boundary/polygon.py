import numpy as np
import pandas as pd
import cv2
from shapely.geometry import Polygon
from typing import List, Union
import anndata as ad

def assign_polygon_regions_to_spots(
    spot_adata: ad.AnnData,
    polys_sorted: List[Union[Polygon]],
    id:str = 'barcode',
    spatial_key: str = "spatial_um",
    label_prefix: str = "region"
) -> pd.DataFrame:
    """
    根据多边形将 spot 坐标映射为所属区域标签。

    Parameters:
    - spot_adata: AnnData 对象，含有 spot 的空间坐标
    - polys_sorted: 多边形列表，shapely.geometry.Polygon 对象
    - spatial_key: obsm 中空间坐标的键（默认是 'spatial_um'）
    - label_prefix: 区域标签的前缀，例如 'region' -> 'region_0'

    Returns:
    - pd.DataFrame，包含 'barcode' 和 'vessel_region_label'
    """
    
    coords = spot_adata.obsm[spatial_key]
    min_x, min_y = coords.min(axis=0).astype(int)
    max_x, max_y = coords.max(axis=0).astype(int)
    width, height = max_x - min_x + 1, max_y - min_y + 1

    # 创建空掩膜图
    label_map = np.full((height, width), fill_value=-1, dtype=np.int16)  

    # 将每个多边形绘制到 label_map 上
    for idx, poly in enumerate(polys_sorted, start=0):
        if poly.is_empty:
            continue
        contour = np.array(poly.exterior.coords, dtype=np.int32)
        contour -= np.array([min_x, min_y])  # 平移坐标到图像范围
        cv2.fillPoly(label_map, [contour], idx)

    # 映射坐标点到掩膜图像坐标
    xy = np.round(coords).astype(int)
    xy_shifted = xy - np.array([min_x, min_y])

    # 防止越界
    in_bounds = (xy_shifted[:, 0] >= 0) & (xy_shifted[:, 0] < width) & \
                (xy_shifted[:, 1] >= 0) & (xy_shifted[:, 1] < height)

    # 初始化结果
    region_ids = np.full(coords.shape[0], fill_value=None, dtype=object)

    # 读取落在掩膜图范围内的点的区域标签
    ys, xs = xy_shifted[in_bounds, 1], xy_shifted[in_bounds, 0]
    region_idx = label_map[ys, xs]
    region_ids[in_bounds] = [
        f"{label_prefix}_{i}" if i >= 0 else None for i in region_idx
    ]

    return pd.DataFrame({
        "id": spot_adata.obs[id].values,
        "vessel_region_label": region_ids
    })

