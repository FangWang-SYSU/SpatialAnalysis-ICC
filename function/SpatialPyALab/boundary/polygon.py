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
    
    coords = spot_adata.obsm[spatial_key]
    min_x, min_y = coords.min(axis=0).astype(int)
    max_x, max_y = coords.max(axis=0).astype(int)
    width, height = max_x - min_x + 1, max_y - min_y + 1
    label_map = np.full((height, width), fill_value=-1, dtype=np.int16)  
    for idx, poly in enumerate(polys_sorted, start=0):
        if poly.is_empty:
            continue
        contour = np.array(poly.exterior.coords, dtype=np.int32)
        contour -= np.array([min_x, min_y])  
        cv2.fillPoly(label_map, [contour], idx)

    xy = np.round(coords).astype(int)
    xy_shifted = xy - np.array([min_x, min_y])

    in_bounds = (xy_shifted[:, 0] >= 0) & (xy_shifted[:, 0] < width) & \
                (xy_shifted[:, 1] >= 0) & (xy_shifted[:, 1] < height)

    region_ids = np.full(coords.shape[0], fill_value=None, dtype=object)

    ys, xs = xy_shifted[in_bounds, 1], xy_shifted[in_bounds, 0]
    region_idx = label_map[ys, xs]
    region_ids[in_bounds] = [
        f"{label_prefix}_{i}" if i >= 0 else None for i in region_idx
    ]

    return pd.DataFrame({
        "id": spot_adata.obs[id].values,
        "vessel_region_label": region_ids
    })

