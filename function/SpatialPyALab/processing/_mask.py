from __future__ import annotations

from typing import Literal, TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.ndimage import binary_dilation
from skimage.filters import threshold_otsu
from skimage.filters import gaussian
from skimage.color import rgb2hsv


if TYPE_CHECKING:
    from anndata import AnnData
    from numpy import ndarray


def binary_fill_holes(
        input, 
        structure=None, 
        output=None, 
        origin=0
    ):

    mask = np.logical_not(input)
    tmp = np.zeros(mask.shape, bool)
    output = binary_dilation(tmp, structure, -1, mask, None, 1,
                            origin)
    
    np.logical_not(output, output)
    
    return output

def mask_tissue(
    image: ndarray,
    mask_gaussian_blur: float = 5,
    fill_holes: bool = True
):
    hsv_image = rgb2hsv(image)
    s_image_gaussian = gaussian(hsv_image[..., 1], sigma=mask_gaussian_blur)
    thresh = threshold_otsu(s_image_gaussian)
    s_image_gaussian_binary = s_image_gaussian > thresh
    if fill_holes:
        mask = binary_fill_holes(s_image_gaussian_binary).astype(int)
    else:
        mask = s_image_gaussian_binary

    return mask


def dilated_mask(
        input: ndarray,
        dilation_kernel: Literal['square', 'cross', 'diamond'] | ndarray | None = None,
        dilation_step: int = 1,
        mask: ndarray = None,
        output_contain: ndarray = None,
        border_value: Literal[0, 1] = 0,
        origin: int | tuple[int] = 0,
        brute_force: bool = False
    ):

    if isinstance(dilation_kernel, str):
        if dilation_kernel == 'square':
            dilation_kernel = np.ones((3, 3))
        elif dilation_kernel == 'cross':
            dilation_kernel = np.array([[0, 1, 0], 
                                        [1, 1, 1], 
                                        [0, 1, 0]])
        elif dilation_kernel == 'diamond':
            dilation_kernel = np.array([[0, 0, 0],
                                        [0, 1, 0],
                                        [0, 0, 0]])
        else:
            raise ValueError("Invalid dilation_kernel value. Please provide 'square', 'cross', 'diamond' or a custom kernel.")

    return binary_dilation(
        input, 
        structure = dilation_kernel,
        iterations = dilation_step,
        mask = mask,
        output = output_contain,
        border_value = border_value,
        origin = origin,
        brute_force = brute_force)


def filter_mask(
        mask: ndarray,
        top_n: int = 1
    ):
    labeled_array, num_features = ndimage.label(mask)
    sizes = ndimage.sum(mask, labeled_array, range(1, num_features + 1))
    labels = [l+1 for l in np.argsort(sizes)[::-1]]
    save_label = labels[:top_n]
    filter_mask = np.zeros_like(mask)
    filter_mask[np.isin(labeled_array, save_label)] = 1

    return filter_mask


import numpy as np
import pandas as pd
from anndata import AnnData
from numpy import ndarray

def mask_adata(
    adata: AnnData,
    mask: ndarray,
    spatial_key: str = "spatial_manual_fine",
):
    height, width = mask.shape

    coords = adata.obsm[spatial_key]
    x = coords[:, 0]
    y = coords[:, 1]

    x = np.asarray(np.round(x), dtype=int)
    y = np.asarray(np.round(y), dtype=int)

    in_bounds = (x >= 0) & (x < width) & (y >= 0) & (y < height)

    in_tissue = np.zeros(adata.n_obs, dtype=bool)
    in_tissue[in_bounds] = mask[y[in_bounds], x[in_bounds]].astype(bool)

    out_tissue = ~in_tissue 

    adata_in  = adata[in_tissue, :].copy()
    adata_out = adata[out_tissue, :].copy()
    return adata_in, adata_out
  
