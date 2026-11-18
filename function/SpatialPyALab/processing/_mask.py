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
    """
    Fill the holes in binary objects, compatible with dask.

    Args:
        input : array_like
            N-D binary array with holes to be filled
        structure : array_like, optional
            Structuring element used in the computation; large-size elements
            make computations faster but may miss holes separated from the
            background by thin regions. The default element (with a square
            connectivity equal to one) yields the intuitive result where all
            holes in the input have been filled.
        output : ndarray, optional
            Array of the same shape as input, into which the output is placed.
            By default, a new array is created.
        origin : int, tuple of ints, optional
            Position of the structuring element.

    Returns:
        out : ndarray
            Transformation of the initial image `input` where holes have been
            filled.

    Notes:
        This was adapted from skimage
    """

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
    """
    Prepare an image for feature matching by applying a series of image processing steps.

    Args:
        image (np.ndarray): Input RGB image to be prepared.
        mask_gaussian_blur (float, optional): Standard deviation for Gaussian blurring of the saturation channel.

    Returns:
        mask: Binary mask of the image after applying Gaussian blurring and Otsu's thresholding.
    """

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
    """
    Dilate a binary mask using a specified kernel.

    Args:
        input (ndarray): Input binary mask.
        dilation_kernel (str or ndarray, optional): Dilation kernel. Defaults to None.
        dilation_step (int, optional): Number of dilation steps. Defaults to 1.
        mask (ndarray, optional): Mask to apply dilation to. Defaults to None.
        output_contain (ndarray, optional): Output array to contain the dilated mask. Defaults to None.
        border_value (int, optional): Value to use for the border. Defaults to 0.
        origin (int or tuple, optional): Position of the structuring element applied for the pixel. Defaults to 0.
        brute_force (bool, optional): Whether to use brute force method. Defaults to False.
    
    Returns:
        output_contain (ndarray): Output array containing the dilated mask.
    """

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
    """
    Filter out small objects from a binary mask.
    
    Args:
        mask (ndarray): Binary mask to filter.
        top_n (int, optional): Number of top objects to keep. Defaults to 1.
    
    Returns:
        filter_mask: Filtered binary mask.
    """

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

    # 如果你的坐标可能是 float，建议转 int（否则索引 mask 会出问题）
    x = np.asarray(np.round(x), dtype=int)
    y = np.asarray(np.round(y), dtype=int)

    in_bounds = (x >= 0) & (x < width) & (y >= 0) & (y < height)

    in_tissue = np.zeros(adata.n_obs, dtype=bool)
    in_tissue[in_bounds] = mask[y[in_bounds], x[in_bounds]].astype(bool)

    out_tissue = ~in_tissue  # 组织外（包含越界点）

    adata_in  = adata[in_tissue, :].copy()
    adata_out = adata[out_tissue, :].copy()
    return adata_in, adata_out



# def mask_adata(
#         adata: AnnData,
#         mask: ndarray,
#         spatial_key: str = 'spatial_fine'
#     ):

#     # Get the mask from the image within adata object 
#     height, width = mask.shape

#     # Get the spots from the adata object
#     spots = pd.DataFrame(adata.obsm[spatial_key], columns = ['X', 'Y'], index = adata.obs_names)
#     x_coords = adata.obsm[spatial_key][:, 0]
#     y_coords = adata.obsm[spatial_key][:, 1]

#     # Crop the spots within the image bounds
#     in_bounds = (
#         (x_coords >= 0) & 
#         (x_coords < width) & 
#         (y_coords >= 0) & 
#         (y_coords < height)
#     )

#     # Crop the spots within the mask
#     valid_mask = np.zeros(x_coords.shape[0], dtype=bool)
#     valid_mask[in_bounds] = mask[
#         y_coords[in_bounds],  # 行索引对应Y坐标
#         x_coords[in_bounds]   # 列索引对应X坐标
#     ]

#     spots_filtered = spots[valid_mask].copy()
#     tissue_adata = adata[spots_filtered.index, :].copy()

#     return tissue_adata

    