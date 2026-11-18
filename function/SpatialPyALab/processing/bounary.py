def find_large_holes(init_mask: np.ndarray, min_hole_area: int = 10) -> np.ndarray:
    """
    Identify holes in a binary mask and return only those with area > min_hole_area.

    Args:
        init_mask (np.ndarray): Binary mask (True = tissue, False = background).
        min_hole_area (int): Minimum area (in pixels) for a hole to be retained.

    Returns:
        large_hole_mask (np.ndarray): Boolean mask with only large holes retained.
    """

    # Step 1: 找到所有的“洞”区域：也就是 init_mask == False
    holes = ~init_mask

    # Step 2: 连通区域标记
    labeled_holes, num_holes = label(holes)

    # Step 3: 使用 regionprops 来分析每个洞的面积
    props = regionprops(labeled_holes)

    # Step 4: 初始化一个空掩膜，只保留面积 > min_hole_area 的洞
    large_hole_mask = np.zeros_like(init_mask, dtype=bool)

    for prop in props:
        if prop.area > min_hole_area:
            large_hole_mask[labeled_holes == prop.label] = True

    return large_hole_mask
    
