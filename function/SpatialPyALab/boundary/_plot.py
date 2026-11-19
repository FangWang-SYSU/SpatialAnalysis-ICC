
def plot_layers_boundaries(layer_adata,
                           layers_dict,
                           ax=None, 
                           edgecolor='black', 
                           linewidth=1, 
                           alpha=1,
                           linestyle='-',
                           boundary=None):
 
    import matplotlib.pyplot as plt
    from shapely.geometry import MultiPolygon,Polygon

    x_min, x_max = layer_adata.obsm['spatial_um'][:, 0].min(), layer_adata.obsm['spatial_um'][:, 0].max()
    y_min, y_max = layer_adata.obsm['spatial_um'][:, 1].min(), layer_adata.obsm['spatial_um'][:, 1].max()
    boundary = Polygon([(x_min, y_min), (x_min, y_max), (x_max, y_max), (x_max, y_min)])
    new_layers_dict = {k: [layer[0] for layer in layers] for k, layers in layers_dict.items()}



    if ax is None:
        fig, ax = plt.subplots()

    for layers in new_layers_dict.values():
        for layer in layers:
            if boundary:
                layer = layer.intersection(boundary)

            if isinstance(layer, MultiPolygon):
                for geom in layer.geoms:
                    xs, ys = geom.exterior.xy
                    ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha, linestyle=linestyle)
                    for interior in geom.interiors:
                        xs, ys = interior.xy
                        ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha, linestyle=linestyle)
            else:
                xs, ys = layer.exterior.xy
                ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha, linestyle=linestyle)
                for interior in layer.interiors:
                    xs, ys = interior.xy
                    ax.plot(xs, ys, color=edgecolor, linewidth=linewidth, alpha=alpha, linestyle=linestyle)

    return ax

