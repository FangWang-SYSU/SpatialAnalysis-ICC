import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_areas(areas, bins = 100):

    fig, ax = plt.subplots(figsize=(7, 4))

    min_area = max(1, np.min(areas))  
    max_area = np.max(areas)
    bins = np.logspace(np.log10(min_area), np.log10(max_area), bins)

    counts, bins, patches = ax.hist(areas, 
                                bins=bins,
                                edgecolor='black',
                                alpha=0.7)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Area (log scale)', fontsize=10)
    ax.set_ylabel('Count (log scale)', fontsize=10)
    ax.set_title('Connected Components Area Distribution', fontsize=10)

    ax.tick_params(axis='both', which='both', labelsize=10)

    ax.set_xlim(left=min_area*0.9, right=max_area*1.1)
    ax.set_ylim(bottom=0.9, top=np.max(counts)*10)

    cmap = plt.cm.viridis
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    for c, p in zip(bin_centers, patches):
        plt.setp(p, 'facecolor', cmap(np.log10(c)/np.log10(max_area)))

    sm = plt.cm.ScalarMappable(cmap=cmap, 
                            norm=plt.Normalize(vmin=np.log10(min_area), 
                                            vmax=np.log10(max_area)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Log10(Area)', rotation=270, labelpad=15)

    plt.tight_layout()
    plt.show()
