import pandas as pd
import os, sys, gzip
from scipy.spatial import cKDTree
import argparse


def obtain_barcodes(_args):
    # Parse command-line arguments
    parser = argparse.ArgumentParser(prog="obtain_barcodes",
                                     description="Match spatial cluster coordinates with barcodes and split by cluster or cell type.")
    parser.add_argument('--coor_barcode', type=str, required=True, 
                            help='Path to the gzipped TSV file containing spatial barcode coordinates. Expected columns: X, Y, barcode, X_spatial, Y_spatial, etc.')
    parser.add_argument('--ficture_cluster', type=str,  required=True, 
                        help='Ficture Result')
    parser.add_argument('--output_dir', type=str, required=True, default='', 
                       help='Maximum Euclidean distance to match clusters with barcodes. Default is 10.')
    parser.add_argument('--distance_threshold', type=float, default=10, 
                        help='Scaling factor for spatial coordinates')
    parser.add_argument('--anno_path', type=str, required=False, 
                        help='Annotation file path')
    args = parser.parse_args(_args)

    # If no arguments are provided, print help
    if len(_args) == 0:
        parser.print_help()
        return
    
    output_dir = args.output_dir
    distance_threshold = args.distance_threshold

    
    os.makedirs(output_dir, exist_ok=True)
    
    barcode = pd.read_csv(args.coor_barcode, sep='\t', compression='gzip')
    cluster = pd.read_csv(args.ficture_cluster, sep='\t', compression='gzip', skiprows=3)
    
    plot_cor = cluster[['X', 'Y', 'K1']]
    
    tree = cKDTree(barcode[['X', 'Y']])
    distances, indices = tree.query(plot_cor[['X', 'Y']], distance_upper_bound=distance_threshold)
    
    valid = distances != float('inf')
    result = pd.concat(
        [plot_cor[valid].reset_index(drop=True), barcode.iloc[indices[valid]].reset_index(drop=True)],
        axis=1
    )
    
    result.columns = ['X', 'Y', 'K1', 'barcode', 'X_spatial', 'Y_spatial', 'X_scale', 'Y_scale', 'X_self', 'Y_self']
    
    result_path = os.path.join(output_dir, 'match_barcode_cluster.tsv.gz')
    with gzip.open(result_path, 'wt', encoding='utf-8') as f:
        result.to_csv(f, sep='\t', index=False, header=True)

    df = result[['barcode', 'K1', 'X', 'Y']]
    
    if args.anno_path:
        anno = pd.read_csv(args.anno_path, sep='\t', header=None)
        anno.columns = ['K1', 'celltypes']
        merged_data = pd.merge(df, anno, on='K1', how='left')
        split_col = 'celltypes'
    else:
        split_col = 'K1'
        merged_data = df
    
    grouped = merged_data.groupby(split_col)

    for col, group in grouped: 
        barcodes = group['barcode']
        output_path = os.path.join(output_dir, f'{col}_barcodes.txt')
        barcodes.to_csv(output_path, sep='\t', index=False, header=False)
    
    return merged_data


if __name__ == "__main__":
    obtain_barcodes(sys.argv[1:])