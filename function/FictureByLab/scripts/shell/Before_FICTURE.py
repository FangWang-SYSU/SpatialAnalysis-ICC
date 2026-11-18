import os, sys
import pandas as pd
import re
import anndata as ad
import numpy as np
import argparse
 
def process_h5ad_to_transcripts(_args):
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(prog="process_h5ad_to_transcripts")
    parser.add_argument('--input_h5ad', type=str, required=True, help='Input h5ad file')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save output files')
    parser.add_argument('--coord_key', type=str, default='spatial', help='Key in the h5ad file for spatial coordinates')
    parser.add_argument('--dummy_genes', type=str, default='RPS|MT-|mt-|RPL', help='A regex describing the names of genes to exclude')
    #parser.add_argument('--mu_scale', type=float, default=1, help='Scaling factor for spatial coordinates')
    args = parser.parse_args(_args)

    # If no arguments are provided, print help
    if len(_args) == 0:
        parser.print_help()
        return

    # Extract arguments from parsed input
    input_h5ad = args.input_h5ad
    output_dir = args.output_dir
    coord_key = args.coord_key
    dummy_genes = args.dummy_genes
    #mu_scale = args.mu_scale


    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Read the h5ad file
    adata = ad.read_h5ad(input_h5ad)

    # Check if the coordinate key exists
    if coord_key not in adata.obsm:
        raise ValueError(f"The coordinate key '{coord_key}' was not found in the h5ad file. Available keys: {list(adata.obsm.keys())}")


    # Get the coordinate and gene expression matrix
    row_indices, column_indices = adata.X.tocoo().nonzero()
    xcoord, ycoord = np.split(adata.obsm[coord_key][row_indices, :].transpose(), 2, axis=0)

    #spatial_coords = adata.obsm[coord_key][row_indices]
    #scaled_coords = np.round(spatial_coords / mu_scale, 2).T
    #xcoord, ycoord = np.split(spatial_coords, 2, axis = 0)

    # Calculate the coordinate range
    xmin, xmax = xcoord.min(), xcoord.max()
    ymin, ymax = ycoord.min(), ycoord.max()
    

    # Save the coordinate range to a file
    coord_file = os.path.join(output_dir, 'coordinate_minmax.tsv')
    df_range = pd.DataFrame([xmin, xmax, ymin, ymax], index=['xmin', 'xmax', 'ymin', 'ymax'])
    df_range.to_csv(coord_file, sep='\t', header=False)


    # Create the transcripts.tsv.gz file
    gene_names = adata.var_names[column_indices].to_numpy().reshape(1, -1)
    counts = adata.X[row_indices, column_indices].astype(int)
    data = np.stack((xcoord, ycoord, gene_names, counts), axis=0).transpose()
    df_transcripts = pd.DataFrame(data, columns=['X', 'Y', 'gene', 'Count']).sort_values(by='Y')
    transcripts_file = os.path.join(output_dir, 'transcripts.tsv.gz')
    df_transcripts.to_csv(transcripts_file, sep='\t', index=False, compression='gzip')

    # 计算 feature.clean.tsv.gz 文件
    feature = df_transcripts.groupby('gene').agg({'Count': 'sum'}).reset_index()

    # Filter out specified dummy genes
    if dummy_genes:
        feature = feature[~feature['gene'].str.contains(dummy_genes, flags=re.IGNORECASE, regex=True)]

    feature_clean_file = os.path.join(output_dir, 'feature.clean.tsv.gz')
    feature.to_csv(feature_clean_file, sep='\t', index=False, compression='gzip')

    return transcripts_file, coord_file, feature_clean_file

# 示例调用
if __name__ == '__main__':
  
    try:
        transcripts, coord_file, feature_clean_file = process_h5ad_to_transcripts(sys.argv[1:])
        print(f"Transcripts: {transcripts}")
        print(f"Coordinate range: {coord_file}")
        print(f"Feature clean: {feature_clean_file}")
    except ValueError as e:
        print(f"Error: {e}")
