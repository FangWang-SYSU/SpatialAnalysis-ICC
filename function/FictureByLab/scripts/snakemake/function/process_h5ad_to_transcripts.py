import os, sys
import pandas as pd
import re
import anndata as ad
import numpy as np
import argparse

def _parse_genes_arg(genes_str: str | None):
    if not genes_str:
        return None
    parts = re.split(r"[,\s;]+", genes_str.strip())
    parts = [g for g in parts if g]
    return parts if parts else None

def _load_genes_file(path: str | None):
    if not path:
        return None
    if not os.path.exists(path):
        raise ValueError(f"The genes_file does not exist.: {path}")
    genes = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.append(g)
    return genes if genes else None

def process_h5ad_to_transcripts(_args):

    parser = argparse.ArgumentParser(prog="process_h5ad_to_transcripts")
    parser.add_argument('--input_h5ad', type=str, required=True, help='Input h5ad file')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to save output files')
    parser.add_argument('--coord_key', type=str, default='spatial', help='Key in the h5ad file for spatial coordinates')
    parser.add_argument('--dummy_genes', type=str, default='RPS|MT-|mt-|RPL', help='A regex describing the names of genes to exclude')
    parser.add_argument('--mu_scale', type=float, default=1, help='Scaling factor for spatial coordinates')

    parser.add_argument("--genes_file", nargs="?",  default=None, const=None,  
        help="Optional. Path to a gene whitelist file (one gene per line). If not provided, no gene filtering is applied.")
    
    args = parser.parse_args(_args)
    if len(_args) == 0:
        parser.print_help()
        return

    input_h5ad = args.input_h5ad
    output_dir = args.output_dir
    coord_key = args.coord_key
    dummy_genes = args.dummy_genes
    mu_scale = args.mu_scale

    os.makedirs(output_dir, exist_ok=True)

   
    adata = ad.read_h5ad(input_h5ad)

    if coord_key not in adata.obsm:
        raise ValueError(f"The coordinate key '{coord_key}' was not found. Available keys: {list(adata.obsm.keys())}")

    
    adata = adata[adata.obs['total_counts'] > 0, :].copy()

      
    genes_file = args.genes_file
    if genes_file:
        if not os.path.exists(genes_file):
            raise ValueError(f"genes_file not found: {genes_file}")

        keep_genes = []
        with open(genes_file, "r", encoding="utf-8") as f:
            for line in f:
                g = line.strip()
                if g and (not g.startswith("#")):
                    keep_genes.append(g)

        keep_set = set(keep_genes)
        present = [g for g in adata.var_names if g in keep_set]
        missing = sorted(list(keep_set - set(present)))

        if len(present) == 0:
            raise ValueError(
                "None of the provided genes were found in adata.var_names "
                "(gene naming may be inconsistent, e.g., ENSG IDs vs gene symbols).")

        if missing:
            print(
                f"[WARN] {len(missing)} provided genes were not found in adata.var_names and will be ignored. "
                f"Examples: {missing[:10]}")

        adata = adata[:, present].copy()


    spatial_coords = adata.obsm[coord_key]
    scaled_coords = np.round(spatial_coords / mu_scale, 2).T
    xcoord, ycoord = np.split(scaled_coords, 2, axis=0)
    xcoord, ycoord = xcoord.flatten(), ycoord.flatten()

    xmin, xmax = xcoord.min(), xcoord.max()
    ymin, ymax = ycoord.min(), ycoord.max()

    barcodes = adata.obs.index
    _X = np.round((xcoord - xmin) * 100).astype(int)
    _Y = np.round((ycoord - ymin) * 100).astype(int)

    bc_coor = pd.DataFrame({
        'barcode': barcodes.str.split(':').str[0],
        'X_spatial': spatial_coords[:, 0],
        'Y_spatial': spatial_coords[:, 1],
        'X_scale': xcoord,
        'Y_scale': ycoord,
        'X': _X,
        'Y': _Y
    })

    bc_coor_file = os.path.join(output_dir, 'barcode.coor.match.tsv.gz')
    bc_coor.to_csv(bc_coor_file, sep='\t', index=False, compression='gzip')

    coord_file = os.path.join(output_dir, 'coordinate_minmax.tsv')
    df_range = pd.DataFrame([xmin, xmax, ymin, ymax], index=['xmin', 'xmax', 'ymin', 'ymax'])
    df_range.to_csv(coord_file, sep='\t', header=False)


    X = adata.X
    if not hasattr(X, "tocoo"):
        from scipy import sparse
        X = sparse.csr_matrix(X)

    Xcoo = X.tocoo()
    row_indices, column_indices = Xcoo.row, Xcoo.col
    gene_names = adata.var_names[column_indices].to_numpy()
    counts = Xcoo.data.astype(int)

    x_used = xcoord[row_indices]
    y_used = ycoord[row_indices]

    df_transcripts = pd.DataFrame({
        "X": x_used,
        "Y": y_used,
        "gene": gene_names,
        "Count": counts
    }).sort_values(by="Y")

    transcripts_file = os.path.join(output_dir, 'transcripts.tsv.gz')
    df_transcripts.to_csv(transcripts_file, sep='\t', index=False, compression='gzip')

    feature = df_transcripts.groupby('gene', as_index=False)['Count'].sum()

    if dummy_genes:
        feature = feature[~feature['gene'].str.contains(dummy_genes, flags=re.IGNORECASE, regex=True)]

    feature_clean_file = os.path.join(output_dir, 'feature.clean.tsv.gz')
    feature.to_csv(feature_clean_file, sep='\t', index=False, compression='gzip')

    return transcripts_file, coord_file, feature_clean_file, bc_coor_file


if __name__ == '__main__':
    try:
        transcripts, coord_file, feature_clean_file, bc_coor_file = process_h5ad_to_transcripts(sys.argv[1:])
        print(f"Transcripts: {transcripts}")
        print(f"Coordinate range: {coord_file}")
        print(f"Feature clean: {feature_clean_file}")
        print(f"Barcode file: {bc_coor_file}")
    except ValueError as e:
        print(f"Error: {e}")
