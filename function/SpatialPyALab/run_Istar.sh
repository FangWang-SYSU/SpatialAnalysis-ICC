#!/bin/bash
set -e

prefix=$1   
istar_path=/Volumes/B/Projects/2025Project/SpatialAnalysis/istar
device="cpu"  
pixel_size=0.5
n_genes=1000 
gene_path=/Volumes/B/Projects/2025Project/SpatialAnalysis/0_myanalysis/scRNA_ref/marker
# preprocess histology image
echo $pixel_size > ${prefix}/pixel-size.txt
python ${istar_path}/rescale.py ${prefix} --image
python  ${istar_path}/preprocess.py ${prefix} --image
# extract histology features
python  ${istar_path}/extract_features.py ${prefix} --device=${device}

# auto detect tissue mask
# If you have a user-defined tissue mask, put it at `${prefix}mask-raw.png` and comment out the line below
python  ${istar_path}/get_mask.py ${prefix}/embeddings-hist.pickle ${prefix}/mask-small.png

# select most highly variable genes to predict
# If you have a user-defined list of genes, put it at `${prefix}gene-names.txt` and comment out the line below
python  ${istar_path}/select_genes.py --n-top=${n_genes} "${prefix}/cnts.tsv"  "${prefix}/gene-names.txt"
python ${istar_path}/merge_hvg_deg.py "${prefix}/cnts.tsv" "${prefix}/gene-names.txt" "${gene_path}/DEG_genes.txt"

# rescale coordinates and spot radius
python ${istar_path}/rescale.py ${prefix} --locs --radius

# train gene expression prediction model and predict at super-resolution
python ${istar_path}/impute.py ${prefix} --epochs=400 --device=${device}
# # If you want to retrain model, you need to delete the existing model:
# rm -r ${prefix}states
# visualize imputed gene expression
python ${istar_path}/plot_imputed.py ${prefix}
# segment image by gene features
python ${istar_path}/cluster.py --filter-size=8 --min-cluster-size=20 --n-clusters=10 --mask=${prefix}mask-small.png ${prefix}embeddings-gene.pickle ${prefix}clusters-gene/

# differential analysis by clusters
python ${istar_path}/aggregate_imputed.py ${prefix}
python ${istar_path}/reorganize_imputed.py ${prefix}
python ${istar_path}/differential.py ${prefix}

# visualize spot-level gene expression data
python ${istar_path}/plot_spots.py ${prefix}

# python marker_score.py ${prefix} ${prefix}data/markers/M1.txt ${prefix}data/markers/M1
for marker_file in "${gene_path}/sub_markers/"*.txt; do
    module_name=$(basename "$marker_file" .txt)
    echo $marker_file    
    python ${istar_path}/marker_score.py \
        "${prefix}" \
        "$marker_file" \
        "$prefix/data/markers/$module_name"
done