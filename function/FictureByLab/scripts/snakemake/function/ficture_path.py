# dir structure

# reference file
ref_path = config.get("reference_path", "")

pseudo_file = os.path.join(ref_path, "pseudobulk.tsv")
ct_marker = os.path.join(ref_path, "MarkerGene.tsv")
ct_anno = os.path.join(ref_path, "annotation.tsv")
palette_file = os.path.join(ref_path, "palette.tsv")
 
root_dir = config['root_dir']

# Generate raw data
raw_path = root_dir+"/raw_data"
transcripts = raw_path + "/transcripts.tsv.gz"
coord_file = raw_path + "/coordinate_minmax.tsv"
feature_file = raw_path + "/feature.clean.tsv.gz"
bc_coor_file = raw_path + "/barcode.coor.match.tsv.gz"

# Generate process data
process_path = root_dir + "/process_data"
batched_matrix = process_path + "/batched.matrix.tsv"
batched_matrix_gz = batched_matrix + ".gz"
hexagon_matrix = process_path + "/" + f"hexagon.d_{config['hex_width']}.tsv"
hexagon_matrix_gz = hexagon_matrix + ".gz"


# output dir
output_path = root_dir + "/analysis/"  
figure_path = output_path + "figure/"
DE_path = output_path + "DE/"


# model fit
fit_prefix = output_path + f"{config['model_id']}"

if bool(config['pseudo']):
    model = fit_prefix + ".model_matrix.tsv.gz"
else:
    model = fit_prefix + ".model.p"

fit_posterior = fit_prefix + ".posterior.count.tsv.gz"
fit_result = fit_prefix + ".fit_result.tsv.gz"


# output file
## figure
color_prefix = figure_path + f"{config['model_id']}"
cmap = color_prefix + ".rgb.tsv"
coarse_prefix = color_prefix + ".coarse"
coarse_plot = coarse_prefix + ".png"


# anchor file
anchor_prefix = output_path + f"{config['model_id']}.{config['anchor_info']}"
anchor_projection = anchor_prefix + ".fit_result.tsv.gz"

outid = f"{config['model_id']}.decode.{config['anchor_info']}_{config['radius']}"

# pixel file
final_prefix = output_path + outid
pixel_decode = final_prefix + ".pixel.tsv.gz"
pixel_posterior = final_prefix + ".posterior.count.tsv.gz"
pixel_sort_decode = final_prefix + ".pixel.sorted.tsv.gz"
pixel_sort_decode_unique = final_prefix + ".pixel.sorted.unique.tsv.gz"
pixel_sort_decode_filter = final_prefix + ".pixel.sorted.filter.tsv.gz"

## debulk
debulk_file = DE_path + f"{config['model_id']}.decode.{config['anchor_info']}_{config['radius']}.bulk_chisq.tsv"

## report
report_file = final_prefix + ".factor.info.tsv"
reportHtml = final_prefix + ".factor.info.html"


pixel_plot = figure_path + outid + ".pixel.png"
plot_tsv_gz = figure_path + outid + ".plot.tsv.gz"
align_tsv_gz = final_prefix + ".align.tsv.gz"
align_coord_log = final_prefix + ".align.txt"


single_prefix = figure_path + "sub/" + outid
single_plots = [single_prefix + f".F_{k}.png" for k in range(config['nFactor'])]

barcode_path = output_path+'barcode'
final_barcode_coor = barcode_path + '/match_barcode_cluster.tsv.gz'
