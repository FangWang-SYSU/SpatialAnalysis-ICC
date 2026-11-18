# Snakefile for FICTURE analysis
# FICTURE analysis.smk


configfile: "scripts/snakemake/config.yaml"

REQUIRED = ["h5ad_file", "root_dir"]
missing = [k for k in REQUIRED if k not in config]
if missing:
    raise ValueError(f"缺少必需参数: {missing}，请用 --config k=v 传入")


# === 可选项：这里才给默认 ===
config.setdefault("coord_key",  "spatial")
config.setdefault("figure_bg",  "000000")
config.setdefault("mu_scale",   1.0)
config.setdefault("pseudo",     False)
config.setdefault("column_name","TopGene_pval")
config.setdefault("nFactor",6)
config.setdefault("dummy_genes", r"RP|MT-|ENSG|LINC|HLA-|IGK|IGH")
config.setdefault("pixels_cut", 0.6)
config.setdefault("genes_file", "")

  

config['model_id'] = f"nF{config['nFactor']}.d_{config['hex_width']}"
config['anchor_info'] = f"prj_{config['fit_width']}.r_{config['anchor_res']}"
config['radius'] = str(config['anchor_res'] + 1)
config['fit_nmove'] = int(config['fit_width']/config['anchor_res'])

include: "function/ficture_path.py"
include: "function/smk_function.py"
include: "function/factor_report_change.py"
include: "function/remove_duplicates.py"

# ignore warnings
os.environ["PYTHONWARNINGS"] = "ignore"


rule run_ficture:
    input:
        expand(fit_result,root_dir= config['root_dir']),
        expand(coarse_plot,root_dir = config['root_dir']),
        expand(anchor_projection,root_dir = config['root_dir']),
        expand(debulk_file,root_dir = config['root_dir']),
        expand(reportHtml,root_dir = config['root_dir']),
        expand(pixel_sort_decode_unique,root_dir = config['root_dir']),
        expand(pixel_plot,root_dir = config['root_dir']),
        expand(single_plots,root_dir = config['root_dir']),
        expand(final_barcode_coor,root_dir = config['root_dir'])
            
## coor/mu_scale
rule process_h5ad:
    input:
        h5ad = config['h5ad_file']
    output:
        transcripts = transcripts,   
        coord_file = coord_file,   
        feature_file = feature_file,
        bc_coor_file = bc_coor_file
    params:
        outpath = raw_path,
        coord_key = config['coord_key'],   
        dummy_genes = config["dummy_genes"],  
        mu_scale = config["mu_scale"],
        genes_file = config["genes_file"]

    shell:
        "python ./scripts/snakemake/function/process_h5ad_to_transcripts.py "
        "--input_h5ad {input.h5ad} "
        "--output_dir {params.outpath} "
        "--coord_key {params.coord_key} "
        "--dummy_genes '{params.dummy_genes}' "
        "--mu_scale {params.mu_scale} "
        "--genes_file {params.genes_file}"


## create batched_matrix
rule make_spatial_minibatch:
    input:
        transcripts
    output:
        temp(batched_matrix)
    params:
        batch_buff = config['batch_buff'],
        batch_size = config['batch_size'],
        major_axis = config['major_axis']
    shell:
        "ficture make_spatial_minibatch "
        "--input {input} "
        "--output {output} "
        "--mu_scale 1 "
        "--batch_size {params.batch_size} "
        "--batch_buff {params.batch_buff} "
        "--major_axis {params.major_axis}"


rule sort_spatial_minibatch:
    input:
        batched_matrix
    output:
        batched_matrix_gz
    shell:
        "sort -S 4G -k2,2n -k1,1g {input} | gzip -c > {output}"

## create hexagon_matrix
rule make_dge:
    input:
        transcripts
    output:
        temp(hexagon_matrix)
    params:
        key = config['key'],
        hex_width = config['hex_width'],
        sliding_step = config['sliding_step'],
        min_ct_per_unit = config['min_ct_per_unit'],
        major_axis = config['major_axis']
    shell:
        "ficture make_dge "
        "--input {input} "
        "--output {output} "
        "--key {params.key} "
        "--count_header {params.key} "
        "--hex_width {params.hex_width} "
        "--n_move {params.sliding_step} "
        "--min_ct_per_unit {params.min_ct_per_unit} "
        "--mu_scale 1 "
        "--precision 2 "
        "--major_axis {params.major_axis}"
#


rule sort_dge:
    input:
        hexagon_matrix
    output:
        hexagon_matrix_gz
    shell:
        "sort -S 4G -k1,1n {input} | gzip -c > {output}"
#

 

## train model
rule fit_model:
    input:
        hexagons = hexagon_matrix_gz,
        feature = feature_file
    output:
        model = model,
        fit_posterior = fit_posterior,
        fit_result = fit_result
    params:
        pseudo = bool(config['pseudo']),
        pseudo_matrix = pseudo_file,
        model_path = output_path,
        model_prefix = fit_prefix,
        nFactor = config["nFactor"],
        train_nEpoch = config["train_nEpoch"],
        key = config["key"],
        min_ct_per_feature = config["min_ct_per_feature"],
        epoch_id_length = config['epoch_id_length'],
        test_split = config['test_split'],
        R = config["R"],
    threads: 32
    run:
        shell("[ -d {params.model_path} ] || mkdir -p {params.model_path}")
        if params.pseudo:
            cmd = "cp " + f"{params.model_prefix}.posterior.count.tsv.gz " + f"{params.model_prefix}.model_matrix.tsv.gz"
            shell(
                "ficture init_model_from_pseudobulk "
                "--input {input.hexagons} "
                "--output {params.model_prefix} "
                "--feature {input.feature} "
                "--model {params.pseudo_matrix} "
                "--scale_model_rel -1 "
                "--epoch 0 "
                "--key {params.key} "
                "--min_ct_per_feature {params.min_ct_per_feature} "
                "--thread {threads};"
                f"{cmd};"
            )
        else:
            shell(
                "ficture fit_model "
                "--input {input.hexagons} "
                "--output {params.model_prefix} "
                "--feature {input.feature} "
                "--nFactor {params.nFactor} "
                "--epoch {params.train_nEpoch} "
                "--epoch_id_length {params.epoch_id_length} "
                "--unit_attr X Y "
                "--key {params.key} "
                "--min_ct_per_feature {params.min_ct_per_feature} "
                "--test_split {params.test_split} "
                "--R {params.R} "
                "--thread {threads}"
            )


rule choose_color:
    input:
        model_fit = fit_result
    output:
        cmap = cmap
    params:
        pseudo = bool(config['pseudo']),
        custom_palette = bool(config['palette']),
        palette = palette_file,
        ct_anno = ct_anno,
        color_prefix = color_prefix,
        cmap_name = config['cmap_name']
    run:
        cmd = ""
        if params.pseudo and os.path.exists(f"{params.ct_anno}"):
            cmd += " --annotation {params.ct_anno}"
            if params.custom_palette and os.path.exists(f"{params.palette}"):
                cmd += " --palette {params.palette}"
        shell(
            "ficture choose_color "
            "--input {input.model_fit} "
            "--output {params.color_prefix} "
            "--cmap_name {params.cmap_name} "
            f"{cmd}"
            )



rule plot_base:
    input:
        cmap = cmap,
        fit_result = fit_result
    output:
        coarse_plot
    params:
        coarse_prefix = coarse_prefix,
        fillr = config['hex_width']/2 + 1, # um
    shell:
        "ficture plot_base "
        "--input {input.fit_result} "
        "--output {params.coarse_prefix} "
        "--fill_range {params.fillr} "
        "--color_table {input.cmap} "
        "--plot_um_per_pixel 1 "
        "--plot_discretized"



rule transform:
    input:
        model = model,
        pixel = transcripts
    output:
        anchor_projection
    params:
        anchor_prefix = anchor_prefix,
        key = config["key"],
        major_axis = config['major_axis'],
        fit_width = config['fit_width'],
        fit_nmove = config['fit_nmove'],
        min_ct_per_fit = config['min_ct_per_fit'],
    threads: 4
    shell:
        "ficture transform "
        "--input {input.pixel} "
        "--output_pref {params.anchor_prefix} "
        "--model {input.model} "
        "--key {params.key} "
        "--major_axis {params.major_axis} "
        "--hex_width {params.fit_width} "
        "--n_move {params.fit_nmove} "
        "--min_ct_per_unit {params.min_ct_per_fit} "
        "--mu_scale 1 "
        "--thread {threads} "
        "--precision 2"


rule slda_decode:
    input:
        model = model,
        batched_matrix = batched_matrix_gz,
        anchor = anchor_projection
    output:
        pixel_decode,
        pixel_posterior
    params:
        final_prefix = final_prefix,
        topk = config['topk'],
        radius = config['radius'],
        key = config['key']
    threads: 12
    shell:
        "ficture slda_decode "
        "--input {input.batched_matrix} "
        "--output {params.final_prefix} "
        "--model {input.model} "
        "--anchor {input.anchor} "
        "--anchor_in_um "
        "--neighbor_radius {params.radius} "
        "--mu_scale 1 "
        "--key {params.key} "
        "--precision 0.01 "
        "--lite_topk_output_pixel {params.topk} "
        "--lite_topk_output_anchor {params.topk} "
        "--thread {threads}"



rule sort_pixel:
    input:
        pixel_decode = pixel_decode,
        coor = coord_file
    output:
        pixel_sort_decode
    params:
        K = config['nFactor'],
        topk = config['topk'],
        bsize = config['bsize'],
        scale = config['scale']
    run:
        sort_pixel_fun(coor = input.coor, pixel_decode = input.pixel_decode, K = params.K, 
                       topk = params.topk, bsize = params.bsize, scale = params.scale, output = output)

        shell("tabix -f -s1 -b3 -e3 {output}")


rule de_bulk:
    input:
        pixel = pixel_posterior
    output:
        debulk_file
    params:
        min_ct_per_feature = config['min_ct_per_feature'],
        max_pval_output = config['max_pval_output'],
        min_fold_output = config['min_fold_output'],
    threads: 4
    shell:
        "ficture de_bulk "
        "--input {input.pixel} "
        "--output {output} "
        "--min_ct_per_feature {params.min_ct_per_feature} "
        "--max_pval_output {params.max_pval_output} "
        "--min_fold_output {params.min_fold_output} "
        "--thread {threads}"




rule factor_report:
    input:
        cmap = cmap,
        posterior = pixel_posterior,
        debulk = debulk_file
    output:
        reportHtml = reportHtml
    params:
        pixel_path = output_path,
        prefix = outid,
        report_file = report_file,
        pseudo = bool(config['pseudo']),
        ct_anno = ct_anno,
        ct_marker = ct_marker,
        column_name=config['column_name']
    run:
        cmd = ""
        if os.path.exists(f"{params.ct_anno}"):
            cmd += " --annotation {params.ct_anno}"
        shell(
            "ficture factor_report "
            "--path {params.pixel_path} "
            "--pref {params.prefix} "
            "--color_table {input.cmap} "
            "--n_top_gene 30 "
            f"{cmd}"
            )
        if params.pseudo and os.path.exists(f"{params.ct_marker}"):
            factor_report_change(['--factor_info', params.report_file, 
                                  '--ref_gene', params.ct_marker, 
                                  '--outpath', params.pixel_path, 
                                  '--pref', params.prefix,
                                  '--column_name', params.column_name])

rule remove_duplicates:
    input:
        pixel = pixel_sort_decode
    output:
        pixel_sort_decode_unique
    run:
        remove_duplicates(input = input.pixel, output = output)


rule filter_pixels_by_posterior:
    input:
        pixel_sort_decode_unique
    output:
        pixel_sort_decode_filter
    params:
        pixels_cut = config['pixels_cut']
    run:
        cmd = ""
        shell(
            "python ./scripts/snakemake/function/filter_pixels_by_posterior.py "
            "--input {input} "
            "--output {output} "
            "--pcut {params.pixels_cut}" 
            f"{cmd}"
            )
        

rule plot_pixel_full:
    input:
        cmap = cmap,
        pixel = pixel_sort_decode_filter
    output:
        full = pixel_plot,
        plot_data = plot_tsv_gz,
    params:
        figure_bg=config['figure_bg']
    shell:
        "ficture plot_pixel_full "
        "--input {input.pixel} "
        "--color_table {input.cmap} "
        "--output {output.full} "
        "--plot_um_per_pixel 1.3 "
        "--background {params.figure_bg} "
        "--plot_top"

rule plot_pixel_single:
    input:
        cmap = cmap,
        pixel = pixel_sort_decode_filter,
    output:
        single_plots
    params:
        output_prefix = single_prefix,
        figure_bg=config['figure_bg']
    shell:
        "ficture plot_pixel_single "
        "--input {input.pixel} "
        "--color_table {input.cmap} "
        "--output {params.output_prefix} "
        "--background {params.figure_bg} "
        "--plot_um_per_pixel 1.5 "
        "--all" 

rule obtain_barcodes:
    input:
        coor_barcode = bc_coor_file,
        cluster = pixel_sort_decode_filter,
    output:
        final_barcode_coor
    params:
        ct_anno = ct_anno,
        barcode_path = barcode_path,
    run:
        cmd = ""
        if os.path.exists(f"{params.ct_anno}"):
            cmd += " --anno_path {params.ct_anno}"
        shell(
            "python ./scripts/snakemake/function/obtain_barcodes.py "
            "--coor_barcode {input.coor_barcode} "
            "--ficture_cluster {input.cluster} "
            "--output_dir {params.barcode_path} "
            f"{cmd}"
            )