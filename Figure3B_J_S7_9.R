library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(ggsci)
library(patchwork)
library(stringr)
library(future)

setwd("./")
patient_analysis_dir <- file.path(project_dir, "Data/ST/inferCNV/Patient_Cluster")
patient_plot_dir <- file.path(project_dir, "Data/ST/inferCNV/Patient_Cluster")

counts <- Matrix::readMM('Data/ST/inferCNV/RawData/create_Seurat_data/counts.mtx.gz')
counts = t(counts)
genes <- fread('Data/ST/inferCNV/RawData/create_Seurat_data/genes.tsv', header = FALSE)$V1
barcodes <- fread('Data/ST/inferCNV/RawData/create_Seurat_data/barcodes.tsv', header = FALSE)$V1
rownames(counts) <- genes
colnames(counts) <- barcodes

meta <- fread('Data/ST/inferCNV/RawData/create_Seurat_data/metadata.csv')
meta <- as.data.frame(meta)
rownames(meta) = meta$barcode
meta <- meta[colnames(counts), , drop = FALSE]

seu <- CreateSeuratObject(counts = counts, meta.data = meta)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
saveRDS(object = seu,file = 'Data/ST/inferCNV/RawData/All_hexagon_in_seurat_original.RDS')

seu_Clean = subset(seu,CellType!='conflict')
seu_Clean$PatientID <- factor(seu_Clean$PatientID,levels = paste0('Pt',1:17))

pats <- unique(as.character(seu_Clean$SampleID))
main_id <- as.integer(sub("^Pt(\\d+)-.*$", "\\1", pats))
sub_num <- as.integer(sub("^Pt\\d+-(\\d+)[a-z]?$", "\\1", pats))

suf <- ifelse(grepl("[a-z]$", pats), sub(".*([a-z])$", "\\1", pats), "")
suf_rank <- match(suf, letters)
suf_rank[is.na(suf_rank)] <- 0L
ord <- order(main_id, sub_num, suf_rank)
sample_order = pats[ord]
seu_Clean$SampleID <- factor(seu_Clean$SampleID, levels = sample_order)

saveRDS(object = seu_Clean, file = 'Data/ST/inferCNV/All_Hexagons_in_seurat_Clean.RDS')

seu_Clean = subset(seu_Clean, nCount_RNA >= 1000 & nCount_RNA <= 25000 & nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 40)
seu_sub = subset(seu_Clean,Type == 'Tumor'&CellType =='Epithelial')

options(future.globals.maxSize = 500 * 1024^3)
plan("multisession", workers = 1)

seu_half <- SCTransform(seu_sub, vars.to.regress = "percent.mt", vst.flavor = "v2", variable.features.n = 3000,ncells = 5000,
                        return.only.var.genes = TRUE,
                        conserve.memory = TRUE)
seu_half <- RunPCA(seu_half, assay = "SCT", features = VariableFeatures(seu_half, assay = "SCT"), npcs = 50)
seu_half <- FindNeighbors(seu_half, reduction = "pca", dims = 1:30, k.param = 40, graph.name = "pca_snn")
seu_half <- FindClusters(seu_half, graph.name = "pca_snn", resolution = 0.3)

seu_half <- RunUMAP(seu_half,
                    reduction = "pca",
                    dims = 1:30,
                    n.neighbors = 30,
                    min.dist = 0.2,
                    reduction.name = "umap.pca",
                    reduction.key = "UMAPPCA_")

p1 = DimPlot(seu_half,
             reduction = 'umap.pca',
             group.by = c('seurat_clusters','PatientID'), ncol = 2,
             cols = pal_igv()(30),
             label = F)
ggsave('Figure/ST/InferCNV/UMAP_Plot/1_allPatients_UMAP.pdf',p1,width = 12,height = 5)

#
CNV_info <- read.csv("Data/ST/inferCNV/inferCNV_result/anno_CNV_subcluster_smoothed.csv",header = T)
rownames(CNV_info) <- CNV_info$Hex_ID
seu$Clone <- NA_character_
seu$Clone_correct <- NA_character_
common_cells <- intersect(colnames(seu), rownames(CNV_info))
seu$Clone[common_cells] <- CNV_info[common_cells, "CNV_subcluster"]
seu$Clone_correct[common_cells] <- CNV_info[common_cells, "CNV_subcluster_smooth"]

saveRDS(seu,file = 'Data/ST/inferCNV/hexagon_36um_add_inferCNV.RDS')

patients = unique(seu$PatientID)

for(pat in patients) {
    
    sub_clean <- subset(seu, PatientID == pat & CellType == "Epithelial" & !is.na(Clone_correct))
    keep_levels <- names(which(table(sub_clean$Clone_correct) >= 10))
    sub_clean <- subset(sub_clean, subset = Clone_correct %in% keep_levels)
    samp_levels <- unique(sub_clean$SampleID)
    samp_part <- sub("^Pt\\d+-", "", samp_levels)
    samp_levels <- samp_levels[order(as.integer(sub("^([0-9]+).*$", "\\1", samp_part)), sub("^[0-9]+", "", samp_part))]
    sub_clean$SampleID <- factor(sub_clean$SampleID, levels = samp_levels)
    
    options(future.globals.maxSize = 24 * 1024^3)
    sub_clean <- SCTransform(sub_clean, vars.to.regress = "percent.mt",vst.flavor = "v2", variable.features.n = 3000, return.only.var.genes = FALSE )
    
    genes <- rownames(sub_clean)
    hvg <- VariableFeatures(sub_clean, assay = "SCT")
    rm_genes <- unique(c(genes[str_detect(genes, "^RPL|^RPS")], genes[str_detect(genes, "^IGH|^IGK|^IGL")], genes[str_detect(genes, "^MT-")]))
    hvg_use <- setdiff(hvg, rm_genes)
    VariableFeatures(sub_clean, assay = "SCT") <- hvg_use
    
    sub_clean <- RunPCA(sub_clean, assay = "SCT", npcs = 30)
    sub_clean <- FindNeighbors(sub_clean, reduction = "pca", dims = 1:20, k.param = 40)
    sub_clean <- FindClusters(sub_clean, resolution = 0.3)
    sub_clean <- RunUMAP(sub_clean, reduction = "pca", dims = 1:15, n.neighbors = 10, min.dist = 0.2)
    
    saveRDS(sub_clean,file = paste0(paste0(patient_analysis_dir, "/"),pat,'-hexagon-SCT-HighVar.Rds'))
    
}

## plotS7
plot_Sample = list()
plot_Clone = list()
cb_clone = c('A' = '#f3b0b3',"B" = '#5e8fac',"C" = "#f5d05a","D" = "#f08b1a","E" = "#238444","F" = "#82539e","G" = "#ac7447")
for(pat in patients) {
    sub_clean <- readRDS(file = paste0(paste0(patient_plot_dir, "/"), pat, "-hexagon-SCT-HighVar.Rds"))
    plot_Sample[[pat]] = DimPlot(sub_clean, group.by = c('SampleID'),size= 1, cols = pal_igv()(30), label = F)+ ggtitle(pat)
    plot_Clone[[pat]] = DimPlot(sub_clean, group.by = c('Clone_correct'), cols =  cb_clone, label = F)+ ggtitle(pat)
    
}

p_all1 <- wrap_plots(plotlist = plot_Sample, ncol = 5)

ggsave(filename = 'Figure/ST/InferCNV/UMAP_Plot/1_UMAP_sampleID-HighVar.pdf',plot = p_all1, width = 20,height = 12)
p_all2 <- wrap_plots(plotlist = c(plot_Clone[5],plot_Clone[7:14],plot_Clone[17]), ncol = 5)
ggsave(filename = 'Figure/ST/InferCNV/UMAP_Plot/2_UMAP_Clone_correct-HighVar.pdf',plot = p_all2,width = 20,height = 7.5)
 
## Figure3B
seu <- subset(seu, CellType == "Epithelial" & !is.na(Clone_correct))
meta <- seu@meta.data

idx <- meta$PatientID == "Pt4"
meta$Clone_correct[idx] <- ifelse(meta$Clone_correct[idx] == "A", "B", ifelse(meta$Clone_correct[idx] == "B", "A", meta$Clone_correct[idx]))

df_tab <- table(PatientID = meta$PatientID, Clone = meta$Clone_correct)
prop_mat <- prop.table(df_tab, margin = 1)

order_df <- data.frame(
    PatientID = rownames(df_tab),
    n_clone = rowSums(df_tab > 0),
    main_clone = colnames(prop_mat)[max.col(prop_mat, ties.method = "first")],
    main_prop = apply(prop_mat, 1, max)) %>%
    arrange(n_clone, desc(main_prop), main_clone)
 
patient_order <- order_df$PatientID
plot_df <- as.data.frame(df_tab) %>%
    mutate(PatientID = as.character(PatientID)) %>%
    mutate(PatientID = factor(PatientID, levels = patient_order))

cb_pattele = c("A" = "#FFB6B9", "B" = "#6096B4", "C" = "#FFD758", "D" =  "#FF9100", "E" = "#238b45", "F" = "#B648F1", "G" = "#B17648")

p1 <- ggplot(plot_df, aes(PatientID, Freq, fill = Clone)) +
    geom_col(position = "fill", width = 0.8) +
    #facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = cb_pattele) +
    theme_few() +
    labs(x = "Patient", y = "Proportion", fill = "Clone") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"))

p1
ggsave("Figure/ST/InferCNV/1_Patient-clone-proportion.pdf", p1, width = 5, height = 3.5)


## Figure3D
norm_shannon <- function(p, K_fixed) {
    p <- as.numeric(p)
    p <- p[!is.na(p)]
    p <- p[p > 0]
    if (length(p) == 0 || K_fixed <= 1) return(0)
    H <- -sum(p * log(p, base = 2))
    H / log(K_fixed, base = 2)
}
meta2 <- meta %>%
    dplyr::select(PatientID, SampleID,Sample_Leiden, Clone_correct) %>%
    filter(!is.na(PatientID), !is.na(SampleID), !is.na(Clone_correct)) %>%
    #filter( !PatientID =='Pt15') %>%
    mutate(across(c(PatientID, SampleID, Clone_correct), as.character)) %>%
    group_by(PatientID) %>%
    filter(n_distinct(SampleID) > 1, n_distinct(Clone_correct) > 1) %>%
    ungroup()

pat_list <- split(meta2, meta2$PatientID)

clone_disp_res <- lapply(pat_list, function(df_pat) {
    
    tab_pat <- table(df_pat$SampleID, df_pat$Clone_correct)    
    frac_pat <- prop.table(tab_pat, margin = 2)   
    
    R_patient <- nrow(tab_pat)
    K_patient <- ncol(tab_pat)
    
    clone_entropy <- apply(frac_pat, 2, norm_shannon, K_fixed = R_patient)
    clone_df <- data.frame(
        PatientID = unique(df_pat$PatientID),
        Clone = colnames(tab_pat),
        entropy = as.numeric(clone_entropy),
        weight = colSums(tab_pat),
        R_patient = R_patient,
        K_patient = K_patient,
        stringsAsFactors = FALSE
    )
    patient_df <- data.frame(
        PatientID = unique(df_pat$PatientID),
        R_patient = R_patient,
        K_patient = K_patient,
        n_clone = ncol(tab_pat),
        entropy_mean = mean(clone_df$entropy, na.rm = TRUE),
        entropy_weighted_mean = weighted.mean(clone_df$entropy, w = clone_df$weight, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    
    list(detail = clone_df, patient = patient_df)
})

patient_clone_disp_df <- bind_rows(lapply(clone_disp_res, `[[`, "patient"))
patient_cell_count_df <- meta2 %>% group_by(PatientID) %>% summarise(n_cells = n(),.groups = "drop")
patient_bubble_df <- patient_clone_disp_df %>%
    dplyr::select(PatientID, clone_dispersion = entropy_mean, clone_dispersion_weighted = entropy_weighted_mean, n_clone, R_patient,K_patient) %>%
    left_join(patient_region_mix_df %>% dplyr::select( PatientID,
                                                       clonal_mixing = entropy_mean,
                                                       clonal_mixing_weighted = entropy_weighted_mean,
                                                       n_region),by = "PatientID") %>%
    left_join(patient_cell_count_df, by = "PatientID") %>%
    mutate(pie_radius = ifelse(max(n_cells) == min(n_cells), 0.035, scales::rescale(sqrt(n_cells), to = c(0.018, 0.05))))
 
p_bubble <- ggplot(patient_bubble_df, aes(x = clone_dispersion, y = clonal_mixing)) +
    geom_vline(xintercept = median(patient_bubble_df$clone_dispersion, na.rm = TRUE), 
               linetype = "dashed", linewidth = 0.4,color = "grey60") +
    geom_hline(yintercept = 0.35, linetype = "dashed", linewidth = 0.4, color = "grey60") +
    geom_point(shape = 21,fill = '#56B4CC',size = 4)+
    geom_text(aes(label = PatientID), size = 2.7, vjust = -1.2, color = "grey20") +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_few() +
    labs(x = "Clone dispersion across regions", y = "Clonal mixing within regions", title = "") +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(10, 20, 10, 10)
    )

p_bubble
ggsave(filename = "Figure/ST/InferCNV/Statistics_New/5_Patient_level_clone_region_RNR_pie_bubble_plot.pdf",plot = p_bubble, width = 5, height = 4)



## Figure3C
umap_fig_dir <- file.path(project_dir, "Figure/ST/InferCNV/UMAP_Plot")
stacked_output_dir <- file.path(umap_fig_dir, "Region_Clone_stacked")
stacked_patients <- c("Pt5", "Pt7", "Pt8", "Pt9", "Pt10", "Pt11", "Pt12", "Pt13", "Pt14", "Pt17")
stacked_settings <- list(reduction_name = "umap", region_col = "SampleID", clone_col = "Clone_correct", plane_width = 1,
    depth_val = 0.48, shear_val = -0.35, tilt_val = 0.12, plane_gap = 0.76, upper_offset_x = -0.04, plane_u_limits = c(-0.5,
        0.5), plane_v_limits = c(-0.5, 0.5), frame_width_factor = 1.1, frame_height_factor = 1.1, frame_pad_x = 0.015,
    frame_pad_y = 0.035, view_padding = 0.05, frame_color = "#8D9297", frame_linewidth = 0.18, point_size = 1,
    max_links = 80L, link_seed = 123L, link_color = "#E2E5E9", link_linewidth = 0.08)
stacked_pdf_width <- 6
stacked_pdf_height <- 5

merge_parallel_links <- function(links) {
    if (nrow(links) == 0L) {
        return(links)
    }
    dx <- links$xend[1] - links$x[1]
    dy <- links$yend[1] - links$y[1]
    line_length <- sqrt(dx^2 + dy^2)
    if (line_length <= 1e-12) {
        return(links[0, , drop = FALSE])
    }
    ux <- dx/line_length
    uy <- dy/line_length
    nx <- -uy
    ny <- ux
    links$normal <- links$x * nx + links$y * ny
    t_start <- links$x * ux + links$y * uy
    t_end <- links$xend * ux + links$yend * uy
    links$lower <- pmin(t_start, t_end)
    links$upper <- pmax(t_start, t_end)
    groups <- split(seq_len(nrow(links)), round(links$normal, 10))
    merged <- lapply(groups, function(indices) {
        d <- links[indices, , drop = FALSE]
        d <- d[order(d$lower, d$upper), , drop = FALSE]
        normal_value <- mean(d$normal)
        current_lower <- d$lower[1]
        current_upper <- d$upper[1]
        pieces <- list()
        piece_index <- 1L
        make_segment <- function(lower, upper) {
            data.frame(x = normal_value * nx + lower * ux, y = normal_value * ny + lower * uy, xend = normal_value *
                nx + upper * ux, yend = normal_value * ny + upper * uy)
        }
        if (nrow(d) > 1L) {
            for (j in seq.int(2L, nrow(d))) {
                if (d$lower[j] <= current_upper + 1e-10) {
                  current_upper <- max(current_upper, d$upper[j])
                } else {
                  pieces[[piece_index]] <- make_segment(current_lower, current_upper)
                  piece_index <- piece_index + 1L
                  current_lower <- d$lower[j]
                  current_upper <- d$upper[j]
                }
            }
        }
        pieces[[piece_index]] <- make_segment(current_lower, current_upper)
        do.call(rbind, pieces)
    })
    result <- do.call(rbind, merged)
    rownames(result) <- NULL
    result
}

make_stacked_plot <- function(obj, patient, settings = stacked_settings) {
    umap <- Embeddings(obj, reduction = settings$reduction_name)
    df <- data.frame(cell = rownames(umap), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], Region = obj@meta.data[rownames(umap),
        settings$region_col], Clone = obj@meta.data[rownames(umap), settings$clone_col], stringsAsFactors = FALSE)
    valid <- is.finite(df$UMAP_1) & is.finite(df$UMAP_2)
    df <- df[valid, , drop = FALSE]
 
    df$Region <- as.character(df$Region)
    df$Clone <- as.character(df$Clone)
    df$Region[is.na(df$Region) | df$Region == ""] <- "Unknown"
    df$Clone[is.na(df$Clone) | df$Clone == ""] <- "Unknown"
    cx <- mean(range(df$UMAP_1))
    cy <- mean(range(df$UMAP_2))
    span <- max(diff(range(df$UMAP_1)), diff(range(df$UMAP_2)))
 
    df$u <- (df$UMAP_1 - cx)/span
    df$v <- (df$UMAP_2 - cy)/span
    project_umap <- function(u, v, offset_x = 0, offset_y = 0) {
        data.frame(x = settings$plane_width * (u + settings$shear_val * v) + offset_x, y = settings$tilt_val *
            u + settings$depth_val * v + offset_y)
    }
    make_plane <- function(offset_x = 0, offset_y = 0) {
        u_center <- mean(settings$plane_u_limits)
        v_center <- mean(settings$plane_v_limits)
        u_half <- diff(settings$plane_u_limits)/2 * settings$frame_width_factor + settings$frame_pad_x
        v_half <- diff(settings$plane_v_limits)/2 * settings$frame_height_factor + settings$frame_pad_y
        corners <- data.frame(u = c(u_center - u_half, u_center + u_half, u_center + u_half, u_center - u_half),
            v = c(v_center - v_half, v_center - v_half, v_center + v_half, v_center + v_half))
        project_umap(corners$u, corners$v, offset_x, offset_y)
    }
    region_offset_x <- settings$upper_offset_x
    region_offset_y <- settings$plane_gap/2
    clone_offset_x <- 0
    clone_offset_y <- -settings$plane_gap/2
    region_df <- cbind(df, project_umap(df$u, df$v, region_offset_x, region_offset_y))
    clone_df <- cbind(df, project_umap(df$u, df$v, clone_offset_x, clone_offset_y))
    plane_region <- make_plane(region_offset_x, region_offset_y)
    plane_clone <- make_plane(clone_offset_x, clone_offset_y)
    cell_links <- data.frame(x = clone_df$x, y = clone_df$y, xend = region_df$x, yend = region_df$y)
    set.seed(settings$link_seed)
    link_indices <- sample.int(nrow(cell_links), min(settings$max_links, nrow(cell_links)))
    cell_links <- cell_links[link_indices, , drop = FALSE]
    cell_links <- merge_parallel_links(cell_links)
    region_levels <- sort(unique(region_df$Region))
    base_region_cols <- ggsci::pal_igv()(51)
    if (length(region_levels) <= length(base_region_cols)) {
        region_colors <- base_region_cols[seq_along(region_levels)]
    } else {
        region_colors <- grDevices::colorRampPalette(base_region_cols)(length(region_levels))
    }
    region_cols <- setNames(region_colors, region_levels)
    region_cols["Unknown"] <- "#B8BDC3"
    clone_levels <- sort(unique(clone_df$Clone))
    clone_cols <- cb_clone
    extra_clones <- setdiff(clone_levels, names(clone_cols))
    if (length(extra_clones)) {
        extra_palette <- grDevices::hcl.colors(360, palette = "Dark 3")
        color_index <- vapply(extra_clones, function(label) {
            code <- utf8ToInt(label)
            as.integer(sum(code * seq_along(code))%%360) + 1L
        }, integer(1))
        clone_cols <- c(clone_cols, setNames(extra_palette[color_index], extra_clones))
    }
    clone_cols["Unknown"] <- "#B8BDC3"
    region_df$color_key <- paste0("Region::", region_df$Region)
    clone_df$color_key <- paste0("Clone::", clone_df$Clone)
    point_df <- rbind(region_df, clone_df)
    point_colors <- c(setNames(unname(region_cols), paste0("Region::", names(region_cols))), setNames(unname(clone_cols),
        paste0("Clone::", names(clone_cols))))
    if (!all(point_df$color_key %in% names(point_colors))) {
        stop("Some cell categories have no assigned color.")
    }
    all_planes <- rbind(plane_region, plane_clone)
    view_xlim <- range(all_planes$x) + c(-1, 1) * settings$view_padding
    view_ylim <- range(all_planes$y) + c(-1, 1) * settings$view_padding
    p <- ggplot() + geom_polygon(data = plane_clone, aes(x = x, y = y), fill = "white", color = NA) + geom_segment(data = cell_links,
        aes(x = x, y = y, xend = xend, yend = yend), color = settings$link_color, linewidth = settings$link_linewidth,
        alpha = 1, show.legend = FALSE) + geom_point(data = point_df, aes(x = x, y = y, color = color_key), shape = 16,
        size = settings$point_size, alpha = 1, show.legend = FALSE) + scale_color_manual(values = point_colors,
        guide = "none") + geom_polygon(data = plane_clone, aes(x = x, y = y), fill = NA, color = settings$frame_color,
        linewidth = settings$frame_linewidth) + coord_equal(xlim = view_xlim, ylim = view_ylim, expand = FALSE,
        clip = "off") + labs(title = NULL, subtitle = NULL, caption = NULL, x = NULL, y = NULL) + theme_void() +
        theme(plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white",
            color = NA), plot.margin = margin(3, 3, 3, 3, unit = "mm"), legend.position = "none")
    return(p)
}

 
patient_files <- list.files(patient_plot_dir, pattern = "-hexagon-SCT-HighVar\\.Rds$", full.names = TRUE, ignore.case = TRUE)
patient_ids <- sub("-hexagon-SCT-HighVar\\.Rds$", "", basename(patient_files), ignore.case = TRUE)
if (!is.null(stacked_patients)) {
    keep <- patient_ids %in% stacked_patients
    patient_files <- patient_files[keep]
    patient_ids <- patient_ids[keep]
}
 
patient_numbers <- suppressWarnings(as.integer(sub("^[^0-9]*", "", patient_ids)))
patient_order <- order(patient_numbers, patient_ids, na.last = TRUE)
patient_files <- patient_files[patient_order]
patient_ids <- patient_ids[patient_order]
 

export_stacked_patient <- function(input_file, patient) {
    pdf_file <- file.path(stacked_output_dir, paste0(patient, "_UMAP_no_labels_no_arrow.pdf"))
    patient_warnings <- character()
    error_message <- ""
    n_cells <- NA_integer_
    temporary_pdf <- tempfile(pattern = paste0(patient, "_"), tmpdir = stacked_output_dir, fileext = ".pdf")
    on.exit(unlink(temporary_pdf), add = TRUE)
    success <- tryCatch(withCallingHandlers({
        obj <- readRDS(input_file)
        n_cells <- ncol(obj)
        p <- make_stacked_plot(obj, patient)
        ggsave(filename = temporary_pdf, plot = p, device = grDevices::cairo_pdf, width = stacked_pdf_width, height = stacked_pdf_height,
            units = "in", bg = "white")
        TRUE
    }, warning = function(w) {
        patient_warnings <<- c(patient_warnings, conditionMessage(w))
        message("  WARNING: ", conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) {
        error_message <<- conditionMessage(e)
        FALSE
    })
    data.frame(patient = patient, cells_in_object = n_cells, status = if (success)
        "success" else "failed", pdf_file = if (success)
        pdf_file else NA_character_, warnings = paste(unique(patient_warnings), collapse = " | "), error = error_message, stringsAsFactors = FALSE)
}

results <- vector("list", length(patient_files))
for (i in seq_along(patient_files)) {
    results[[i]] <- export_stacked_patient(patient_files[i], patient_ids[i])
    invisible(gc())
}

#-----------------S9-------------------
 
library(ggplot2)
library(ggsci)
library(dplyr)
library(infercnv)
library(ggpubr)
library(Seurat)
library(tibble)
library(CytoTRACE)
library(stringr)

setwd('/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC_spatial_V2')
 
seu <- readRDS('Data/SpatialData/1_Integration/inferCNV/hexagon_36um_add_inferCNV.RDS')
info <- seu@meta.data
info$Hex_ID <- rownames(info)
cnv_all <- setNames(rep(NA_real_, ncol(seu)), colnames(seu))

out_dir <- "Figure/SpatialFigure/1_Integration/3_InferCNV/inferCNV_score"

patients <- unique(info$PatientID)
Clone_order <- list(Pt1 = c('A','B'),
                    Pt2 = c('A','B'),
                    Pt3 = c('A','B'),
                    Pt4 = c('A','B'),
                    Pt5 = c('D','F','B','E','A','C'),
                    Pt6 = c('A','B'),
                    Pt7 = c('E','A','C','G','B','F','D'),
                    Pt8 = c('B','C','A'),
                    Pt9 = c('A','B'),
                    Pt10 = c('A','B'),
                    Pt11 = c('A','B'),
                    Pt12 = c('B','A','E','C','D'),
                    Pt13 =c('A','B'),
                    Pt14 = c('B','A','C'),
                    Pt15 = c('A'),
                    Pt16 = c('A','B'),
                    Pt17 = c('E','D','C','B','A'))

cnv_all <- setNames(rep(NA_real_, ncol(seu)), colnames(seu))

for (pat in patients) {
    
    cat("Processing:", pat, "\n")
    obj_path <- file.path("/Project_ICC/Data/ST/inferCNV/patient_inferCNV/",paste0(pat, "_hex_36um"),"run.final.infercnv_obj")
    obj <- readRDS(obj_path)
    
    sub_info <- subset(info,PatientID == pat  & !is.na(info$Clone_correct))
    
    if (nrow(sub_info) == 0) {
        cat("  -> no sub_info, skip\n")
        next
    }
    
    expr_mat <- as.matrix(obj@expr.data[,sub_info$Hex_ID])
    expr_shift <- expr_mat - 1
    inferCNV_score <- colMeans(abs(expr_shift), na.rm = TRUE)
    
    cnv_all[names(inferCNV_score)] <- inferCNV_score
    
    plot_df <- data.frame(HexID = names(inferCNV_score), CNV_score = inferCNV_score, stringsAsFactors = FALSE)
    plot_df$Clone <- sub_info[plot_df$HexID, "Clone_correct"]
    plot_df <- plot_df %>% filter(!is.na(Clone))
    
    clone_count <- plot_df %>% count(Clone, name = "n")
    
    plot_df <- plot_df %>% semi_join(clone_count %>% filter(n >= 10), by = "Clone")
    n_clone <- length(unique(plot_df$Clone))
    
    plot_df2 <- plot_df %>% group_by(Clone) %>% mutate(mean_CNV = mean(CNV_score, na.rm = TRUE)) %>% ungroup()
    clone_l <- Clone_order[[pat]]
    plot_df2$Clone <- factor(plot_df2$Clone, levels = rev(clone_l))
    
    library(ggthemes)
    if (n_clone >= 2) {
        p <- ggplot(plot_df2, aes(x = Clone, y = CNV_score, fill = mean_CNV)) +
            geom_violin(trim = FALSE, scale = "width", width = 0.7) +
            geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
            scale_fill_gradientn(colours = c("#fff7f3","#fde0dd","#fcc5c0","#fa9fb5","#f768a1","#dd3497")) +
            coord_flip() +
            theme_few() +
            labs(x = "", y = "inferCNV", fill = "Mean CNV")
        
        ggsave(filename = file.path(out_dir, paste0(pat, "_inferCNV_score_by_clone.pdf")),plot = p,width = 2.8,height = n_clone*0.3 + 1.2)
    }
    
}

cnv_all = as.data.frame(cnv_all)
cnv_all_df <- data.frame(Hex = rownames(cnv_all),cnv_score = cnv_all$cnv_all)
write.table(cnv_all_df,file = 'Data/SpatialData/1_Integration/inferCNV/inferCNV_result/cnv_score.tsv',sep = '\t')


#-------------Cytotrace------------
base_outdir <- "Data/SpatialData/1_Integration/CytoTRACE/Sample/"
cyto_all <- setNames(rep(NA_real_, ncol(seu)), colnames(seu))

for (pat in patients) {

    message(">>> Processing patient: ", pat)
    cells_use <- rownames(seu@meta.data)[seu$PatientID %in% pat & !is.na(seu$Clone_correct)]
    
    if (length(cells_use) == 0) {
        message("No cells found for ", pat, ", skip.")
        next
    }
    obj_pat <- subset(seu, cells = cells_use)
    options(future.globals.maxSize = 24 * 1024^3)   
    obj_pat = SCTransform(obj_pat, vars.to.regress = "percent.mt",vst.flavor = "v2", variable.features.n = 2000, return.only.var.genes = FALSE )
    genes <- rownames(obj_pat)
    hvg <- VariableFeatures(obj_pat, assay = "SCT")
    rm_genes <- unique(c(
        genes[str_detect(genes, "^RPL|^RPS")],   
        genes[str_detect(genes, "^IGH|^IGK|^IGL")],   
        genes[str_detect(genes, "^MT-")]
    ))
    hvg_use <- setdiff(hvg, rm_genes)
    
    VariableFeatures(obj_pat, assay = "SCT") <- hvg_use
    obj_pat <- RunPCA(obj_pat, assay = "SCT", npcs = 30)
    obj_pat <- FindNeighbors(obj_pat, reduction = "pca", dims = 1:20, k.param = 40)
    obj_pat <- FindClusters(obj_pat, resolution = 0.3)
    obj_pat <- RunUMAP(obj_pat, reduction = "pca", dims = 1:15, n.neighbors = 10, min.dist = 0.1)
    
    outdir <- file.path(base_outdir, pat)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    count_mat <- as.matrix(GetAssayData(obj_pat, layer = "counts", assay = "RNA"))
    ct_res <- CytoTRACE(count_mat, ncores = 16, subsamplesize = 1000)
    
    cyto_vec <- ct_res$CytoTRACE
    names(cyto_vec) <- gsub("\\.", "-", names(cyto_vec))
    
    cyto_df <- data.frame(CytoTRACE = cyto_vec)
    
    common_cells <- intersect(rownames(cyto_df), colnames(obj_pat))
    cyto_df <- cyto_df[common_cells, , drop = FALSE]
    obj_pat <- AddMetaData(obj_pat, cyto_df)
    
    cyto_all[common_cells] <- cyto_df[common_cells, "CytoTRACE"]
    
    meta_sub <- obj_pat@meta.data 
    n_niche <- meta_sub %>% filter(!is.na(Clone_correct)) %>% pull(Clone_correct) %>% unique() %>% length()
    
    
    phenotype_all <- obj_pat@meta.data$Clone_correct
    names(phenotype_all) <- colnames(obj_pat)
    plotCytoTRACE(ct_res, phenotype = phenotype_all, emb = obj_pat@reductions$umap@cell.embeddings, outputDir = paste0(outdir,'/'))
    
    saveRDS(obj_pat, file = paste0(outdir, "/with_CytoTRACE.rds")) 
}

cyto_all = as.data.frame(cyto_all)
cyto_all_df <- data.frame(Hex = rownames(cyto_all),cyto_score = cyto_all$cyto_all)
head(cyto_all_df)

write.table(cyto_all_df,file = 'Data/SpatialData/1_Integration/inferCNV/inferCNV_result/cytoTrace_score.tsv',sep = '\t')
 
  
for (pat in patients) {
    #pat = patients[1]
    print(pat)
    outdir <- file.path(base_outdir, pat)
    obj_pat <- readRDS(paste0(outdir, "/with_CytoTRACE.rds"))
    meta_sub <- obj_pat@meta.data
    meta_sub <- meta_sub %>% filter(!is.na(Clone_correct), !is.na(CytoTRACE))
    
    head(meta_sub)
    
    plot_df <- meta_sub[,c('barcode','Clone_correct','CytoTRACE')]
    colnames(plot_df) <- c('Hex_ID','Clone','CytoTRACE')
    
    plot_df <- plot_df %>% filter(!is.na(Clone))
    clone_count <- plot_df %>% count(Clone, name = "n")
    plot_df <- plot_df %>% semi_join(clone_count %>% filter(n >= 10), by = "Clone")
    
    n_clone <- length(unique(plot_df$Clone))
    n_clone
    
    plot_df2 <- plot_df %>% group_by(Clone) %>% mutate(mean_CytoTRACE = mean(CytoTRACE, na.rm = TRUE)) %>% ungroup()
    clone_l <- Clone_order[[pat]]
    plot_df2$Clone <- factor(plot_df2$Clone, levels = rev(clone_l))
    
    
    library(ggthemes)
    if (n_clone >= 2) {
        p <- ggplot(plot_df2, aes(x = Clone, y = CytoTRACE, fill = mean_CytoTRACE)) +
            geom_violin(trim = FALSE, scale = "width", width = 0.7) +
            geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
            scale_fill_gradientn(colours = c("#edf8b1","#7fcdbb","#1d91c0")) +
            coord_flip() +
            theme_few() +
            labs(x = "", y = "CytoTRACE", fill = "CytoTRACE")
        
        p
        ggsave(filename = paste0("Figure/SpatialFigure/1_Integration/3_InferCNV/CytoTRACE_Score/", pat, "_CytoTRACE.pdf"), plot = p, width = 2.8, height = n_clone*0.3 + 1.2)
        
    }
}
 
coord_path <- "Data/SpatialData/1_Integration/hex_36um/merge_sample_files/coordinates.tsv"
cyto_path  <- "Data/SpatialData/1_Integration/inferCNV/inferCNV_result/cytoTrace_score.tsv"
cnv_path  <- "Data/SpatialData/1_Integration/inferCNV/inferCNV_result/cnv_score.tsv"

coords <- read.table(coord_path, sep = "\t",header = T)
colnames(coords) <- c('Hex','x','y','SampleID','PatientID')
cyto <- read.table(cyto_path, sep = "\t",header = T)
cnv <- read.table(cnv_path, sep = "\t",header = T)

plot_df <- coords %>% left_join(cyto, by = "Hex")  %>% left_join(cnv, by = "Hex")
make_hex_df <- function(df, r = 16, angle0 = pi / 6) {
    angles <- angle0 + seq(0, 5) * pi / 3
    idx <- rep(seq_len(nrow(df)), each = 6)
    hex_df <- df[idx, , drop = FALSE]
    hex_df$.hex_id <- idx
    hex_df$.vertex <- rep(1:6, times = nrow(df))
    hex_df$x_hex <- hex_df$x + rep(r * cos(angles), times = nrow(df))
    hex_df$y_hex <- hex_df$y + rep(r * sin(angles), times = nrow(df))
    hex_df
}

cyto_list <- list()
cnv_list <- list()
for (pid in patients) {
    
    df_sub <- plot_df %>% filter(PatientID == pid)
    hex_df <- make_hex_df(df_sub, r = 30)
    q_low  <- quantile(hex_df$cnv_score, 0.05, na.rm = TRUE)
    q_high <- quantile(hex_df$cnv_score, 0.95, na.rm = TRUE)
    hex_df <- hex_df %>%mutate(cnv_clip = pmin(pmax(cnv_score, q_low), q_high),cnv_scaled = rescale(cnv_clip, to = c(0, 1)))
    cyto_list[[pid]] <- ggplot(hex_df, aes(x = x_hex, y = y_hex, group = .hex_id, fill = cyto_score)) +
        geom_polygon(color = NA, alpha = 0.95) +
        scale_fill_viridis_c(option = "mako", na.value = "#f0f0f0") +
        facet_wrap(~ SampleID, nrow = 1) +
        coord_fixed() +
        theme_void() +
        scale_y_reverse() +
        theme(
            strip.text = element_text(size = 10, face = "bold"),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
        labs(fill = "cyto_score")  
    
    cnv_list[[pid]] <- ggplot(hex_df, aes(x = x_hex, y = y_hex, group = .hex_id, fill = cnv_score)) +
        geom_polygon(color = NA, alpha =1) +
        scale_fill_viridis_c(option = "turbo", na.value = "#f0f0f0") +
        facet_wrap(~ SampleID, nrow = 1) +
        coord_fixed() +
        theme_void() +
        scale_y_reverse() +
        theme(
            strip.text = element_text(size = 10, face = "bold"),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        ) +
        labs(fill = "CNV score")  
}
for (pid in patients) {
    p = cyto_list[[pid]]
    ggsave(paste0('Figure/SpatialFigure/1_Integration/3_InferCNV/CytoTRACE_Score/',pid,'_cyto_spatial.pdf'),plot = p)
    
}

for (pid in patients) {
    p = cnv_list[[pid]]
    ggsave(paste0('Figure/SpatialFigure/1_Integration/3_InferCNV/inferCNV_score/',pid,'_cnv_spatial.pdf'),plot = p)
    
}

##-Figure-Heatmap
max_match_distance_um <- 36
clone_order <- list(
    Pt1 = c('A', 'B'),
    Pt2 = c('A', 'B'),
    Pt3 = c('A', 'B'),
    Pt4 = c('A', 'B'),
    Pt5 = c('D', 'F', 'B', 'E', 'A', 'C'),
    Pt6 = c('A', 'B'),
    Pt7 = c('E', 'A', 'C', 'G', 'B', 'F', 'D'),
    Pt8 = c('B', 'C', 'A'),
    Pt9 = c('A', 'B'),
    Pt10 = c('A', 'B'),
    Pt11 = c('A', 'B'),
    Pt12 = c('B', 'A', 'E', 'C', 'D'),
    Pt13 = c('A', 'B'),
    Pt14 = c('B', 'A', 'C'),
    Pt15 = 'A',
    Pt16 = c('A', 'B'),
    Pt17 = c('E', 'D', 'C', 'B', 'A')
)

clone_colors <- c(A = '#FFB6B9', B = '#6096B4', C = '#FFD758', D = '#FF9100', E = '#238b45', F = '#B648F1', G = '#B17648')
sep_colors <- c(M_SEP1 = '#e41a1c', M_SEP2 = '#2373F4', M_SEP3 = '#FCAD38', M_SEP4 = '#b3de69', M_SEP5 = '#A290B7')
cnv_colors <- c('0' = '#2166ac', '1' = '#92c5de', '2' = '#f0f0f0', '3' = '#fddbc7', '4' = '#f46d43', '6' = '#b2182b')
status_colors <- c(Mixed = '#999999')

Clone_obj <- readRDS('Data/ST/inferCNV/hexagon_36um_add_inferCNV.RDS')
Clone_meta <- Clone_obj@meta.data

SEP_obj <- readRDS('Data/ST/SEP/All_hex_SEP_CelltypeConstrained.RDS')
SEP_meta <- SEP_obj@meta.data

clean_label <- function(x) {
    x <- trimws(as.character(x))
    x[tolower(x) %in% c('', 'na', 'nan', 'none', 'unassigned')] <- NA_character_
    x
}
prepare_metadata <- function(metadata, columns, label) {
    if (is.data.table(metadata) && 'HexID' %in% names(metadata)) {
        metadata <- copy(metadata)
    } else {
        metadata <- as.data.frame(metadata)
        metadata$HexID <- rownames(metadata)
        metadata <- as.data.table(metadata)
    }
    metadata[, (columns) := lapply(.SD, clean_label), .SDcols = columns]
    metadata
}

parse_accepted <- function(x) {
    x <- tolower(trimws(as.character(x)))
    allowed <- c('true', 'false', '1', '0', 't', 'f', '', 'na', 'nan')
    x %in% c('true', '1', 't')
}

clone_for_display <- function(patient_id, clone_label) {
    labels <- as.character(clone_label)
    swap <- which(patient_id %in% 'Pt4' & labels %in% c('A', 'B'))
    labels[swap] <- unname(c(A = 'B', B = 'A')[labels[swap]])
    labels
}

save_heatmap <- function(heatmap, filename, width, height, extra_legends = list()) {
    render <- function() {
        draw(
            heatmap, heatmap_legend_side = 'right', annotation_legend_side = 'right',
            annotation_legend_list = extra_legends, merge_legends = TRUE
        )
    }
    
    pdf(paste0('Figure/ST/InferCNV/Heatmap_domSEP/',filename, '.pdf'), width = width, height = height)
    tryCatch(render(), finally = dev.off())
}
Clone_meta <- prepare_metadata(Clone_meta, c('HexID', 'PatientID', 'SampleID', 'Clone', 'Clone_correct'), 'Clone_meta')
SEP_meta <- prepare_metadata(SEP_meta, c('HexID', 'SampleID', 'sample_id', 'CellType', 'SEP', 'Clone'), 'SEP_meta')
SEP_meta[, SEP := gsub('-', '_', SEP, fixed = TRUE)]
observed_msep <- unique(SEP_meta[grepl('^M_SEP', SEP), SEP])
SEP_meta[, Include_in_MSEP_composition := SEP %in% names(sep_colors)]

SEP_meta[, `:=`(
    upstream_match_accepted = parse_accepted(hex36_match_accepted),
    Hex36_candidate = clean_label(matched_hex36_id),
    Clone_inherited = clean_label(Clone),
    distance_um = suppressWarnings(as.numeric(as.character(dist_to_hex36)))
)]
SEP_meta[, distance_accepted := is.finite(distance_um) & distance_um >= 0 &
             distance_um <= max_match_distance_um]
SEP_meta[, match_accepted := upstream_match_accepted & distance_accepted]

sample_prefix <- paste0(SEP_meta$sample_id, '__')
target_index <- match(SEP_meta$Hex36_candidate, Clone_meta$HexID)
strip_prefix <- which(is.na(target_index) & !is.na(SEP_meta$Hex36_candidate) & startsWith(SEP_meta$Hex36_candidate, sample_prefix))
SEP_meta[strip_prefix, Hex36_candidate := substring(Hex36_candidate, nchar(sample_prefix[strip_prefix]) + 1L)]
target_index <- match(SEP_meta$Hex36_candidate, Clone_meta$HexID)
resolved <- which(SEP_meta$upstream_match_accepted & !is.na(target_index))

SEP_meta[, `:=`(
    TargetHexID = Clone_meta$HexID[target_index],
    TargetClone = Clone_meta$Clone_correct[target_index],
    TargetClone_display = clone_for_display(
        Clone_meta$PatientID[target_index], Clone_meta$Clone_correct[target_index]
    )
)]

SEP_meta[, Clone_agrees := fifelse(
    !is.na(Clone_inherited) & !is.na(TargetClone), Clone_inherited == TargetClone, NA
)]

SEP_meta[, Mapping_status := fcase(
    !upstream_match_accepted, 'Rejected upstream',
    !is.finite(distance_um) | distance_um < 0, 'Invalid distance',
    !distance_accepted, 'Rejected by stricter distance',
    is.na(Hex36_candidate), 'Missing target ID',
    is.na(TargetHexID), 'Target absent from clone file',
    is.na(TargetClone), 'Target has no corrected clone',
    default = 'Mapped to corrected clone'
)]


mapping_status <- SEP_meta[, .N, by = .(SampleID, Mapping_status)]
distance_qc <- SEP_meta[, .(
    Max_match_distance_um = max_match_distance_um,
    N_accepted_upstream = sum(upstream_match_accepted),
    N_accepted_current = sum(match_accepted),
    N_MSEP_removed_by_distance = sum(
        upstream_match_accepted & !distance_accepted & Include_in_MSEP_composition
    )
), by = SampleID]

prepare_patient_cnv <- function(patient, clone_metadata) {
    cell_meta <- copy(clone_metadata[PatientID == patient & !is.na(Clone_correct)])
    if (!nrow(cell_meta)) return(NULL)
    
    hmm_file <- file.path('Data/ST/inferCNV/patient_inferCNV', paste0(patient, '_hex_36um'), '20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj')
    
    hmm_object <- readRDS(hmm_file)
    cnv_matrix <- hmm_object@expr.data
    hex_ids <- cell_meta$HexID[cell_meta$HexID %in% colnames(cnv_matrix)]
    
    cell_meta <- cell_meta[match(hex_ids, HexID)]
    cnv_matrix <- cnv_matrix[, hex_ids, drop = FALSE]
    
    clone_groups <- split(seq_along(hex_ids), cell_meta$Clone)
    for (label in names(clone_groups)) {
        indices <- clone_groups[[label]]
        differs <- vapply(indices, function(i) {
            any(cnv_matrix[, i] != cnv_matrix[, indices[1]])
        }, logical(1))
    }
    
    representatives <- vapply(clone_groups, function(indices) indices[1], integer(1))
    reference_matrix <- cnv_matrix[, representatives, drop = FALSE]
    colnames(reference_matrix) <- names(representatives)
    
    corrected_matrix <- reference_matrix[, cell_meta$Clone_correct, drop = FALSE]
    colnames(corrected_matrix) <- hex_ids
    
    cell_meta[, Clone_correct_original := Clone_correct]
    cell_meta[, Clone_correct := clone_for_display(PatientID, Clone_correct_original)]
    clone_levels <- clone_order[[patient]]
    
    sample_levels <- unique(cell_meta$SampleID)
    sample_part <- sub('^Pt[0-9]+-', '', sample_levels)
    sample_number <- suppressWarnings(as.integer(sub('^([0-9]+).*$', '\\1', sample_part)))
    sample_levels <- sample_levels[order(sample_number, sample_part)]
    row_order <- order(
        match(cell_meta$Clone_correct, clone_levels),
        match(cell_meta$SampleID, sample_levels), cell_meta$HexID
    )
    cell_meta <- cell_meta[row_order]
    
    gene_order <- hmm_object@gene_order
    genes <- rownames(gene_order)[rownames(gene_order) %in% rownames(corrected_matrix)]
    
    chromosome_labels <- sub('^chr', '', as.character(gene_order[genes, 'chr']))
    chromosome <- factor(chromosome_labels, levels = unique(chromosome_labels))
    heatmap_matrix <- t(corrected_matrix[genes, cell_meta$HexID, drop = FALSE] * 2)
    
    
    list(metadata = cell_meta, matrix = heatmap_matrix, chromosome = chromosome,
         clone_levels = clone_levels, sample_levels = sample_levels)
}


summarize_patient_sep <- function(patient, cell_meta, sep_metadata) {
    cell_meta <- copy(cell_meta)
    hex_ids <- cell_meta$HexID
    mapped_upstream <- sep_metadata[upstream_match_accepted & TargetHexID %in% hex_ids]
    
    mapped <- copy(mapped_upstream[match_accepted == TRUE])
    
    sep_counts <- matrix(
        0L, length(hex_ids), length(sep_colors),
        dimnames = list(hex_ids, names(sep_colors))
    )
    counts_by_sep <- mapped[Include_in_MSEP_composition == TRUE, .(n = .N), by = .(TargetHexID, SEP)]
    if (nrow(counts_by_sep)) {
        count_index <- cbind(
            match(counts_by_sep$TargetHexID, hex_ids), match(counts_by_sep$SEP, names(sep_colors))
        )
        sep_counts[count_index] <- counts_by_sep$n
    }
    n_msep <- rowSums(sep_counts)
    sep_fraction <- sep_counts / pmax(n_msep, 1)
    
    support <- mapped[, {
        msep_distance <- distance_um[Include_in_MSEP_composition]
        list(
            N_hex15_matched_in_input = .N,
            N_hex15_NSEP_excluded = sum(grepl('^N_SEP', SEP)),
            Median_MSEP_distance_um = if (length(msep_distance)) median(msep_distance) else NA_real_,
            Max_MSEP_distance_um = if (length(msep_distance)) max(msep_distance) else NA_real_,
            N_clone_disagreements = sum(Clone_agrees %in% FALSE)
        )
    }, by = TargetHexID]
    distance_removed <- mapped_upstream[
        !distance_accepted & Include_in_MSEP_composition,
        .(N_MSEP_removed_by_distance = .N), by = TargetHexID
    ]
    
    support <- support[match(hex_ids, TargetHexID)]
    distance_removed <- distance_removed[match(hex_ids, TargetHexID)]
    count_columns <- c('N_hex15_matched_in_input', 'N_hex15_NSEP_excluded', 'N_clone_disagreements')
    for (column in count_columns) {
        set(support, which(is.na(support[[column]])), column, 0L)
    }
    distance_removed[is.na(N_MSEP_removed_by_distance), N_MSEP_removed_by_distance := 0L]
    
    top_category <- max.col(sep_fraction, ties.method = 'first')
    dominant_fraction <- sep_fraction[cbind(seq_along(hex_ids), top_category)]
    tied <- rowSums(abs(sweep(sep_fraction, 1, dominant_fraction, '-')) < 1e-8) > 1
    dominant_sep <- colnames(sep_fraction)[top_category]
    dominant_sep[tied | n_msep == 0] <- 'Mixed'
    assignment_reason <- fcase(n_msep == 0, 'No M-SEP', tied, 'Tied M-SEP proportions', n_msep == 1, 'Dominant from one M-SEP hex', default = 'Dominant M-SEP')
    dominant_fraction[n_msep == 0] <- NA_real_
    n_labelled <- n_msep + support$N_hex15_NSEP_excluded
    
    cell_meta[, `:=`(
        Dom_SEP = dominant_sep,
        Dominant_fraction = dominant_fraction,
        SEP_assignment_reason = assignment_reason,
        N_hex15_matched_in_input = support$N_hex15_matched_in_input,
        N_hex15_MSEP = n_msep,
        N_hex15_NSEP_excluded = support$N_hex15_NSEP_excluded,
        MSEP_fraction_of_labelled_hex15 = n_msep / ifelse(n_labelled > 0, n_labelled, NA_real_),
        Median_MSEP_distance_um = support$Median_MSEP_distance_um,
        Max_MSEP_distance_um = support$Max_MSEP_distance_um,
        N_clone_disagreements = support$N_clone_disagreements,
        N_MSEP_removed_by_distance = distance_removed$N_MSEP_removed_by_distance
    )]
    mapped[, Heatmap_row := match(TargetHexID, hex_ids)]
    setorder(mapped, Heatmap_row, HexID)
    
    list(metadata = cell_meta, counts = sep_counts, fractions = sep_fraction, mapping = mapped)
}


make_row_annotation <- function(cell_meta, clone_levels, sample_levels) {
    clone <- factor(cell_meta$Clone_correct, levels = clone_levels)
    sample <- factor(cell_meta$SampleID, levels = sample_levels)
    dominant_levels <- c(names(sep_colors), names(status_colors))
    dominant <- factor(cell_meta$Dom_SEP, levels = dominant_levels[dominant_levels %in% cell_meta$Dom_SEP])
    sample_colors <- setNames(hcl.colors(length(sample_levels), 'Dark 3'), sample_levels)
    dominant_colors <- c(sep_colors, status_colors)
    
    rowAnnotation(
        Clone = clone, Sample = sample, `Dom SEP` = dominant,
        col = list(Clone = clone_colors[clone_levels], Sample = sample_colors,
                   `Dom SEP` = dominant_colors[levels(dominant)]),
        simple_anno_size = unit(4, 'mm'), annotation_name_gp = gpar(fontsize = 9),
        annotation_legend_param = list(
            Clone = list(title = 'Clone'), Sample = list(title = 'Samples'),
            `Dom SEP` = list(title = 'Dominant M-SEP\n(15 um-derived)')
        )
    )
}

plot_patient_heatmaps <- function(patient, cnv, sep) {
    cell_meta <- sep$metadata
    clone_levels <- cnv$clone_levels[cnv$clone_levels %in% cell_meta$Clone_correct]
    clone <- factor(cell_meta$Clone_correct, levels = clone_levels)
    row_annotation <- function() make_row_annotation(cell_meta, clone_levels, cnv$sample_levels)
    chromosome <- cnv$chromosome
    chromosome_text <- levels(chromosome)
    chromosome_text[as.integer(table(chromosome)) < 200] <- ''
    chromosome_annotation <- HeatmapAnnotation(
        chr = anno_block(gp = gpar(fill = rep(c('#f0f0f0', '#bababa'), length.out = nlevels(chromosome))), labels = chromosome_text, labels_gp = gpar(fontsize = 10)),
        show_annotation_name = FALSE)
    cnv_heatmap <- Heatmap(
        cnv$matrix, name = 'HMM', col = cnv_colors,
        cluster_rows = FALSE, cluster_row_slices = FALSE,
        cluster_columns = FALSE, cluster_column_slices = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        row_split = clone, row_gap = unit(0.8, 'mm'), row_title = NULL,
        column_split = chromosome, column_gap = unit(1.1, 'mm'),
        top_annotation = chromosome_annotation, left_annotation = row_annotation(),
        column_title = patient, use_raster = TRUE,
        heatmap_legend_param = list(at = c(0, 1, 2, 3, 4, 6))
    )
    save_heatmap(cnv_heatmap, paste0(patient, '_CNV_domSEP'), width = 13, height = 5)
    
    composition <- sep$fractions
    present <- colSums(composition) > 0
    if (!any(present)) present[] <- TRUE
    composition <- composition[, present, drop = FALSE]
    composition_colors <- sep_colors[colnames(composition)]
    composition_annotation <- rowAnnotation(
        `M-SEP composition` = anno_barplot(
            composition, beside = FALSE, border = FALSE, bar_width = 1,
            gp = gpar(fill = unname(composition_colors), col = NA),
            ylim = c(0, 1), extend = 0, width = unit(85, 'mm'),
            axis_param = list(side = 'top', at = c(0, .5, 1), labels = c('0%', '50%', '100%'))
        ),
        annotation_name_gp = gpar(fontsize = 9)
    )
    
    anchor <- matrix(0, nrow(cell_meta), 1, dimnames = list(cell_meta$HexID, ''))
    composition_heatmap <- Heatmap(
        anchor, name = 'anchor', col = c('0' = 'white'),
        rect_gp = gpar(type = 'none'), width = unit(0.1, 'mm'), show_heatmap_legend = FALSE,
        cluster_rows = FALSE, cluster_row_slices = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, show_column_names = FALSE,
        row_split = clone, row_gap = unit(0.8, 'mm'), row_title = NULL,
        left_annotation = row_annotation(), right_annotation = composition_annotation,
        column_title = paste0(patient, ' | M-SEP composition per 36 um hex')
    )
    composition_legend <- Legend(title = 'SEP', labels = names(composition_colors), legend_gp = gpar(fill = unname(composition_colors)), ncol = 2)
    save_heatmap(composition_heatmap, paste0(patient, '_SEP_composition_per_hex'), width = 12, height = 8, extra_legends = list(composition_legend))
}


for (patient in patients) {
    message('Patient: ', patient)
    cnv <- prepare_patient_cnv(patient, Clone_meta)
    sep <- summarize_patient_sep(patient, cnv$metadata, SEP_meta)
    plot_patient_heatmaps(patient, cnv, sep)
}


## Figure_Tree
library(ggalluvial)   
source('./function/tree.R')
cb_pattle <- pal_igv()(50)
S_colors <- setNames(cb_pattle, LETTERS[seq_along(cb_pattle)])
S_colors["unassigned"] <- "#BDBDBD"
 
for (patientId in patients) {
    
    sub_obj = subset(Clone_obj,PatientID == patientId)
    final_clone <- names(which(table(sub_obj$Clone_correct) >= 10))
    min_cells <- 10
    datapath <- file.path(root_dir, "Data/ST/inferCNV/patient_inferCNV/")
    fig_dir <- file.path(datapath, paste0(patientId, "_hex_36um")) 
    
    obj <- readRDS(paste0(datapath, patientId,'_hex_36um/20_HMM_pred.repr_intensitiesHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.infercnv_obj'))
    ord_list <- run_infercnv_subcluster_tree(obj = obj,
                                             min_cells = min_cells,
                                             sample_id = patientId, 
                                             make_tree_min_k = 1,
                                             final_clone = final_clone,
                                             fig_dir = "Figure/ST/InferCNV/Tree-Plot/",
                                             width = 6,
                                             height = 4,
                                             cb_pattle = S_colors)
}

 
#  
Clone_obj <- subset(Clone_obj, subset = !is.na(Clone_correct) & Clone_correct != "NA")
for(patientId in patients){
    datapath <- file.path(root_dir, "Data/SpatialData/1_Integration/inferCNV/subclusters_V2/")
    seu_f <- subset(Clone_obj,PatientID == patientId)
    Clone <- as.character(seu_f$Clone_correct)
    Sample <- as.character(seu_f$SampleID)

    samp_levels <- unique(Sample)
    samp_part <- sub("^Pt\\d+-", "", samp_levels)
    samp_levels <- samp_levels[order(
        as.integer(sub("^([0-9]+).*$", "\\1", samp_part)),
        sub("^[0-9]+", "", samp_part)
    )]
    Sample <- factor(Sample, levels = samp_levels)
    Clone <- factor(Clone, levels = Clone_order[[patientId]])
    ord <- order(Clone)
    Clone <- Clone[ord]
    Sample <- Sample[ord]
 
    df = data.frame(Clone = Clone,Sample = Sample)
    tab <- df %>% count(Clone, Sample, name = "n") %>% group_by(Clone) %>%
        mutate(prop = n / sum(n),lab  = paste0(Sample, "\n", scales::percent(prop, accuracy = 0.1))) %>%
        ungroup()
    
    p_pie <- tab %>%
        ggplot(aes(x = 1, y = prop, fill = Sample)) +
        geom_col(width = 1, color = "white", linewidth = 0.2) +
        coord_polar(theta = "y") +
        facet_wrap(~ Clone,ncol = 1) +
        scale_fill_manual(values = pal_sample) +
        theme_void() +
        theme(strip.text = element_text(size = 12),legend.position = "right")
    
    ggsave(filename = paste0('Figure/SpatialFigure/1_Integration/3_InferCNV/Pie-correct_Patient/',patientId,'_pie_subcluster_sample.pdf'),
           plot = p_pie ,width = 2,height = 6 )
}



library(data.table)
library(networkD3)
library(htmlwidgets)
library(jsonlite)


root <- '/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC'
input_dir <- file.path(root, 'Figure/ST/InferCNV/Heatmap_domSEP')
out_dir <- file.path(root, 'Figure/ST/InferCNV/Patient_Clone_DomSEP_Sankey_NoMixed')

patients <- paste0('Pt', 1:17) 
show_plots <- interactive()
save_pdf <- TRUE   
save_svg <- TRUE  
plot_width <- 700    
plot_height <- 550
font_size <- 13
font_family <- 'Arial'  
link_color <- '#BDBDBD'
node_width <- 28
node_padding <- 14

read_heatmap_rows <- function(patient, path) {
    rows <- fread(path)
    required <- c('HexID', 'PatientID', 'SampleID', 'Clone_correct', 'Dom_SEP',
                  'SEP_assignment_reason', 'Heatmap_row')
    
    label_columns <- setdiff(required, 'Heatmap_row')
    rows[, (label_columns) := lapply(.SD, as.character), .SDcols = label_columns]
    levels_clone <- clone_order[[patient]]
    setorder(rows, Heatmap_row)
    rows[, ..required]
}

stopifnot(length(patients) > 0L, !anyNA(patients), !anyDuplicated(patients))
input_files <- file.path(input_dir, paste0(patients, '_SEP_composition_per_hex.csv'))

heatmap_rows <- rbindlist(Map(read_heatmap_rows, patients, input_files))

plot_rows <- heatmap_rows[Dom_SEP != 'Mixed']
links <- plot_rows[, .(N_hex36 = .N), by = .(PatientID, Clone_correct, Dom_SEP)]
links[, Fraction_within_patient := N_hex36 / sum(N_hex36), by = PatientID]
links[, Fraction_within_clone := N_hex36 / sum(N_hex36), by = .(PatientID, Clone_correct)]
links[, Fraction_within_Dom_SEP := N_hex36 / sum(N_hex36), by = .(PatientID, Dom_SEP)]


row_counts <- heatmap_rows[, .(
    N_heatmap_rows = .N,
    N_Mixed_excluded = sum(Dom_SEP == 'Mixed'),
    N_heatmap_rows_retained = sum(Dom_SEP != 'Mixed')
), by = .(PatientID, Clone_correct)]
link_counts <- links[, .(N_sankey_hex36 = sum(N_hex36)), by = .(PatientID, Clone_correct)]
count_check <- merge(row_counts, link_counts, by = c('PatientID', 'Clone_correct'), all = TRUE)
count_check[is.na(N_sankey_hex36), N_sankey_hex36 := 0L]


make_patient_sankey <- function(patient, patient_links) {
    levels_clone <- clone_order[[patient]]
    if (is.null(levels_clone)) levels_clone <- names(clone_colors)
    present_clones <- levels_clone[levels_clone %in% patient_links$Clone_correct]
    present_seps <- names(sep_colors)[names(sep_colors) %in% patient_links$Dom_SEP]
    
    nodes <- rbindlist(list(
        data.table(id = paste0('CL_', present_clones), label = present_clones),
        data.table(id = paste0('SEP_', present_seps),
                   label = gsub('_', '-', present_seps, fixed = TRUE))
    ))
    nodes[, color_key := id]
    
    plot_links <- patient_links[, .(
        source = match(paste0('CL_', Clone_correct), nodes$id) - 1L,
        target = match(paste0('SEP_', Dom_SEP), nodes$id) - 1L,
        value = N_hex36,
        link_group = 'link_grey'
    )]
    
    setorder(plot_links, source, target)
    
    node_colors <- c(
        setNames(clone_colors, paste0('CL_', names(clone_colors))),
        setNames(sep_colors, paste0('SEP_', names(sep_colors))),
        link_grey = link_color
    )
    colour_scale <- htmlwidgets::JS(paste0(
        'd3.scaleOrdinal().domain(',
        jsonlite::toJSON(names(node_colors), auto_unbox = FALSE),
        ').range(',
        jsonlite::toJSON(unname(node_colors), auto_unbox = FALSE), ')'
    ))
    
    networkD3::sankeyNetwork(
        Links = as.data.frame(plot_links), Nodes = as.data.frame(nodes),
        Source = 'source', Target = 'target', Value = 'value',
        NodeID = 'label', NodeGroup = 'color_key', LinkGroup = 'link_group',
        colourScale = colour_scale,
        fontSize = font_size, fontFamily = font_family,
        nodeWidth = node_width, nodePadding = node_padding,
        width = plot_width, height = plot_height,
        iterations = 0, sinksRight = FALSE
    )
}

export_editable_svg <- function(html_file, svg_file) {
    browser <- chromote::ChromoteSession$new(width = plot_width, height = plot_height)
    on.exit(browser$close(), add = TRUE)
    page_loaded <- browser$Page$loadEventFired(wait_ = FALSE)
    local_path <- normalizePath(html_file, winslash = '/', mustWork = TRUE)
    browser$Page$navigate(paste0('file://', utils::URLencode(local_path)), wait_ = FALSE)
    browser$wait_for(page_loaded)
    
    svg_result <- browser$Runtime$evaluate(
        expression = r"(
      (async function() {
        const started = Date.now();
        let svg;
        while (true) {
          svg = document.querySelector('.sankeyNetwork svg');
          if (svg && svg.querySelector('.node text') && svg.querySelector('.link')) break;
          if (Date.now() - started > 10000) throw new Error('Sankey SVG did not render');
          await new Promise(resolve => setTimeout(resolve, 100));
        }
        await document.fonts.ready;
        await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));

        const clone = svg.cloneNode(true);
        const sourceElements = [svg, ...svg.querySelectorAll('*')];
        const targetElements = [clone, ...clone.querySelectorAll('*')];
        const properties = [
          'fill', 'fill-opacity', 'stroke', 'stroke-width', 'stroke-opacity',
          'stroke-linecap', 'stroke-linejoin', 'opacity', 'font-family',
          'font-size', 'font-weight', 'font-style', 'text-anchor',
          'dominant-baseline', 'letter-spacing'
        ];
        sourceElements.forEach((element, index) => {
          const style = getComputedStyle(element);
          properties.forEach(property => {
            const value = style.getPropertyValue(property);
            if (value) {
              targetElements[index].setAttribute(property, value);
              targetElements[index].style.setProperty(property, value);
            }
          });
        }); 
 
        const bounds = svg.getBBox();
        const left = Math.min(0, bounds.x - 10);
        const top = Math.min(0, bounds.y - 10);
        const width = Math.max(svg.clientWidth, bounds.x + bounds.width + 10) - left;
        const height = Math.max(svg.clientHeight, bounds.y + bounds.height + 10) - top;
        clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
        clone.setAttribute('width', width);
        clone.setAttribute('height', height);
        clone.setAttribute('viewBox', [left, top, width, height].join(' '));
        clone.querySelectorAll('script, style, title').forEach(element => element.remove());
        if (!clone.querySelector('text') || clone.querySelector('image')) {
          throw new Error('SVG must contain editable text and no raster images');
        }
        return new XMLSerializer().serializeToString(clone);
      })()
    )",
    returnByValue = TRUE, awaitPromise = TRUE
    )
    
    svg_text <- svg_result$result$value
    writeLines(enc2utf8(svg_text), svg_file, useBytes = TRUE)
}


can_save_pdf <- save_pdf && requireNamespace('webshot2', quietly = TRUE) 
if ((save_svg || can_save_pdf) && !nzchar(Sys.getenv('CHROMOTE_CHROME'))) {
    chrome_candidates <- c(
        '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
        '/Applications/Google Chrome 2.app/Contents/MacOS/Google Chrome'
    )
    available_chrome <- chrome_candidates[file.exists(chrome_candidates)]
    if (length(available_chrome)) Sys.setenv(CHROMOTE_CHROME = available_chrome[1])
}

for (patient in patients) {
    patient_links <- links[PatientID == patient]
    
    plot <- make_patient_sankey(patient, patient_links)
    
    filename <- file.path(out_dir, paste0(patient, '_Clone_DomSEP_Sankey'))
    out_html <- paste0(filename, '.html')
    htmlwidgets::saveWidget(plot, out_html, selfcontained = FALSE, title = paste0(patient, ' | Clone - Dominant SEP'))
    if (save_svg) {
        export_editable_svg(out_html, paste0(filename, '.svg'))
    }
    if (can_save_pdf) {
        tryCatch(
            webshot2::webshot(
                url = out_html, file = paste0(filename, '.pdf'),
                vwidth = plot_width, vheight = plot_height, zoom = 1.5, delay = 1
            ),
            error = function(e) warning(patient, 'error', conditionMessage(e))
        )
    }
}
