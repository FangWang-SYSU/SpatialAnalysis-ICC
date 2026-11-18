library(cytosignal)
library(dplyr)
library(ggplot2) 
library(ggrepel)

setwd('')
info_df <- read.csv("Table/Patient_Info.csv", stringsAsFactors = FALSE)
samples <- unique(info_df$SampleName)
data_dir <- "./Data/ST/Samples"
paths <- file.path(data_dir, samples, "cytoSignal", "cytosignal.rds")
keep <- file.exists(paths)

samples <- samples[keep]
paths <- paths[keep]

cs_list <- lapply(seq_along(paths), function(i) {
    message(sprintf(i, length(paths), paths[i]))
    readRDS(paths[i])
})
names(cs_list) <- samples

dataset.meta <- as.data.frame(info_df[,c("SampleName","SampleID","Type")])
rownames(dataset.meta) <- dataset.meta$SampleName
dataset.meta$SampleName <- factor(dataset.meta$SampleName,levels = unique(dataset.meta$SampleName))

multics <- mergeCytoSignal(cs_list, metadata = dataset.meta,name.by = 'SampleName')
multics$clusters <- gsub('_','-',multics$clusters)
multics
saveRDS(multics,file = 'Data/ST/CytoSignal/Final/All_Cyto_Sample_mergeData.rds')


sample_class <- read.csv("Table/ST/SEP_res/two_component_90_classification.csv") 
sample_class <- sample_class[,c("SampleID","Group")]
multics = readRDS('Data/ST/CytoSignal/Final/All_Cyto_Sample_mergeData.rds')

sample_idx <- match(multics@metadata$SampleID, sample_class$SampleID)
multics@metadata$Sample_class <- sample_class$Group[sample_idx]

multics_pt <- subset(x = multics, multics$Sample_class %in% c("R-like", "NR-like") & multics$clusters == "Epithelial", drop = TRUE)

multics_pt <- runNEBULA(multics_pt, covariates = c('Sample_class'), cpc_thresh = 0.001, ncore = 12)
saveRDS(multics_pt,file = 'Data/ST/CytoSignal/Final/multics_epithelial_two_component.RDS')
multics_pt = readRDS('Data/ST/CytoSignal/Final/multics_epithelial_two_component.RDS')

diff_Niche <- multics_pt@diff.results$`Sample_class_R-like`
df <- diff_Niche
df2 <- df %>%
    mutate(padj_plot = ifelse(is.na(padj) | padj <= 0, 1e-300, padj),
           neglog10_padj = -log10(padj_plot),
           sig = case_when(padj < 0.05 & logFC >= 0.5 ~ "R-like",
                           padj < 0.05 & logFC <=  -0.5 ~ "NR-like", 
                           TRUE ~ "Not significant"))

df2$neglog10_padj = ifelse(df2$neglog10_padj>=15,15,df2$neglog10_padj)

df2$neglog10_padj <- pmin(df2$neglog10_padj, 15)
group_levels <- c("NR-like", "R-like")

diff_use <- df2 %>%
    transmute(interaction = trimws(gsub("_HUMAN", "", interaction)),
              selected_group = trimws(sig), logFC, p, padj) %>%
    filter(selected_group %in% group_levels, !is.na(interaction),
           interaction != "", is.finite(padj), is.finite(logFC)) %>%
    mutate(pair_key = toupper(interaction)) %>%
    arrange(padj, interaction) %>%
    distinct(selected_group, pair_key, .keep_all = TRUE)

top_lr <- diff_use %>%
    group_by(selected_group) %>%
    arrange(padj, interaction, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(ligand = sub("-.*$", "", interaction),
           receptor = sub("^[^-]+-", "", interaction))

write.csv(top_lr, "Data/ST/CytoSignal/Final/top10_LR_by_group.csv", row.names = FALSE)

label_df <- top_lr %>%
    mutate(sig = selected_group,
           neglog10_padj = pmin(-log10(pmax(padj, 1e-300)), 15))

p <- ggplot(df2, aes(x = logFC, y = neglog10_padj)) +
    geom_point(size = 1.5, alpha = 0.8,aes(color = sig)) +
    geom_vline(xintercept = c(-0.5, 0.5),linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05),linetype = "dashed", color = "black") +
    geom_text_repel(data = label_df,
                    aes(label = interaction),
                    size = 3,box.padding = 0.35,
                    point.padding = 0.2, segment.color = "black", max.overlaps = Inf) +
    scale_color_manual(values = c("R-like" = "#fdb462", "NR-like" = "#80b1d3","Not significant" = "grey70")) +
    labs(title = "", 
         x = "log2FC",
         y = expression(-log[10](FDR)),
         color = "") +
    theme_classic() +
    xlim(-6,6)+
    theme(plot.title = element_text(hjust = 0.5, size = 16))
ggsave(filename = 'Figure/ST/CytoSignal/Final/Figure6A_DEGs_LR.pdf',width = 7,height = 5,plot = p)

top_lr <- read.csv("Data/ST/CytoSignal/Final/top10_LR_by_group.csv", stringsAsFactors = FALSE)
top_ligand <- top_lr %>%
    count(selected_group, ligand, name = "n_LR") %>%
    group_by(selected_group) %>%
    arrange(desc(n_LR), ligand, .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()

selected_lr <- top_lr %>%
    inner_join(select(top_ligand, selected_group, ligand),
               by = c("selected_group", "ligand")) %>%
    select(any_of(c("pair_key", "interaction", "ligand", "receptor", "selected_group", "logFC", "p", "padj"))) %>%
    distinct(pair_key, .keep_all = TRUE)

write.csv(selected_lr, "Data/ST/CytoSignal/Final/selected_LR.csv", row.names = FALSE)
sample_map <- read.csv( "Table/Patient_Info.csv") %>%
    transmute(SampleName = trimws(SampleName), SampleID = trimws(SampleID)) %>%
    distinct()
class_map <- read.csv("Table/ST/SEP_res/two_component_90_classification.csv", stringsAsFactors = FALSE) %>%
    transmute(SampleID = trimws(SampleID), Sample_class = trimws(Group)) %>%
    distinct()
group_levels <- c("NR-like", "R-like")
sample_use <- sample_map %>%
    left_join(class_map, by = "SampleID") %>%
    filter(Sample_class %in% group_levels) %>%
    mutate(rds_file = file.path("Data/ST/Samples", SampleName, "cytoSignal", "cytosignal_add_LR_score.rds"))

group_colors <- c("NR-like" = "#80B1D3", "R-like" = "#FDB462")
slot_use <- "diffusion-contact_smooth"
epithelial_label <- "Epithelial"
high_fraction <- 0.5
min_epithelial_positions <- 5
min_samples_per_group <- 2

expected_keys <- crossing(SampleName = sample_use$SampleName, pair_key = selected_lr$pair_key)

score_cache <- tibble(
    SampleName = character(), pair_key = character(), slot_use = character(),
    mean_score = double(), status = character()
)

missing_keys <- expected_keys %>% anti_join(score_cache, by = c("SampleName", "pair_key"))
missing_samples <- unique(missing_keys$SampleName)
new_scores <- vector("list", length(missing_samples))

for (i in seq_along(missing_samples)) {
    sample_name <- missing_samples[i]
    pair_needed <- missing_keys %>%
        filter(SampleName == sample_name) %>%
        pull(pair_key)
    rds_file <- sample_use$rds_file[match(sample_name, sample_use$SampleName)]
    message("[", i, "/", length(missing_samples), "] ", sample_name)
    
    out_i <- tibble(
        SampleName = sample_name, pair_key = pair_needed, slot_use = slot_use,
        mean_score = NA_real_, status = "not_processed"
    )
    if (!file.exists(rds_file)) {
        out_i$status <- "missing_rds"
        new_scores[[i]] <- out_i
        next
    }
    
    cs <- readRDS(rds_file)
    
    score <- cs@lrscore[[slot_use]]@score
    epithelial_ids <- intersect(
        names(cs@clusters)[!is.na(cs@clusters) & cs@clusters == epithelial_label],
        rownames(score)
    )
    if (length(epithelial_ids) < min_epithelial_positions) {
        out_i$status <- "too_few_epithelial_positions"
        new_scores[[i]] <- out_i
        rm(cs, score)
        next
    }
    
    intr_index <- as.data.frame(cs@intr.valid$intr.index, stringsAsFactors = FALSE)
    needed_cols <- c("id_cp_interaction", "protein_name_a", "protein_name_b",
                     "partner_a", "partner_b")
    stopifnot(all(needed_cols %in% names(intr_index)))
    
    pair_map <- intr_index %>%
        transmute(
            intr_id = as.character(id_cp_interaction),
            ligand = coalesce(na_if(trimws(as.character(protein_name_a)), ""),
                              as.character(partner_a)),
            receptor = coalesce(na_if(trimws(as.character(protein_name_b)), ""),
                                as.character(partner_b))
        ) %>%
        mutate(
            ligand = trimws(gsub("_HUMAN", "", ligand)),
            receptor = trimws(gsub("_HUMAN", "", receptor)),
            pair_key = toupper(paste(ligand, receptor, sep = "-"))
        ) %>%
        filter(pair_key %in% pair_needed, intr_id %in% colnames(score)) %>%
        distinct(pair_key, intr_id)
    
    out_i$status <- "lr_not_matched"
    if (nrow(pair_map) > 0) {
        score_epi <- score[epithelial_ids, pair_map$intr_id, drop = FALSE]
        invalid_n <- as.integer(Matrix::colSums(
            !is.finite(score_epi) | score_epi < 0, na.rm = TRUE
        ))
        hit <- match(pair_map$pair_key, out_i$pair_key)
        out_i$mean_score[hit] <- as.numeric(Matrix::colMeans(score_epi))
        out_i$status[hit] <- ifelse(invalid_n == 0, "ok", "invalid_score")
    }
    
    new_scores[[i]] <- out_i
    rm(cs, score)
}

score_cache <- bind_rows(score_cache, bind_rows(new_scores)) %>%
    semi_join(expected_keys, by = c("SampleName", "pair_key")) %>%
    distinct(SampleName, pair_key, .keep_all = TRUE)

write.csv(score_cache, "Data/ST/CytoSignal/Final/sample_selected_LR_mean_score.csv", row.names = FALSE)
valid_scores <- score_cache %>%
    filter(status == "ok", is.finite(mean_score)) %>%
    left_join(select(sample_use, SampleName, Sample_class), by = "SampleName") %>%
    left_join(select(selected_lr, pair_key, interaction), by = "pair_key")

ranked_scores <- valid_scores %>%
    group_by(pair_key) %>%
    arrange(desc(mean_score), SampleName, .by_group = TRUE) %>%
    mutate(
        n_NR = sum(Sample_class == "NR-like"),
        n_R = sum(Sample_class == "R-like"),
        cutoff_score = mean_score[ceiling(n() * high_fraction)],
        is_high = n_NR >= min_samples_per_group &
            n_R >= min_samples_per_group &
            n_distinct(mean_score) > 1 &
            mean_score >= cutoff_score & mean_score > 0
    ) %>%
    ungroup()

plot_df <- ranked_scores %>%
    filter(is_high) %>%
    count(pair_key, interaction, Sample_class, name = "n_high") %>%
    complete(nesting(pair_key, interaction), Sample_class = group_levels,
             fill = list(n_high = 0L)) %>%
    group_by(pair_key) %>%
    mutate(fraction = n_high / sum(n_high)) %>%
    ungroup()

lr_order <- plot_df %>%
    filter(Sample_class == "NR-like") %>%
    arrange(fraction, interaction) %>%
    pull(pair_key)

plot_df <- plot_df %>%
    mutate(
        pair_key = factor(pair_key, levels = lr_order),
        Sample_class = factor(Sample_class, levels = group_levels)
    )

lr_labels <- setNames(selected_lr$interaction, selected_lr$pair_key)
p_fraction <- ggplot(plot_df, aes(x = fraction, y = pair_key, fill = Sample_class)) +
    geom_col(width = 0.7, color = "black", position = position_stack(reverse = TRUE)) +
    scale_y_discrete(labels = lr_labels) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(x = "Fraction (%)", y = "Ligand-receptor pair", fill = "Group") +
    theme_classic(base_size = 11)

ggsave("/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC/Figure/ST/CytoSignal/Final/Figure6B_LR_pair.pdf", p_fraction, width = 5, height = 4)


celltype_levels <- c("Bcell", "Tcell", "Myeloid", "Hepatocyte", "Endothelial", "Fibroblast", "Epithelial", "Unknown")
celltype_colors <- c(Bcell = "#BA55D3", Tcell = "#FF0000", Myeloid = "#FFD500", 
                     Hepatocyte = "#639D63", Endothelial = "#1E90FF",
                     Fibroblast = "#BD8F8F", Epithelial = "#9FA6AB", Unknown = "grey80")
sample_dirs <- list.dirs("Data/ST/Samples", recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[file.exists(file.path(sample_dirs, "cytoSignal", "cytosignal_add_LR_score.rds"))]
result_list <- vector("list", length(sample_dirs))
normalize_lr <- function(x) toupper(trimws(gsub("_HUMAN", "", x)))

selected_lr <- read.csv("Data/ST/CytoSignal/Final/selected_LR.csv", stringsAsFactors = FALSE) %>%
    transmute(pair_key = normalize_lr(pair_key), interaction) %>%
    distinct(pair_key, .keep_all = TRUE)
for (i in seq_along(sample_dirs)) {
    sample_path <- sample_dirs[i]
    sample_id <- basename(sample_path)
    message("[", i, "/", length(sample_dirs), "] ", sample_id)
    
    cs <- readRDS(file.path(sample_path, "cytoSignal", "cytosignal_add_LR_score.rds"))
    
    all_intr <- showIntr(
        cs, slot.use = "diffusion-contact_smooth", signif.use = "result.spx", return.name = TRUE
    )
    lr_map <- tibble(
        intr_id = names(all_intr),
        pair_key = normalize_lr(as.character(all_intr))
    ) %>%
        filter(pair_key %in% selected_lr$pair_key) %>%
        distinct(pair_key, intr_id)
    
    if (nrow(lr_map) == 0) next
    
    res_spx <- cs@lrscore[["diffusion-contact_smooth"]]@res.list[["result.spx"]]
    lr_map <- filter(lr_map, intr_id %in% names(res_spx))
    if (nrow(lr_map) == 0) next
    
    sample_result <- lapply(seq_len(nrow(lr_map)), function(j) {
        receivers <- unique(as.character(res_spx[[lr_map$intr_id[j]]]))
        receivers <- intersect(receivers[!is.na(receivers)], names(cs@clusters))
        if (length(receivers) == 0) return(NULL)
        
        receiver_ct <- as.character(cs@clusters[receivers])
        receiver_ct[is.na(receiver_ct)] <- "Unknown"
        tibble(pair_key = lr_map$pair_key[j], receiver_ct = receiver_ct) %>%
            count(pair_key, receiver_ct, name = "n_receiver_cells")
    })
    
    result_list[[i]] <- bind_rows(sample_result)
    rm(cs)
}



receiver_counts <- bind_rows(result_list)

receiver_comp_by_lr <- receiver_counts %>%
    group_by(pair_key, receiver_ct) %>%
    summarise(n_receiver_cells = sum(n_receiver_cells), .groups = "drop") %>%
    complete(
        pair_key = selected_lr$pair_key,
        receiver_ct = celltype_levels,
        fill = list(n_receiver_cells = 0L)
    ) %>%
    group_by(pair_key) %>%
    mutate(frac_receiver_cells = n_receiver_cells / sum(n_receiver_cells)) %>%
    ungroup() %>%
    left_join(selected_lr, by = "pair_key") %>%
    select(pair_key, interaction, receiver_ct,
           n_receiver_cells, frac_receiver_cells)

write.csv(receiver_comp_by_lr, "Data/ST/CytoSignal/Final/selected_LR_receiver_celltype_fraction.csv", row.names = FALSE)


receiver_comp_by_lr <- receiver_comp_by_lr %>%
    mutate(
        receiver_ct = factor(receiver_ct, levels = celltype_levels),
        interaction = factor(interaction, levels = selected_lr$interaction)
    )

p_receiver <- ggplot(receiver_comp_by_lr,aes(x = frac_receiver_cells, y = interaction, fill = receiver_ct)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.25) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = celltype_colors, drop = FALSE) +
    labs(x = "Fraction (%)", y = "Ligand-receptor pair", fill = "Cell type") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave("Figure/ST/CytoSignal/Final/Figure6C_receiver_celltype_fraction.pdf", p_receiver, width = 4, height = 5)

