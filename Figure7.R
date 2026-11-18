library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggthemes)
library(ggsci)
library(ggrepel)
library(data.table)
library(circlize)
library(ComplexHeatmap)

Tcell_sub <- readRDS('Data/scRNA/Tcell/1_Tcell_obj.Rds')
CD8_Tcell <- subset(Tcell_sub, T_main == 'CD8_Tcell') 
CD8_Tcell <- seurat_process1(CD8_Tcell, harmony = T, vars = c("Sample"), npcs = 15, res=0.5) 

CD8_Tcell$Tsubcell2 <- case_when(
    
    CD8_Tcell$celltype %in% c("CD8_Tn") ~ "CD8_Tn",
    CD8_Tcell$celltype %in% c( "CD8_Tc17") ~ "CD8_Tc17",
    CD8_Tcell$celltype %in% c("CD8_GZMK_early_Tem") ~ "CD8_Tm",
    CD8_Tcell$celltype %in% c("CD8_IL7R_Tm") ~ "CD8_Tm",
    CD8_Tcell$celltype %in% c("CD8_Tex","CD8_GZMK_Tex", "CD8_OXPHOS_Tex", "CD8_TCF7+Tex") ~ "CD8_Tex",
    CD8_Tcell$celltype %in% c("CD8_Temra") ~ "CD8_Temra",
    CD8_Tcell$celltype %in% c("CD8_MKI67") ~ "CD8_MKI67",
    CD8_Tcell$celltype %in% c("CD8_ISG") ~ "CD8_ISG",
    CD8_Tcell$celltype %in% c("CD8_KIR+TXK+NK-like") ~ "CD8_NK_like",
    TRUE ~ "Other"
)

p2 <- DimPlot(CD8_Tcell,reduction='umap', group.by ="Tsubcell2", pt.size = 0.1, label = T) + scale_color_igv()
ggsave(filename = 'Figure/scRNA/Tcell/UMAP_Plot_CD8_Tcell.pdf', p2, width = 5, height = 4)


cd8_marker_panel_short <- c("CD8A", "CD8B",
                            "IL7R", "CCR7", "SELL", "TCF7", "LEF1",
                            "GZMK", "CCL5", "CXCR3", "NKG7",
                            "GNLY", "PRF1", "GZMB", "GZMH", "FGFBP2", "CX3CR1", "KLRG1",
                            "PDCD1", "TOX", "TIGIT", "LAG3", "HAVCR2", "CXCL13", "ENTPD1",
                            "SLAMF6", "CXCR5", "ID3",
                            "IL17A", "IL17F", "RORC", "CCR6",
                            "KIR2DL3", "KIR3DL1", "TXK", "KLRC1", "KLRD1",
                            "MKI67", "TOP2A", "STMN1",
                            "ISG15", "IFIT1", "MX1", "OAS1", "IFI6", "STAT1"
)


p2 <- DotPlot(CD8_Tcell, features = cd8_marker_panel_short,group.by = 'Tsubcell2') +  RotatedAxis()
p2
ggsave(filename = 'Figure/scRNA/Tcell/4_marker_Plot_CD8_Tcell.pdf', p2)


# ------------------- CD8 T cells: non-Tex / Tex UMAP -------------------

CD8_Tcell$Tex_status <- case_when(
    CD8_Tcell$celltype %in% c("CD8_Tex", "CD8_GZMK_Tex", "CD8_OXPHOS_Tex", "CD8_TCF7+Tex") ~ "Tex",
    CD8_Tcell$celltype %in% c("CD8_GZMK_early_Tem", "CD8_IL7R_Tm", "CD8_ISG", "CD8_MKI67", "CD8_Tc17", "CD8_Temra", "CD8_Tn", "CD8_KIR+TXK+NK-like") ~ "non-Tex",
    TRUE ~ "Other"
)

CD8_Tcell$Tex_status <- factor(CD8_Tcell$Tex_status, levels = c("non-Tex", "Tex", "Other"))

p_tex <- DimPlot(CD8_Tcell, reduction = "umap",group.by = "Tex_status",
                 cols = c("non-Tex" = "#66A3BF", "Tex" = "#FFA259", "Other" = "#D9D9D9"),pt.size = 0.1, label = FALSE) +
    labs(title = NULL, color = NULL) +
    theme_classic() +
    theme(legend.position = "right", axis.title = element_text(size = 11), axis.text = element_text(size = 10))
ggsave(filename = "Figure/scRNA/Tcell/UMAP_Plot_CD8_Tcell_nonTex_Tex.pdf",plot = p_tex, width = 5, height = 4)



  
cb_pattle2= c("E" = "#b2abd2", "NE" = "#fee0b6")

scRNA_obj = readRDS('Data/scRNA/rawData/icc_scRNA_clean.RDS')
sampleClone <- read.csv(file = 'Data/scRNA/Tcell/expansion_info.csv',row.names = 1)
sampleClone2 <- sampleClone %>% select(Sample,Expan) %>% distinct()

meta_df <- scRNA_obj@meta.data %>% tibble::rownames_to_column("cell") %>% left_join(sampleClone2, by = "Sample") %>% tibble::column_to_rownames("cell")
scRNA_obj@meta.data <- meta_df

table(scRNA_obj$Expan,scRNA_obj$orig.ident)

scRNA_pre = subset(scRNA_obj, Treatment == 'pre' & !is.na(Expan))
expr_df <- FetchData(object = scRNA_pre, vars = c("Sample", "Expan", c('LAMB3')), layer = "data")

sample_df <- expr_df %>%
    group_by(Sample, Expan) %>%
    summarise(mean_expr = mean(LAMB3), pos_frac  = mean(LAMB3 > 0), .groups = "drop")

p1 <- ggplot(sample_df, aes(x = Expan, y = mean_expr, fill = Expan)) +
    geom_boxplot(width = 0.45, outlier.shape = NA,  linewidth = 0.7) +
    geom_jitter(aes(fill = Expan), åcolor = "black", shape = 21, width = 0.12, size = 3.5, alpha = 0.95) +
    stat_compare_means(comparisons = list(c("E", "NE")), method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = cb_pattle2) +
    scale_color_manual(values = cb_pattle2) +
    labs(x = NULL, y = "LAMB3 expression in Cells") +
    theme_classic( ) +
    ylim(0,1)+
    theme(legend.position = "none",
          axis.title.x = element_blank())
p1

p2 <- ggplot(sample_df, aes(x = Expan, y = pos_frac, fill = Expan)) +
    geom_boxplot(width = 0.45, outlier.shape = NA, color = "#4D4D4D", linewidth = 0.7) +
    geom_jitter(aes(fill = Expan), color = "black", shape = 21, width = 0.12, size = 3.5, alpha = 0.95) +
    stat_compare_means(comparisons = list(c("E", "NE")), 
                       method = "wilcox.test", label = "p.format") +
    scale_fill_manual(values = cb_pattle2) +
    scale_color_manual(values = cb_pattle2) +
    labs(x = NULL, y = "% LAMB3+ cells") +
    theme_classic() +
    ylim(0,0.52)+
    theme(legend.position = "none",axis.title.x = element_blank())
p2

p = (p1 + p2) + plot_layout(ncol = 2)
p
ggsave('Figure/scRNA/Figure7/Figure7A_LAMB3_perSample_Expan.pdf',p,width = 4.5,height = 2.6)

 
T_obj <- readRDS('Data/scRNA/Tcell/1_Tcell_obj.Rds')

sampleClone <- read.csv(file = 'Data/scRNA/Tcell/expansion_info.csv')
sampleClone2 <- sampleClone %>% dplyr::select(Sample,Expan,num) %>% distinct()
colnames(sampleClone2) <- c('Sample','Expan','TCR_ratio')
head(sampleClone2)

meta_df <- T_obj@meta.data %>% tibble::rownames_to_column("cell") %>% left_join(sampleClone2, by = "Sample") %>% tibble::column_to_rownames("cell")
T_obj@meta.data <- meta_df
 
CD8_Tcell <- subset(T_obj, T_main == 'CD8_Tcell'&Treatment =='pre' )  
CD8_Tcell = subset(CD8_Tcell,!is.na(Expan))

CD8_Tcell$ER_group <- paste(CD8_Tcell$Expan, CD8_Tcell$response, sep = "-")
CD8_Tcell <- subset(CD8_Tcell, subset = ER_group %in% c("E-R", "E-NR", "NE-R", "NE-NR"))
CD8_Tcell$ER_group <- factor(CD8_Tcell$ER_group, levels = c("E-R", "E-NR", "NE-R", "NE-NR"))

Idents(CD8_Tcell) <- "ER_group"
markers_all <- FindAllMarkers(object = CD8_Tcell, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25, test.use = "wilcox")

 
top10_markers <- markers_all %>% 
    dplyr::filter(p_val_adj < 0.05,
                  !grepl("^MT-", gene, ignore.case = TRUE),
                  !grepl("^RPS", gene, ignore.case = TRUE),
                  !grepl("^RPL", gene, ignore.case = TRUE),
                  !grepl("^HSP", gene, ignore.case = TRUE),
                  !grepl("^LINC", gene, ignore.case = TRUE), 
                  !grepl("^LINC", gene, ignore.case = TRUE),
                  !gene %in% c("MALAT1", "NEAT1","AC020916.1")) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC),p_val_adj) %>%
    slice_head(n = 10) %>%
    ungroup()

top10_genes <- unique(top10_markers$gene)
avg_exp <- AverageExpression(CD8_Tcell, group.by = "ER_group", features = top10_genes, assays = "RNA",slot = "data")$RNA

mat_z <- t(scale(t(avg_exp)))
mat_z[is.na(mat_z)] <- 0

gene_order_df <- top10_markers %>% dplyr::select(gene, cluster, avg_log2FC) %>% distinct() %>%
    arrange(factor(cluster, levels = c("E-R", "E-NR", "NE-R", "NE-NR")),desc(avg_log2FC))
mat_z <- mat_z[gene_order_df$gene, c("E-R", "E-NR", "NE-R", "NE-NR")]

col_anno <- data.frame(Expan = c("E", "E", "NE", "NE"),Response = c("R", "NR", "R", "NR"))
rownames(col_anno) <- c("E-R", "E-NR", "NE-R", "NE-NR")
anno_colors <- list( Expan = c("E" = "#bf812d","NE" = "#f6e8c3"),
                     Response = c("R" = "#529386","NR" = "#BDE3C3"))

top_ha <- HeatmapAnnotation(df = col_anno, col = anno_colors,
                            annotation_name_side = "left",
                            annotation_name_gp = gpar(fontsize = 8),
                            simple_anno_size = unit(0.3, "cm"))
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#01665e", "white", "#762a83"))

ht <- Heatmap(mat_z, 
              name = "Z-score",
              col = col_fun,
              top_annotation = top_ha,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = TRUE,
              show_column_names = F,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 10),
              rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(at = c(-1.5, 0, 1.5), title = "Expression", legend_height = unit(4, "cm")))

pdf('Figure/scRNA/Tcell/CD8_Tcell/1_CD8_Tcell_top10_Gene_expr_expansion_response.pdf',width = 4,height = 6)
draw(ht)
dev.off()



# -------------------------------------------
cytotoxic_genes <- c("NKG7", "PRF1", "GNLY", "GZMB", "CTSW","CCL5", "GZMH", "CST7", "IFNG", "CX3CR1")
proliferation_genes <- c( "MKI67", "TOP2A", "TYMS", "PCNA", "HMGB2",
                          "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
                          "FOXM1", "MYBL2", "PLK1", "AURKA", "BUB1",
                          "UBE2C", "BIRC5", "TK1", "CENPF", "ZWINT")
cytotoxic_genes_use <- intersect(cytotoxic_genes, rownames(CD8_Tcell))
proliferation_genes_use <- intersect(proliferation_genes, rownames(CD8_Tcell))

CD8_Tcell <- AddModuleScore(object = CD8_Tcell, features = list(cytotoxic_genes_use), name = "Cytotoxic_Score")
CD8_Tcell <- AddModuleScore(object = CD8_Tcell, features = list(proliferation_genes_use), name = "Proliferation_Score")

head(CD8_Tcell@meta.data)
CD8_Tcell$Cytotoxic <- CD8_Tcell$Cytotoxic_Score1
CD8_Tcell$Proliferation <- CD8_Tcell$Proliferation_Score1

saveRDS(CD8_Tcell,file = 'Data/scRNA/Tcell/2_CD8_Tcell_seu_Obj.Rds')

top10_genes = top10_markers%>% dplyr::filter(cluster == 'NE-NR')%>%pull(gene)

expr_df <- FetchData(CD8_Tcell, vars = c(top10_genes, "Cytotoxic", "Proliferation", "Sample"))
expr_df$TCR_ratio <- CD8_Tcell$TCR_ratio

sample_df <- expr_df %>%
    group_by(Sample) %>%
    summarise(
        Cytotoxic = mean(Cytotoxic, na.rm = TRUE),
        Proliferation = mean(Proliferation, na.rm = TRUE),
        TCR_ratio = dplyr::first(TCR_ratio),
        n_cells = n(),
        across(all_of(top10_genes), ~mean(.x, na.rm = TRUE))
    ) %>%
    ungroup()

sample_df_use <- sample_df %>% filter(n_cells >= 10)
sample_df_use

make_cor_plot <- function(df, gene, score_var) {
    tmp_cor <- cor.test(x = df[[score_var]], y = df[[gene]], method = "pearson")
    
    ggplot(df, aes_string(x = score_var, y = gene)) +
        geom_point(size = 2, color = "#01665e") +
        geom_smooth(method = "lm", se = TRUE, color = "#bf812d") +
        annotate("text", x = Inf, y = Inf, 
                 label = paste0("R = ", round(tmp_cor$estimate, 3), "\nP = ", signif(tmp_cor$p.value, 3)),
                 hjust = 1.1, vjust = 1.5, size = 4) +
        theme_few() +
        labs(title = NULL, x = score_var, y = paste0(gene, " expression \nin CD8 T cell")) +
        theme(
            plot.title = element_blank(),
            axis.title = element_text(size = 11)
        )
}

gene_plot_list <- list()
for (gene in top10_genes) {
    p1 <- make_cor_plot(sample_df_use, gene, "Cytotoxic")
    p2 <- make_cor_plot(sample_df_use, gene, "Proliferation")
    p3 <- make_cor_plot(sample_df_use, gene, "TCR_ratio")
    p_gene <- p1 | p2 | p3
    gene_plot_list[[gene]] <- p_gene
}
p = gene_plot_list[['IL7R']] /gene_plot_list[['NCR3']]
p <- wrap_plots(gene_plot_list, ncol = 2)
ggsave(filename = 'Figure/scRNA/Tcell/CD8_Tcell/2_corrlation_gene_with_cytotoxix.pdf',p,width = 16,height = 12)


CD8_Tcell <- readRDS(file = 'Data/scRNA/Tcell/2_CD8_Tcell_seu_Obj.Rds')
plot_gene_expr_pos <- function(object, genes, sample_var = "Sample",group_var, group_levels, group_colors, compare_groups, y_label_suffix = "CD8+T Cells", layer = "data", test_method = "t.test") {
    
    meta_vars <- c(sample_var, group_var, genes)
    expr_df <- FetchData(object = object, vars = meta_vars, layer = layer)
    
    expr_long <- expr_df %>%
        pivot_longer(
            cols = all_of(genes),
            names_to = "Gene",
            values_to = "Expression"
        ) %>%
        mutate(
            Gene = factor(Gene, levels = genes),
            !!group_var := factor(.data[[group_var]], levels = group_levels)
        )
    
    sample_df <- expr_long %>%
        group_by(.data[[sample_var]], .data[[group_var]], Gene) %>%
        summarise(
            mean_expr = mean(Expression, na.rm = TRUE),
            pos_frac  = mean(Expression > 0, na.rm = TRUE),
            .groups = "drop"
        )
    
    plot_list <- lapply(genes, function(g) {
        
        df_g <- sample_df %>% filter(Gene == g)
        
        p_mean <- ggplot(df_g, aes(x = .data[[group_var]], y = mean_expr, fill = .data[[group_var]])) +
            geom_boxplot(width = 0.45, outlier.shape = NA, linewidth = 0.7) +
            geom_jitter(color = "black", shape = 21, width = 0.12, size = 3, alpha = 0.95) +
            stat_compare_means(
                comparisons = list(compare_groups),
                method = test_method,
                label = "p.format",
                size = 4,
                label.y.npc = 0.98
            ) +
            scale_fill_manual(values = group_colors) +
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
            coord_cartesian(clip = "off") +
            labs(x = NULL, y = paste0(g, " expression in\n", y_label_suffix)) +
            theme_classic() +
            theme(
                legend.position = "none",
                axis.title.x = element_blank(),
                plot.margin = margin(8, 12, 8, 8)
            )
        
        p_pos <- ggplot(df_g, aes(x = .data[[group_var]], y = pos_frac, fill = .data[[group_var]])) +
            geom_boxplot(width = 0.45, outlier.shape = NA, color = "#4D4D4D", linewidth = 0.7) +
            geom_jitter(color = "black", shape = 21, width = 0.12, size = 3, alpha = 0.95) +
            stat_compare_means(
                comparisons = list(compare_groups),
                method = test_method,
                label = "p.format",
                size = 4,
                label.y.npc = 0.98
            ) +
            scale_fill_manual(values = group_colors) +
            scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
            coord_cartesian(clip = "off") +
            labs(x = NULL, y = paste0("% ", g, "+ ", y_label_suffix)) +
            theme_classic() +
            theme(
                legend.position = "none",
                axis.title.x = element_blank(),
                plot.margin = margin(8, 12, 8, 8)
            )
        
        p_mean + p_pos
    })
    
    wrap_plots(plot_list, ncol = 1)
}

p1 <- plot_gene_expr_pos(object = CD8_Tcell, genes = 'IL7R', group_var = "Expan", group_levels = c("E", "NE"), group_colors = c("E" = "#B2ABD2", "NE" = "#FEE0B6"),
                         compare_groups = c("E", "NE"),
                         y_label_suffix = "CD8+T Cells",
                         test_method = "t.test")
p2 <- plot_gene_expr_pos(object = CD8_Tcell,
                         genes = 'NCR3',
                         group_var = "Expan",
                         group_levels = c("E", "NE"),
                         group_colors = c("E" = "#B2ABD2", "NE" = "#FEE0B6"),
                         compare_groups = c("E", "NE"),
                         y_label_suffix = "CD8+T Cells",
                         test_method = "t.test")
ggsave(filename = 'Figure/scRNA/Tcell/CD8_Tcell/3_IL7R_NCR3_exprssion_per_sample.pdf',plot = p1+p2)

obj = readRDS('Data/scRNA/rawData/icc_scRNA_clean.RDS')
meta <- obj@meta.data
sample_col <- "Sample"
celltype_col <- "cellsubtype"
genes_use <- c("LAMB3", "NCR3")
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
expr_mat <- counts[genes_use, , drop = FALSE]
df <- data.frame(cell_id = colnames(obj), Sample = meta$Sample, 
                 cellsubtype = meta$cellsubtype,
                 LAMB3_count = as.numeric(expr_mat["LAMB3", ]),
                 NCR3_count  = as.numeric(expr_mat["NCR3", ]))

df <- df %>% mutate(LAMB3_pos = LAMB3_count > 0,NCR3_pos  = NCR3_count > 0, is_CD8T = grepl("^CD8_", cellsubtype), is_NCR3_pos_CD8T = is_CD8T & NCR3_pos)
sample_prop_df <- df %>% group_by(Sample) %>%
    summarise(n_all_cells = n(),
              # LAMB3:All cells
              n_LAMB3_pos = sum(LAMB3_pos),
              LAMB3_pos_prop_all = n_LAMB3_pos / n_all_cells,
              # NCR3 :all CD8 T cells
              n_CD8T = sum(is_CD8T),
              n_NCR3_pos_CD8T = sum(is_NCR3_pos_CD8T),
              NCR3_pos_prop_in_CD8T = ifelse(n_CD8T > 0, n_NCR3_pos_CD8T / n_CD8T, NA_real_),
              .groups = "drop") %>%
    mutate(LAMB3_pos_pct_all = LAMB3_pos_prop_all * 100, NCR3_pos_pct_in_CD8T = NCR3_pos_prop_in_CD8T * 100)
 
cor_df <- sample_prop_df %>% filter(!is.na(LAMB3_pos_prop_all), !is.na(NCR3_pos_prop_in_CD8T), n_CD8T >= 20)
 
p <- ggplot(cor_df, aes( x = LAMB3_pos_pct_all, y = NCR3_pos_pct_in_CD8T)) +
    geom_smooth(
        method = "lm",
        se = TRUE,
        color = "#1C4D8D",fill = '#92c5de',
        linewidth = 0.8) +
    geom_point(aes(size = n_CD8T),alpha = 0.85,shape = 21,fill = '#2166ac') +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 4) +
    labs(
        x = "LAMB3+ cells / all cells (%)",
        y = "NCR3+ cells / CD8 T cells (%)",
        size = "CD8 T cells") +
    theme_classic()  

p
ggsave('Figure/scRNA/Tcell/NCR3-LAMB3_corrlation_plot.pdf',width = 4,height = 3)


 

CD8_Tcell <- readRDS(file = 'Data/scRNA/Tcell/2_CD8_Tcell_seu_Obj.Rds')
CD8_Tcell$NCR3_expr <- FetchData(CD8_Tcell, vars = "NCR3")[, 1]
ncr3_cutoff <- 0
CD8_Tcell$NCR3_group <- ifelse(CD8_Tcell$NCR3_expr > ncr3_cutoff, "pos", "neg")
CD8_Tcell$NCR3_group <- factor(CD8_Tcell$NCR3_group, levels = c("pos", "neg"))

sig_list <- list(Effector_Memory = unique(c("IL7R","CCL4L1","CCL4L2","HLA-DRA","HLA-DRB5","NKG7","HLA-DRB6","HLA-DRB1","GZMH","CCL4",
                                            "GPR183","TNFAIP3","HLA-DPA1","HSPA1B","CD55","CD74","HLA-DPB1","HSPA1A","TOB1","CCR7",
                                            "NFKBIA","FOS","KLRB1","PIM2","PFKFB3","GZMB","SATB1","SC5DL","NFKBIZ","SELL",
                                            "PDE4D","IFNG","CCR5","GNLY","SNHG16","TCF7","PYHIN1","LMNA","PTGER2","RORA",
                                            "TMEM123","TNFRSF25","TOX","AREG","SERINC5","CRTAM","IDH2","HLA-DQB1","NELL2","CARD16",
                                            "KIAA0101","DPP4","CCR6","ICAM2","CDKN1A","RRM2","LOC541471","ITM2C","OASL","STMN1")),
                 
                 Chemokine_IFNg = unique(c("RGS1","GNLY","CCL4","GZMK","CCL4L1","CCL4L2","IL7R","CD74","ZFP36","PTGER2",
                                           "CD55","HSPA1A","HSPA1B","KLRD1","DUSP4","HLA-DRB1","ZNF683","PLAC8","HLA-DRB5","CRTAM",
                                           "S100A4","GZMB","KLRB1","LGALS1","HLA-DRA","HLA-DPA1","HOPX","TMEM2","HLA-DQB1","FLT3LG",
                                           "S1PR1","PIM1","FOSB","DSTN","ADAM19","SERINC5","CCL3","FKBP11","TMEM123","HLA-DQA1",
                                           "LAIR1","TOX","TNFRSF9","RRM2","HLA-DMA","TNFSF14","IFNG","BTG3","SERTAD1")),
                 
                 Cytotoxicity = c("NKG7", "PRF1", "GNLY", "GZMB", "CTSW","CCL5", "GZMH", "CST7", "IFNG", "CX3CR1"),
                 
                 Interferon = unique(c("CLIC1","CD74","PSME2","APOBEC3G","HLA-DRB1","HLA-DRA","MT2A","IL7R","HLA-DRB5","IFI6",
                                       "CRIP1","APOBEC3C","COTL1","COX8A","TSC22D3","IDH2","STMN1","CCR5","UBE2L6","PLAC8",
                                       "APOBEC3H","STAT1","TXNIP","PTTG1","ANXA5","CALM3","OASL","ISG15","GBP1","ICOS",
                                       "ITM2C","CD38","HMGN2","HMGB2","COX5A","CARD16","PCNA","NDUFA12","DNPH1","MCM5",
                                       "PYCARD","CCDC167","EZH2","FABP5","SIVA1","NDUFB6","RASA3","KIAA0101","HSPB11","MX1",
                                       "MCM7","SKA2","PLSCR1","LOC541471","TNFSF10","NDUFB3","AIM1","OAS3","COX7B","STRA13",
                                       "KLRD1","SNRPG","CD55","RRM2","ATHL1","SAMD9L","SMC4","MT1X","PSMD14","LGALS1",
                                       "IFI44","SF3B14","CTLA4","PPP5C","HSPA1B","NT5C3","GLUL","IFI35","FOSB","RRM1",
                                       "EPSTI1","TUBB","KIF11","CENPN","DUSP1","TRAFD1","CLIC3","LSM7","IRF7","RSAD2",
                                       "CCR6","LOC284454","KLRB1","TUBA1B")))


sig_list_use <- lapply(sig_list, function(x) intersect(x, rownames(CD8_Tcell)))

for (nm in names(sig_list_use)) {
    CD8_Tcell <- AddModuleScore(object = CD8_Tcell, features = list(sig_list_use[[nm]]), name = paste0(nm, "_Score"),assay = "RNA")
    CD8_Tcell[[nm]] <- CD8_Tcell[[paste0(nm, "_Score1")]]
}


plot_df <- CD8_Tcell@meta.data %>%
    dplyr::select(NCR3_group,  Effector_Memory,Chemokine_IFNg,Cytotoxicity, Interferon) %>%
    pivot_longer(cols = -NCR3_group, names_to = "Signature", values_to = "Score")

plot_df$NCR3_group <- factor(plot_df$NCR3_group, levels = c("pos", "neg"))

plot_df2 <- plot_df %>% group_by(Signature, NCR3_group) %>% mutate(mean_score = mean(Score, na.rm = TRUE)) %>% ungroup()

plot_one <- function(sig_name) {
    df_sub <- plot_df2 %>% filter(Signature == sig_name)
    
    ggplot(df_sub, aes(x = NCR3_group, y = Score, fill = mean_score)) +
        geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = 0.3) +
        geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
        stat_compare_means(comparisons = list(c("pos", "neg")), method = "wilcox.test", label = "p.format", size = 3.5) +
        scale_fill_gradientn(colours = c("#c7eae5", "#80cdc1", "#35978f")) +
        labs(title = sig_name, x = NULL, y = "Module score") +
        theme_bw() +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.1)))+
        theme(plot.title = element_text(hjust = 0.5,size = 11),
              strip.background = element_rect(fill = "white", color = "black"),
              axis.text.x = element_text(angle = 0, hjust = 0.5),
              panel.grid = element_blank(),
              legend.position = 'none'
        )
}
plot_list <- lapply(names(sig_list_use), plot_one)
p <- wrap_plots(plotlist = plot_list, ncol = 4)
ggsave('Figure/scRNA/Tcell/CD8_Tcell/5_NCR3_vs_multi_score_VlnPlot.pdf',p,width = 7,height = 2.5)

#Figure7I
CD8_Tcell <-  readRDS(file = 'Data/scRNA_Data/2_Tcell/2_CD8_Tcell_seu_Obj.Rds')
meta <- CD8_Tcell@meta.data
meta$NCR3_expr <- FetchData(CD8_Tcell, vars = "NCR3")[, 1]
plot_df1 <- meta %>%
    mutate(
        Patient = sub("_(Pre|Post|pre|post)$", "", Sample),
        Treatment2 = case_when(
            grepl("post", tolower(Treatment)) ~ "Post",
            grepl("pre",  tolower(Treatment)) ~ "Pre",
            grepl("_Post$", Sample) ~ "Post",
            grepl("_Pre$",  Sample) ~ "Pre",
            TRUE ~ NA_character_
        ),
        IL7R_Tm_flag = Tsubcell2 == "CD8_Tn"
    ) %>%
    filter(!is.na(Treatment2)) %>%
    group_by(Patient, Treatment2) %>%
    summarise(
        n_T = n(),
        n_IL7R_Tm = sum(IL7R_Tm_flag, na.rm = TRUE),
        prop_IL7R_Tm = n_IL7R_Tm / n_T,
        .groups = "drop"
    )

plot_df1_pair <- plot_df1 %>%
    group_by(Patient) %>%
    filter(all(c("Pre", "Post") %in% Treatment2)) %>%
    ungroup()

plot_df1_pair$Treatment2 <- factor(plot_df1_pair$Treatment2, levels = c("Pre", "Post"))

test_df1 <- plot_df1_pair %>%
    select(Patient, Treatment2, prop_IL7R_Tm) %>%
    pivot_wider(names_from = Treatment2, values_from = prop_IL7R_Tm)

wilcox_res1 <- t.test(test_df1$Pre, test_df1$Post, paired = TRUE, alternative = "less")
p_label1 <- paste0("p = ", signif(wilcox_res1$p.value, 3))

p1 <- ggplot(plot_df1_pair, aes(x = Treatment2, y = prop_IL7R_Tm, group = Patient)) +
    geom_line(color = "#45B5D8", linewidth = 1) +
    geom_point(size = 3, shape = 21, fill = "#45B5D8", color = "black", stroke = 0.3) +
    annotate("text", x = 1.5, y = max(plot_df1_pair$prop_IL7R_Tm, na.rm = TRUE) * 1.08, label = p_label1, size = 4) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.15))) +
    labs(x = NULL, y = "% IL7R+ Tm in CD8+T") +
    theme_classic( )  
 
plot_df2 <- meta %>%
    mutate(
        Patient = sub("_(Pre|Post|pre|post)$", "", Sample),
        Treatment2 = case_when(
            grepl("post", tolower(Treatment)) ~ "Post",
            grepl("pre",  tolower(Treatment)) ~ "Pre",
            grepl("_Post$", Sample) ~ "Post",
            grepl("_Pre$",  Sample) ~ "Pre",
            TRUE ~ NA_character_
        ),
        NCR3_pos_flag = NCR3_expr > 0
    ) %>%
    filter(!is.na(Treatment2), celltype == "CD8_IL7R_Tm") %>%
    group_by(Patient, Treatment2) %>%
    summarise(
        n_IL7R_Tm = n(),
        n_NCR3_pos = sum(NCR3_pos_flag, na.rm = TRUE),
        prop_NCR3_in_IL7R_Tm = n_NCR3_pos / n_IL7R_Tm,
        .groups = "drop"
    )

plot_df2_pair <- plot_df2 %>%
    group_by(Patient) %>%
    filter(all(c("Pre", "Post") %in% Treatment2)) %>%
    ungroup()
plot_df2_pair$Treatment2 <- factor(plot_df2_pair$Treatment2, levels = c("Pre", "Post"))
test_df2 <- plot_df2_pair %>%
    select(Patient, Treatment2, prop_NCR3_in_IL7R_Tm) %>%
    pivot_wider(names_from = Treatment2, values_from = prop_NCR3_in_IL7R_Tm)
wilcox_res2 <- t.test(test_df2$Pre, test_df2$Post, paired = TRUE, alternative = "greater")
p_label2 <- paste0("p = ", signif(wilcox_res2$p.value, 3))
p2 <- ggplot(plot_df2_pair, aes(x = Treatment2, y = prop_NCR3_in_IL7R_Tm, group = Patient)) +
    geom_line(color = "#7BAF3A", linewidth = 1) +
    geom_point(size = 3, shape = 21, fill = "#7BAF3A", color = "black", stroke = 0.3) +
    annotate("text",x = 1.5, y = max(plot_df2_pair$prop_NCR3_in_IL7R_Tm, na.rm = TRUE) * 1.08, label = p_label2, size = 4) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.15))) +
    labs(x = NULL, y = "% NCR3+ IL7R+ Tm") +
    theme_classic() 
ggsave('Figure/scRNA_Figure/4_Tcell/5_CD8_Tcell/6_IL7R_expr_pre_post.pdf',plot = p1+p2)



 


library(CellChat)

scRNA_obj <- readRDS(file = 'Data/scRNA/rawData/icc_scRNA_clean.RDS')
sampleClone <- read.csv(file = 'Data/scRNA_Data/2_Tcell/expansion_info.csv')
sampleClone2 <- sampleClone %>% dplyr::select(Sample,Expan) %>%distinct()
meta_df <- scRNA_obj@meta.data %>% tibble::rownames_to_column("cell") %>% left_join(sampleClone2, by = "Sample") %>% tibble::column_to_rownames("cell")
scRNA_obj@meta.data <- meta_df
scRNA_obj$cellchat_label <- case_when(scRNA_obj$celltype %in% c("C0_ADH1B", "C1_TIMP1", "C2_APOA1", "C3_APOB") ~ "Hepatocyte",
                                      scRNA_obj$celltype %in% c("Neutrophil") ~ "Neutrophil",
                                      scRNA_obj$celltype %in% c("gamma_delt_T") ~ "γδ T",
                                      scRNA_obj$celltype %in% c("CD4_IL21_Tfh", "CD4_IL7R_Tn", "CD4_Temra","CD4_TIMP1_Tm", "CD4_TNFRSF9_Treg") ~ "CD4",
                                      scRNA_obj$celltype %in% c("B_cell", "Naive_B", "Centroblast", "Plasma") ~ "Bcell",
                                      scRNA_obj$celltype %in% c("cDC_CD1C", "cDC_CLEC9A") ~ "cDC",
                                      scRNA_obj$celltype %in% c("Fibor_apCAF", "Fibro_eCAF", "Fibro_iCAF",
                                                                "Fibro_mCAF", "Fibro_MKI", "Fibro_vCAF") ~ "Fibroblast",
                                      scRNA_obj$celltype %in% c("Endo_C0_COL4A1", "Endo_C1_ACKR1", "Endo_C2_FCN3",
                                                                "Endo_C3_FBLN5", "EndoC4_CCL21") ~ "Endothelial",
                                      scRNA_obj$celltype %in% c("Tumor") ~ "Tumor",
                                      scRNA_obj$celltype %in% c("Mono_CD14") ~ "Monocyte",
                                      scRNA_obj$celltype %in% c("Mast") ~ "Mast",
                                      scRNA_obj$celltype %in% c("Macro_C1QC") ~ "Macrophage",
                                      scRNA_obj$celltype %in% c("CD8_GZMK_Tex", "CD8_OXPHOS_Tex","CD8_TCF7+Tex", "CD8_Tex") ~ "CD8_Tex",
                                      scRNA_obj$celltype %in% c("CD8_GZMK_early_Tem", "CD8_IL7R_Tm","CD8_ISG", "CD8_MKI67", "CD8_Tc17",
                                                                "CD8_Temra", "CD8_Tn", "CD8_KIR+TXK+NK-like") ~ "CD8_non_Tex",
                                      scRNA_obj$celltype %in% c("NK") ~ "NK",TRUE ~ NA_character_)

scRNA_obj$cellchat_label <- factor(scRNA_obj$cellchat_label,
                                   levels = c("Hepatocyte", "Neutrophil", "γδ T", "CD4", "Bcell", "cDC", "Fibroblast", "Endothelial", "Tumor", "Monocyte",
                                              "Mast", "Macrophage", "CD8_non_Tex", "CD8_Tex", "NK"))

scRNA_obj <- subset(scRNA_obj,!is.na(Expan)&Treatment == 'pre')
 
obj_s2 <- subset(scRNA_obj, Expan == "E")
obj_s3 <- subset(scRNA_obj, Expan == "NE")

data_s2 <- GetAssayData(obj_s2, assay = "RNA", layer = "data")
meta_s2 <- obj_s2@meta.data

data_s3 <- GetAssayData(obj_s3, assay = "RNA", layer = "data")
meta_s3 <- obj_s3@meta.data

cellchat_s2 <- createCellChat(object = data_s2, meta = meta_s2, group.by = 'cellchat_label')
cellchat_s3 <- createCellChat(object = data_s3, meta = meta_s3, group.by = 'cellchat_label')

CellChatDB <- CellChatDB.human
cellchat_s2@DB <- CellChatDB
cellchat_s3@DB <- CellChatDB
cellchat_s2 <- subsetData(cellchat_s2)
cellchat_s3 <- subsetData(cellchat_s3)

cellchat_s2 <- identifyOverExpressedGenes(cellchat_s2)
cellchat_s2 <- identifyOverExpressedInteractions(cellchat_s2)

cellchat_s3 <- identifyOverExpressedGenes(cellchat_s3)
cellchat_s3 <- identifyOverExpressedInteractions(cellchat_s3)
cellchat_s2 <- computeCommunProb(cellchat_s2, raw.use = TRUE)
cellchat_s3 <- computeCommunProb(cellchat_s3, raw.use = TRUE)

cellchat_s2 <- filterCommunication(cellchat_s2, min.cells = 10)
cellchat_s3 <- filterCommunication(cellchat_s3, min.cells = 10)
cellchat_s2 <- computeCommunProbPathway(cellchat_s2)
cellchat_s3 <- computeCommunProbPathway(cellchat_s3)
cellchat_s2 <- aggregateNet(cellchat_s2)
cellchat_s3 <- aggregateNet(cellchat_s3)
cellchat_s2 <- netAnalysis_computeCentrality(cellchat_s2, slot.name = "netP")
cellchat_s3 <- netAnalysis_computeCentrality(cellchat_s3, slot.name = "netP")
object.list <- list(E = cellchat_s2, NE = cellchat_s3)
cellchat_merge <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat_merge, file = "Data/scRNA_Data/2_Tcell_analysis//cellchat_E_vs_NE_merged.rds")
 
cellchat_res <- netVisual_bubble(cellchat_merge, targets.use = c('CD8_Tex', 'CD8_non_Tex','NK'), 
                                 sources.use  = c('Hepatocyte','Neutrophil','γδ T','CD4','Bcell', 'cDC', 'Fibroblast', 'Endothelial', 'Tumor','Monocyte','Mast','Macrophage'),
                                 return.data = T,
                                 comparison = c(1, 2), angle.x = 45,thresh = 0.05)

df <- cellchat_res[["communication"]]
df <- subset(df, interaction_name_2 == 'BAG6 - NCR3-PS')

pairLR.use <- data.frame(interaction_name = unique(df$interaction_name), stringsAsFactors = FALSE)
 

source_order <- c("Hepatocyte", "Neutrophil", "γ δ T", "CD4", "Bcell", "cDC", "Fibroblast", "Endothelial", "Tumor", "Monocyte", "Mast", "Macrophage")
target_order  <- c('CD8_non_Tex', 'CD8_Tex','NK')
dataset_order <- c("E", "NE")

plot_df <- df %>% mutate(source  = as.character(source),
                         target  = as.character(target),
                         dataset = as.character(dataset))
plot_df <- expand_grid(source  = source_order,target  = target_order, dataset = dataset_order) %>%
    left_join(plot_df %>% dplyr::select(source, target, dataset, prob, pval), by = c("source", "target", "dataset")) %>%
    mutate(
        prob = ifelse(is.na(prob), 0, prob),
        pval = ifelse(is.na(pval), 1, pval),
        p_value = ifelse(pval > 1,  "p < 0.05","p > 0.05"),
        source  = factor(source, levels = rev(source_order)),
        target  = factor(target, levels = target_order),
        dataset = factor(dataset, levels = dataset_order),
        p_value = factor(p_value, levels = c("p > 0.05", "p < 0.05"))
    )
p <- ggplot(plot_df, aes(x = dataset, y = source)) +
    geom_point(aes(size = p_value, fill = prob), shape = 21, color = "black", stroke = 0.35) +
    facet_grid(. ~ target, scales = "free_x", space = "free_x") +
    scale_fill_gradientn(colors = c('#dfc27d','#f5f5f5','#80cdc1'))+
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        axis.ticks.length = unit(0.15, "cm"),
        panel.spacing.x = unit(0.15, "cm"))


p
ggsave('Figure/scRNA_Figure/5_Tcell/CD8_Tcell/4_BAG6_NCR3_LR.pdf',p,width = 3.8,height = 3.4)