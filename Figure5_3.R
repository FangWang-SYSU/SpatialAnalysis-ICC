rm(list = ls())
gc()

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(clusterProfiler)
library(stringr)
scRNA_obj <- readRDS(file = 'Data/scRNA/rawData/icc_scRNA_clean.RDS')
Fibro_sub <- subset(scRNA_obj, ct == 'Fibroblast') 
 
marker.gene <- FindMarkers(Fibro_sub, ident.1 = "NR", ident.2 = "R", 
                           group.by = "response",logfc.threshold = 0, min.pct = 0.05)


deg_gene <- marker.gene
deg_gene$gene <- rownames(deg_gene) 
deg_gene <- deg_gene %>% mutate(change = case_when(p_val_adj < 0.05 & avg_log2FC >= 1  ~  "NR", 
                                                   p_val_adj < 0.05 & avg_log2FC <= -1 ~  "R", TRUE ~ "NS"),
                                neglog10_padj = -log10(p_val_adj + 1e-300))
 
exclude_pattern <- "^(MT-|RPL|RPS|MALAT1$|HBA1$|HBA2$|HBB$|LINC)"
deg_plot <- deg_gene %>% filter(!grepl(exclude_pattern, gene))  

top_up_p <- deg_plot %>% filter(change == "NR") %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice_head(n = 10) %>% pull(gene)
top_up_fc <- deg_plot %>% filter(change == "NR") %>% arrange(desc(avg_log2FC), p_val_adj) %>% slice_head(n = 10)  %>% pull(gene)
top_down_p <- deg_plot %>% filter(change == "R") %>% arrange(p_val_adj, avg_log2FC) %>% slice_head(n = 10)  %>% pull(gene)
top_down_fc <- deg_plot %>% filter(change == "R") %>% arrange(avg_log2FC, p_val_adj) %>%slice_head(n = 10) %>% pull(gene)

deg_plot %>% filter(change == "NR") %>%   pull(gene)
deg_plot %>% filter(change == "NR") %>%  pull(gene)

label_genes <- unique(c(top_up_p,top_down_p))

label_df <- deg_plot %>% filter(gene %in% label_genes)
  
format_df <- function(df, direction = c("NR", "R")) {
    df <- as.data.frame(df)
    direction <- match.arg(direction)
    
    if (nrow(df) == 0) return(NULL)
    
    df$plot_q <- if ("qvalue" %in% colnames(df)) df$qvalue else df$p.adjust
    df$plot_q[is.na(df$plot_q)] <- df$p.adjust[is.na(df$plot_q)]
    
    df <- df %>%
        filter(!is.na(plot_q), plot_q > 0) %>%
        arrange(plot_q) %>% 
        mutate(Description = str_wrap(Description, width = 100),
               score = -log10(plot_q),
               direction = direction)
    
    if (direction == "R") {
        df$score <- -df$score
    }
    df
}
 
deg_plot <- deg_plot %>% mutate(change = case_when(p_val_adj < 0.05 & avg_log2FC >= 1  ~  "NR", 
                                                   p_val_adj < 0.05 & avg_log2FC <= -1 ~  "R", TRUE ~ "NS"),
                                neglog10_padj = -log10(p_val_adj + 1e-300))

bg_genes <- rownames(marker.gene )
up_genes <- deg_plot %>% filter(change == "NR") %>% pull(gene) %>% unique()
up_genes

down_genes <- deg_plot %>% filter(change == "R") %>% pull(gene) %>% unique()
down_genes

gene_bg_df <- bitr(bg_genes,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_up_df <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_down_df <- bitr(down_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

up_entrez <- unique(gene_up_df$ENTREZID)
down_entrez <- unique(gene_down_df$ENTREZID)
bg_entrez <- unique(gene_bg_df$ENTREZID)

ego_up <- enrichGO(gene = up_entrez,universe = bg_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,readable = TRUE)

ego_down <- enrichGO(gene = down_entrez,universe = bg_entrez, OrgDb = org.Hs.eg.db,keyType = "ENTREZID",
                     ont =  "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)

ratio_to_numeric <- function(x) {
    vapply(
        strsplit(as.character(x), "/", fixed = TRUE),
        function(z) as.numeric(z[1]) / as.numeric(z[2]),
        numeric(1)
    )
}

go_all <- bind_rows(
    as.data.frame(ego_down) %>% mutate(direction = "R"),
    as.data.frame(ego_up) %>% mutate(direction = "NR")
)

go_sig <- go_all %>%
    filter(!is.na(p.adjust), p.adjust <= 0.01) %>%
    mutate(
        FoldEnrichment = ratio_to_numeric(GeneRatio) / ratio_to_numeric(BgRatio),
        log10padj = -log10(p.adjust + 1e-300)) %>%
    arrange(direction, p.adjust, desc(Count)) %>%
    relocate(direction, ID, Description)
 
go_sig %>% filter(direction == "R") %>% pull(Description)
go_sig %>% filter(direction == "NR") %>% pull(Description)
write.csv(go_sig, file = "Data/scRNA/Fibroblast/GO_BP_NR_significant.csv") 
 

term_map <- tibble::tribble(
    ~direction, ~Description,
    
    "R", "myeloid leukocyte migration",
    "R", "leukocyte chemotaxis",
    "R", "monocyte chemotaxis",
    "R", "mononuclear cell migration",
    "R", "regulation of leukocyte migration",
    "R", "extracellular matrix organization",
    "NR", "phospholipase C-activating G protein-coupled receptor signaling pathway",
    "NR", "calcium ion transmembrane transport",
    "NR", "vasoconstriction",
    "NR", "endothelial cell migration",
    "NR", "regulation of blood circulation",
    "NR", "cytokine-mediated signaling pathway",
    "NR", "response to type I interferon"
)
 
ratio_num <- function(x) {
    sapply(strsplit(as.character(x), "/"), function(z) {
        as.numeric(z[1]) / as.numeric(z[2])
    })
}

plot_df <- go_sig %>%
    filter(p.adjust < 0.01) %>%
    inner_join(term_map, by = c("direction", "Description")) %>%
    mutate(
        FoldEnrichment = ratio_num(GeneRatio) / ratio_num(BgRatio),
        log10padj = -log10(p.adjust + 1e-300),
        score = if_else(direction == "R", -FoldEnrichment, FoldEnrichment)
    )
 
plot_df <- bind_rows(
    plot_df %>% filter(direction == "R") %>% arrange(desc(FoldEnrichment)),
    plot_df %>% filter(direction == "NR") %>% arrange(FoldEnrichment)
)

plot_df$Description <- factor(
    plot_df$Description,
    levels = rev(plot_df$Description)
)

 
p_go_ora <- ggplot(plot_df, aes(x = score, y = Description)) +
    geom_col(
        aes(fill = log10padj),
        width = 0.75, color = "black", linewidth = 0.25
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_x_continuous(
        labels = abs,
        expand = expansion(mult = c(0.08, 0.08))
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
    scale_fill_gradientn(
        colours = c("#fee6ce", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
        name = expression(-log[10]("adjusted P value"))
    ) +
    labs(
        title = "Representative GO Biological Processes in Fibroblasts",
        subtitle = "Left: R; Right: NR",
        x = "Fold enrichment",
        y = NULL,
        caption = "Selected representative terms; BH-adjusted P < 0.01."
    ) +
    theme_classic(base_size = 11) +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        plot.caption = element_text(hjust = 0)
    )

p_go_ora
ggsave("Figure/scRNA/Fibroblast/GO_BP_0928.pdf", plot = p_go_ora, width = 12, height = 9)

write.csv(plot_df, "Table/scRNA/Fibroblast/GO_BP_selected_padj001.csv", row.names = FALSE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggrepel)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(clusterProfiler)
library(stringr)

scRNA_obj <- readRDS(file = 'Data/scRNA/rawData/icc_scRNA_clean.RDS')
Endo_sub <- subset(scRNA_obj, ct == 'Endothelia') 
marker.gene <- FindMarkers(Endo_sub, ident.1 = "NR", ident.2 = "R", group.by = "response",logfc.threshold = 0, min.pct = 0.05)
deg_gene <- marker.gene
deg_gene$gene <- rownames(deg_gene) 
deg_gene <- deg_gene %>% mutate(change = case_when(p_val_adj < 0.05 & avg_log2FC >= 1  ~  "NR", 
                                                   p_val_adj < 0.05 & avg_log2FC <= -1 ~  "R", TRUE ~ "NS"),
                                neglog10_padj = -log10(p_val_adj + 1e-300))
 
exclude_pattern <- "^(MT-|RPL|RPS|MALAT1$|HBA1$|HBA2$|HBB$|LINC)"
deg_plot <- deg_gene %>% filter(!grepl(exclude_pattern, gene))  

top_up_p <- deg_plot %>% filter(change == "NR") %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice_head(n = 10) %>% pull(gene)
top_up_fc <- deg_plot %>% filter(change == "NR") %>% arrange(desc(avg_log2FC), p_val_adj) %>% slice_head(n = 10)  %>% pull(gene)
top_down_p <- deg_plot %>% filter(change == "R") %>% arrange(p_val_adj, avg_log2FC) %>% slice_head(n = 10)  %>% pull(gene)
top_down_fc <- deg_plot %>% filter(change == "R") %>% arrange(avg_log2FC, p_val_adj) %>%slice_head(n = 10) %>% pull(gene)

deg_plot %>% filter(change == "NR") %>%   pull(gene)
deg_plot %>% filter(change == "NR") %>%  pull(gene)

label_genes <- unique(c(#top_up_label_genes, top_down_label_genes,
    top_up_p,top_down_p))

label_df <- deg_plot %>% filter(gene %in% label_genes)

p_volcano <- ggplot(deg_plot, aes(x = avg_log2FC, y = neglog10_padj, color = change)) +
    geom_point(size = 1.2, alpha = 1) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
    geom_text_repel(data = label_df,
                    aes(label = gene),
                    color = 'black',
                    size = 3.2,
                    max.overlaps = 100,
                    box.padding = 0.35,
                    point.padding = 0.2,
                    segment.color = "black",
                    show.legend = FALSE) + 
    scale_color_manual(values = c("R" = "#fdb462", "NR" = "#80b1d3", "NS" = "grey75")) +
    labs(title = "scRNA",x = "Log2(Fold change)",y = "-log10(Padj)",color = '') +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = 'none',
          panel.grid = element_blank())

p_volcano
ggsave('Figure/scRNA/Endothelial/volcano_cell-level.pdf',plot = p_volcano, width = 5.5,height = 5)
 
format_df <- function(df, direction = c("NR", "R")) {
    df <- as.data.frame(df)
    direction <- match.arg(direction)
    
    if (nrow(df) == 0) return(NULL)
    
    df$plot_q <- if ("qvalue" %in% colnames(df)) df$qvalue else df$p.adjust
    df$plot_q[is.na(df$plot_q)] <- df$p.adjust[is.na(df$plot_q)]
    
    df <- df %>%
        filter(!is.na(plot_q), plot_q > 0) %>%
        arrange(plot_q) %>% 
        mutate(Description = str_wrap(Description, width = 100),
               score = -log10(plot_q),
               direction = direction)
    
    if (direction == "R") {
        df$score <- -df$score
    }
    
    df
}
 
deg_plot <- deg_plot %>% mutate(change = case_when(p_val_adj < 0.05 & avg_log2FC >= 0.5  ~  "NR", 
                                                   p_val_adj < 0.05 & avg_log2FC <= -0.5 ~  "R", TRUE ~ "NS"),
                                neglog10_padj = -log10(p_val_adj + 1e-300))

bg_genes <- rownames(marker.gene )
up_genes <- deg_plot %>% filter(change == "NR") %>% pull(gene) %>% unique()
up_genes

down_genes <- deg_plot %>% filter(change == "R") %>% pull(gene) %>% unique()
down_genes

gene_bg_df <- bitr(bg_genes,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_up_df <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_down_df <- bitr(down_genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

up_entrez <- unique(gene_up_df$ENTREZID)
down_entrez <- unique(gene_down_df$ENTREZID)
bg_entrez <- unique(gene_bg_df$ENTREZID)

ego_up <- enrichGO(gene = up_entrez,universe = bg_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,readable = TRUE)

ego_down <- enrichGO(gene = down_entrez,universe = bg_entrez, OrgDb = org.Hs.eg.db,keyType = "ENTREZID",
                     ont =  "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
 

ratio_to_numeric <- function(x) {
    vapply(
        strsplit(as.character(x), "/", fixed = TRUE),
        function(z) as.numeric(z[1]) / as.numeric(z[2]),
        numeric(1)
    )
}

go_all <- bind_rows(as.data.frame(ego_down) %>% mutate(direction = "R"), as.data.frame(ego_up) %>% mutate(direction = "NR"))
go_sig <- go_all %>%
    filter(!is.na(p.adjust), p.adjust <= 0.01) %>%
    mutate(
        FoldEnrichment = ratio_to_numeric(GeneRatio) / ratio_to_numeric(BgRatio),
        log10padj = -log10(p.adjust + 1e-300)) %>%
    arrange(direction, p.adjust, desc(Count)) %>%
    relocate(direction, ID, Description)
 
go_sig %>% filter(direction == "R") %>% pull(Description)
go_sig %>% filter(direction == "NR") %>% pull(Description)

write.csv(go_sig, file = "Data/scRNA/Endothelial/GO_BP_NR_significant.csv") 

library(dplyr)
library(ggplot2)
library(stringr)

 
term_map <- tibble::tribble(
    ~direction, ~Description,
    "R", "MHC class II protein complex assembly",
    "R", "antigen processing and presentation of peptide antigen via MHC class II",
    "R", "positive regulation of T cell activation",
    "R", "positive regulation of leukocyte cell-cell adhesion",
    "R", "integrin-mediated signaling pathway",
    "R", "actin cytoskeleton organization",
    "R", "regulation of nitric oxide biosynthetic process",

    "NR", "angiogenesis",
    "NR", "endothelial cell migration",
    "NR", "sprouting angiogenesis",
    "NR", "branching involved in blood vessel morphogenesis",
    "NR", "artery development",
    "NR", "cell junction organization",
    "NR", "regulation of cellular response to growth factor stimulus"
)

ratio_num <- function(x) {
    sapply(strsplit(as.character(x), "/"), function(z) {
        as.numeric(z[1]) / as.numeric(z[2])
    })
}

plot_df <- go_sig %>%
    filter(p.adjust < 0.01) %>%
    inner_join(term_map, by = c("direction", "Description")) %>%
    mutate(
        FoldEnrichment = ratio_num(GeneRatio) / ratio_num(BgRatio),
        log10padj = -log10(p.adjust + 1e-300),
        score = if_else(direction == "R", -FoldEnrichment, FoldEnrichment)
    )
 
plot_df <- bind_rows(plot_df %>% filter(direction == "R") %>% arrange(desc(FoldEnrichment)), plot_df %>% filter(direction == "NR") %>% arrange(FoldEnrichment))
plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
p_go_ora <- ggplot(plot_df, aes(x = score, y = Description)) +
    geom_col(
        aes(fill = log10padj),
        width = 0.75, color = "black", linewidth = 0.25
    ) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_x_continuous(
        labels = abs,
        expand = expansion(mult = c(0.08, 0.08))
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
    scale_fill_gradientn(
        colours = c("#fee6ce", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
        name = expression(-log[10]("adjusted P value"))
    ) +
    labs(
        title = "Representative GO Biological Processes in Fibroblasts",
        subtitle = "Left: R; Right: NR",
        x = "Fold enrichment",
        y = NULL,
        caption = "Selected representative terms; BH-adjusted P < 0.01."
    ) +
    theme_classic(base_size = 11) +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"),
        plot.caption = element_text(hjust = 0))
 
ggsave("Figure/scRNA/Endothelial/GO_BP_0928.pdf", plot = p_go_ora)
write.csv(plot_df, "Table/scRNA/Endothelial/GO_BP_selected_padj001.csv", row.names = FALSE)
