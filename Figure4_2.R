 
 
source('Code/function/seurat_function.R')

library(ggplot2)
library(ggsci)
library(Seurat)
library(ggpubr)
library(dplyr)
library(harmony)
library(magrittr)
library(cowplot)
library(reshape2)
library(scales)  
library(ggrepel) 

scRNA_obj = readRDS('Data/scRNA/rawData/icc_scRNA_clean.RDS')

ct_df <- as.data.frame(table(scRNA_obj$ct))
colnames(ct_df) <- c("Celltype", "Number")
ct_df <- ct_df %>% arrange(Number) %>% mutate(Celltype = factor(Celltype, levels = Celltype))
p_ct <- ggplot(ct_df, aes(x = Celltype, y = Number, fill = Celltype)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = Number), hjust = -0.1, size = 4) +
    coord_flip() +
    scale_fill_igv() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Cell Number") +
    theme_few() +
    theme(legend.position = "none",
          panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))
 
ggsave('Figure/scRNA/AllCells/1_CellNumbers.pdf',plot = p_ct)
 

Fibro_sub <- subset(scRNA_obj, ct == 'Fibroblast')
Fibro_sub <- seurat_process1(Fibro_sub, harmony = T, vars = c("Sample"), npcs = 15, res=0.5) 
p1 <- DimPlot(Fibro_sub,reduction='umap', group.by ="celltype", pt.size = 0.2, label = T) + scale_color_igv()
p1
ggsave(filename = 'Figure/scRNA/Fibroblast/CellTyping/1_UMAP_Plot_Fibroblast.pdf', p1,width = 5,height = 4)

saveRDS(Fibro_sub,file = 'Data/scRNA/Fibroblast/1_Fibroblast_obj.Rds')

pan_CAF <- c("ACTA2","COL1A2","PDGFRB","IL6","MCAM","RGS5",'POSTN', 'COL6A3','FN1',
             'FBLN1', 'C3', 'IGF1','HLA-DRB1', 'CD74', 'HLA-DRA',"KRT19","KRT8","SAA1")

p2 <- DotPlot(Fibro_sub,features = pan_CAF,group.by = 'celltype') + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p2
ggsave(filename = 'Figure/scRNA/Fibroblast/CellTyping/2_Dot_Plot_Fibroblast.pdf', p2, width = 7.5,height = 3.5)


T_obj <- subset(scRNA_obj, ct == 'Tcell') 
T_obj@meta.data$T_main <- as.character(T_obj@meta.data$celltype)
T_obj@meta.data$T_main <- dplyr::case_when(grepl("^CD4", T_obj@meta.data$celltype) ~ "CD4_Tcell",
                                               grepl("^CD8", T_obj@meta.data$celltype) ~ "CD8_Tcell",
                                               TRUE ~ T_obj@meta.data$celltype)

T_obj$T_main <- factor(T_obj$T_main, levels = c("CD4_Tcell", "CD8_Tcell", "gamma_delt_T", "NK"))
T_obj <- seurat_process1(T_obj, harmony = T, vars = c("Sample"), npcs = 15, res=0.5) 
saveRDS(T_obj,file = 'Data/scRNA/Tcell/1_Tcell_obj.Rds')

p = DimPlot(T_obj,group.by = 'T_main',pt.size = 0.01) + scale_color_manual(values = c('#80b1d3','#fdb462','#8dd3c7','#bebada'))
ggsave('Figure/scRNA/Tcell/1_UMAP_Plot_Tcell.pdf',plot = p,width = 5.5,height = 4.5)

p2 = FeaturePlot(T_obj,features = c("CD3G","CD3D","CD8A","CD8B","CD4","TRDC","NKG7","GNLY") ,ncol = 4)
ggsave(filename = 'Figure/scRNA/Tcell/2_Tmain_markerGenes.pdf', p2, width = 11, height = 5)



library(dplyr)
T_obj$TCRFreq_bin <- case_when(is.na(T_obj$TCRFreq) ~ "[0,1)",
                               T_obj$TCRFreq < 1 ~ "[0,1)",
                               T_obj$TCRFreq >= 1 & T_obj$TCRFreq < 2 ~ "[1,2)",
                               T_obj$TCRFreq >= 2 & T_obj$TCRFreq < 3 ~ "[2,3)",
                               T_obj$TCRFreq >= 3 & T_obj$TCRFreq < 4 ~ "[3,4)",
                               T_obj$TCRFreq >= 4 & T_obj$TCRFreq < 5 ~ "[4,5)",
                               T_obj$TCRFreq >= 5 & T_obj$TCRFreq < 6 ~ "[5,6)",
                               T_obj$TCRFreq >= 6 & T_obj$TCRFreq < 10 ~ "[6,10)",
                               T_obj$TCRFreq >= 10 ~ "[10,Inf)")

T_obj$TCRFreq_bin <- factor(T_obj$TCRFreq_bin, levels = c("[0,1)", "[1,2)", "[2,3)", "[3,4)", "[4,5)", "[5,6)", "[6,10)", "[10,Inf)"))

tcr_bin_cols <- c("[0,1)" = "#e0e0e0", "[1,2)" = "#fee5d9", "[2,3)" = "#fcbba1", "[3,4)" = "#fc9272", "[4,5)" = "#fb6a4a",
                  "[5,6)" = "#ef3b2c", "[6,10)" = "#cb181d", "[10,Inf)" = "#99000d")

T_obj_pre <- subset(T_obj,Treatment =='pre')
p1 <- DimPlot(T_obj_pre, group.by = "TCRFreq_bin", pt.size = 0.001) + scale_color_manual(values = tcr_bin_cols) + theme_void() +
    labs(title = 'TCR Frequency')+
    theme(plot.title = element_text(hjust = 0.5))
p1
ggsave('Figure/scRNA/Tcell/1_Clone_expr_Tcell.pdf',width = 5,height = 4.5)

sub_Tcell <- subset(T_obj, T_main != 'NK' & T_main != 'gamma_delt_T')
sub_Tcell <- subset(sub_Tcell,Treatment =='pre')
sub_Tcell <- subset(sub_Tcell,!is.na(TCRFreq))

metaInfo <- sub_Tcell@meta.data
cb_response = c(NR = '#9ecae1', R = '#238b45')

sampleClone <- metaInfo %>% group_by(Sample, response) %>% 
    summarise(num = sum(TCRFreq >= 2, na.rm = TRUE) / n(), .groups = "drop")%>%
    as.data.frame() %>% 
    arrange(desc(num)) %>%
    mutate(sample = factor(Sample, levels = Sample))

mean_num <- mean(sampleClone$num, na.rm = TRUE)
p1_r <- ggplot(data = sampleClone, mapping = aes(x = sample,y = num, fill = response)) +
    geom_hline(yintercept = mean_num,
               linetype = "dashed",
               color = "black",
               linewidth = 0.4) +
    annotate("text", x = Inf, y = mean_num, label = paste0("Mean = ", round(mean_num,3)), hjust = 1.1, vjust = -0.4, size = 4) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = cb_response) +
    theme_bw() +
    labs(x = NULL, y = "% cells in expanded clonotypes", title = "") +
    theme(axis.text.x = element_text(angle = 45, hjust =1, vjust = 1)) 
p1_r
ggsave('Figure/scRNA/Tcell/Expan/Expan_TCR_response.pdf',plot = p1_r,width = 7,height = 6)
 
 


Myeloid_sub <- subset(scRNA_obj, ct == 'Myeloid')
Myeloid_sub <- seurat_process1(Myeloid_sub, harmony = T, vars = c("Sample"), npcs = 20, res=0.5) 

ct_map <- c("cDC_CD1C" = "cDC",
            "cDC_CLEC9A" = "cDC",
            "Macro_C1QC" = "Macrophage",
            "Mono_CD14" = "Monocyte",
            "Mast" = "Mast",
            "Neutrophil" = "Neutrophil")

meta_df <- Myeloid_sub@meta.data

meta_df$ct_raw <- as.character(meta_df$celltype )
meta_df$ct_major <- unname(ct_map[meta_df$ct_raw])
Myeloid_sub$ct_major <- meta_df$ct_major

p1 <- DimPlot(Myeloid_sub,reduction='umap', group.by ="ct_major", pt.size = 0.01, label = T) + 
    theme(legend.position = 'none') + 
    scale_color_jama()
ggsave(filename = 'Figure/Myeloid/1_UMAP_Plot_Myeloid.pdf', p1, width = 5, height = 5)

Mast=c("KIT","TPSAB1","CPA3","MS4A2") 
cDC=c("CLEC9A","XCR1", "BATF3","FLT3","CD1C")
Monotypes=c("CD14","FCN1","LYZ", "VCAN","S100A9","S100A8") 
Macrophage=c("C1QC","C1QA","APOE")
Neutrophil=c("S100A8","CXCR2","FCGR3B","CSF3R")

totalGene = unique(c(Mast,cDC,Macrophage,Monotypes,Neutrophil))
p2 <- DotPlot(Myeloid_sub,features = totalGene, group.by = 'ct_major') + RotatedAxis()
ggsave(filename = 'Figure/scRNA/Myeloid/2_DotPlot_markergenes.pdf', p2, width = 9, height = 3)
saveRDS(Myeloid_sub, file = 'Data/scRNA/Myeloid/1_Myeloid_obj.Rds')


# ------------------------------ Endothelial 
Endothelial_sub <- subset(scRNA_obj, ct == 'Endothelia') 
Endothelial_sub <- seurat_process1(Endothelial_sub, harmony = T, vars = c("Sample"), npcs = 20, res=0.5) 
p1 <- DimPlot(Endothelial_sub,reduction='umap', group.by ="celltype", pt.size = 0.01, label = T) + theme(legend.position = 'none') + scale_color_igv()
ggsave(filename = 'Figure/scRNA/Endothelial/1_UMAP_Plot_Endothelial.pdf', p1, width = 5, height = 5)
saveRDS(Endothelial_sub, file = 'Data/scRNA/Endothelial/1_Endothelial_obj.Rds')
