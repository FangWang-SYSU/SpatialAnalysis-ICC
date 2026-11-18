rm(list = ls())
gc()

library(UCell)
library(Seurat)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(patchwork)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(fgsea)
library(corrplot)


setwd("./")
obj <- readRDS("Data/scRNA/rawData/icc_scRNA_clean")
celltype_labels <- c(B_cell="B cells", Plasma="Plasma cells", 
                     CD4_Tcell="CD4 T", CD8_Tcell="CD8 T", gamma_delta_T="γδ T", NK="NK",
                     DC="DC", Macrophages="Macrophages", Mast="Mast cells", Monotype="Monocytes", Neutrophil="Neutrophils",
                     Endo_C0_COL4A1="COL4A1+ EC", Endo_C1_ACKR1="ACKR1+ EC", Endo_C2_FCN3="FCN3+ EC", Endo_C3_FBLN5="FBLN5+ EC", 
                     Endo_C4_CCL21="CCL21+ EC",
                     Fibor_apCAF="apCAF", Fibro_eCAF="eCAF", Fibro_iCAF="iCAF", Fibro_mCAF="mCAF", Fibro_MKI="Cycling CAF", Fibro_vCAF="vCAF",
                     Hepatocytes="Hepatocytes", Tumor="Tumor")

original_celltype <- as.character(obj$celltype_s1)
adjusted_celltype <- unname(celltype_labels[original_celltype])
adjusted_celltype[is.na(adjusted_celltype)] <- original_celltype[is.na(adjusted_celltype)]
obj$celltype_plot <- adjusted_celltype

cb_celltype = c('ACKR1+ EC' = '#35978f','CCL21+ EC' = '#2196F3',apCAF = '#FFD6E0',
                'mCAF' = '#A6E3E9','COL4A1+ EC' = '#FFDC91FF' ,'Cycling CAF' = '#FFA6FB',
                'eCAF' = '#80cdc1', 'iCAF'  = '#FF80C7','DC' = '#CAAA98', '
                FBLN5+ EC' = '#E37222FF', 'FCN3+ EC' = '#F0AB00FF', 'Hepatocytes' = '#34A99D',
                'Macrophages' = '#525EA7', 'Mast cells' = '#8FA28A', 'Monocytes' = '#FF9494',
                'Neutrophils' = '#A45FE6',
                'B cells' = '#F4D35E','CD4 T' = '#b3de69','NK' = '#66BB6A', 
                'Plasma cells' = '#C06C84','Tumor' = '#808180FF',
                'vCAF' = '#7D5CC6FF', 'γδ T' = '#AF719D',
                'CD8 T' = '#fbb4ae')

p = DimPlot(obj,group.by = 'celltype_plot') + scale_color_manual(values = cb_celltype) + ggtitle("")+
    theme_few()
ggsave('Figure/scRNA/SEP/Figure4_scRNA_all_celltypes.pdf',p,width = 8,height = 5.5)
p = DimPlot(obj,group.by = 'response') + scale_color_manual(values = c('#FFD758','#218DAE')) + 
    theme_few()
ggsave('Figure/scRNA/SEP/Figure4_scRNA_all_celltypes_response.pdf',p,width = 6.5,height = 5.5)

p = DimPlot(obj,group.by = 'Sample') + ggtitle("") + 
    scale_color_igv()+
    theme_few()
ggsave('Figure/scRNA/SEP/Figure4_scRNA_all_celltypes_Sample.pdf',p,width = 8,height = 5.5)

marker_list <- list(
    `B cells`       = c("MS4A1", "CD79A"),
    `Plasma cells`  = c("MZB1", "JCHAIN"),
    `T cells`       = c("CD3D", "CD3E",'CD3G', "CD4", "IL7R",'TRGC1','TRGC2', "CD8A", "CD8B", "TRDC"),
    `NK`            = c("NKG7", "GNLY", "KLRD1"),
    `Myeloid`       = c("LYZ", "LST1", "FCN1", "C1QA", "CD68"),
    `DC`            = c("CD1C", "FCER1A", "CLEC9A"),
    `Mast cells`    = c("TPSAB1", "KIT"),
    `Neutrophils`   = c("FCGR3B", "CSF3R"),
    `Endothelial`   = c("PECAM1", "VWF", "COL4A1", "ACKR1", "FCN3", "FBLN5", "CCL21"),
    `Fibroblast`    = c("DCN", "COL1A1", "COL1A2", "FAP", "ACTA2", "TAGLN", "CXCL12"),
    `Antigen presentation` = c("CD74", "HLA-DRA","HLA-DRB1"),
    `Proliferation` = c("MKI67", "TOP2A"),
    `Hepatocytes`   = c("ALB", "APOA1"),
    `Epithelial`    = c("EPCAM", "KRT19", "KRT8")
)

genes_available <- rownames(obj[["RNA"]])
marker_use <- lapply(marker_list, function(x) intersect(x, genes_available))
marker_use <- marker_use[lengths(marker_use) > 0]

celltype_present <- unique(as.character(obj$celltype_plot))
celltype_present <- celltype_present[!is.na(celltype_present) & nzchar(celltype_present)]
celltype_order <- c(unname(celltype_labels)[unname(celltype_labels) %in% celltype_present], setdiff(celltype_present, unname(celltype_labels)))
celltype_order <- unique(celltype_order)

p_dot <-  DotPlot(object = obj, assay = "RNA", features =  unique(unlist(marker_use, use.names = FALSE)), 
                  group.by = "celltype_plot", scale    = TRUE, dot.min  = 0, dot.scale = 5, col.min = -1,col.max = 2) +
    scale_y_discrete(limits = rev(celltype_order)) +
    scale_color_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), limits = c(-1, 2), oob = scales::squish, name = "Scaled mean\nexpression") +
    theme_few() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 9),
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 90, size = 8, face = "bold"),
        panel.spacing.x = grid::unit(0.12, "cm"),
        legend.position = "right")

ggsave("Figure/scRNA/SEP/FigureS6_celltype_marker_DotPlot.pdf", plot = p_dot, width = 12, height = 5.63 )


MP_list <- readRDS("Data/ST/SEP/Module1/03_metaprograms/all_MP_list_New.RDS")
MP_use <- lapply(MP_list, function(x) {
    intersect(unique(x), rownames(obj))
})

gene_check <- data.table(MP = names(MP_use), n_input = lengths(MP_list), n_found = lengths(MP_use))
obj_MPscore <- AddModuleScore_UCell(obj, features = MP_use)
#saveRDS(obj_MPscore,file = 'Data/scRNA/SEP/scRNA_add_MPscore.RDS')
obj_use =  readRDS('Data/scRNA/SEP/scRNA_add_MPscore.RDS')

meta <- as.data.table(obj_use@meta.data, keep.rownames="Cell")
score_cols <- intersect(paste0(names(MP_use), "_UCell"), colnames(meta))

plot_names <- gsub("_UCell$", "", score_cols)
plot_names <- gsub("_", "-", plot_names)

heatmap_wide <- meta[, lapply(.SD, mean, na.rm=TRUE), by=.(CellGroup=celltype_plot), .SDcols=score_cols]
heatmap_long <- melt(heatmap_wide, id.vars="CellGroup", measure.vars=score_cols, variable.name="MP", value.name="Score")
heatmap_long[, MP:=sub("_UCell$", "", MP)]
heatmap_long <- heatmap_long[is.finite(Score)]

heatmap_long[, z_score:={
    score_sd <- sd(Score, na.rm=TRUE)
    if(!is.finite(score_sd) || score_sd==0) rep(0, .N) else (Score-mean(Score, na.rm=TRUE))/score_sd
}, by=MP]

heatmap_mat_dt <- dcast(heatmap_long, MP~CellGroup, value.var="z_score", fill=0)
heatmap_mat <- as.matrix(heatmap_mat_dt[, -1])
rownames(heatmap_mat) <- heatmap_mat_dt$MP

row_sd <- apply(heatmap_mat, 1, sd, na.rm=TRUE)
heatmap_mat <- heatmap_mat[is.finite(row_sd) & row_sd>0, , drop=FALSE]
heatmap_mat <- pmax(pmin(heatmap_mat, 2.5), -2.5)

celltype_order <- c("B cells", "Plasma cells", "CD4 T", "CD8 T", "γδ T", "NK",
                    "DC", "Macrophages", "Mast cells", "Monocytes", "Neutrophils",
                    "COL4A1+ EC", "ACKR1+ EC", "FCN3+ EC", "FBLN5+ EC", "CCL21+ EC",
                    "apCAF", "eCAF", "iCAF", "mCAF", "Cycling CAF", "vCAF",
                    "Hepatocytes", "Tumor")

celltype_order <- c(intersect(celltype_order, colnames(heatmap_mat)), setdiff(colnames(heatmap_mat), celltype_order))
heatmap_mat <- heatmap_mat[, celltype_order, drop=FALSE]

cell_class_map <- c(
    "B cells"="Lymphoid", "Plasma cells"="Lymphoid", "CD4 T"="Lymphoid", "CD8 T"="Lymphoid", "γδ T"="Lymphoid", NK="Lymphoid",
    DC="Myeloid", Macrophages="Myeloid", "Mast cells"="Myeloid", Monocytes="Myeloid", Neutrophils="Myeloid",
    "COL4A1+ EC"="Endothelial", "ACKR1+ EC"="Endothelial", "FCN3+ EC"="Endothelial", "FBLN5+ EC"="Endothelial", "CCL21+ EC"="Endothelial",
    apCAF="Fibroblast", eCAF="Fibroblast", iCAF="Fibroblast", mCAF="Fibroblast", "Cycling CAF"="Fibroblast", vCAF="Fibroblast",
    Hepatocytes="Epithelial", Tumor="Epithelial"
)

cell_classes <- unname(cell_class_map[colnames(heatmap_mat)])
cell_classes[is.na(cell_classes)] <- "Other"
cell_classes <- factor(cell_classes, levels=c("Lymphoid", "Myeloid", "Endothelial", "Fibroblast", "Epithelial", "Other"))

cell_class_colors <- c(Lymphoid="#DF301C", Myeloid="#FFCC4D", Endothelial="#30AFFF", Fibroblast="#FFABAB",Epithelial="#999999")

mp_names <- rownames(heatmap_mat)
is_tumor_mp <- grepl("^M[_-]?SEP|^Tumor[_-]?MP", mp_names, ignore.case=TRUE)
mp_groups <- factor(ifelse(is_tumor_mp, "Malignant SEPs", "Normal SEPs"), levels=c("Normal SEPs", "Malignant SEPs"))

mp_group_colors <- c("Normal SEPs"="#91C8E4", "Malignant SEPs"="#BBBFCA")

sep_short_names <- c(N_SEP1="POSTN+CAF",
                     N_SEP2="HLA-II+Immune",
                     N_SEP3="Hep-Secretory",
                     N_SEP4="Vascular EC",
                     N_SEP5="TLS Plasma-Stromal",
                     N_SEP6="Contractile VSMC",
                     N_SEP7="Xeno-Hep",
                     N_SEP8="Angiogenic EC",
                     N_SEP9="TLS B-T",
                     N_SEP10="Sinusoidal Niche",
                     N_SEP11="Pericyte-Mural",
                     N_SEP12="COL11A1+CAF",
                     M_SEP1="Abs-Met",
                     M_SEP2="LAMC2/LAMB3 Basal-Inv",
                     M_SEP3="S100P+Muc-Sec",
                     M_SEP4="SPP1+Trans-Sec",
                     M_SEP5="Cycling-WNT" )

display_rownames <- unname(sep_short_names[rownames(heatmap_mat)])
display_rownames[is.na(display_rownames)] <- rownames(heatmap_mat)[is.na(display_rownames)]

standard_mp_names <- gsub("-", "_", rownames(heatmap_mat))
matched_labels <- unname(sep_short_names[standard_mp_names])
display_rownames[!is.na(matched_labels)] <- matched_labels[!is.na(matched_labels)]

top_annotation <- HeatmapAnnotation(CellType=cell_classes, 
                                    col=list(CellType=cell_class_colors),
                                    annotation_name_side="left",
                                    annotation_name_gp=gpar(fontsize=9, fontface="bold"),
                                    simple_anno_size=unit(4, "mm"),
                                    show_legend=TRUE)

left_annotation <- rowAnnotation(SEPs=mp_groups, col=list(SEPs=mp_group_colors),
                                 annotation_name_rot=90,
                                 annotation_name_gp=gpar(fontsize=9, fontface="bold"),
                                 simple_anno_size=unit(4, "mm"), show_legend=TRUE)

heatmap_col_fun <- colorRamp2(c(-2.5, -1.25, 0, 1.25, 2.5),c("#244A7C", "#8BB3CF", "#FFFFFF", "#E59A91", "#9E3038"))

ht <- Heatmap(heatmap_mat,
              name="Z-score",
              col=heatmap_col_fun,
              na_col="#EFEFEF",
              
              row_split=mp_groups,
              cluster_rows=TRUE,
              cluster_row_slices=FALSE,
              clustering_distance_rows="pearson",
              clustering_method_rows="average",
              row_gap=unit(1, "mm"),
              
              column_split=cell_classes,
              cluster_columns=FALSE,
              cluster_column_slices=FALSE,
              column_gap=unit(1, "mm"),
              
              top_annotation=top_annotation,
              left_annotation=left_annotation,
              
              #row_labels=display_rownames,
              show_row_names=TRUE,
              row_names_side="right",
              row_names_gp=gpar(fontsize=9),
              row_names_max_width=unit(48, "mm"),
              
              show_column_names=TRUE,
              column_names_rot=90,
              column_names_gp=gpar(fontsize=9),
              column_names_centered=FALSE,
              
              row_title_gp=gpar(fontsize=11, fontface="bold"),
              row_title_rot=0,
              column_title_gp=gpar(fontsize=10, fontface="bold"),
              
              rect_gp=gpar(col="white", lwd=0.8),
              
              cell_fun=function(j, i, x, y, width, height, fill){
                  if(heatmap_mat[i, j] >= 0.5){
                      grid.text(
                          "*",
                          x, y,
                          gp=gpar(fontsize=10, fontface="bold")
                      )
                  }
              },
              
              heatmap_legend_param=list(
                  title="Z-score",
                  at=c(-2.5, -1, 0, 1, 2.5),
                  labels=c("-2.5", "-1", "0", "1", "2.5"),
                  legend_height=unit(35, "mm"),
                  title_gp=gpar(fontsize=10, fontface="bold"),
                  labels_gp=gpar(fontsize=9)
              ),
              
              column_title="Metaprogram activity across cell types"
)

pdf(file = 'Figure/scRNA/SEP/Figure2_scRNA_heatmap_SEPs_celltypes.pdf',width = 10,height = 4.5)
draw(ht, heatmap_legend_side="right", annotation_legend_side="right",  merge_legends=FALSE, padding=unit(c(5, 5, 5, 5), "mm"))
dev.off() 

 
obj_use <- subset(obj_use,Treatment=="pre")
response_levels <- c("NR","R")
min_cells_per_sample <- 10
min_samples_per_group <- 3

mp_order <- c(paste0("N_SEP",1:12),paste0("M_SEP",1:5))
n_sep <- paste0("N_SEP",1:12)
m_sep <- paste0("M_SEP",1:5)
mp_cols <- paste0(mp_order,"_UCell")

panel_specs <- list(Nontumor=n_sep,Tumor=m_sep)

obj_use$response <- toupper(trimws(as.character(obj_use$response)))
obj_use$Sample <- trimws(as.character(obj_use$Sample))
obj_use$Panel <- ifelse(obj_use$celltype_s1=="Tumor","Tumor","Nontumor")

meta <- as.data.table(obj_use@meta.data,keep.rownames="CellID")
valid_cells <- meta[!is.na(Sample) & Sample!="" & !is.na(response) & response %chin% response_levels,CellID]
obj_use <- subset(obj_use,cells=valid_cells)

meta <- as.data.table(obj_use@meta.data,keep.rownames="CellID")

# Prepare panel
prepare_panel <- function(panel){
    panel_meta <- meta[Panel==panel]
    sample_n <- panel_meta[,.N,by=.(Sample,response)]
    keep_samples <- sample_n[N>=min_cells_per_sample,Sample]
    panel_meta <- panel_meta[Sample %chin% keep_samples]
    
    sample_info <- unique(panel_meta[,.(Sample,Response=response)])[order(Sample)]
    sample_info[,Response:=factor(Response,levels=response_levels)]
    
    n_NR <- sample_info[Response=="NR",.N]
    n_R <- sample_info[Response=="R",.N]
    
    if(n_NR<min_samples_per_group || n_R<min_samples_per_group)
        stop(panel,": insufficient patients after filtering; NR = ",n_NR,", R = ",n_R)
    
    list(meta=panel_meta,sample_info=sample_info,sample_n=sample_n,n_NR=n_NR,n_R=n_R)
}

make_pathways <- function(target_MPs,genes,panel){
    missing_target <- setdiff(target_MPs,names(MP_list))
    if(length(missing_target)>0) stop(panel,": missing MPs in MP_list: ",paste(missing_target,collapse=", "))
    
    pathways <- lapply(MP_list[target_MPs],function(x) intersect(unique(x),genes))
    gene_check <- data.table(Panel=panel,MP=names(pathways),n_original=lengths(MP_list[target_MPs]),n_tested=lengths(pathways))
    gene_check[,`:=`(fraction_tested=n_tested/n_original,retained_for_fgsea=n_tested>=10)]
    
    pathways <- pathways[lengths(pathways)>=10]
    if(length(pathways)==0) stop(panel,": no MP retains at least 10 genes")
    
    list(pathways=pathways,gene_check=gene_check)
}

run_fgsea <- function(pathways,rank_stat,panel,method,panel_data){
    set.seed(42)
    
    fg <- as.data.table(fgseaMultilevel(pathways=pathways,stats=rank_stat,minSize=10,maxSize=500,eps=0))
    setnames(fg,"pathway","MP")
    
    fg[,`:=`(
        Panel=panel,Method=method,direction=fifelse(NES>0,"R enriched","NR enriched"),
        n_NR_samples=panel_data$n_NR,n_R_samples=panel_data$n_R,
        n_NR_cells=panel_data$meta[response=="NR",.N],n_R_cells=panel_data$meta[response=="R",.N],
        leadingEdge_genes=vapply(leadingEdge,paste,collapse=",",FUN.VALUE=character(1))
    )]
    
    fg[]
}

run_cell_gsea <- function(panel,target_MPs){
    panel_data <- prepare_panel(panel)
    panel_obj <- subset(obj_use,cells=panel_data$meta$CellID)
    Idents(panel_obj) <- "response"
    
    de <- FindMarkers(panel_obj,ident.1="R",ident.2="NR",assay="RNA",test.use="wilcox",
                      min.pct=0,logfc.threshold=0,only.pos=FALSE,verbose=FALSE)
    
    de <- as.data.table(de,keep.rownames="gene")
    fc_col <- intersect(c("avg_log2FC","avg_logFC"),names(de))[1]
    if(length(fc_col)==0 || is.na(fc_col)) stop(panel,": cannot find avg_log2FC/avg_logFC")
    
    de[,rank_stat:=get(fc_col)]
    de <- de[is.finite(rank_stat) & !duplicated(gene)]
    rank_stat <- sort(setNames(de$rank_stat,de$gene),decreasing=TRUE)
    
    pathway_data <- make_pathways(target_MPs,names(rank_stat),panel)
    
    list(gsea=run_fgsea(pathway_data$pathways,rank_stat,panel,"Cell_Wilcoxon",panel_data),
         DE=de,gene_check=pathway_data$gene_check,panel_data=panel_data)
}


# Run GSEA
cell_results <- lapply(names(panel_specs),function(panel) run_cell_gsea(panel,panel_specs[[panel]]))
names(cell_results) <- names(panel_specs)

cell_summary <- rbindlist(lapply(cell_results,`[[`,"gsea"),fill=TRUE)
cell_summary[,padj_global:=p.adjust(pval,method="BH")]

# GSEA plot
plot_gsea <- function(result,title){
    d <- copy(result)
    
    d[,SEP_type:=fifelse(grepl("^N_SEP",MP),"N-SEPs","M-SEPs")]
    d[,SEP_type:=factor(SEP_type,levels=c("N-SEPs","M-SEPs"))]
    d[,shape_type:=SEP_type]
    
    d[,p_group:=fcase(
        is.na(padj_global) | padj_global>=0.05,"p.adj >= 0.05",
        padj_global>=0.01,"p.adj < 0.05",
        padj_global>=0.001,"p.adj < 0.01",
        padj_global>=0.0001,"p.adj < 0.001",
        default="p.adj < 0.0001"
    )]
    
    d[,p_group:=factor(p_group,levels=c("p.adj > 0.05","p.adj < 0.05","p.adj < 0.01","p.adj < 0.001","p.adj < 0.0001"))]
    d[,padj_plot:=pmax(fifelse(is.na(padj_global),1,padj_global),1e-4)]
    d[,significance:=-log10(padj_plot)]
    
    setorder(d,SEP_type,NES)
    d[,MP:=factor(MP,levels=rev(MP))]
    
    ggplot(d,aes(y=MP,x=NES)) +
        geom_vline(xintercept=0,linetype=1,linewidth=0.45,color="grey45") +
        geom_vline(xintercept=1.5,linetype=2,linewidth=0.45,color="grey45") +
        geom_vline(xintercept=-1.5,linetype=2,linewidth=0.45,color="grey45") +
        
        geom_segment(aes(x=0,xend=NES,yend=MP),linewidth=0.65,color="grey35") +
        
        geom_point(aes(fill=p_group,size=significance,shape=shape_type),color="grey15",stroke=0.55) +
        scale_shape_manual(values=c("N-SEPs"=21,"M-SEPs"=23),name=NULL) +
        scale_fill_manual(values=c(
            "p.adj > 0.05"="#92c5de",
            "p.adj < 0.05"="#fddbc7",
            "p.adj < 0.01"="#f4a582",
            "p.adj < 0.001"="#d6604d",
            "p.adj < 0.0001"="#b2182b"
        ),name="Significance",drop=FALSE) +
        scale_size_continuous(range=c(2.5,7),name="Adjusted p value",
                              breaks=c(-log10(0.05),2,3,4),labels=c("0.05","0.01","0.001","<=0.0001")) +
        guides(fill=guide_legend(order=1,override.aes=list(shape=21,size=5,color="grey15")),
               size=guide_legend(order=2,override.aes=list(shape=21,fill="grey75",color="grey15")),
               shape=guide_legend(order=3,override.aes=list(size=5,fill="white"))) +
        scale_x_continuous(expand=expansion(mult=c(0.12,0.12))) +
        labs(x="GSEA NES",y=NULL,title=title) +
        theme_classic(base_size=12) +
        theme(legend.position="right",legend.box="vertical",
              axis.line.y=element_blank(),axis.ticks.y=element_blank(),
              axis.text.y=element_text(color="black"),
              plot.title=element_text(hjust=0.5,face="bold"))
}

p_cell <- plot_gsea(cell_summary,"R versus NR GSEA Enrichment")
ggsave("Figure/scRNA/SEP/response/Figure4_scRNA_response_GSEA.pdf",p_cell,width=6,height=5.5)


obj_use <- readRDS("Data/scRNA/SEP/scRNA_add_MPscore.RDS")

score_cols <- intersect(paste0(names(MP_list),"_UCell"),colnames(obj_use@meta.data))
plot_names <- gsub("_","-",sub("_UCell$","",score_cols))

m_score_cols <- c("M_SEP1_UCell","M_SEP2_UCell", "M_SEP3_UCell","M_SEP5_UCell")
m_score_names <- gsub("_","-",sub("_UCell$","",m_score_cols))

combined_score_cols <- c("M_sensitive_score","M_resistance_score","M_response_score")
combined_score_names <- c("M-sensitive score","M-resistance score","M-response score")

sep_colors <- c("#E8E8E8","#F6D58A","#E98B58","#C7473A","#7F1D1D")

response_colors <- c(low="#2166AC",mid="#F2F2F2",high="#B2182B")
raw_score_range <- c(0,0.5)
score_theme <- theme_void() +
    theme(plot.title=element_text(hjust=0.5,size=10,face="bold"),
          legend.title=element_text(size=9),
          legend.text=element_text(size=8))

p_list <- FeaturePlot(obj_use,features=score_cols,combine=FALSE,pt.size=0.1)
for(i in seq_along(p_list)){
    p_list[[i]] <- p_list[[i]] +
        scale_color_gradientn(colors=sep_colors,limits=raw_score_range,
                              oob=scales::squish,name="UCell score") +
        ggtitle(plot_names[i]) + score_theme +
        theme(plot.title=element_text(size=11))
}

p <- wrap_plots(p_list,ncol=6,guides="collect") & theme(legend.position="right")
ggsave("Figure/scRNA/Figure4B_scRNA_allCells_SEP_score.pdf",p,width=10.5,height=5.5)


tumor_obj <- subset(obj_use,celltype_s1=="Tumor")

for(i in seq_along(m_score_cols)){
    score_name <- sub("_UCell$","",m_score_cols[i])
    tumor_obj[[paste0(score_name,"_zscore")]] <-
        as.numeric(scale(tumor_obj@meta.data[[m_score_cols[i]]]))
}

tumor_obj$M_sensitive_score <- (tumor_obj$M_SEP3_zscore + tumor_obj$M_SEP5_zscore) / 2
tumor_obj$M_resistance_score <- (tumor_obj$M_SEP1_zscore + tumor_obj$M_SEP2_zscore) / 2

tumor_obj$M_response_score <- tumor_obj$M_resistance_score - tumor_obj$M_sensitive_score

run_umap_cluster <- function(seu_obj){
    seu_obj <- NormalizeData(seu_obj,verbose=FALSE)
    seu_obj <- FindVariableFeatures(seu_obj,nfeatures=2000,verbose=FALSE)
    features <- VariableFeatures(seu_obj)
    seu_obj <- ScaleData(seu_obj,features=features,verbose=FALSE)
    seu_obj <- RunPCA(seu_obj,features=features,verbose=FALSE)
    seu_obj <- FindNeighbors(seu_obj,reduction="pca",dims=1:30,verbose=FALSE)
    seu_obj <- FindClusters(seu_obj,resolution=0.5,verbose=FALSE)
    seu_obj <- RunUMAP(seu_obj,reduction="pca",dims=1:30,seed.use=42,verbose=FALSE)
    return(seu_obj)
}

all_tumor_obj <- run_umap_cluster(tumor_obj) 
p1 <- DimPlot(all_tumor_obj,group.by="response") +
    scale_color_manual(values=c(R="#FFD758",NR="#218DAE")) + theme_few()
p2 <- DimPlot(all_tumor_obj,group.by="Sample") + scale_color_igv() + theme_few()

ggsave('Figure/scRNA/Figure4E_tumor_scRNA_res_sample.pdf',p1, width = 4.7,height = 4)
ggsave('Figure/scRNA/FigureS11F_tumor_scRNA_sample.pdf',p2, width = 4.7,height = 4)


tumor_obj_list <- list(All=all_tumor_obj, 
                       R=subset(all_tumor_obj,response=="R"),
                       NR=subset(all_tumor_obj,response=="NR"))

sr_score_cols <- c("M_resistance_score", "M_sensitive_score")
sr_score_names <- c("M-resistance score", "M-sensitive score")

sr_limit <- unname(quantile(abs(as.numeric(as.matrix(tumor_obj@meta.data[,sr_score_cols,drop=FALSE]))),
                            probs=0.99,
                            na.rm=TRUE))

sr_score_range <- c(-sr_limit,sr_limit)

panel_score_cols <- c(m_score_cols,sr_score_cols)
panel_score_names <- c(m_score_names,sr_score_names)
panel_groups <- c("R","NR") 

tumor_umap <- Embeddings(all_tumor_obj,reduction="umap")[,1:2,drop=FALSE]
panel_xlim <- range(tumor_umap[,1],finite=TRUE)
panel_ylim <- range(tumor_umap[,2],finite=TRUE)

panel_theme <- score_theme +
    theme(panel.border=element_rect(colour="grey50",fill=NA,linewidth=0.35),
          plot.title=element_text(hjust=0.5,size=10,face="plain",
                                  margin=margin(b=6)),
          plot.margin=margin(0,0,0,0))

p_MSEP_list <- list()

for(j in panel_groups){
    
    p_list <- FeaturePlot(tumor_obj_list[[j]],features=panel_score_cols,
                          combine=FALSE,pt.size=0.05,reduction="umap",order=TRUE)
    
    for(i in seq_along(p_list)){
        is_sep <- i <= length(m_score_cols)
        p_list[[i]] <- p_list[[i]] +
            scale_color_gradientn(colors=sep_colors,
                                  limits=if(is_sep) raw_score_range else sr_score_range,
                                  oob=scales::squish,
                                  name=if(is_sep) "SEP score (UCell)" else "Mean z-score") +
            coord_fixed(xlim=panel_xlim,ylim=panel_ylim) +
            labs(title=if(j==panel_groups[1]) panel_score_names[i] else NULL,
                 x=NULL,y=if(i==1) j else NULL) + panel_theme
        if(i==1){
            p_list[[i]] <- p_list[[i]] +
                theme(axis.title.y=element_text(size=10,angle=90,hjust=0.5,
                                                margin=margin(r=5)))
        }
        p_MSEP_list[[paste(j,panel_score_cols[i],sep="_")]] <- p_list[[i]]
    }
}

p_MSEP <- wrap_plots(p_MSEP_list,ncol=length(panel_score_cols),byrow=TRUE,guides="collect") &
    theme(legend.position="right")

ggsave("Figure/scRNA/Figure4F_scRNA_Tumor_All_R_NR_M-SEP_score.pdf", p_MSEP,width=12,height=7)

p_response_list <- list()

for(j in panel_groups){
    
    p_list <- FeaturePlot(tumor_obj_list[[j]], features="M_response_score", reduction="umap", combine=FALSE,pt.size=0.05, order=TRUE)
    p_list[[1]] <- p_list[[1]] + 
        scale_color_gradient2(low=response_colors["low"], mid=response_colors["mid"], high=response_colors["high"],
                              midpoint=0, limits=sr_score_range, oob=scales::squish, name="Resistance - sensitivity") +
        coord_fixed(xlim=panel_xlim,ylim=panel_ylim) +
        labs(title=if(j==panel_groups[1]) "M-response score" else NULL,
             x=NULL,y=j) + panel_theme +
        theme(axis.title.y=element_text(size=10,angle=90,hjust=0.5,margin=margin(r=5)))
    
    p_response_list[[j]] <- p_list[[1]]
    
}

p_response <- wrap_plots(p_response_list, ncol=1, guides="collect") &
    theme(legend.position="right")
ggsave("Figure/scRNA/Figure4F_scRNA_Tumor_All_R_NR_response_score.pdf", p_response, width=4.7, height=5.5)

library(tibble)
library(ggpubr)

obj_tumor <- subset(obj_use, Treatment == 'pre'&celltype_s1=="Tumor")

sep_score <- obj_tumor@meta.data %>%
    rownames_to_column("Cell") %>%
    transmute(
        Sample=as.character(Sample),
        Response=factor(response,levels=c("NR","R")),
        M_SEP1=M_SEP1_UCell,
        M_SEP2=M_SEP2_UCell,
        M_SEP3=M_SEP3_UCell,
        M_SEP5=M_SEP5_UCell
    ) %>%
    group_by(Sample,Response) %>%
    filter(n()>=10) %>%
    summarise(across(starts_with("M_SEP"), ~mean(.x,na.rm=TRUE)),.groups="drop") %>%
    mutate(
        across(c(M_SEP1,M_SEP2,M_SEP3,M_SEP5), ~as.numeric(scale(.x))),
        M_sensitive_score=(M_SEP3+M_SEP5)/2,
        M_resistance_score=(M_SEP1+M_SEP2)/2,
        M_response_score=M_resistance_score-M_sensitive_score
    )

p3 <- ggplot(sep_score,aes(Response,M_response_score)) +
    geom_boxplot(width=0.6,outlier.shape=NA,alpha=0.7) +
    geom_jitter(aes(fill=Response),width=0.12,size=3,shape=21,alpha=0.8) +
    stat_compare_means(
        comparisons=list(c("NR","R")),
        method="t.test",
        method.args=list(alternative="greater"),
        label="p.format"
    ) +
    scale_fill_manual(values=c(R="#FFB900", NR="#4DBBD5")) +
    scale_y_continuous(expand=expansion(mult=c(0.02,0.15))) +
    labs(x=NULL,y="resistance-sensitive") +
    theme_classic() +
    theme(legend.position="none")
ggsave("Figure/scRNA/Figure4G_scRNA_response_Predicted_Score_Pre.pdf",p3,width=2,height=2.7)

library(tidyr)
expan <- read.csv("Data/scRNA/Tcell/expansion_info.csv")
expan$Sample <- as.character(expan$Sample)
expan$Treatment = substr(expan$sample,5,10)
expan = subset(expan,Treatment == 'Pre')

tumor_mps <- c(paste0("M_SEP", 1:5))
tumor_mp_cols <- paste0(tumor_mps, "_UCell")

expansion_info <- expan %>%
    transmute(Sample = as.character(Sample), Expansion_score = as.numeric(num)) %>%
    distinct(Sample, .keep_all = TRUE) %>%
    filter(is.finite(Expansion_score))


tumor_mp_patient <- meta_expansion %>%
    filter(CellType == "Tumor") %>%
    select(Sample, all_of(tumor_mp_cols)) %>%
    pivot_longer( cols = all_of(tumor_mp_cols), names_to = "MP", values_to = "UCell") %>%
    mutate(MP = sub("_UCell$", "", MP)) %>%
    group_by(Sample, MP) %>%
    summarise(n_Tumor = sum(is.finite(UCell)),
              MP_mean = mean(UCell, na.rm = TRUE), 
              MP_median = median(UCell, na.rm = TRUE), .groups = "drop") %>%
    filter(n_Tumor >= 10, is.finite(MP_mean)) %>%
    left_join(expansion_info, by = "Sample") %>%
    filter(is.finite(Expansion_score)) %>%
    mutate(MP = factor(MP, levels = tumor_mps))

safe_cor_test <- function(df) {
    x <- df$Expansion_score
    y <- df$MP_mean
    
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    
    if(length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
        return(tibble(n_patients = length(x), rho = NA_real_, p_value = NA_real_))
    }
    
    result <- cor.test(x, y, method = "spearman", exact = FALSE)
    tibble(n_patients = length(x), rho = unname(result$estimate), p_value = result$p.value)
}

tumor_cor_results <- tumor_mp_patient %>%group_by(MP) %>% group_modify(~safe_cor_test(.x)) %>% ungroup() %>%
    mutate(FDR = p.adjust(p_value, method = "BH"),
           Direction = case_when(rho > 0 ~ "Positive", rho < 0 ~ "Negative", TRUE ~ "Undetermined"),
           label = paste0("Spearman r = ", sprintf("%.3f", rho),"\n P = ",
                          format.pval(p_value,digits = 3,eps = 0.001))) %>% arrange(FDR)

plot_label <- tumor_cor_results %>% select(MP, label)

p_tumor_cor <- ggplot(tumor_mp_patient, aes(x = Expansion_score, y = MP_mean)) +
    geom_point(size = 2.8, shape = 21, fill = "#B4E1EB") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#2166AC", fill = "#92C5DE") +
    geom_text(data = plot_label, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 3.2) +
    facet_wrap(~MP, scales = "free_y", ncol = 6) +
    labs(x = "TCR expansion ratio", y = "Patient-level mean UCell score") +
    theme_few() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(filename = 'Figure/scRNA/Figure4I_TCR_clone_expansion_ratio_cor_MSEP.pdf',p_tumor_cor,width = 13, height = 3 )
 
sep_celltype <- list(
    N_SEP1=c("Fibor_apCAF","Fibro_iCAF","Fibro_mCAF","Fibro_MKI","Fibro_vCAF","Fibro_eCAF"),
    N_SEP2=c("B_cell", "CD4_Tcell","CD8_Tcell","DC","Macrophages","Monotype"),
    N_SEP3=c("Hepatocytes"),
    N_SEP4=c("Endo_C0_COL4A1","Endo_C1_ACKR1","Endo_C2_FCN3","Endo_C3_FBLN5","Endo_C4_CCL21","Fibro_vCAF"),
    N_SEP5=c("Plasma","Fibor_apCAF","Fibro_iCAF","Fibro_mCAF"),
    N_SEP6=c("Endo_C1_ACKR1","Endo_C2_FCN3","Endo_C3_FBLN5","Endo_C4_CCL21","Fibro_eCAF", "Fibro_MKI", "Fibro_vCAF"),
    N_SEP7=c("Hepatocytes"),
    
    N_SEP8=c("Endo_C0_COL4A1","Endo_C1_ACKR1","Endo_C2_FCN3","Endo_C3_FBLN5","Endo_C4_CCL21"),
    
    N_SEP9=c("B_cell", "CD4_Tcell","CD8_Tcell","gamma_delta_T","NK"),
    
    N_SEP10=c("Endo_C1_ACKR1","Endo_C2_FCN3","Endo_C3_FBLN5", "Fibor_apCAF"),
    
    N_SEP11=c("Fibor_apCAF","Fibro_iCAF","Fibro_mCAF","Fibro_MKI","Fibro_vCAF","Fibro_eCAF","Endo_C3_FBLN5"),
    N_SEP12=c("Fibor_apCAF","Fibro_iCAF","Fibro_mCAF","Fibro_MKI","Fibro_vCAF","Fibro_eCAF"),
    
    M_SEP1="Tumor",
    M_SEP2="Tumor",
    M_SEP3="Tumor",
    M_SEP4="Tumor",
    M_SEP5="Tumor"
)

calc_sep_patient <- function(sep,cells){
    meta_all_expansion %>%
        filter(celltype_s1 %in% cells) %>%
        select(Sample,all_of(paste0(sep,"_UCell"))) %>%
        rename(score=paste0(sep,"_UCell")) %>%
        group_by(Sample) %>%
        summarise(SEP_score= mean(score,na.rm=TRUE),n_cells=n(),.groups="drop") %>%
        filter(n_cells>=10) %>%
        mutate(SEP=sep)
}

sep_patient <- bind_rows(lapply(names(sep_celltype),function(x){calc_sep_patient(x,sep_celltype[[x]])})) %>%
    left_join(expansion_info,by="Sample") %>%
    filter(is.finite(Expansion_score))

safe_cor <- function(df){
    df <- df %>% filter(is.finite(SEP_score), is.finite(Expansion_score))
    if(nrow(df)<3 || length(unique(df$SEP_score))<3 || length(unique(df$Expansion_score))<3){
        return(tibble(rho=NA, p=NA))
    }
    test <- cor.test(df$SEP_score, df$Expansion_score, method="spearman")
    tibble(rho=unname(test$estimate),p=test$p.value)
}

plot_sep_expansion <- function(df){
    cor_res <- safe_cor(df)
    label <- paste0("Spearman rho = ", sprintf("%.2f",cor_res$rho),"\nP = ",format.pval(cor_res$p,digits=2,eps=0.001))
    
    ggplot(df,aes(Expansion_score, y=SEP_score))+
        geom_point(size=2,shape=21,fill="#B4E1EB")+
        geom_smooth(method="lm",formula=y~x,se=TRUE,color="#2166AC",fill="#92C5DE")+
        annotate("text",x=-Inf,y=Inf,label=label,hjust=-0.05,vjust=1.2,size=3)+
        theme_few()+
        theme(plot.title=element_text(hjust=0.5,size = 10))+
        labs(x="TCR clone expansion ratio",y="SEP score")
}

expansion_cor <- sep_patient %>% group_by(SEP) %>% group_modify(~safe_cor(.x)) %>%
    ungroup() %>% mutate(FDR=p.adjust(p,"BH"))

p_list <- lapply(c(paste0('N_SEP',1:12)),function(x){
    plot_sep_expansion(subset(sep_patient,SEP==x))+
        ggtitle(paste0(sep_short_names[x],"\n",x))
})
p_sep_celltype_expansion = wrap_plots(p_list,ncol=6)  
ggsave('Figure/scRNA/Figure4JK_scRNA_TCR_ratio_vs_N_SEP_resistance_celltype.pdf',p_sep_celltype_expansion,width = 16,height = 6)


sep_score2 = merge(sep_score,expan,by = "Sample" )
df_plot <- sep_score2  %>%
    filter(Expan %in% c("E", "NE"), !is.na(M_response_score)) %>%
    mutate(Expan = factor(Expan, levels = c("NE", "E")))



p <- ggplot(df_plot, aes(x = Expan, y = M_response_score, fill = Expan)) +
    geom_boxplot(width = 0.65,  fill = 'white', outlier.shape = NA) +
    geom_jitter(width = 0.15,shape = 21, size = 4) +
    stat_compare_means(
        comparisons = list(c("E", "NE")),
        method = "t.test",
        method.args = list(alternative = "less", exact = TRUE),
        label = "p.format"
    ) +
    scale_fill_manual(values = c("E" = "#bcbddc", "NE" = "#f3e9bb")) +
    labs(x = NULL,y = "Responsed predicted Score", title = "") +
    ylim(-3,2)+
    theme_few() +
    theme(legend.position = "none")
ggsave(filename = 'Figure/scRNA/Figure4H_boxplor.pdf',width = 3,height = 4)

