rm(list = ls())
gc()

#' ----------------------------------
#' inferCNV-patient
#' Author zhangyanru 
#' Time: 26.01.06
#' 

library(dplyr)
library(infercnv)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(dendextend)


data_path <- "/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC_spatial_V2/"  
setwd(data_path)

## patient Info
info_df = read.csv("info_file/Patient_Info.csv")
meta <- fread('Data/SpatialData/1_Integration/hex_36um/merge_sample_files/metadata.csv')
meta <- as.data.frame(meta)
rownames(meta) = meta$barcode


cell_type_obs = subset(meta, Type == 'Tumor'&CellType =='Epithelial')
cell_type_obs$Hex_ID = rownames(cell_type_obs)
cell_type_obs = cell_type_obs[,c('Hex_ID','sample_id','CellType')]
  
cell_type_ref = subset(meta, Type == 'Normal'&CellType =='Hepatocyte'&sample_id !='SR250811M4')
cell_type_ref$Hex_ID = rownames(cell_type_ref)
cell_type_ref = cell_type_ref[,c('Hex_ID','sample_id','CellType')]
 
# merge
cell_type_anno = rbind(cell_type_ref,cell_type_obs)
cell_info = merge(cell_type_anno,info_df,by.x = 'sample_id',by.y = 'SampleName')
patients = unique(cell_info$PatientID)
patients

tumor = subset(cell_info, Type == 'Tumor')
normal_cells = subset(cell_info,Type != 'Tumor')

#---------
count.data <- Matrix::readMM('Data/SpatialData/1_Integration/hex_36um/merge_sample_files/counts.mtx.gz')
count.data = Matrix::t(count.data)
genes <- fread('Data/SpatialData/1_Integration/hex_36um/merge_sample_files/genes.tsv', header = FALSE)$V1
barcodes <- fread('Data/SpatialData/1_Integration/hex_36um/merge_sample_files/barcodes.tsv', header = FALSE)$V1
rownames(count.data) <- genes
colnames(count.data) <- barcodes
 
 
gene.loc <- read.table('/Volumes/B/Projects/2025Project/SpatialAnalysis/Code/inferCNV/ref/gene.loc.hg38.txt',row.names = 1)

for (pati in patients) {
    pat_tumor = subset(tumor,PatientID == pati)
    sub_normal <- normal_cells %>% slice_sample(n = 2000, replace = FALSE)
    sub_normal = as.data.frame(sub_normal)
    sub_normal$CellType = 'reference'

    cells = rbind(pat_tumor,sub_normal)
    cellL = as.character(names(table(cells$SampleID)))[table(cells$SampleID)>100]
    cells = cells[cells$SampleID %in%cellL,]
    rownames(cells) = cells$Hex_ID
    head(cells)
    table(cells$CellType)

    cells = cells[,c('CellType','Hex_ID')]

    infer.count = count.data[,cells$Hex_ID]
    infer.count = as.matrix(infer.count)

    save(infer.count,gene.loc,cells,file = paste0('Data/SpatialData/1_Integration/inferCNV/run_inferCNV/',pati,'_inferCNV.Rdata'))

    load(paste0('Data/SpatialData/1_Integration/2_inferCNV/run_inferCNV/',pati,'_inferCNV.Rdata'))

    infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = infer.count,
                                         gene_order_file = gene.loc,
                                         annotations_file = cells,
                                         ref_group_names = c('reference'),
                                         chr_exclude = c('chrY', 'chrM', 'chrMT'))
    
    suppressWarnings(infercnv::run(infercnv_obj, 
                                   cutoff=0.1,
                                   out_dir=paste0('Data/ST/inferCNV/patient_inferCNV/',pati,'_hex_36um'), 
                                   cluster_by_groups = F, 
                                   analysis_mode = 'subclusters', 
                                   HMM=T,   
                                   write_expr_matrix=T,
                                   num_threads=24, 
                                   tumor_subcluster_pval = 0.05,
                                   no_plot = F,
                                   k_nn = 50,
                                   leiden_resolution = 0.005,
                                   output_format = 'png'))
    
}


datapath <- file.path(root_dir, "Data/ST/inferCNV/patient_inferCNV")
suffix_dir <- "_hex_36um"
anno_list <- lapply(patients, function(pid) {
    
    anno_path <- file.path(datapath, paste0(pid, suffix_dir), "infercnv.observation_groupings.txt")
    infer_anno <- read.table(anno_path,header = TRUE,sep = " ",stringsAsFactors = FALSE)
    
    infer_anno$cell_id <- rownames(infer_anno)
    group_stat <- infer_anno %>% dplyr::count(Dendrogram.Group, name = "group_total")
    keep_groups <- group_stat %>% filter(group_total >= 80) %>% pull(Dendrogram.Group)
    
    infer_anno_keep <- infer_anno %>%
        inner_join(group_stat, by = "Dendrogram.Group") %>%
        filter(Dendrogram.Group %in% keep_groups) %>%
        mutate(PatientID = pid) %>%
        relocate(PatientID, cell_id, Dendrogram.Group, group_total)
    
    group_counts <- table(infer_anno_keep$Dendrogram.Group)  
    sorted_groups <- names(sort(group_counts, decreasing = TRUE)) 
 
    new_labels <- LETTERS[seq_along(sorted_groups)]
    name_map <- setNames(new_labels, sorted_groups)
    infer_anno_keep$CNV_subcluster <- name_map[infer_anno_keep$Dendrogram.Group]
 
    return(infer_anno_keep)
})

anno_all <- bind_rows(anno_list)
anno_all <- anno_all[,c('cell_id','CNV_subcluster')]
colnames(anno_all) <- c('Hex_ID','CNV_subcluster')
write.table(anno_all,file = 'Data/ST/inferCNV_result/subcluster_anno.tsv',sep = "\t",quote = FALSE,row.names = FALSE)


