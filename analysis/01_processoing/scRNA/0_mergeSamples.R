#'
#'Time:36.03.05
#' Zhangyanru
rm(list = ls())
gc()

library(Seurat)
library(dplyr)
root_dir <- "./"
setwd(root_dir)

sample_dirs <- list.dirs('Data/scRNA_Data/0.rawData/0_Samples/', full.names = TRUE, recursive = FALSE)
sample_names <- basename(sample_dirs)
sample_names

obj_list <- vector("list", length(sample_dirs))
names(obj_list) <- sample_names

for (i in seq_along(sample_dirs)) {
    sdir <- sample_dirs[i]
    sname <- sample_names[i]
    h5 <- file.path(sdir, "filtered_feature_bc_matrix.h5")
    x <- Read10X_h5(h5)
    obj <- CreateSeuratObject(counts = x,project = sname,min.cells = 0, min.features = 3)
    obj$sample <- sname
    
    obj_list[[sname]] <- obj
    message("✅ loaded: ", sname, "  cells=", ncol(obj), "  genes=", nrow(obj))
}

Merge_seurat <- merge(x = obj_list[[1]],y = obj_list[-1],add.cell.ids = names(obj_list),project = "ICC_scRNA")
Merge_seurat
 
patient_Info <- read.csv('./info_file/scPatient_info.csv',row.names = 1)
patient_Info <- patient_Info[Merge_seurat$sample,]
Merge_seurat@meta.data[,c('Patient','Treatment', 'Respond','SampleID')] <- patient_Info[,c('Patient','Treatment', 'Respond','SampleID')]
  

Merge_seurat[["RNA"]] <- JoinLayers(Merge_seurat[["RNA"]])
Merge_seurat$SampleName <- paste0(Merge_seurat$Patient,'_',Merge_seurat$Treatment) 
Merge_seurat[["percent.mt"]] <- PercentageFeatureSet(Merge_seurat, pattern = "^MT-")
 
saveRDS(Merge_seurat, file = "Data/scRNA_Data/0.rawData/1_Integration/ICC_scRNA_merged_seurat_original.rds")


Merge_seurat <- subset(Merge_seurat,nCount_RNA <= 15000 & nCount_RNA >= 500 & nFeature_RNA <= 7500 & nFeature_RNA >= 300 & percent.mt <= 50)

Merge_seurat = readRDS("Data/scRNA_Data/0.rawData/1_Integration/ICC_scRNA_merged_seurat_Clean.rds") 
# LogNormalize
Merge_seurat <- NormalizeData(Merge_seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
Merge_seurat <- FindVariableFeatures(Merge_seurat, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
Merge_seurat <- ScaleData(Merge_seurat, features = VariableFeatures(Merge_seurat), vars.to.regress = "percent.mt", verbose = FALSE)
Merge_seurat <- RunPCA(Merge_seurat, features = VariableFeatures(Merge_seurat), npcs = 50, seed.use = 42)
Merge_seurat <- FindNeighbors(Merge_seurat, reduction = "pca", dims = 1:20, k.param = 20)

Merge_seurat <- FindClusters(Merge_seurat,resolution = 0.5, random.seed = 42)
Merge_seurat <- RunUMAP(Merge_seurat, reduction = "pca", dims = 1:20, n.neighbors = 30, seed.use = 42)
 
DimPlot(Merge_seurat, reduction = "umap", group.by = c('seurat_clusters','SampleName'),label = TRUE, repel = TRUE, raster = TRUE)
 