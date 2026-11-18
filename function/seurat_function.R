
seurat_process1<-function(obj,harmony=T,vars="orig.ident",npcs=10,res=0.3,max.iter.harmony = 50){
    library(harmony)
    
    all.genes <- rownames(obj)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst",nfeatures = 2000, verbose = FALSE)
    
    VariableGenes = VariableFeatures(object = obj)
    IGgenes = VariableGenes[grep("^IG",VariableGenes)]
    TRgenes = VariableGenes[grep("^TR",VariableGenes)]
    ACgenes = VariableGenes[grep("^AC",VariableGenes)]
    ALgenes = VariableGenes[grep("^AL",VariableGenes)]
    HISTgenes = VariableGenes[grep("^HIST",VariableGenes)]
    LINgenes = VariableGenes[grep("^LIN",VariableGenes)]
    MTgenes = VariableGenes[grep("^MT",VariableGenes)]
    RPSgenes = VariableGenes[grep("^RPS",VariableGenes)]
    
    featuregene = VariableGenes[!VariableGenes %in% c(RPSgenes,TRgenes, IGgenes,ACgenes, ALgenes,HISTgenes,LINgenes,MTgenes)]
    
    obj@assays$RNA@var.features = featuregene
    
    obj <- ScaleData(obj, verbose = FALSE,features = all.genes)
    obj <- RunPCA(obj,features = featuregene,verbose = FALSE)
    
    if(harmony==T){
        #before dim=20
        obj <- RunHarmony(obj, group.by.vars = vars, max_iter = max.iter.harmony)
        obj <- RunUMAP(obj, reduction = "harmony", dims = 1:npcs)
        obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:npcs) %>% FindClusters(resolution = res)
    }else{
        obj <- RunUMAP(obj, dims = 1:npcs)
        obj <- FindNeighbors(obj, dims = 1:npcs) %>% FindClusters(resolution = res)
    }
    
    return(obj)
}

