
library(Seurat)
library(Matrix)
library(NMF)
library(data.table)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(viridis)
library(scales)
library(patchwork)
library(ggsci)
library(grid)
library(ggplot2)
library(reshape2)
library(viridis)
library(scales)
library(patchwork)
setwd("./")
root_dir <- "Data/ST/SEP/Module1/"
nmf_input_dir <- file.path(root_dir, "01_nmf_inputs")
nmf_res_dir <- file.path(root_dir, "02_nmf_results")
mp_dir <- file.path(root_dir, "03_metaprograms")
figure_dir = file.path('Figure/ST/SEP/SEP_Heatmap')

obj_hex <- readRDS('Data/ST/Spatial_Niche/All_hex_in_seurat_clean_final.RDS')
sample_names <- as.character(unique(obj_hex$SampleID))

is_bad_gene <- function(gene){
    technical <- grepl("^MT-|^RPL|^RPS|^LINC|^LOC|^MIR", gene)
    housekeeping <- gene %in% c("MALAT1", "ACTB", "ACTG1", "GAPDH", "EEF1A1")
    vdj <- grepl("^(IGHV|IGKV|IGLV|IGHJ|IGKJ|IGLJ|TRAV|TRBV|TRDV|TRGV|TRAJ|TRBJ|TRDJ|TRGJ)", gene)
    technical | housekeeping | vdj
}
pad_program <- function(x, n = 50){
    x <- unique(x[!is.na(x) & x != ""])
    if(length(x) > n) x <- x[seq_len(n)]
    length(x) <- n
    x
}
overlap_n <- function(x, y){
    x <- unique(x[!is.na(x) & x != ""])
    y <- unique(y[!is.na(y) & y != ""])
    length(intersect(x, y))
}

leiden_clustering_list <- setNames(
    lapply(sample_names, function(s_id){
        sample_obj <- subset(obj_hex, SampleID == s_id)
        sample_obj <- NormalizeData(sample_obj, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
        Idents(sample_obj) <- "Niche"
        markers <- FindAllMarkers(sample_obj, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.25, random.seed = 42, return.thresh = 0.05,verbose = F)
        markers$Sample <- s_id
        markers
    }),
    sample_names
)
saveRDS(leiden_clustering_list, file.path(root_dir, "all_leiden_programs.rds"))


leiden_programs_by_sample <- setNames(
    lapply(sample_names, function(sample_id){

        sample_table <- leiden_clustering_list[[sample_id]]
        
        cluster_ids <- unique(as.character(sample_table$cluster))
        cluster_ids <- cluster_ids[order(suppressWarnings(as.numeric(cluster_ids)), cluster_ids)]
        
        sample_programs <- setNames(
            lapply(cluster_ids, function(cluster_id){
                cluster_table <- sample_table[as.character(sample_table$cluster) == cluster_id & sample_table$p_val_adj <= 0.05 & !is_bad_gene(sample_table$gene),, drop = FALSE]
                genes <- cluster_table$gene[order(cluster_table[["avg_log2FC"]], decreasing = TRUE)]
                head(unique(genes), 50)
            }),
            paste0("Leiden.", cluster_ids)
        ) 
        sample_programs[lengths(sample_programs) >= 10]
    }),
    sample_names
)

all_leiden_programs <- unlist(leiden_programs_by_sample, recursive = FALSE, use.names = TRUE)
all_leiden_programs <- lapply(all_leiden_programs, pad_program, n = 50)
leiden_programs_mat <- do.call(cbind, all_leiden_programs)
rownames(leiden_programs_mat) <- paste0("Gene_", seq_len(nrow(leiden_programs_mat)))

saveRDS(leiden_programs_mat,file =file.path(root_dir,'leiden_programs_mat.RDS') )
  

obj_pseudo <- readRDS('./Project_ICC/Data/ST/SEP/Module1/00_pseudo/obj_pseudo.RDS')

for(s_id in sample_names){
    
    cat("\nPreparing NMF input:", s_id, "\n")
    
    sample_obj <- subset(obj_pseudo, SampleID == as.character(s_id))
    
    m <- LayerData(sample_obj, assay = "RNA", layer = "counts")
    m <- m[!is_bad_gene(rownames(m)), , drop = FALSE]
    m <- m[, colSums(m) > 0, drop = FALSE]
    
    min_detect_n <- max(10, ceiling(0.01 * ncol(m)))
    detected_n <- Matrix::rowSums(m > 0)
    keep_detected <- detected_n >= min_detect_n
    m <- m[keep_detected, ,drop = FALSE]
    
    cat("Genes after detection filtering:", nrow(m), "\n")
    
    m <- as.matrix(m)
    scaling_factor <- 1000000 / colSums(m)
    m_cpm <- sweep(m, 2, scaling_factor, "*")
    m_log <- log2(1 + m_cpm / 10)
    
    gene_mean <- rowMeans(m_log)
    gene_var <- matrixStats::rowVars(m_log)
    keep_basic <- (gene_mean > 0.4 & gene_var > 0)
    m_proc <- m_log[keep_basic, ,drop = FALSE]
    
    gene_var_use <- gene_var[keep_basic]
    
    if(nrow(m_proc) > 4000){
        
        top_gene_index <- order(gene_var_use, decreasing = TRUE)[seq_len(4000)]
        m_proc <- m_proc[top_gene_index, ,drop = FALSE]
    }
    m_nmf <- m_proc - rowMeans(m_proc)
    m_nmf[m_nmf < 0] <- 0
    m_nmf <- m_nmf[rowSums(m_nmf) > 0, , drop = FALSE]
    
    saveRDS(m_nmf, file.path(nmf_input_dir, paste0(s_id, "_nmf_input.rds")))
    
    rm(sample_obj, m, m_cpm, m_log, m_proc, m_nmf)
    gc()
}

 
set.seed(42)
for(s_id in sample_names){

    cat("\nRunning test NMF:", s_id, "\n")

    m <- readRDS(file.path(nmf_input_dir, paste0(s_id, "_nmf_input.rds")))
    m <- m[Matrix::rowSums(is.na(m)) == 0, , drop = FALSE]
    m <- m[Matrix::rowSums(m) > 0, , drop = FALSE]
    m <- m[, Matrix::colSums(m) > 0, drop = FALSE]

    cat("Genes:", nrow(m), "\n")
    cat("Hexagons:", ncol(m), "\n")

    run_time <- system.time({
        nmf_res <- NMF::nmf(x = m, rank = 2:11, nrun = 10, method = "snmf/r",
                            .opt = list(debug = FALSE, parallel = FALSE, shared.memory = FALSE, verbose = TRUE))
    })

    saveRDS(nmf_res, file.path(nmf_res_dir, paste0(s_id, "_nmf_res.rds")))

    cat("Elapsed time:", run_time["elapsed"], "seconds\n")
    rm(m, nmf_res)
    gc()
}

prog_genes_ls <- list()
for(s_id in sample_names){
    
    cat("\nExtracting NMF programs:", s_id, "\n")
    nmf_obj <- readRDS(file.path(nmf_res_dir, paste0(s_id, "_nmf_res.rds")))
    nmf_mats <- c()
    
    for(j in names(nmf_obj$fit)){
        basis_mat <- NMF::basis(nmf_obj$fit[[j]])
        rank_value <- ncol(basis_mat)
        colnames(basis_mat) <- paste0(s_id, ".NMF.r", rank_value, ".f", seq_len(rank_value)) 
        nmf_mats <- cbind(nmf_mats, basis_mat)
    }
    prog_genes_ls[[s_id]] <- nmf_mats
}
saveRDS(prog_genes_ls,file =file.path(root_dir,'prog_genes_ls_nmf.RDS'))

prog_genes_ls = readRDS(file =file.path(root_dir,'prog_genes_ls_nmf.RDS'))
Genes_nmf_w_basis <- prog_genes_ls
nmf_programs_original <- lapply(prog_genes_ls, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
all_nmf_mat <- do.call(cbind, nmf_programs_original)
dim(all_nmf_mat)
 



Robust_NMF_Programs <- function(nmf_programs, intra_min = 30, intra_max = 10, inter_min = 10){
    
    intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x, y)))))
    intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = TRUE)[2]))
    
    nmf_sel <- lapply(names(nmf_programs), function(s) nmf_programs[[s]][, intra_intersect_max[[s]] >= intra_min, drop = FALSE])
    names(nmf_sel) <- names(nmf_programs)
    nmf_sel <- nmf_sel[sapply(nmf_sel, ncol) > 0]
    nmf_all <- do.call(cbind, nmf_sel)
    inter_intersect <- apply(nmf_all, 2, function(x) apply(nmf_all, 2, function(y) length(intersect(x, y))))
    final_filter <- character()
    
    for(s_id in names(nmf_sel)){
        current_ids <- colnames(nmf_sel[[s_id]])
        other_ids <- setdiff(colnames(nmf_all), current_ids)
        if(length(current_ids) == 0 || length(other_ids) == 0) next
        
        cross_overlap <- inter_intersect[other_ids, current_ids, drop = FALSE]
        best_overlap <- sort(apply(cross_overlap, 2, max), decreasing = TRUE)
        best_overlap <- best_overlap[best_overlap >= inter_min]
        
        if(length(best_overlap) == 0) next
        
        selected <- names(best_overlap)[1]
        
        if(length(best_overlap) > 1){
            for(i in 2:length(best_overlap)){
                candidate <- names(best_overlap)[i]
                if(max(inter_intersect[selected, candidate, drop = FALSE]) <= intra_max) selected <- c(selected, candidate)
            }
        }
        
        final_filter <- c(final_filter, selected)
    }
    
    unique(final_filter)
}

nmf_robust_id <- Robust_NMF_Programs(nmf_programs_original, 
                                     intra_min = 35,
                                     intra_max = 10, 
                                     inter_min = 10)

robust_nmf_mat <- all_nmf_mat[, nmf_robust_id, drop = FALSE]
robust_nmf_mat[1:4,1:4]
dim(robust_nmf_mat) 
saveRDS(robust_nmf_mat,file = file.path(root_dir,'robust_nmf_programs_mat.RDS') )


robust_nmf_mat <- readRDS(file = file.path(root_dir,'robust_nmf_programs_mat.RDS'))
leiden_programs_mat <- readRDS(file = file.path(root_dir,'leiden_programs_mat.RDS'))

all_programs <- cbind(robust_nmf_mat, leiden_programs_mat)
all_programs <- all_programs[,colSums(is.na(all_programs)) == 0, drop = FALSE]

# calculate similarity between programs
nmf_intersect <- apply(all_programs, 2, function(x) apply(all_programs, 2,  function(y)  length(intersect(na.omit(x), na.omit(y)))))
 
# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_all <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc_all <- reorder(as.dendrogram(nmf_intersect_hc_all), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc_all), order.dendrogram(nmf_intersect_hc_all)]


save(nmf_intersect, all_programs,file = file.path(root_dir,'nmf_intersect_and_all_programs.rdata'))

load( file.path(root_dir,'nmf_intersect_and_all_programs.rdata'))

nmf_intersect_KEEP <- nmf_intersect
work_programs = all_programs

Min_intersect_initial <- 10
Min_intersect_cluster <- 10
Min_group_size <- 4
Min_program_per_MP <- 5

Sorted_intersection <- sort(apply(nmf_intersect_KEEP , 2, function(x) (length(which( x >= Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list <- list()   ### Every entry contains the NMFs of a chosec cluster
k <- 1
Curr_cluster <- c()
MP_list <- list()

while (Sorted_intersection[1] > Min_group_size) {   
    
    Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
    
    Genes_MP <- work_programs[,names(Sorted_intersection[1])] 
    work_programs  <- work_programs[,-match(names(Sorted_intersection[1]) , colnames(work_programs))]  
    Intersection_with_Genes_MP <- sort(apply(work_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)  
    NMF_history <- Genes_MP  
    
    while (Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {
        
        new_program_name <- names(Intersection_with_Genes_MP)[1]
        Curr_cluster <- c(Curr_cluster, new_program_name)
        new_program <- work_programs[, new_program_name]
        
        all_history_genes <- c(NMF_history, new_program)
        all_history_genes <- all_history_genes[!is.na(all_history_genes) & all_history_genes != ""]
        
        gene_frequency <- sort(table(all_history_genes), decreasing = TRUE)
  
        cutoff_frequency <- unname(gene_frequency[50])
        genes_above_cutoff <- names(gene_frequency[gene_frequency > cutoff_frequency])
        genes_at_border <- names(gene_frequency[gene_frequency == cutoff_frequency])
        n_border_needed <- 50 - length(genes_above_cutoff)
        
        border_score <- setNames(rep(-Inf, length(genes_at_border)), genes_at_border)
        
        for (program_id in Curr_cluster) {
            
            curr_study <- strsplit(program_id, ".", fixed = TRUE)[[1]][1]
            
            if (!curr_study %in% names(Genes_nmf_w_basis)) next
            
            curr_basis <- Genes_nmf_w_basis[[curr_study]]
            
            if (!program_id %in% colnames(curr_basis)) next
            
            gene_index <- match(toupper(genes_at_border), toupper(rownames(curr_basis)))
            valid_gene <- !is.na(gene_index)
            
            if (!any(valid_gene)) next
            
            curr_score <- as.numeric(curr_basis[gene_index[valid_gene], program_id, drop = TRUE])
            names(curr_score) <- genes_at_border[valid_gene]
            curr_score <- curr_score[is.finite(curr_score)]
            
            if (length(curr_score) == 0) next
            
            for (gene in names(curr_score)) {
                border_score[gene] <- max(border_score[gene], curr_score[gene], na.rm = TRUE)
            }
        }
        
        ranked_border_genes <- names(sort(border_score, decreasing = TRUE))
        Genes_MP <- c(genes_above_cutoff, head(ranked_border_genes, n_border_needed))
        Genes_MP <- unique(Genes_MP[!is.na(Genes_MP) & Genes_MP != ""])
    
        NMF_history <- all_history_genes
        
        work_programs <- work_programs[, -match(new_program_name, colnames(work_programs)), drop = FALSE]
        
        if (ncol(work_programs) == 0) break
        
        Intersection_with_Genes_MP <- sort(apply(work_programs, 2, function(x) length(intersect(Genes_MP, x))), decreasing = TRUE)
    }
    
    Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
    MP_list[[paste0("MP_",k)]] <- Genes_MP
    k <- k+1
    
    nmf_intersect_KEEP <- nmf_intersect_KEEP[-match(Curr_cluster,rownames(nmf_intersect_KEEP) ) , -match(Curr_cluster,colnames(nmf_intersect_KEEP) ) ]  # remove current chosen cluster
    Sorted_intersection <- sort(apply(nmf_intersect_KEEP , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   ### Sort intersection of remaining NMFs not included in any of the previous clusters
    
    Curr_cluster <- c()
    print(dim(nmf_intersect_KEEP)[2])
}

MP_program_number <- lengths(Cluster_list)
keep_MP <- names(MP_program_number[MP_program_number >= Min_program_per_MP])
Cluster_list <- Cluster_list[keep_MP]
MP_list <- MP_list[sub("Cluster_", "MP_", keep_MP)]

score_one_sample <- function(counts_mat, MP_list, sample_id, celltype_vec, type_vec,min_mean_expression = 0.1, conserved_genes = 0.3, min_celltype_hex = 20L, seed = 42) {
 
   
    counts_mat <- as(counts_mat, "dgCMatrix")
    
    all_spots <- colnames(counts_mat)
    sep_names <- names(MP_list)
    keep_genes <- !grepl("^MT-|^RPL|^RPS", rownames(counts_mat))
    m <- counts_mat[keep_genes, , drop = FALSE]
    lib_size <- Matrix::colSums(m)
    good <- is.finite(lib_size) & lib_size > 0
    m <- m[, good, drop = FALSE]
    lib_size <- lib_size[good]
    
    # Same normalization as the original: log2(1 + CPM/10).
    m_log <- m %*% Matrix::Diagonal(x = 1e5 / lib_size)
    m_log <- as(m_log, "dgCMatrix")
    dimnames(m_log) <- dimnames(m)
    m_log@x <- log2(1 + m_log@x)

    if (ncol(m_log)) {
        gene_mean <- Matrix::rowMeans(m_log)
        m_sq <- m_log
        m_sq@x <- m_sq@x ^ 2
        gene_var <- pmax(0, Matrix::rowMeans(m_sq) - gene_mean ^ 2)
        rm(m_sq)
    } else {
        gene_mean <- setNames(rep(0, nrow(m_log)), rownames(m_log))
        gene_var <- gene_mean
    }
    
    
    variable_genes <- rownames(m_log)[gene_var > 0]
    legacy_genes <- rownames(m_log)[gene_var > 0 & gene_mean > min_mean_expression]

    ct <- as.character(celltype_vec[match(colnames(m_log), names(celltype_vec))])
    md_good <- data.table(SpotID = colnames(m_log), CellType = ct)
    celltype_qc <- md_good[, .(n_hex = .N), by = CellType]
    celltype_qc[, fraction_of_nonzero_library_hex := n_hex / max(1L, ncol(m_log))]
    
    celltype_qc[, eligible_for_support :=
        !is.na(CellType) & nzchar(trimws(CellType)) &
        !(tolower(trimws(CellType)) %in% c("unassigned", "unknown", "na")) &
        n_hex >= min_celltype_hex]
    
    support_types <- sort(celltype_qc[eligible_for_support == TRUE, CellType])
    celltype_genes <- setNames(vector("list", length(support_types)), support_types)
    celltype_means <- setNames(vector("list", length(support_types)), support_types)
    
    for (ct_name in support_types) {
        ct_mean <- Matrix::rowMeans(m_log[, which(ct == ct_name), drop = FALSE])
        celltype_means[[ct_name]] <- ct_mean
        celltype_genes[[ct_name]] <- rownames(m_log)[
            rownames(m_log) %in% variable_genes & ct_mean > min_mean_expression]
    }

    
    coverage_by_celltype <- rbindlist(lapply(sep_names, function(mp) {
        z <- copy(celltype_qc)
        z[, `:=`(MP = mp, n_total = length(MP_list[[mp]]), n_present = NA_integer_)]
        for (ct_name in support_types) {
            z[CellType == ct_name, n_present :=
                length(intersect(MP_list[[mp]], celltype_genes[[ct_name]]))]
        }
        z[, fraction_present := n_present / n_total]
        z
    }), use.names = TRUE, fill = TRUE)
    
    selected_signatures <- setNames(vector("list", length(sep_names)), sep_names)
    coverage_parts <- gene_parts <- vector("list", length(sep_names))

    for (j in seq_along(sep_names)) {
        mp <- sep_names[j]
        candidates <- copy(coverage_by_celltype[MP == mp & eligible_for_support == TRUE])
        # Deterministic ties: more hex, then CellType name. No SEP score is used.
        setorderv(candidates, c("n_present", "n_hex", "CellType"), c(-1L, -1L, 1L))
        
        support <- if (nrow(candidates)) candidates$CellType[1L] else NA_character_
        
        retained <- if (!is.na(support)) {
            MP_list[[mp]][MP_list[[mp]] %in% celltype_genes[[support]]]
        } else character()
        
        selected_signatures[[mp]] <- retained
        n_total <- length(MP_list[[mp]])
        n_present <- length(retained)
        n_legacy <- sum(MP_list[[mp]] %in% legacy_genes)
        
        coverage_parts[[j]] <- data.table(
            MP = mp, n_total = n_total, n_present = n_present,
            fraction_present = n_present / n_total,
            analysis_eligible = n_present > 0L & n_present / n_total >= conserved_genes,
            selected_CellType = support,
            selected_CellType_n_hex = if (nrow(candidates)) candidates$n_hex[1L] else 0L,
            selected_CellType_fraction = if (nrow(candidates))
                candidates$fraction_of_nonzero_library_hex[1L] else NA_real_,
            n_tied_best_CellTypes = if (nrow(candidates))
                sum(candidates$n_present == candidates$n_present[1L]) else 0L,
            legacy_n_present = n_legacy,
            legacy_fraction_present = n_legacy / n_total,
            legacy_analysis_eligible = n_legacy / n_total >= conserved_genes,
            restored_eligibility = n_present / n_total >= conserved_genes &
                n_legacy / n_total < conserved_genes,
            status = if (is.na(support)) "no_eligible_CellType" else if (!n_present)
                "no_retained_signature_genes" else if (n_present / n_total < conserved_genes)
                "low_coverage" else "eligible",
            selected_genes = paste(retained, collapse = ";")
        )
    
        
        gene_parts[[j]] <- data.table(
            MP = mp, gene = MP_list[[mp]], selected_CellType = support,
            in_retained_count_universe = MP_list[[mp]] %in% rownames(m_log),
            nonzero_sample_variance = MP_list[[mp]] %in% variable_genes,
            sample_mean = unname(gene_mean[MP_list[[mp]]]),
            selected_CellType_mean = if (!is.na(support))
                unname(celltype_means[[support]][MP_list[[mp]]]) else NA_real_,
            retained_for_score = MP_list[[mp]] %in% retained,
            retained_by_legacy_filter = MP_list[[mp]] %in% legacy_genes
        )
    }
    
    MP_coverage <- rbindlist(coverage_parts)
    selected_genes <- rbindlist(gene_parts)
    MP_list_use <- selected_signatures[lengths(selected_signatures) > 0L]

    score_universe <- rownames(m_log)[rownames(m_log) %in%
        union(legacy_genes, unlist(MP_list_use, use.names = FALSE))]
    
    expr_nbin_use <- min(30L, max(1L, floor(length(score_universe) / 120)))
    
    signatures_raw <- matrix(NA_real_, nrow = length(all_spots), ncol = length(sep_names),
                             dimnames = list(all_spots, sep_names))
    
    sample_seed <- as.integer((as.double(seed) +
        sum(as.double(utf8ToInt(sample_id)) * seq_along(utf8ToInt(sample_id)))) %%
        (.Machine$integer.max - 1L) + 1L)
    
    if (length(MP_list_use) && ncol(m_log)) {
        m_score <- as.matrix(m_log[score_universe, , drop = FALSE])
        stopifnot(all(unlist(MP_list_use, use.names = FALSE) %in% rownames(m_score)))
        set.seed(sample_seed)
        computed <- tryCatch(scalop::sigScores(
            m_score, MP_list_use,
            center.rows = TRUE, center = TRUE, expr.center = TRUE,
            conserved.genes = 1, expr.nbin = expr_nbin_use,
            expr.binsize = 100, replace = FALSE
        ), error = function(e) stop("SigScore failed in ", sample_id, ": ",
                                    conditionMessage(e), call. = FALSE))
        computed <- as.matrix(computed)
        stopifnot(setequal(rownames(computed), colnames(m_score)),
                  setequal(colnames(computed), names(MP_list_use)),
                  all(is.finite(computed)))
        signatures_raw[rownames(computed), colnames(computed)] <- computed
        rm(m_score, computed)
    }
    
    weight_vec <- setNames(pmin(MP_coverage$fraction_present / conserved_genes, 1),MP_coverage$MP)
    
    signatures_weighted <- as.data.frame(sweep(signatures_raw, 2L, weight_vec[colnames(signatures_raw)], "*"))
    signatures_raw <- as.data.frame(signatures_raw)
    
    sample_qc <- data.table(
        n_input_hex = length(all_spots), n_nonzero_library_hex = sum(good),
        n_zero_library_hex = sum(!good), n_support_CellTypes = length(support_types),
        n_legacy_genes = length(legacy_genes), n_score_universe_genes = length(score_universe),
        n_restored_signature_genes = length(setdiff(score_universe, legacy_genes)),
        n_scored_SEPs = length(MP_list_use),
        n_eligible_SEPs = sum(MP_coverage$analysis_eligible),
        expr_nbin = expr_nbin_use, seed = sample_seed
    )
    print(MP_coverage[, .(MP, selected_CellType, n_present, fraction_present,
                         analysis_eligible, legacy_fraction_present)])
    
    return(list(
        raw_scores = signatures_raw, weighted_scores = signatures_weighted,
        coverage = MP_coverage, coverage_by_celltype = coverage_by_celltype,
        celltype_qc = celltype_qc, selected_genes = selected_genes,
        selected_signatures = selected_signatures, sample_qc = sample_qc
    ))
}

obj_hex <- readRDS('Data/ST/Spatial_Niche/All_hex_in_seurat_original_niche.RDS')
MP_list <- readRDS('Data/ST/SEP/Module1/03_metaprograms/MP_List_original.RDS')

min_mean_expression <- 0.1
conserved_genes <- 0.3
min_celltype_hex <- 20L

meta <- obj_hex@meta.data
obj_hex$CellType_original <- as.character(meta$CellType)

sample_names <- unique(as.character(meta$SampleID))
sample_names

counts_all <- LayerData(obj_hex, assay = "RNA", layer = "counts")

score_result <- setNames(vector("list", length(sample_names)), sample_names)
score_result

for (sample_id in sample_names) {
    
    cat("Sample:", sample_id,"\n")
    spot_ids <- rownames(meta)[as.character(meta[["SampleID"]]) == sample_id]
    counts_sample <- counts_all[,spot_ids, drop = FALSE]
    
    celltype_sample <- as.character(meta[spot_ids,"CellType"])
    names(celltype_sample) <- spot_ids
    type_sample <- as.character(meta[spot_ids, "Type"])
    names(type_sample) <- spot_ids
    
    score_result[[sample_id]] <- score_one_sample(counts_mat = counts_sample,
                                                  MP_list = MP_list,
                                                  sample_id = sample_id,
                                                  celltype_vec = celltype_sample,
                                                  type_vec = type_sample,
                                                  min_mean_expression = min_mean_expression,
                                                  conserved_genes = conserved_genes,
                                                  min_celltype_hex = min_celltype_hex,
                                                  seed = 42)
    invisible(gc())
}



# ------------------------- result -----------------------------
coverage_all <- rbindlist(
    lapply(names(score_result), function(sample_id) {
        x <- as.data.table(score_result[[sample_id]]$coverage)
        x[, SampleID := sample_id]
        x
    }),
    fill = TRUE)

collect_audit <- function(field) rbindlist(lapply(names(score_result), function(sid) {
    x <- copy(as.data.table(score_result[[sid]][[field]]))
    x[, SampleID := rep(sid, .N)]
    x
}), use.names = TRUE, fill = TRUE)
 
raw_score_all <- rbindlist(lapply(score_result,function(x){
    score_df <- as.data.frame(x$raw_scores)
    score_df$SpotID <- rownames(score_df)
    score_df}), fill = TRUE)

raw_score_all <- as.data.frame(raw_score_all)
rownames(raw_score_all) <- raw_score_all$SpotID
raw_score_all$SpotID <- NULL

score_idx <- match(colnames(obj_hex), rownames(raw_score_all))

for (mp in colnames(raw_score_all)) {
    obj_hex[[paste0(mp,"_score")]] <- raw_score_all[score_idx, mp]
}

weighted_score_all <- rbindlist(lapply(score_result,function(x){
    score_df <- as.data.frame(x$weighted_scores)
    score_df$SpotID <- rownames(score_df)
    score_df
}), fill = TRUE)

weighted_score_all <- as.data.frame(weighted_score_all)
rownames(weighted_score_all) <- weighted_score_all$SpotID
weighted_score_all$SpotID <- NULL

weighted_idx <- match(colnames(obj_hex), rownames(weighted_score_all))

for (mp in colnames(weighted_score_all)) {
    obj_hex[[paste0(mp, "_weighted_score")]] <- weighted_score_all[weighted_idx, mp]
}

meta <- obj_hex@meta.data
idx_normal_epi <- meta$Type %in% "Normal" & meta$CellType %in% "Epithelial"
celltype_updated <- as.character(meta$CellType)
celltype_updated[idx_normal_epi] <- "unassigned"
obj_hex$CellType <- celltype_updated
 

# ------------
meta_export <- as.data.table(obj_hex@meta.data,keep.rownames = "SpotID")
head(meta_export)

raw_score_cols <- paste0(names(MP_list), "_score")
raw_score_cols <- intersect(raw_score_cols, colnames(meta_export))

weighted_score_cols <- paste0(names(MP_list), "_weighted_score")
weighted_score_cols <- intersect(weighted_score_cols, colnames(meta_export))

export_cols <- c(
    "SpotID", "sample_id", "SampleID", "PatientID", "Type", "Region", "CellType", "CellType_original",
    raw_score_cols,
    weighted_score_cols
)

export_cols <- intersect(export_cols, colnames(meta_export))
mp_export <- meta_export[, ..export_cols]
head(mp_export)

# ---------------------- Plot -----------------------------
meta_st <- as.data.table(obj_hex@meta.data, keep.rownames="SpotID")
meta_st <- meta_st[!is.na(CellType) & CellType!="" & CellType!="unassigned"]
score_cols <- grep("MP_[0-9]+_score$", colnames(meta_st), value=TRUE)
score_cols

heatmap_full <- meta_st[, lapply(.SD, function(x) {
    if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
}), by = CellType, .SDcols = score_cols]

 
for (sample_id in unique(meta_st$SampleID)) {
    
    low_cov_mps <- coverage_all[SampleID == sample_id & !analysis_eligible,MP]
    idx_sample <- which(meta_st$SampleID == sample_id)
    
    for (mp in low_cov_mps) {
        score_col <- paste0(mp, "_score")
        if (score_col %in% colnames(meta_st)) {
            set(meta_st, i = idx_sample, j = score_col, value = NA_real_)
        }
    }
}

heatmap_celltype <- meta_st[,lapply(.SD,function(x) {
    if (all(!is.finite(x))) {
        NA_real_
    } else {
        mean(x[is.finite(x)])
    }}),
    by = CellType, .SDcols = score_cols]


heatmap_long <- melt(heatmap_celltype, id.vars="CellType", measure.vars=score_cols, variable.name="SEP", value.name="Score")
heatmap_long[, SEP:=sub("_score$", "", SEP)]
heatmap_long[!is.finite(Score), Score:=NA_real_]

heatmap_long[, z_score:={
    valid <- is.finite(Score)
    ans <- rep(NA_real_, .N)
    if(sum(valid)>=2){
        s <- sd(Score[valid])
        if(is.finite(s) && s>0) ans[valid] <- (Score[valid]-mean(Score[valid]))/s else ans[valid] <- 0
    }
    ans
}, by=SEP]


heatmap_dt <- dcast(heatmap_long, SEP~CellType, value.var="z_score", fill=NA_real_)
heatmap_mat <- as.matrix(heatmap_dt[, -1])
rownames(heatmap_mat) <- heatmap_dt$SEP

heatmap_mat <- pmax(pmin(heatmap_mat, 2.5), -2.5)

celltype_order <- c("Bcell", "Tcell", "Myeloid", "Endothelial", "Fibroblast", "Hepatocyte", "Epithelial")

for (ct_missing in setdiff(celltype_order, colnames(heatmap_mat))) {
    extra <- matrix(NA_real_, nrow(heatmap_mat), 1L,
                    dimnames = list(rownames(heatmap_mat), ct_missing))
    heatmap_mat <- cbind(heatmap_mat, extra)
}
heatmap_mat <- heatmap_mat[, celltype_order, drop=FALSE]
colnames(heatmap_mat) <- c("B cells","T cells","Myeloid","Endothelial","Fibroblast","Hepatocytes","Epithelial")

sep_names <- rownames(heatmap_mat)
sep_groups <- factor(ifelse(grepl("^N_SEP", sep_names), "N-SEPs", "M-SEPs"),levels = c("N-SEPs", "M-SEPs"))

plot_mat <- t(heatmap_mat)
colnames(plot_mat) <-  gsub('_','-',colnames(plot_mat))

sep_group_colors <- c("N-SEPs"="#7BB6D8", "M-SEPs"="#B8A7C8")

top_annotation <- HeatmapAnnotation(SEPs=sep_groups,
                                    col=list(SEPs=sep_group_colors), 
                                    simple_anno_size=unit(4, "mm"),
                                    annotation_name_side="left", 
                                    annotation_name_gp=gpar(fontsize=9, fontface="bold"), show_legend=TRUE)


heatmap_col_fun <- colorRamp2(c(-2.5, -1.25, 0, 1.25, 2.5), c("#244A7C", "#8BB3CF", "#FFFFFF", "#E59A91", "#9E3038"))
ht_st <- Heatmap(plot_mat, 
                 name="Z-score", 
                 col=heatmap_col_fun, 
                 na_col="#EFEFEF", 
                 cluster_rows=FALSE,
                 column_split=sep_groups, 
                 # Correlation clustering is undefined for all-NA/constant SEP columns.
                 cluster_columns=all(is.finite(plot_mat)) &&
                     all(apply(plot_mat, 2, function(x) sd(x) > 0)), 
                 cluster_column_slices=FALSE,
                 clustering_distance_columns="pearson", 
                 clustering_method_columns="average", column_gap=unit(2, "mm"),
                 top_annotation=top_annotation, 
                 show_row_names=TRUE,
                 row_names_side="left", row_names_gp=gpar(fontsize=9),
                 show_column_names=TRUE, column_names_rot=90,
                 column_names_gp=gpar(fontsize=9), column_names_max_height=unit(40, "mm"),
                 column_title_gp=gpar(fontsize=10, fontface="bold"), rect_gp=gpar(col="white", lwd=0.8),
                 heatmap_legend_param=list(
                     title="Z-score", at=c(-2.5, -1, 0, 1, 2.5), legend_height=unit(35, "mm"),
                     title_gp=gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=9)
                 ))


pdf(file = "Figure/ST/SEP/SEP_Heatmap/ST_heatmap_N_SEPs_condidata.pdf", width=7, height=3.5, useDingbats=FALSE)
draw(ht_st, heatmap_legend_side="right", annotation_legend_side="right", merge_legends=FALSE, padding=unit(c(5, 5, 5, 5), "mm"))
dev.off()


Cluster_list = readRDS(file.path(mp_dir,"Cluster_list_ori.RDS"))
MP_list = readRDS(file.path(mp_dir,"MP_List_original.RDS"))
sample_names <- as.character(unique(obj_hex$SampleID))
patient_names <- as.character(unique(obj_hex$PatientID))

sample_levels <- levels(obj_hex$SampleID)
patient_levels <- paste0('Pt',1:17)

sample_cols <- c("#5050FFFF", "#CE3D32FF", "#749B58FF", "#F0E685FF", "#466983FF", "#BA6338FF", "#802268FF", "#6BD76BFF",
                "#924822FF", "#837B8DFF", "#C75127FF", "#D58F5CFF", "#E4AF69FF", "#3B1B53FF", "#CDDEB7FF",
                "#E7C76FFF", "#5A655EFF", "#CC9900FF", "#99CC00FF", "#A9A9A9FF", "#AE1F63FF", "#D595A7FF",
                "#612A79FF", "#0099CCFF", "#0A47FFFF", "#FFC20AFF", "#991A00FF", "#996600FF", "#809900FF",
                "#79AF97FF", "#374E55FF", "#339900FF", "#008099FF", "#003399FF", "#990080FF", "#00D68FFF",
                "#33CC00FF", "#B24745FF", "#DF8F44FF", "#4DBBD5FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
                "#91D1C2FF", "#7E6148FF", "#B09C85FF", "#3B4992FF", "#008B45FF", "#631879FF", "#008280FF",
                "#BB0021FF", "#5F559BFF", "#A20056FF", "#808180FF")
names(sample_cols) = sample_names

patient_cols <- c("#fb8500", "#8dd3c7", "#d4a373", "#fdb462", "#bebada", "#fb8072", "#80b1d3",
                "#ffc9b9", "#ccebc5", "#bc80bd", "#f58549",  "#2ec4b6", "#d9d9d9",
                "#ffed6f", "#ffffb3", "#fccde5", "#b3de69")

names(patient_cols) = patient_names
custom_magma <- c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18)))

annotation_bar <- function(x, title, color, legend_ncol = 1, height_ratio = 0.025){
    d <- data.frame(x_pos = seq_along(x), value = x)
    n <- nrow(d)
    
    ggplot(d, aes(x = x_pos, y = 1, fill = value)) +
        geom_tile(width = 1, height = 1) +
        scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
        scale_fill_manual(values = color, name = title, drop = FALSE,
                          guide = guide_legend(ncol = legend_ncol, byrow = TRUE,
                                               keywidth = grid::unit(3, "mm"),
                                               keyheight = grid::unit(3, "mm"))) +
        coord_fixed(ratio = height_ratio * n) +
        theme_void() +
        theme(plot.margin = margin(0, 0, 0, 0))
}

plot_program_jaccard <- function(all_programs, overlap_mat, Cluster_list, clustered_only = FALSE, similarity_limits = c(2, 25)){
    clustered_programs <- unique(unlist(Cluster_list, use.names = FALSE))
    clustered_programs <- clustered_programs[clustered_programs %in% colnames(overlap_mat)]
    
    program_order <- if(clustered_only) clustered_programs else c(clustered_programs, setdiff(colnames(overlap_mat), clustered_programs))
 
    n_prog <- length(program_order)
    x_map <- setNames(seq_len(n_prog), program_order)
    
    program_mat <- all_programs[, program_order, drop = FALSE]
    n_genes <- vapply(seq_len(ncol(program_mat)), function(i){
        x <- program_mat[, i]
        length(unique(x[!is.na(x) & x != ""]))
    }, integer(1))
    names(n_genes) <- program_order
    
    overlap <- overlap_mat[program_order, program_order, drop = FALSE]
    union_size <- outer(n_genes, n_genes, "+") - overlap
    jaccard_mat <- 100 * overlap / union_size
    jaccard_mat[union_size == 0] <- NA_real_
    
    jaccard_df <- reshape2::melt(jaccard_mat, varnames = c("Program_y", "Program_x"), value.name = "Similarity")
    jaccard_df$x <- unname(x_map[as.character(jaccard_df$Program_x)])
    jaccard_df$y <- unname(x_map[as.character(jaccard_df$Program_y)])
    
    mp_order <- names(Cluster_list)
    mp_ranges <- data.frame()
    
    for(mp in mp_order){
        
        curr_programs <- Cluster_list[[mp]]
        curr_programs <- curr_programs[curr_programs %in% program_order]
        
        if(length(curr_programs) == 0) next
        idx <- x_map[curr_programs]
        mp_ranges <- rbind(mp_ranges, data.frame(MP = mp, xmin = min(idx)-0.5, xmax = max(idx)+0.5, ymin = min(idx)-0.5, ymax = max(idx)+0.5))
    }

    p_heatmap <- ggplot(jaccard_df, aes(x = x, y = y, fill = Similarity)) +
        geom_tile(width = 1, height = 1) +
        geom_rect(data = mp_ranges, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.7) +
        scale_x_continuous(limits = c(0.5, n_prog + 0.5), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0.5, n_prog + 0.5), expand = c(0, 0)) +
        scale_fill_gradientn(colours = custom_magma, limits = similarity_limits, oob = scales::squish,
                             name = "Similarity\n(Jaccard × 100)",
                             guide = guide_colourbar(barwidth = grid::unit(4, "mm"),
                                                     barheight = grid::unit(20, "mm"))) +
        coord_fixed() +
        theme_void() +
        theme(plot.margin = margin(0, 0, 0, 0))
    
    cluster_map <- setNames(rep(paste0("MP_", seq_along(Cluster_list)), lengths(Cluster_list)),
                            unlist(Cluster_list, use.names = FALSE))
    
    program_annotation <- data.frame(Program = program_order, stringsAsFactors = FALSE) 
    program_annotation$Sample <- sub("\\.(NMF|Leiden)\\..*$", "", program_annotation$Program)
    program_annotation$Sample <- sub("Tumor.", "", program_annotation$Sample)
    program_annotation$Patient <- sub("-.*$", "", program_annotation$Sample)
    program_annotation$MP <- unname(cluster_map[program_annotation$Program]) 
    
    mp_levels_all <- paste0("MP_", seq_along(Cluster_list))
    mp_levels <- if(clustered_only) mp_levels_all else c(mp_levels_all, "Unassigned")
    
    program_annotation$Sample <- factor(program_annotation$Sample, levels = sample_levels) 
    program_annotation$Patient <- factor(program_annotation$Patient, levels = patient_levels)
    program_annotation$MP <- factor(program_annotation$MP, levels = mp_levels)
     
    #sample_cols <- setNames(cb_palette[seq_along(sample_levels)], sample_levels)
    #patient_cols <- setNames(cb_patient[seq_along(patient_levels)], patient_levels)
    #mp_cols <- setNames(c(ggsci::pal_igv()(length(mp_levels_all)), "#D9D9D9"), c(mp_levels_all, "Unassigned"))
    
    p_sample <- annotation_bar(program_annotation$Sample, "Sample", sample_cols, legend_ncol = 2)
    p_patient <- annotation_bar(program_annotation$Patient, "Patient", patient_cols, legend_ncol = 2)
    #p_mp <- annotation_bar(program_annotation$MP, "MP", mp_cols[mp_levels], legend_ncol = 1)
 
   (p_patient / p_sample / p_heatmap) + patchwork::plot_layout(guides = "collect")
}

options(repr.plot.width = 12, repr.plot.height = 10)

p_jaccard_clustered <- plot_program_jaccard(all_programs = all_programs,
                                            overlap_mat = nmf_intersect,
                                            Cluster_list = Cluster_list,
                                            clustered_only = TRUE)
p_jaccard_clustered
ggsave(file.path(figure_dir,"all_MP_pragrams_original.pdf"),p_jaccard_clustered,width = 12,height = 10)

mp_program_annotation <- rbindlist(lapply(seq_along(Cluster_list), function(i){
    programs <- Cluster_list[[i]]
    SampleID = sub("\\.(NMF|Leiden)\\..*$", "", programs)
    data.table(
        MP = paste0("MP_", i),
        Program = programs,
        Method = ifelse(grepl("\\.Leiden\\.", programs), "Leiden", "NMF"),
        Sample = SampleID,
        Patient = sub("-.*$", "", SampleID)  
    )
}))

mp_summary <- mp_program_annotation[, {
    sample_count <- sort(table(Sample), decreasing = TRUE)
    nmf_programs <- Program[Method == "NMF"]
    leiden_programs <- Program[Method == "Leiden"]
    .(n_programs = .N,
      n_NMF = length(nmf_programs),
      n_Leiden = length(leiden_programs),
      n_samples = uniqueN(Sample),
      n_patient = uniqueN(Patient))
}, by = MP]

mp_summary[, MP_num := as.integer(sub("MP_", "", MP))]
setorder(mp_summary, MP_num)
mp_summary[, MP_num := NULL]
mp_summary


library(UCell)
h_list = list(MP_list[[10]],MP_list[[16]])
MP_use <- lapply(h_list, function(x) {
    intersect(unique(x), rownames(obj_use))
})  

obj_use <-  NormalizeData(obj_hex, normalization.method = "LogNormalize", scale.factor = 10000)
obj_use <- AddModuleScore_UCell(obj_use, features = MP_use)
cor(obj_use$signature_1_UCell,obj_use$signature_2_UCell, method="spearman")
MP10_programs <- Cluster_list[["Cluster_10"]]
MP16_programs <- Cluster_list[["Cluster_16"]]
Merged_hepatocyte_programs <- c(MP10_programs, MP16_programs)

all_history_genes <- c()
for(program_id in Merged_hepatocyte_programs){
    if(program_id %in% colnames(all_programs)){
        genes <- all_programs[, program_id]
        all_history_genes <- c(all_history_genes,genes)
    }
}
all_history_genes <- all_history_genes[!is.na(all_history_genes) & all_history_genes != ""]
gene_frequency <- sort(table(all_history_genes),decreasing = TRUE)
Merged_MP10_MP16_genes <- names(head(gene_frequency,50))
Merged_MP10_MP16_genes

Final_MP_list <- MP_list
 
Final_MP_list$MP_2 <- NULL
Final_MP_list$MP_7 <- NULL
Final_MP_list$MP_9 <- NULL
Final_MP_list$MP_17 <- NULL
Final_MP_list$MP_18 <- NULL
Final_MP_list$MP_10 <- Merged_MP10_MP16_genes
Final_MP_list$MP_16 <- NULL

names(Final_MP_list) =   paste0("MP", seq_along(Final_MP_list))
Final_MP_list
write.csv(Final_MP_list,file.path(mp_dir,"N_SEP_List.csv"))
saveRDS(Final_MP_list,file.path(mp_dir,"N_SEP_List.RDS"))
p_jaccard_clustered <- plot_program_jaccard(all_programs = all_programs,overlap_mat = nmf_intersect,Cluster_list = Final_Cluster_list,clustered_only = TRUE)
marker_list <- list(Epithelial = unique(c("EPCAM", "TACSTD2","KRT7", "KRT8", "KRT18", "KRT19", "KRT17", "KRT23", "KRT6A", "KRT6B",
                                          "CLDN1", "CLDN3", "CLDN4", "CLDN7", "CLDN10","MUC1", "MUC13","ELF3", "MAL2", "DSG2", "DSP",
                                          "SPINT1", "SPINT2", "LAD1","ANXA4", "AGR2", "CFTR", "MMP7", "SPINK1", "CEACAM5", "CEACAM6", "SPP1",
                                          "LAMC2", "LAMB3", "LAMA3", "ITGB6", "S100P", "S100A14", "WFDC2", "TFF1", "TFF3", 
                                          "CXCL1", "CXCL5", "CXCL6", "CXCL8")),                                     
                    Immune = unique(c("PTPRC", "LAPTM5", "LCP1", "SRGN", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                                      "LYZ", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "CD68", "CD163", "MSR1", "LGMN",
                                      "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "IL7R", "CCL5", "MS4A1", "CD79A", "CD79B", "CD19", "CD37",
                                      "JCHAIN", "MZB1", "TNFRSF17", "NKG7", "GNLY", "KLRD1")),
                    Stromal = unique(c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A3", "COL11A1", "COL14A1", 
                                       "DCN", "LUM", "BGN", "AEBP1", "THBS2", "POSTN", "CTHRC1", "PDGFRA", "PDGFRB",
                                       "ACTA2", "TAGLN", "MYL9", "RGS5", "NOTCH3", "MYH11", "SFRP2", "SFRP4", "COL10A1", "MMP11", "COMP", "INHBA")))
 
epithelial_core <- c("EPCAM", "TACSTD2", "KRT7", "KRT8", "KRT18", "KRT19", "KRT17", "KRT23", "KRT6A", "KRT6B", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "CLDN10",
                     "MUC1", "MUC13", "ELF3", "MAL2", "DSG2", "DSP", "SPINT1", "SPINT2", "LAD1")

 clean_program_genes <- function(x) {    
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    unique(x)
}
 
program_annotation <- do.call(rbind, lapply(colnames(all_nmf_mat),function(program_id) {
    genes <- clean_program_genes(all_nmf_mat[, program_id])
    epithelial_overlap <- intersect(genes,marker_list$Epithelial)
    immune_overlap <- intersect(genes,marker_list$Immune)
    stromal_overlap <- intersect(genes,marker_list$Stromal)
    
    core_overlap <- intersect(genes,epithelial_core)
    
    data.frame(Program = program_id,
               N_Epithelial = length(epithelial_overlap),
               N_Immune = length(immune_overlap),
               N_Stromal = length(stromal_overlap),
               N_Epithelial_core = length(core_overlap),
               
               Epithelial_genes = paste(epithelial_overlap,collapse = ";"),
               Immune_genes = paste(immune_overlap,collapse = ";"),
               Stromal_genes = paste(stromal_overlap,collapse = ";"),
               stringsAsFactors = FALSE
    )}))

rownames(program_annotation) <- NULL

 
score_matrix <- as.matrix(program_annotation[,c("N_Epithelial","N_Immune","N_Stromal")])
colnames(score_matrix) <- c("Epithelial","Immune","Stromal" )

program_annotation$Dominant_class <- apply(score_matrix, 1, function(x) {
        
        if (max(x) == 0) {
            return("Unclassified")
        }
        
        best_class <- names(x)[x == max(x)]
        
        if (length(best_class) > 1) {
            return(paste(best_class, collapse = "_"))
        }
        best_class
    }
)

program_annotation$Tumor_related <- with(program_annotation, N_Epithelial >= 2 & N_Epithelial_core >= 1 & N_Epithelial > N_Immune & N_Epithelial > N_Stromal )

tumor_program_annotation <- program_annotation[program_annotation$Tumor_related, , drop = FALSE]
tumor_program_annotation <- tumor_program_annotation[order(-tumor_program_annotation$N_Epithelial), , drop = FALSE]

tumor_program_ids <- unique(as.character(tumor_program_annotation$Program))
tumor_prog_genes_ls <- lapply(prog_genes_ls,function(nmf_mat) {        
        keep_ids <- intersect(colnames(nmf_mat), tumor_program_ids)
        nmf_mat[, keep_ids, drop = FALSE]
    }
)

tumor_program_n <- sapply(tumor_prog_genes_ls, ncol)
tumor_prog_genes_ls <- tumor_prog_genes_ls[tumor_program_n > 0]

nmf_programs_original <- lapply(tumor_prog_genes_ls, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
tumor_nmf_mat <- do.call(cbind, nmf_programs_original)
dim(tumor_nmf_mat)

nmf_robust_id <- Robust_NMF_Programs(nmf_programs_original, intra_min = 30, intra_max = 15, inter_min = 10)
robust_tumor_nmf_mat <- tumor_nmf_mat[, nmf_robust_id, drop = FALSE]
program_info_nmf <- data.frame(program = colnames(robust_tumor_nmf_mat),
                               sample = sub("\\.NMF\\..*$", "", colnames(robust_tumor_nmf_mat)),
                               stringsAsFactors = FALSE)

program_jaccard <- apply(robust_tumor_nmf_mat, 2, function(x){
    apply(robust_tumor_nmf_mat, 2, function(y){
        length(intersect(x, y)) / length(union(x, y))
    })
})
diag(program_jaccard) <- 1 

library(dynamicTreeCut)
library(pheatmap)
hc <- hclust(as.dist(1 - program_jaccard),method = "average")
cluster_dynamic <- cutreeDynamic(dendro = hc,
                                 distM = as.matrix(1 - program_jaccard),
                                 deepSplit = 2,
                                 pamRespectsDendro = TRUE,
                                 minClusterSize = 4)

names(cluster_dynamic) <- colnames(program_jaccard)
program_info_nmf$cluster <- unname(cluster_dynamic[program_info_nmf$program])
table(program_info_nmf$cluster)

sample_patient <- unique(obj_hex@meta.data[, c("SampleID", "PatientID")])
program_info_nmf$patient <- as.character(sample_patient$PatientID[match(program_info_nmf$sample,as.character(sample_patient$SampleID))])

nmf_family_summary <- program_info_nmf %>%
  filter(!is.na(cluster), cluster > 0) %>%
  group_by(cluster) %>%
  summarise(
    n_programs = n_distinct(program),
    n_samples  = n_distinct(sample),
    n_patients = n_distinct(patient),
    patient_ids = paste(sort(unique(patient)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(cluster) %>%
  mutate(cluster = paste0("NMF_F", cluster)) %>% rename(NMF_family = cluster)


get_family_gene_frequency <- function(cluster_id,cluster_assign,programs_mat){
    
    prog_names <- names(cluster_assign)[cluster_assign == cluster_id]
    genes <- unlist(lapply(prog_names, function(p){programs_mat[, p]}))
    freq <- sort(table(genes), decreasing = TRUE)
    data.frame(gene = names(freq), n_programs = as.integer(freq), fraction = as.integer(freq) / length(prog_names), stringsAsFactors = FALSE)
}

valid_clusters <- sort(unique(cluster_dynamic[cluster_dynamic > 0]))

NMF_family_freq <- lapply(valid_clusters, function(cl){ 
    get_family_gene_frequency(cluster_id = cl,cluster_assign = cluster_dynamic,programs_mat = robust_tumor_nmf_mat)
    }
)
names(NMF_family_freq) <- paste0("NMF_F", valid_clusters)

all_family_genes <- sort(unique(unlist(lapply(NMF_family_freq, function(x){ x$gene}))))
family_freq_mat <- matrix(0, nrow = length(all_family_genes), ncol = length(NMF_family_freq),dimnames = list(all_family_genes,names(NMF_family_freq)))

for(f in names(NMF_family_freq)){
    df <- NMF_family_freq[[f]]
    family_freq_mat[df$gene,f] <- df$fraction
}

n_family <- ncol(family_freq_mat)
gene_family_n <- rowSums(family_freq_mat > 0)
gene_idf <- log((n_family + 1) / (gene_family_n + 1)) + 1
family_freq_weighted <- sweep(family_freq_mat,1,gene_idf, "*")

family_norm <- sqrt(colSums(family_freq_weighted^2))
family_similarity <- crossprod(family_freq_weighted) / outer(family_norm, family_norm)
diag(family_similarity) <- 1
hc_family <- hclust(as.dist(1 - family_similarity), method = "complete")
plot(hc_family, hang = -1, main = "Conservative merging of tumor NMF families", xlab = "NMF families", ylab = "1 - weighted cosine similarity")
pheatmap(family_similarity, cluster_rows = hc_family,
    cluster_cols = hc_family,
    color = colorRampPalette(c("white", "#FDD0A2", "#E6550D", "#7F0000"))(100),
    border_color = NA,
    display_numbers = TRUE,
    number_format = "%.2f",
    main = "Similarity between tumor NMF families")
NMF_family_gene_list <- lapply(NMF_family_freq, function(df) head(df$gene, 50))
saveRDS(NMF_family_gene_list,'./Data/ST/SEP/Module1/epithelial_nmf_original.rds')

MP_list <- readRDS('./Data/ST/SEP/Module1/epithelial_nmf_original.rds')
obj_hex <- readRDS('Data/ST/Spatial_Niche/All_hex_in_seurat_original_niche.RDS')

min_mean_expression <- 0.1
conserved_genes <- 0.3
min_celltype_hex <- 20L
meta <- obj_hex@meta.data
obj_hex$CellType_original <- as.character(meta$CellType)
sample_names <- unique(as.character(meta$SampleID))
counts_all <- LayerData(obj_hex, assay = "RNA", layer = "counts")
score_result <- setNames(vector("list", length(sample_names)), sample_names)


for (sample_id in sample_names) {
    
    cat("Sample:", sample_id,"\n")
    spot_ids <- rownames(meta)[as.character(meta[["SampleID"]]) == sample_id]
    counts_sample <- counts_all[,spot_ids, drop = FALSE]
    
    celltype_sample <- as.character(meta[spot_ids,"CellType"])
    names(celltype_sample) <- spot_ids
    type_sample <- as.character(meta[spot_ids, "Type"])
    names(type_sample) <- spot_ids
    
    score_result[[sample_id]] <- score_one_sample(counts_mat = counts_sample,
                                                  MP_list = MP_list,
                                                  sample_id = sample_id,
                                                  celltype_vec = celltype_sample,
                                                  type_vec = type_sample,
                                                  min_mean_expression = min_mean_expression,
                                                  conserved_genes = conserved_genes,
                                                  min_celltype_hex = min_celltype_hex,
                                                  seed = 42)
    invisible(gc())
}



#result
coverage_all <- rbindlist(
    lapply(names(score_result), function(sample_id) {
        x <- as.data.table(score_result[[sample_id]]$coverage)
        x[, SampleID := sample_id]
        x
    }),
    fill = TRUE)

collect_audit <- function(field) rbindlist(lapply(names(score_result), function(sid) {
    x <- copy(as.data.table(score_result[[sid]][[field]]))
    x[, SampleID := rep(sid, .N)]
    x
}), use.names = TRUE, fill = TRUE)
 
raw_score_all <- rbindlist(lapply(score_result,function(x){
    score_df <- as.data.frame(x$raw_scores)
    score_df$SpotID <- rownames(score_df)
    score_df}), fill = TRUE)

raw_score_all <- as.data.frame(raw_score_all)
rownames(raw_score_all) <- raw_score_all$SpotID
raw_score_all$SpotID <- NULL

score_idx <- match(colnames(obj_hex), rownames(raw_score_all))

for (mp in colnames(raw_score_all)) {
    obj_hex[[paste0(mp,"_score")]] <- raw_score_all[score_idx, mp]
}

weighted_score_all <- rbindlist(lapply(score_result,function(x){
    score_df <- as.data.frame(x$weighted_scores)
    score_df$SpotID <- rownames(score_df)
    score_df
}), fill = TRUE)

weighted_score_all <- as.data.frame(weighted_score_all)
rownames(weighted_score_all) <- weighted_score_all$SpotID
weighted_score_all$SpotID <- NULL

weighted_idx <- match(colnames(obj_hex), rownames(weighted_score_all))

for (mp in colnames(weighted_score_all)) {
    obj_hex[[paste0(mp, "_weighted_score")]] <- weighted_score_all[weighted_idx, mp]
}

meta <- obj_hex@meta.data
idx_normal_epi <- meta$Type %in% "Normal" & meta$CellType %in% "Epithelial"
celltype_updated <- as.character(meta$CellType)
celltype_updated[idx_normal_epi] <- "unassigned"
obj_hex$CellType <- celltype_updated
 

# ------------
meta_export <- as.data.table(obj_hex@meta.data,keep.rownames = "SpotID")
head(meta_export)

raw_score_cols <- paste0(names(MP_list), "_score")
raw_score_cols <- intersect(raw_score_cols, colnames(meta_export))

weighted_score_cols <- paste0(names(MP_list), "_weighted_score")
weighted_score_cols <- intersect(weighted_score_cols, colnames(meta_export))

export_cols <- c(
    "SpotID", "sample_id", "SampleID", "PatientID", "Type", "Region", "CellType", "CellType_original",
    raw_score_cols,
    weighted_score_cols
)

export_cols <- intersect(export_cols, colnames(meta_export))
mp_export <- meta_export[, ..export_cols]
head(mp_export)

 # ---------------------- Plot -----------------------------
meta_st <- as.data.table(obj_hex@meta.data, keep.rownames="SpotID")
meta_st <- meta_st[!is.na(CellType) & CellType!="" & CellType!="unassigned"]
score_cols <- grep("NMF_F[0-9]+_score$", colnames(meta_st), value=TRUE)
score_cols



heatmap_full <- meta_st[, lapply(.SD, function(x) {
    if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
}), by = CellType, .SDcols = score_cols]

 
for (sample_id in unique(meta_st$SampleID)) {
    
    low_cov_mps <- coverage_all[SampleID == sample_id & !analysis_eligible,MP]
    idx_sample <- which(meta_st$SampleID == sample_id)
    
    for (mp in low_cov_mps) {
        score_col <- paste0(mp, "_score")
        if (score_col %in% colnames(meta_st)) {
            set(meta_st, i = idx_sample, j = score_col, value = NA_real_)
        }
    }
}

heatmap_celltype <- meta_st[,lapply(.SD,function(x) {
    if (all(!is.finite(x))) {
        NA_real_
    } else {
        mean(x[is.finite(x)])
    }}),
    by = CellType, .SDcols = score_cols]


heatmap_long <- melt(heatmap_celltype, id.vars="CellType", measure.vars=score_cols, variable.name="SEP", value.name="Score")
heatmap_long[, SEP:=sub("_score$", "", SEP)]
heatmap_long[!is.finite(Score), Score:=NA_real_]

heatmap_long[, z_score:={
    valid <- is.finite(Score)
    ans <- rep(NA_real_, .N)
    if(sum(valid)>=2){
        s <- sd(Score[valid])
        if(is.finite(s) && s>0) ans[valid] <- (Score[valid]-mean(Score[valid]))/s else ans[valid] <- 0
    }
    ans
}, by=SEP]


heatmap_dt <- dcast(heatmap_long, SEP~CellType, value.var="z_score", fill=NA_real_)
heatmap_mat <- as.matrix(heatmap_dt[, -1])
rownames(heatmap_mat) <- heatmap_dt$SEP

heatmap_mat <- pmax(pmin(heatmap_mat, 2.5), -2.5)

celltype_order <- c("Bcell", "Tcell", "Myeloid", "Endothelial", "Fibroblast", "Hepatocyte", "Epithelial")

for (ct_missing in setdiff(celltype_order, colnames(heatmap_mat))) {
    extra <- matrix(NA_real_, nrow(heatmap_mat), 1L,
                    dimnames = list(rownames(heatmap_mat), ct_missing))
    heatmap_mat <- cbind(heatmap_mat, extra)
}
heatmap_mat <- heatmap_mat[, celltype_order, drop=FALSE]
colnames(heatmap_mat) <- c("B cells","T cells","Myeloid","Endothelial","Fibroblast","Hepatocytes","Epithelial")

sep_names <- rownames(heatmap_mat)
sep_groups <- factor(ifelse(grepl("^N_SEP", sep_names), "N-SEPs", "M-SEPs"),levels = c("N-SEPs", "M-SEPs"))

plot_mat <- t(heatmap_mat)
colnames(plot_mat) <-  gsub('_','-',colnames(plot_mat))

sep_group_colors <- c("N-SEPs"="#7BB6D8", "M-SEPs"="#B8A7C8")

top_annotation <- HeatmapAnnotation(SEPs=sep_groups,
                                    col=list(SEPs=sep_group_colors), 
                                    simple_anno_size=unit(4, "mm"),
                                    annotation_name_side="left", 
                                    annotation_name_gp=gpar(fontsize=9, fontface="bold"), show_legend=TRUE)


heatmap_col_fun <- colorRamp2(c(-2.5, -1.25, 0, 1.25, 2.5), c("#244A7C", "#8BB3CF", "#FFFFFF", "#E59A91", "#9E3038"))
ht_st <- Heatmap(plot_mat, 
                 name="Z-score", 
                 col=heatmap_col_fun, 
                 na_col="#EFEFEF", 
                 cluster_rows=FALSE,
                 column_split=sep_groups, 
                 # Correlation clustering is undefined for all-NA/constant SEP columns.
                 cluster_columns=all(is.finite(plot_mat)) &&
                     all(apply(plot_mat, 2, function(x) sd(x) > 0)), 
                 cluster_column_slices=FALSE,
                 clustering_distance_columns="pearson", 
                 clustering_method_columns="average", column_gap=unit(2, "mm"),
                 top_annotation=top_annotation, 
                 show_row_names=TRUE,
                 row_names_side="left", row_names_gp=gpar(fontsize=9),
                 show_column_names=TRUE, column_names_rot=90,
                 column_names_gp=gpar(fontsize=9), column_names_max_height=unit(40, "mm"),
                 column_title_gp=gpar(fontsize=10, fontface="bold"), rect_gp=gpar(col="white", lwd=0.8),
                 heatmap_legend_param=list(
                     title="Z-score", at=c(-2.5, -1, 0, 1, 2.5), legend_height=unit(35, "mm"),
                     title_gp=gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=9)))

final_family_groups <- list(
    Tumor_MP1=c("NMF_F2","NMF_F6","NMF_F10"),
    Tumor_MP2=c("NMF_F3"),
    Tumor_MP3=c("NMF_F4"),
    Tumor_MP4=c("NMF_F5"),
    Tumor_MP5=c("NMF_F7","NMF_F9")
)
make_final_consensus <- function(families, cluster_dynamic, all_programs) {
    family_ids <- as.integer(sub("NMF_F","",families))
    programs <- names(cluster_dynamic)[cluster_dynamic %in% family_ids]
    if(length(programs)==0) stop("No programs found for: ",paste(families,collapse=", "))
    genes <- as.character(unlist(all_programs[,programs,drop=FALSE],use.names=FALSE))
    genes <- genes[!is.na(genes) & genes!=""]
    freq <- sort(table(genes),decreasing=TRUE)
    gene_table <- data.frame(gene=names(freq),n_programs=as.integer(freq),fraction=as.integer(freq)/length(programs),stringsAsFactors=FALSE)
    list(families=families,programs=programs,n_programs=length(programs),table=gene_table)
}

Final_Tumor_res <- lapply(final_family_groups,function(x) make_final_consensus(x,cluster_dynamic,robust_tumor_nmf_mat))
sapply(Final_Tumor_res,function(x) x$n_programs)

mp_names <- names(Final_Tumor_res)
all_tumor_genes <- sort(unique(unlist(lapply(Final_Tumor_res,function(x) x$table$gene))))
count_mat <- matrix(0,nrow=length(all_tumor_genes),ncol=length(mp_names),dimnames=list(all_tumor_genes,mp_names))
fraction_mat <- matrix(0,nrow=length(all_tumor_genes),ncol=length(mp_names),dimnames=list(all_tumor_genes,mp_names))

for(mp in mp_names) {
    df <- Final_Tumor_res[[mp]]$table
    stopifnot(all(c("gene","n_programs","fraction") %in% colnames(df)))
    df <- df[!is.na(df$gene) & df$gene!="",,drop=FALSE]
    df <- df[!duplicated(df$gene),,drop=FALSE]
    count_mat[df$gene,mp] <- df$n_programs
    fraction_mat[df$gene,mp] <- df$fraction
}

stopifnot(identical(dim(count_mat),dim(fraction_mat)))
stopifnot(all(fraction_mat>=0 & fraction_mat<=1))
colSums(count_mat>0)

assignment_table <- do.call(rbind,lapply(mp_names,function(mp) {
    other_mps <- setdiff(mp_names,mp)
    other_max <- apply(fraction_mat[,other_mps,drop=FALSE],1,max)
    data.frame(gene=rownames(fraction_mat),MP=mp,count=count_mat[,mp],fraction=fraction_mat[,mp],other_max=other_max,specificity=fraction_mat[,mp]-other_max,stringsAsFactors=FALSE)
}))

rownames(assignment_table) <- NULL
assignment_table <- assignment_table[assignment_table$count>0,,drop=FALSE] 

top_n <- 50
fraction_threshold <- 0.4
allowed_assignment <- do.call(rbind,lapply(split(assignment_table, assignment_table$gene), function(df) {

            df <- df[order(-df$fraction, -df$count, -df$specificity, df$MP),,drop = FALSE]

            if (nrow(df) == 1) {
                return(df)
            }

            high_fraction <- df[df$fraction > fraction_threshold,,drop = FALSE]

            if (nrow(high_fraction) > 0) {
                return(high_fraction)
            }
            df[1, , drop = FALSE]
        }
    )
)

rownames(allowed_assignment) <- NULL

Final_Tumor_MP_list <- setNames(vector("list", length(mp_names)), mp_names)

for (mp in mp_names) {
    df <- allowed_assignment[allowed_assignment$MP == mp, , drop = FALSE]
    df <- df[order(-df$fraction,-df$count,-df$specificity,df$gene),,drop = FALSE]
    Final_Tumor_MP_list[[mp]] <- head(df$gene, top_n)
}

lengths(Final_Tumor_MP_list)
Final_Tumor_MP_list

saveRDS(Final_Tumor_MP_list,file = file.path(mp_dir,"M_SEP_List.RDS"))
Final_MP_list <- readRDS(file = file.path(mp_dir,"N_SEP_List.RDS"))
allMP = c(Final_MP_list,Final_Tumor_MP_list)
names(allMP) = c(paste0('N_SEP',1:12),paste0('M_SEP',1:5))

write.csv(allMP,file.path(mp_dir,"all_MP_list_New.csv"))
saveRDS(allMP,file.path(mp_dir,"all_MP_list_New.RDS"))



Final_Tumor_Cluster_list <- lapply(final_family_groups,function(families) {
    family_ids <- as.integer(sub("NMF_F","",families))
    names(cluster_dynamic)[cluster_dynamic %in% family_ids]
})

names(Final_Tumor_Cluster_list) <- names(final_family_groups)

tumor_program_overlap <- apply(robust_tumor_nmf_mat,2,function(x) {
    apply(robust_tumor_nmf_mat,2,function(y) length(intersect(x,y)))
})
rownames(tumor_program_overlap) <- colnames(robust_tumor_nmf_mat)
colnames(tumor_program_overlap) <- colnames(robust_tumor_nmf_mat)
Final_Tumor_Cluster_list_ordered <- lapply(
    Final_Tumor_Cluster_list,
    function(programs) {
        if(length(programs)<=2) return(programs)
        sim <- program_jaccard[programs,programs,drop=FALSE]
        hc_local <- hclust(as.dist(1-sim),method="average")
        programs[hc_local$order]
    }
)
p_tumor_program_jaccard <- plot_program_jaccard(
    all_programs=robust_tumor_nmf_mat,
    overlap_mat=tumor_program_overlap,
    Cluster_list=Final_Tumor_Cluster_list_ordered,
    clustered_only=TRUE,
    similarity_limits=c(5,30)
)
ggsave(p_tumor_program_jaccard,file = 'Figure/ST/SEP/SEP_Heatmap/02_M-SEP_heatmap.pdf',width = 6,height = 5)
save(Final_Cluster_list,all_programs,Final_Tumor_Cluster_list,robust_tumor_nmf_mat,file = 'Data/ST/SEP/Plot.rdata')


load('/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC/Data/ST/SEP/Plot.rdata')
obj_hex <- readRDS('Data/ST/Spatial_Niche/All_hex_in_seurat_clean_final.RDS')

N_program_info <- rbindlist(
    lapply(seq_along(Final_Cluster_list), function(i) {
        programs <- Final_Cluster_list[[i]]
        sample_id <- sub("\\.(NMF|Leiden)\\..*$", "", programs)
        data.table(SEP = paste0("N_SEP", i), SEP_type = "N-SEP", Program = programs,Method = ifelse(grepl("\\.Leiden\\.", programs), "Leiden", "NMF"), SampleID = sample_id, PatientID = sub("-.*$", "", sample_id))
    })
)

M_program_info <- rbindlist(
    lapply(seq_along(Final_Tumor_Cluster_list), function(i) {
        programs <- Final_Tumor_Cluster_list[[i]]
        sample_id <- sub("\\.NMF\\..*$", "", programs)
        data.table(SEP = paste0("M_SEP", i), SEP_type = "M-SEP", Program = programs, Method = "NMF", SampleID = sample_id, PatientID = sub("-.*$", "", sample_id))
    })
)

S3_program_info <- rbind(N_program_info,M_program_info, use.names = TRUE, fill = TRUE)

S3A_summary <- S3_program_info[, .(n_candidate_programs = uniqueN(Program),
                                   n_NMF = uniqueN(Program[Method == "NMF"]), 
                                   n_Leiden = uniqueN(Program[Method == "Leiden"])),by = .(SEP_type, SEP)]

sep_levels <- c(paste0("N_SEP", 1:12), paste0("M_SEP", 1:5))
S3A_summary[, SEP := factor(SEP, levels = sep_levels)]

all_sample_info <- unique(as.data.table(obj_hex@meta.data)[,.(SampleID, PatientID)])
all_sample_info <- unique(all_sample_info)
S3C_presence <- unique(S3_program_info[, .(SEP, SampleID, Present = 1L)])
S3C_full <- CJ(SEP = sep_levels, SampleID = as.character(all_sample_info$SampleID), unique = TRUE)
S3C_full <- merge(S3C_full,S3C_presence,by = c("SEP", "SampleID"),all.x = TRUE)
S3C_full[is.na(Present), Present := 0L]
S3C_full <- merge(S3C_full, all_sample_info, by = "SampleID", all.x = TRUE)
sample_order <- all_sample_info[order(PatientID, SampleID), SampleID]
S3C_full[, SampleID := factor(SampleID, levels = sample_order)]
S3C_full[ , SEP := factor(SEP, levels = rev(sep_levels))]

patient_names <- as.character(unique(obj_hex$PatientID))
names(patient_cols) = patient_names

presence_mat <- dcast(S3C_full, SEP ~ SampleID, value.var = "Present")
presence_mat <- as.data.frame(presence_mat)
sep_names <- as.character(presence_mat$SEP)
presence_mat$SEP <- NULL
presence_mat <- as.matrix(presence_mat)
rownames(presence_mat) <- sep_names
presence_mat <- presence_mat[sep_levels, sample_order, drop = FALSE]

sample_patient <- all_sample_info$PatientID[match(colnames(presence_mat), all_sample_info$SampleID)]
names(sample_patient) <- colnames(presence_mat)
ha <- HeatmapAnnotation(Patient = sample_patient, col = list(Patient = patient_cols), annotation_name_side = "left")

ht_S3C <- Heatmap(presence_mat,name = "SEP support", col = c("0" = "#EEEEEE", "1" = "#B2182B"),
    top_annotation = ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_rot = 90,
    column_names_gp = gpar(fontsize = 5),
    row_names_gp = gpar(fontsize = 8),
    rect_gp = gpar(col = "white", lwd = 0.3),
    column_title = "Samples"
)
draw(ht_S3C)


MP_list <- readRDS('Data/ST/SEP/Module1/03_metaprograms/all_MP_list_New.RDS')
obj_hex <- readRDS('Data/ST/Spatial_Niche/All_hex_in_seurat_original_niche.RDS')
DefaultAssay(obj_hex) <- "RNA"
min_mean_expression <- 0.1
conserved_genes <- 0.3
min_celltype_hex <- 20L

meta <- obj_hex@meta.data
obj_hex$CellType_original <- as.character(meta$CellType)

sample_names <- unique(as.character(meta$SampleID))
sample_names

counts_all <- LayerData(obj_hex, assay = "RNA", layer = "counts")

score_result <- setNames(vector("list", length(sample_names)), sample_names)
score_result

for (sample_id in sample_names) {
    
    cat("Sample:", sample_id,"\n")
    spot_ids <- rownames(meta)[as.character(meta[["SampleID"]]) == sample_id]
    counts_sample <- counts_all[,spot_ids, drop = FALSE]
    
    celltype_sample <- as.character(meta[spot_ids,"CellType"])
    names(celltype_sample) <- spot_ids
    type_sample <- as.character(meta[spot_ids, "Type"])
    names(type_sample) <- spot_ids
    
    score_result[[sample_id]] <- score_one_sample(counts_mat = counts_sample,
                                                  MP_list = MP_list,
                                                  sample_id = sample_id,
                                                  celltype_vec = celltype_sample,
                                                  type_vec = type_sample,
                                                  min_mean_expression = min_mean_expression,
                                                  conserved_genes = conserved_genes,
                                                  min_celltype_hex = min_celltype_hex,
                                                  seed = 42)
    invisible(gc())
}

coverage_all <- rbindlist(
    lapply(names(score_result), function(sample_id) {
        x <- as.data.table(score_result[[sample_id]]$coverage)
        x[, SampleID := sample_id]
        x
    }),
    fill = TRUE)

collect_audit <- function(field) rbindlist(lapply(names(score_result), function(sid) {
    x <- copy(as.data.table(score_result[[sid]][[field]]))
    x[, SampleID := rep(sid, .N)]
    x
}), use.names = TRUE, fill = TRUE)

raw_score_all <- rbindlist(lapply(score_result,function(x){
    score_df <- as.data.frame(x$raw_scores)
    score_df$SpotID <- rownames(score_df)
    score_df}), fill = TRUE)

raw_score_all <- as.data.frame(raw_score_all)
rownames(raw_score_all) <- raw_score_all$SpotID
raw_score_all$SpotID <- NULL

score_idx <- match(colnames(obj_hex), rownames(raw_score_all))

for (mp in colnames(raw_score_all)) {
    obj_hex[[paste0(mp,"_score")]] <- raw_score_all[score_idx, mp]
}

weighted_score_all <- rbindlist(lapply(score_result,function(x){
    score_df <- as.data.frame(x$weighted_scores)
    score_df$SpotID <- rownames(score_df)
    score_df
}), fill = TRUE)

weighted_score_all <- as.data.frame(weighted_score_all)
rownames(weighted_score_all) <- weighted_score_all$SpotID
weighted_score_all$SpotID <- NULL

weighted_idx <- match(colnames(obj_hex), rownames(weighted_score_all))

for (mp in colnames(weighted_score_all)) {
    obj_hex[[paste0(mp, "_weighted_score")]] <- weighted_score_all[weighted_idx, mp]
}

meta <- obj_hex@meta.data
idx_normal_epi <- meta$Type %in% "Normal" & meta$CellType %in% "Epithelial"
celltype_updated <- as.character(meta$CellType)
celltype_updated[idx_normal_epi] <- "unassigned"
obj_hex$CellType <- celltype_updated

meta_export <- as.data.table(obj_hex@meta.data,keep.rownames = "SpotID")
head(meta_export)

raw_score_cols <- paste0(names(MP_list), "_score")
raw_score_cols <- intersect(raw_score_cols, colnames(meta_export))

weighted_score_cols <- paste0(names(MP_list), "_weighted_score")
weighted_score_cols <- intersect(weighted_score_cols, colnames(meta_export))

export_cols <- c("SpotID", "sample_id", "SampleID", "PatientID", "Type", "Region", "CellType", "CellType_original", raw_score_cols, weighted_score_cols)
export_cols <- intersect(export_cols, colnames(meta_export))
mp_export <- meta_export[, ..export_cols]
head(mp_export)

fwrite(mp_export, file = "/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC/Data/ST/SEP/All_hex_MP_scores_only_for_h5ad.csv")


# ---------------------- Plot -----------------------------
meta_st <- as.data.table(obj_hex@meta.data, keep.rownames="SpotID")
meta_st <- meta_st[!is.na(CellType) & CellType!="" & CellType!="unassigned"]
score_cols <- grep("^[NM]_SEP[0-9]+_score$", colnames(meta_st), value=TRUE)
score_cols

heatmap_full <- meta_st[, lapply(.SD, function(x) {
    if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
}), by = CellType, .SDcols = score_cols]

for (sample_id in unique(meta_st$SampleID)) {
    
    low_cov_mps <- coverage_all[SampleID == sample_id & !analysis_eligible,MP]
    idx_sample <- which(meta_st$SampleID == sample_id)
    
    for (mp in low_cov_mps) {
        score_col <- paste0(mp, "_score")
        if (score_col %in% colnames(meta_st)) {
            set(meta_st, i = idx_sample, j = score_col, value = NA_real_)
        }
    }
}
heatmap_celltype <- meta_st[,lapply(.SD,function(x) {
    if (all(!is.finite(x))) {
        NA_real_
    } else {
        mean(x[is.finite(x)])
    }}),
    by = CellType, .SDcols = score_cols]

heatmap_long <- melt(heatmap_celltype, id.vars="CellType", measure.vars=score_cols, variable.name="SEP", value.name="Score")
heatmap_long[, SEP:=sub("_score$", "", SEP)]
heatmap_long[!is.finite(Score), Score:=NA_real_]
heatmap_long[, z_score:={
    valid <- is.finite(Score)
    ans <- rep(NA_real_, .N)
    if(sum(valid)>=2){
        s <- sd(Score[valid])
        if(is.finite(s) && s>0) ans[valid] <- (Score[valid]-mean(Score[valid]))/s else ans[valid] <- 0
    }
    ans
}, by=SEP]
heatmap_dt <- dcast(heatmap_long, SEP~CellType, value.var="z_score", fill=NA_real_)
heatmap_mat <- as.matrix(heatmap_dt[, -1])
rownames(heatmap_mat) <- heatmap_dt$SEP
heatmap_mat <- pmax(pmin(heatmap_mat, 2.5), -2.5)
celltype_order <- c("Bcell", "Tcell", "Myeloid", "Endothelial", "Fibroblast", "Hepatocyte", "Epithelial")
for (ct_missing in setdiff(celltype_order, colnames(heatmap_mat))) {
    extra <- matrix(NA_real_, nrow(heatmap_mat), 1L,
                    dimnames = list(rownames(heatmap_mat), ct_missing))
    heatmap_mat <- cbind(heatmap_mat, extra)
}
heatmap_mat <- heatmap_mat[, celltype_order, drop=FALSE]
colnames(heatmap_mat) <- c("B cells","T cells","Myeloid","Endothelial","Fibroblast","Hepatocytes","Epithelial")
sep_names <- rownames(heatmap_mat)
sep_groups <- factor(ifelse(grepl("^N_SEP", sep_names), "N-SEPs", "M-SEPs"),levels = c("N-SEPs", "M-SEPs"))
plot_mat <- t(heatmap_mat)
colnames(plot_mat) <-  gsub('_','-',colnames(plot_mat))
sep_group_colors <- c("N-SEPs"="#7BB6D8", "M-SEPs"="#B8A7C8")
top_annotation <- HeatmapAnnotation(SEPs=sep_groups,
                                    col=list(SEPs=sep_group_colors), 
                                    simple_anno_size=unit(4, "mm"),
                                    annotation_name_side="left", 
                                    annotation_name_gp=gpar(fontsize=9, fontface="bold"), show_legend=TRUE)
heatmap_col_fun <- colorRamp2(c(-2.5, -1.25, 0, 1.25, 2.5), c("#244A7C", "#8BB3CF", "#FFFFFF", "#E59A91", "#9E3038"))
ht_st <- Heatmap(plot_mat, 
                 name="Z-score", 
                 col=heatmap_col_fun, 
                 na_col="#EFEFEF", 
                 cluster_rows=FALSE,
                 column_split=sep_groups, 
                 # Correlation clustering is undefined for all-NA/constant SEP columns.
                 cluster_columns=all(is.finite(plot_mat)) &&
                     all(apply(plot_mat, 2, function(x) sd(x) > 0)), 
                 cluster_column_slices=FALSE,
                 clustering_distance_columns="pearson", 
                 clustering_method_columns="average", column_gap=unit(2, "mm"),
                 top_annotation=top_annotation, 
                 show_row_names=TRUE,
                 row_names_side="left", row_names_gp=gpar(fontsize=9),
                 show_column_names=TRUE, column_names_rot=90,
                 column_names_gp=gpar(fontsize=9), column_names_max_height=unit(40, "mm"),
                 column_title_gp=gpar(fontsize=10, fontface="bold"), rect_gp=gpar(col="white", lwd=0.8),
                 heatmap_legend_param=list(
                     title="Z-score", at=c(-2.5, -1, 0, 1, 2.5), legend_height=unit(35, "mm"),
                     title_gp=gpar(fontsize=10, fontface="bold"), labels_gp=gpar(fontsize=9)
                 ))

pdf(file = "Figure/ST/SEP/SEP_Heatmap/ST_heatmap_SEPs_celltypes.pdf", width=7, height=3.5, useDingbats=FALSE)
draw(ht_st, heatmap_legend_side="right", annotation_legend_side="right", merge_legends=FALSE, padding=unit(c(5, 5, 5, 5), "mm"))
dev.off()


# ---- assign by z-score -0.35 -------------------------
celltype_sep_candidates <- list(
    Bcell = paste0("N_SEP", c(2, 5, 9,10)),
    Tcell = paste0("N_SEP", c(2, 5, 9)),
    Myeloid = paste0("N_SEP", c(2, 10)),
    Fibroblast = paste0("N_SEP", c(1, 11, 12, 4, 6)),
    Endothelial = paste0("N_SEP", c(11, 4, 6, 8)),
    Hepatocyte = paste0("N_SEP", c(3, 7, 10)),
    Epithelial = paste0("M_SEP", 1:5)
)

assign_one_sample <- function(signatures_raw, MP_coverage, sample_id,
                              celltype_vec, type_vec, sep_candidates) {
    stopifnot(is.list(sep_candidates), !is.null(names(sep_candidates)),
              !anyDuplicated(names(sep_candidates)), !anyNA(names(sep_candidates)),
              all(nzchar(names(sep_candidates))),
              all(rownames(signatures_raw) %in% names(celltype_vec)),
              all(rownames(signatures_raw) %in% names(type_vec)),
              all(colnames(signatures_raw) %in% MP_coverage$MP))
    
    score_matrix <- as.matrix(signatures_raw)
    eligible_MPs <- MP_coverage$MP[MP_coverage$analysis_eligible]
    

    celltype_use <- celltype_vec[match(rownames(score_matrix),names(celltype_vec))]
    type_use <- type_vec[match(rownames(score_matrix), names(type_vec))]
    is_normal <- type_use %in% "Normal"
    celltype_use[is_normal & celltype_use %in% "Epithelial"] <- "unassigned"
    
    N_SEPs <- grep("^N_SEP",colnames(score_matrix),value = TRUE)
    
    top1_MP <- rep(NA_character_,nrow(score_matrix))
    top1_score <- rep(NA_real_,nrow(score_matrix))
    top2_MP <- rep(NA_character_,nrow(score_matrix))
    top2_score <- rep(NA_real_,nrow(score_matrix))
    MP_margin <- rep(NA_real_,nrow(score_matrix))
    
    assign_sep <- function(idx,MPs){
        MPs <- colnames(score_matrix)[colnames(score_matrix) %in% MPs & colnames(score_matrix) %in% eligible_MPs]
        
        if(length(idx) == 0 || length(MPs) == 0) return(NULL)
        
        score_sub <- score_matrix[idx,MPs,drop = FALSE]
        
        
        for (i in seq_along(idx)) {
            
            allowed <- if (is_normal[idx[i]]) {
                colnames(score_sub) %in% N_SEPs 
                } else {
                    rep(TRUE, ncol(score_sub))
                }
            
            valid <- which(allowed & is.finite(score_sub[i, ]))
            
            if (length(valid) == 0) next
            
            ranked <- valid[order(-score_sub[i, valid], valid)]
            
            spot_idx <- idx[i]
            
            top1_MP[spot_idx] <<- colnames(score_sub)[ranked[1]]
            top1_score[spot_idx] <<- score_sub[i, ranked[1]]
            
            if (length(ranked) >= 2) {
                top2_MP[spot_idx] <<- colnames(score_sub)[ranked[2]]
                top2_score[spot_idx] <<- score_sub[i, ranked[2]]
                MP_margin[spot_idx] <<- top1_score[spot_idx] - top2_score[spot_idx]
            }
        }
    }
    
    for (celltype in names(sep_candidates)) {
        assign_sep(which(celltype_use == celltype), sep_candidates[[celltype]])
    }

    idx_unassigned <- which(celltype_use == "unassigned")
    assign_sep(idx_unassigned,colnames(score_matrix))

    idx_other <- which(!(celltype_use %in% c(names(sep_candidates), "unassigned")))
    
    if (length(idx_other) > 0) {
        warning("Sample ", sample_id, ": using all N-SEPs for unmapped CellType: ",
                paste(unique(ifelse(is.na(celltype_use[idx_other]), "<NA>",
                                    celltype_use[idx_other])), collapse = ", "))
    }
    assign_sep(idx_other,N_SEPs)
    
    assignment <- data.frame(SpotID = rownames(score_matrix),
                             sample = sample_id,
                             dominant_MP = top1_MP,
                             dominant_MP_score = top1_score,
                             second_MP = top2_MP,
                             second_MP_score = top2_score,
                             MP_margin = MP_margin,
                             stringsAsFactors = FALSE
    ) 
    assignment
}

 
meta <- obj_hex[[]]
score_cols <- paste0(names(MP_list), "_score")
stopifnot(all(score_cols %in% names(meta)), !anyDuplicated(rownames(meta)),
          !anyDuplicated(coverage_all[, .(SampleID, MP)]))

sample_names <- unique(as.character(meta$SampleID))
assignments <- setNames(vector("list", length(sample_names)), sample_names)
    

for (sample_id in sample_names) {
    md <- meta[as.character(meta$SampleID) == sample_id, , drop = FALSE]
    scores <- as.matrix(md[, score_cols, drop = FALSE])
    colnames(scores) <- names(MP_list)
    assignments[[sample_id]] <- assign_one_sample(
        signatures_raw = scores,
        MP_coverage = coverage_all[SampleID == sample_id],
        sample_id = sample_id,
        celltype_vec = setNames(as.character(md$CellType), rownames(md)),
        type_vec = setNames(as.character(md$Type), rownames(md)),
        sep_candidates = celltype_sep_candidates
    )
}
    
assignment_all <- rbindlist(assignments)
stopifnot(nrow(assignment_all) == ncol(obj_hex), !anyDuplicated(assignment_all$SpotID))
idx <- match(colnames(obj_hex), assignment_all$SpotID)
obj_hex$SEP <- assignment_all$dominant_MP[idx]
obj_hex$SEP_score <- assignment_all$dominant_MP_score[idx]
obj_hex$SEP_second <- assignment_all$second_MP[idx]
obj_hex$SEP_second_score <- assignment_all$second_MP_score[idx]
obj_hex$SEP_margin <- assignment_all$MP_margin[idx]

saveRDS(obj_hex,  "Data/ST/SEP/All_hex_SEP_CelltypeConstrained.RDS")

meta_export <- as.data.table(obj_hex[[]], keep.rownames = "SpotID")

export_cols <- intersect(c(
    "SpotID", "sample_id", "SampleID", "PatientID", "Type", "Region",
    "CellType", "CellType_original", "SEP", "SEP_score",
    "SEP_second", "SEP_second_score", "SEP_margin",
    paste0(names(MP_list), "_score"), paste0(names(MP_list), "_weighted_score")
), names(meta_export))

fwrite(meta_export[, ..export_cols], "Data/ST/SEP/All_hex_MP_scores_for_h5ad_0924.csv")
 
candidate_audit <- rbindlist(lapply(names(celltype_sep_candidates), function(ct) {
    genes <- celltype_sep_candidates[[ct]]
    data.table(CellType = ct, SEP = if (length(genes)) genes else NA_character_)
}))

fwrite(candidate_audit, "Data/ST/SEP/SEP_CellType_candidates.csv")

sep_colors <- c("N-SEP1" = "#d7b0b0",
                "N-SEP12" = "#bc5c69",
                "N-SEP11" = "#0496ff",
                "N-SEP4" = "#ff7944",
                "N-SEP8" = "#ff78b3",
                "N-SEP6" = "#a8dadc",
                "N-SEP2" = "#ffbe0b",
                "N-SEP5" = "#ba55d3",
                "N-SEP9" = "#e63946",
                "N-SEP10" = "#2c6e49",
                "N-SEP7" = "#ace894",
                "N-SEP3" = "#78a978",
                "M-SEP1" = "#9fc1da",
                "M-SEP2" = "#FFFC8C",
                "M-SEP3" = "#F99FFF",
                "M-SEP4" = "#F9D2BA",
                "M-SEP5" = "#70FFD2")

sample_levels <- if (is.factor(obj_hex$SampleID)){
    levels(obj_hex$SampleID)
 } else{
    unique(as.character(obj_hex$SampleID))
 }
    
M_order <- paste0("M-SEP", 1:5)
N_order <- paste0("N-SEP", c(2,5,9,1,12,11,4,6,8,3,7,10))
SEP_order_all <- c(M_order, N_order)

meta = obj_hex@meta.data
sep_meta <- as.data.table(copy(meta))
sep_meta[, Category := fifelse(grepl("^N_SEP", SEP), "N-SEPs", fifelse(grepl("^M_SEP", SEP),"M-SEPs",NA_character_))]
sep_meta[, SEP := gsub("_", "-", SEP)]

# Composition plots summarize assigned spots.
sep_meta <- sep_meta[!is.na(SEP)]
sep_meta[, SEP := factor(SEP, levels=SEP_order_all)]

sep_composition <- sep_meta[!is.na(SEP),.N,by = .(SampleID, SEP)]
sep_composition[, Sample_plot := factor(SampleID, levels=sample_levels)]

p_sep_count <- ggplot(sep_composition, aes(Sample_plot, N, fill = SEP)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.15) +
    scale_fill_manual(values = sep_colors, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_few() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),legend.position = "bottom") +
    labs(y = "Number of spots", fill = "SEP")

p_sep_count
ggsave("Figure/ST/SEP/domSEP/01_sep_composition_counts.pdf",p_sep_count,width = 10,height = 4)

# ---- Figure 2: N-SEP and M-SEP composition ---------------------------------

sep_direction <- sep_meta[!is.na(SEP), .N,by = .(SampleID, Category, SEP)]
sep_direction[, SEP := factor(SEP, levels=c(N_order, M_order))]

sep_direction[, Value_abs := fifelse(Category=="M-SEPs", -N, N)]
sep_direction[, Sample_total := sum(N), by=SampleID]
sep_direction[, Proportion := N/Sample_total]
sep_direction[, Value_rel := fifelse(Category=="M-SEPs", -Proportion, Proportion)]
sep_direction[, Sample_plot := factor(SampleID, levels=sample_levels)]

p_sep_direction_abs <- ggplot(sep_direction, aes(Sample_plot, Value_abs, fill = SEP)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = sep_colors, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_few() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "bottom"
    ) +
    labs(y = "Number of spots", fill = "SEP")

p_sep_direction_abs
ggsave("Figure/ST/SEP/domSEP/02_sep_composition_counts_n_M.pdf",p_sep_direction_abs, width = 10,height = 4)


# ---- Figure 3 relative: order samples by malignant proportion --------------

sep_direction_relative <- sep_meta[!is.na(SEP), .N, by=.(SampleID, Type, Category, SEP)]
sep_direction_relative[, SEP := factor(SEP, levels=c(N_order, M_order))]

sep_direction_relative[, Sample_total := sum(N), by=SampleID]
sep_direction_relative[, Proportion := N/Sample_total]
sep_direction_relative[, Value_rel := fifelse(Category=="M-SEPs", -Proportion, Proportion)]

# Tumor first, Normal second; keep original SampleID order within each group
sample_info <- unique(sep_direction_relative[, .(SampleID, Type)])
sample_info[, original_order := match(as.character(SampleID), sample_levels)]

sample_order <- c(
    sample_info[Type=="Tumor"][order(original_order), as.character(SampleID)],
    sample_info[Type=="Normal"][order(original_order), as.character(SampleID)]
)

sep_direction_relative[, Sample_plot_rel := factor(SampleID, levels=sample_order)]

legend_order <- c(paste0("N-SEP", 1:12), paste0("M-SEP", 1:5))

p_sep <- ggplot(sep_direction_relative, aes(Sample_plot_rel, Value_rel, fill = SEP)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    annotate("text", x = Inf, y = 1.05, label = "Non-malignant", size = 4, hjust = 1.1, fontface = "bold") +
    annotate("text", x = Inf, y = -0.5, label = "Malignant", size = 4, hjust = 1.1, fontface = "bold") +
    scale_y_continuous(
        labels = function(x) percent(abs(x), accuracy = 1),
        limits = c(-0.75, 1.1),
        breaks = c(-0.5, 0, 0.5, 1),
        expand = expansion(mult = c(0.01, 0.01))) +
    scale_fill_manual(values = sep_colors, breaks = legend_order, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_few() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank(),
          legend.position = "bottom") +
    labs(y = "Relative proportion of spots", fill = "SEP")

p_sep
ggsave( "Figure/ST/SEP/domSEP/03_sep_direction_relative.pdf",p_sep, width = 10,height = 4)


sep_meta_tumor = subset(sep_meta, Type == 'Tumor')

sep_direction <- sep_meta_tumor[!is.na(SEP), .N,by = .(SampleID, Category, SEP)]
sep_direction[, SEP := factor(SEP, levels=c(N_order, M_order))]

sep_direction[, Value_abs := fifelse(Category=="M-SEPs", -N, N)]
sep_direction[, Sample_total := sum(N), by=SampleID]
sep_direction[, Proportion := N/Sample_total]
sep_direction[, Value_rel := fifelse(Category=="M-SEPs", -Proportion, Proportion)]
sep_direction[, Sample_plot := factor(SampleID, levels=sample_levels)]

p_sep_direction_abs <- ggplot(sep_direction, aes(Sample_plot, Value_abs, fill = SEP)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_y_continuous(labels = abs) +
    scale_fill_manual(values = sep_colors, drop = FALSE) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_few() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "bottom"
    ) +
    labs(y = "Number of spots", fill = "SEP")

p_sep_direction_abs
ggsave("Figure/ST/SEP/domSEP/04_sep_composition_counts_tumor_tissue.pdf",p_sep_direction_abs, width = 9,height = 4.2)

 
coverage_file <- file.path("./SEP_gene_coverage.csv")
score_file <- file.path("./All_hex_MP_scores_only_for_h5ad.csv")

sep_levels <- c(paste0("N_SEP", 1:12), paste0("M_SEP", 1:5))
eligibility_threshold <- 0.3
weight_saturation <- 0.3  
points_per_sample_sep <- 100L   
set.seed(42)
coverage <- fread(coverage_file)
coverage[, SEP_type := fifelse(grepl("^N_", MP), "N-SEP", "M-SEP")]
coverage[, SEP := factor(MP, levels = sep_levels)]
S6A_summary <- coverage[, .(
    n_samples = .N, median_coverage = median(fraction_present),
    min_coverage = min(fraction_present), n_below_threshold = sum(!analysis_eligible),
    fraction_below_threshold = mean(!analysis_eligible)
), by = .(SEP_type, SEP)]

options(repr.plot.width=7,repr.plot.height=4.5)

p_S6A <- ggplot(coverage, aes(SEP, fraction_present, fill = SEP_type)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.45, linewidth = 0.25) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    geom_point(aes(color = SEP_type),position = position_jitter(width = 0.12, seed = 42), size = 0.6, alpha = 1) +
    geom_hline(yintercept = eligibility_threshold, linetype = 2, colour = "#B2182B") +
    scale_fill_manual(values = c("N-SEP" = "#7BB6D8", "M-SEP" = "#B8A7C8")) +
    scale_color_manual(values = c("N-SEP" = "#005f96", "M-SEP" = "#6c488d")) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Signature gene coverage", fill = NULL) +
    theme_few() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

library(data.table)
library(ggplot2)
library(SeuratObject)
library(ggthemes)

setwd("")
metadata_file <- "Data/ST/SEP/All_hex_SEP_CelltypeConstrained.RDS"
sep_sets <- list(`N-SEP` = paste0("N_SEP", 1:12), `M-SEP` = paste0("M_SEP", 1:5))
obj <- readRDS(metadata_file)
meta <- as.data.table(obj@meta.data)[, .(SampleID, PatientID, Type, SEP)]

for (nm in names(meta)) set(meta, j = nm, value = as.character(meta[[nm]]))

meta <- meta[Type == "Tumor" & SEP %in% unlist(sep_sets)]
meta[, Compartment := fifelse(SEP %in% sep_sets[["N-SEP"]], "N-SEP", "M-SEP")]

counts <- meta[, .(n_hex = .N), by = .(SampleID, PatientID, Compartment, SEP)]
totals <- dcast(counts[, .(n_hex = sum(n_hex)), by = .(SampleID, PatientID, Compartment)],
                SampleID + PatientID ~ Compartment, value.var = "n_hex", fill = 0)
samples <- totals[`N-SEP` > 0 & `M-SEP` > 0]
samples <- samples[PatientID %in% samples[, .N, by = PatientID][N >= 2, PatientID]]

composition <- rbindlist(lapply(names(sep_sets), function(comp) {
    d <- merge(CJ(SampleID = samples$SampleID, SEP = sep_sets[[comp]]),
               counts[Compartment == comp, .(SampleID, SEP, n_hex)],
               by = c("SampleID", "SEP"), all.x = TRUE)
    d[is.na(n_hex), n_hex := 0L]
    d[, proportion := n_hex / sum(n_hex), by = SampleID]
    d[, Compartment := comp]
    merge(d, samples[, .(SampleID, PatientID)], by = "SampleID")
}))

js_distance <- function(p, q) {
    m <- (p + q) / 2
    kl <- function(x) { k <- x > 0; sum(x[k] * log2(x[k] / m[k])) }
    sqrt(max(0, (kl(p) + kl(q)) / 2))
}
parts <- list()
for (comp in names(sep_sets)) {
    wide <- dcast(composition[Compartment == comp], SampleID ~ SEP, value.var = "proportion")
    mat <- as.matrix(wide[, sep_sets[[comp]], with = FALSE])
    rownames(mat) <- wide$SampleID
    for (pat in unique(samples$PatientID)) {
        ids <- sort(samples[PatientID == pat, SampleID])
        ij <- combn(ids, 2)
        parts[[length(parts) + 1L]] <- data.table(
            PatientID = pat, Sample1 = ij[1, ], Sample2 = ij[2, ], Compartment = comp,
            Distance = apply(ij, 2, function(z) js_distance(mat[z[1], ], mat[z[2], ])))
    }
}
pairs <- rbindlist(parts)
patient_long <- pairs[, .(Mean_distance = mean(Distance)), by = .(PatientID, Compartment)]
patient <- dcast(patient_long, PatientID ~ Compartment, value.var = "Mean_distance")
patient[, Delta := `N-SEP` - `M-SEP`]
exact_ok <- !any(patient$Delta == 0) && !anyDuplicated(abs(patient$Delta))
p_value <- if (all(patient$Delta == 0)) 1 else wilcox.test(
    patient$`N-SEP`, patient$`M-SEP`, paired = TRUE, alternative = "two.sided", exact = exact_ok)$p.value
statistics <- data.table(n_patients = nrow(patient), n_regions = nrow(samples),
    n_region_pairs = nrow(pairs) / 2, mean_N = mean(patient$`N-SEP`),
    mean_M = mean(patient$`M-SEP`), mean_difference = mean(patient$Delta),
    n_N_greater = sum(patient$Delta > 0), paired_Wilcoxon_P = p_value)


colors <- c(`N-SEP` = "#237C91", `M-SEP` = "#D46A43")
pairs[, Weight := 1 / .N, by = .(PatientID, Compartment)]


grid_x <- seq(0, 1, length.out = 501)
density_data <- pairs[, {
    w <- Weight / sum(Weight)
    y <- vapply(grid_x, function(z) sum(w * (
        dnorm((z - Distance) / .06) + dnorm((z + Distance) / .06) +
        dnorm((z - (2 - Distance)) / .06))) / .06, numeric(1))
    area <- sum(diff(grid_x) * (head(y, -1) + tail(y, -1)) / 2)
    .(Distance = grid_x, Density = y / area)
}, by = Compartment]


p_density <- ggplot(density_data, aes(Distance, Density,color = Compartment, fill = Compartment)) +
    geom_area(position = "identity", alpha = .15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = colors, limits = names(colors)) +
    scale_fill_manual(values = colors, limits = names(colors)) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
        x = "Within-patient Jensen-Shannon distance",
        y = "Density",
        subtitle = sprintf("Paired Wilcoxon P = %.4g", p_value)
    ) +
    theme_few() +
    theme(
        legend.title = element_blank(),
        legend.position = "top"
    )

p_density
