suppressPackageStartupMessages({
    library(phangorn)
    library(ggtree)
    library(ggplot2)
    library(ape)
    library(patchwork)
})

run_infercnv_subcluster_tree <- function(obj,
                                         sample_id = "sample", 
                                         fig_dir = "Figure/SpatialFigure/1_Integration/InferCNV/Tree",
                                         min_cells = 10,
                                         make_tree_min_k = 2,
                                         final_clone = NULL,
                                         cb_pattle,
                                         width = 6, height = 6,
                                         stop_on_inconsistent = TRUE,
                                         verbose = TRUE
) {
    get_tip_order <- function(tr) {
        dd <- ggtree(tr)$data
        tip <- dd[dd$isTip, c("label", "y")]
        tip <- tip[order(-tip$y), , drop = FALSE]
        return(tip$label)
    }
    
    plot_tree <- function(tr, title, pal) {
        ggtree(tr, aes(color = factor(node)),size = 1.1) +
            geom_point2(aes(subset = isTip, fill = label), shape = 21, size = 8, stroke = 0) +
            geom_text2(aes(subset = isTip, label = label), size = 4, color = "black") +
            theme_void() +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14)) +
            ggtitle(title) +
            scale_fill_manual(values = pal, na.value = "grey80") +
            scale_color_manual(values = pal, na.value = "grey80")
    }
    
    mat <- obj@expr.data
    subc <- obj@tumor_subclusters[["subclusters"]][["all_observations"]]
    subc <- subc[lengths(subc) >= min_cells]
    clone <- names(subc)
    
    if (length(clone) == 0) {
        stop(sample_id, ": min_cells >= ", min_cells, " subcluster。")
    }

    rep_cells <- character(0)
    clone_used <- character(0)
    
    for (cl in clone) {
        cells <- intersect(names(subc[[cl]]), colnames(mat))
        if (length(cells) == 0) {
            if (verbose) message("worning", sample_id, " / ", cl, ": no cells to used, next")
            next
        }
        
        ref_cell <- cells[1]
        ref <- mat[, ref_cell]
        
        if (length(cells) >= 2) {
            for (j in 2:length(cells)) {
                cc <- cells[j]
                diff_n <- sum(mat[, cc] != ref)
                if (diff_n > 0) {
                    msg <- paste0(sample_id, " / ", cl,
                                  " differnet ", cc,
                                  "Ref: ", ref_cell, ",different counts", diff_n, "）")
                    if (stop_on_inconsistent) stop(msg) else warning(msg)
                    break
                }
            }
        }
        
        rep_cells <- c(rep_cells, ref_cell)
        clone_used <- c(clone_used, cl)
    }
    
    if (length(rep_cells) == 0) {
        stop("error ", sample_id, ": rep_cells is empty")
    }
    
    # -------------------------
    # 3) new_mat / CNV_mat
    #clone_named <- gsub("all_observations_s", "S", clone_used)
    clone_counts <- lengths(subc)
    
    sorted_clone <- names(sort(clone_counts, decreasing = TRUE))
    new_labels <- LETTERS[seq_along(sorted_clone)]
    name_map <- setNames(new_labels, sorted_clone)

    clone_named <- name_map[clone]
    
    new_mat <- mat[, rep_cells, drop = FALSE]
    colnames(new_mat) <- clone_named
    
    if(is.null(final_clone)){
        CNV_mat = new_mat
    }else{
        CNV_mat = new_mat[, final_clone, drop = FALSE]
    }
   
    CNV_mat <- CNV_mat * 2
    K <- ncol(CNV_mat)
 
    if (K <= make_tree_min_k) {
        if (verbose) message(sample_id, ": K=", K, " <= ", make_tree_min_k, "next")
        return(invisible(list()))
    }
 
    pal <- cb_pattle[clone_named]
    pal_ref <- c(pal, ref = "grey70")
 
    CNV_mat2 <- cbind(CNV_mat, ref = rep(2, nrow(CNV_mat)))
    dm.cell2 <- dist(t(CNV_mat2))
    treeNJ2 <- NJ(dm.cell2)
    
    cnh3_2 <- as.phyDat(t(CNV_mat2), type = "USER", levels = 0:max(CNV_mat2))
    PARtree2 <- pratchet(cnh3_2, minit = 10, maxit = 100, k = 5)
    PARtree2_1 <- acctran(PARtree2, cnh3_2)
    
    reroot_treeNJ <- root(treeNJ2, outgroup = "ref", resolve.root = TRUE)
    reroot_treeMP <- root(PARtree2_1, outgroup = "ref", resolve.root = TRUE)
    
    p3 <- plot_tree(reroot_treeNJ, "rooted_NJ", pal_ref)
    p4 <- plot_tree(reroot_treeMP, "rooted_MP", pal_ref)
    
    if (!is.null(fig_dir) && nzchar(fig_dir)) {
        if (!dir.exists(fig_dir)) {
            dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
        }
        ggsave(
            file.path(fig_dir, paste0(sample_id, "_rooted_tree.pdf")),
            plot =  (p3 + p4),
            width = width,
            height = height
        )
    }
     
    # result_order$NJ_rooted <-  pal[get_tip_order(reroot_treeNJ)]
    # result_order$MP_rooted <-  pal[get_tip_order(reroot_treeMP)]
    # 
    # invisible(result_order)
}
