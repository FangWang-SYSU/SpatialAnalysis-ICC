
library(ggplot2)
library(mixtools)

#------------Figure5A
binwidth <- 0.12
d_all <- read.csv("Table/ST/SEP_res/Sample_R_NR_classification.csv")
d <- d_all[d_all$QC_status == "Pass", , drop = FALSE]
x <- d$M_resistance_balance
fit_equal_variance_gmm <- function(x, n_starts = 40L, seed = 42) {
    set.seed(seed)
    fits <- list()
    for (i in seq_len(n_starts)) {
        mu0 <- if (i == 1L) as.numeric(quantile(x, c(0.25, 0.75))) else sort(sample(x, 2L))
        sigma0 <- sd(x) * runif(1L, 0.4, 1.2)
        w0 <- if (i == 1L) 0.5 else runif(1L, 0.2, 0.8)
        fit <- NULL
        invisible(capture.output(
            fit <- suppressWarnings(tryCatch(
                normalmixEM(
                    x, lambda = c(w0, 1 - w0), mu = mu0, sigma = sigma0,
                    k = 2L, arbvar = FALSE, epsilon = 1e-9,
                    maxit = 5000L, maxrestarts = 5L, verb = FALSE
                ),
                error = function(e) NULL
            ))
        ))
        if (!is.null(fit) && is.finite(fit$loglik) &&
            min(fit$lambda) > 1e-4 && min(fit$sigma) > 1e-4) {
            fits[[length(fits) + 1L]] <- fit
        }
    }
    
    fit <- fits[[which.max(vapply(fits, function(z) z$loglik, numeric(1)))]]
    component_order <- order(fit$mu)
    fit$mu <- fit$mu[component_order]
    fit$lambda <- fit$lambda[component_order]
    fit$sigma <- rep(fit$sigma[1], 2L)
    fit
}

fit <- fit_equal_variance_gmm(x)
mu <- fit$mu
sigma <- fit$sigma[1]
weight <- fit$lambda

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Asia/Shanghai")

plot_xlim <- range(x) + c(-0.55, 0.55)
grid <- seq(plot_xlim[1], plot_xlim[2], length.out = 5000L)
hist_breaks <- seq(floor(min(x) / binwidth) * binwidth,
                   ceiling(max(x) / binwidth) * binwidth, by = binwidth)

weighted_r <- weight[1] * dnorm(grid, mu[1], sigma)
weighted_nr <- weight[2] * dnorm(grid, mu[2], sigma)

curves <- rbind(
    data.frame(score = grid, density = weighted_r, curve = "R-like component"),
    data.frame(score = grid, density = weighted_nr, curve = "NR-like component")
)

total <- data.frame(score = grid, density = weighted_r + weighted_nr)

log_r <- log(weight[1]) + dnorm(x, mu[1], sigma, log = TRUE)
log_nr <- log(weight[2]) + dnorm(x, mu[2], sigma, log = TRUE)
p_nr <- plogis(log_nr - log_r)
p_r <- 1 - p_nr

posterior_boundary <- function(probability_nr) {
    (qlogis(probability_nr) - log(weight[2] / weight[1]) +
         (mu[2]^2 - mu[1]^2) / (2 * sigma^2)) * sigma^2 / (mu[2] - mu[1])
}

classify <- function(confidence) {
    lower <- posterior_boundary(1 - confidence)
    upper <- posterior_boundary(confidence)
    list(
        lower = lower,
        upper = upper,
        class = ifelse(x < lower, "R-like", ifelse(x > upper, "NR-like", "Ambiguous"))
    )
}

result90 <- classify(0.90)
classification <- data.frame(SampleID = d$SampleID, PatientID = d$PatientID, M_resistance_balance = x,
                             p_R = p_r,
                             p_NR = p_nr,
                             Group = result90$class,
                             stringsAsFactors = FALSE)
 

component_colors <- c("R-like component" = "#FDBD3F", 
                      "NR-like component" = "#56B4CC",
                      "Mixture density" = "#222222")

class_colors <- c("R-like" = "#D59627", "Ambiguous" = "#a7a4ac", "NR-like" = "#459BB5")

make_plot <- function(confidence, classified) {
    
    q <- d
    q$class <- factor(classified$class, levels = c("R-like", "Ambiguous", "NR-like"))
    counts <- table(q$class)
    ymax <- 1.25 * max(hist(x, breaks = hist_breaks, plot = FALSE)$density, total$density)
    
    ggplot(q, aes(M_resistance_balance)) +
        annotate("rect", xmin = classified$lower, xmax = classified$upper, ymin = -Inf, ymax = Inf, fill = "#cfcfcf", alpha = 0.12) +
        geom_histogram(aes(y = after_stat(density)), binwidth = binwidth, boundary = 0, fill = "#E8EBED", color = "white", linewidth = 0.35) +
        geom_area(data = curves, aes(score, density, fill = curve), inherit.aes = FALSE, position = "identity", alpha = 0.23) +
        geom_line(data = curves, aes(score, density, color = curve),inherit.aes = FALSE, linewidth = 0.9) +
        geom_line(data = total, aes(score, density, color = "Mixture density"),inherit.aes = FALSE, linewidth = 1) +
        
        geom_vline(xintercept = c(classified$lower, classified$upper),
                   linetype = "dashed", color = "#665A73", linewidth = 0.65) +
        geom_rug(aes(color = class), sides = "b", linewidth = 0.55, show.legend = FALSE) +
        
        annotate("text", x = mu[1], y = ymax * 0.91,
                 label = sprintf("R-like\n(n = %d)", counts[["R-like"]]),
                 color = "#A36B00", size = 3.8) +
        annotate("text", x = mean(c(classified$lower, classified$upper)), y = ymax * 0.77,
                 label = sprintf("Ambiguous\n(n = %d)", counts[["Ambiguous"]]),
                 color = "#665A73", size = 3.8) +
        annotate("text", x = mu[2], y = ymax * 0.91,
                 label = sprintf("NR-like\n(n = %d)", counts[["NR-like"]]),
                 color = "#277C90", size = 3.8) +
        annotate("text", x = classified$lower - 0.03, y = ymax * 0.58,
                 label = sprintf("%.3f", classified$lower), hjust = 1,
                 color = "#665A73", size = 3.2) +
        annotate("text", x = classified$upper + 0.03, y = ymax * 0.58,
                 label = sprintf("%.3f", classified$upper), hjust = 0,
                 color = "#665A73", size = 3.2) +
        scale_color_manual(values = c(component_colors, class_colors),
                           breaks = names(component_colors)) +
        scale_fill_manual(values = component_colors, guide = "none") +
        coord_cartesian(xlim = plot_xlim, ylim = c(0, ymax)) +
        
        labs(
            title = sprintf("Two-component mixture with a %.0f%% posterior-confidence region",
                            confidence * 100),
            subtitle = "Samples in the overlap region are retained as Ambiguous and excluded only from strict two-group comparisons",
            x = "M-resistance balance", y = "Probability density", color = NULL,
            caption = sprintf("Histogram bin width = %.2f. Dashed lines: posterior-confidence boundaries (max posterior = %.2f).",
                              binwidth, confidence)
        ) +
        theme_classic(base_size = 11) +
        theme(
            legend.position = "bottom",
            legend.text = element_text(size = 9),
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 10),
            plot.caption = element_text(hjust = 0, size = 8),
            plot.margin = margin(12, 14, 12, 14)
        )
}
 
plot90 <- make_plot(0.90, result90)
plot90


##-Figure5GH
library(SeuratObject)
library(Matrix)
library(data.table)
library(scalop)
library(ggplot2)
library(ggthemes)
seurat_obj <- readRDS("./Data/ST/Stromal_seu.RDS")
raw_counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")

sep_gene_path <- "./Table/all_MP_list_New.csv"
sep_gene_table <- fread(sep_gene_path, header = TRUE, colClasses = "character")

raw_counts <- drop0(as(raw_counts, "dgCMatrix"))
metadata <- as.data.table(seurat_obj[[]], keep.rownames = "HexID")

metadata[, nCount_RNA := as.numeric(Matrix::colSums(raw_counts))]
metadata[, nFeature_RNA := diff(raw_counts@p)]
metadata <- metadata[nCount_RNA >= 200 & nCount_RNA <= 4000 & nFeature_RNA >= 100 & nFeature_RNA <= 2000 & percent.mt <= 30]

setorder(metadata, SampleID, HexID)
raw_counts <- raw_counts[, metadata$HexID, drop = FALSE]
sample_groups <- unique(metadata[, .(SampleID, PatientID, Group)])
cohort_summary <- metadata[, .(N_hex = .N, N_samples = uniqueN(SampleID), N_patients = uniqueN(PatientID)), by = Group]

downsample_cap <- function(counts, replicate_id,sample_ids, sample_columns,target_umi = 500L) {
    sampled_values <- counts@x
    original_depth <- as.integer(Matrix::colSums(counts))
    for (sample_index in seq_along(sample_ids)) {
        set.seed(42L + 1000L * replicate_id + sample_index)
        for (column_index in sample_columns[[sample_ids[sample_index]]]) {
            if (original_depth[column_index] <= target_umi) next
            value_index <- seq.int(counts@p[column_index] + 1L, counts@p[column_index + 1L])
            molecule_id <- sample.int(original_depth[column_index], target_umi, replace = FALSE)
            gene_id <- findInterval(molecule_id - 1L, cumsum(counts@x[value_index])) + 1L
            sampled_values[value_index] <- as.numeric(tabulate(gene_id, nbins = length(value_index)))
        }
    }
    counts@x <- sampled_values
    drop0(counts)
}

normalize_counts <- function(counts,keep_gene) {
    expression <- counts[keep_gene, , drop = FALSE]
    library_size <- as.numeric(Matrix::colSums(expression))
    stopifnot(all(library_size > 0))
    expression@x <- log2(1 + expression@x * rep.int(1e5 / library_size, diff(expression@p)))
    expression
}

make_reference_weights <- function(mu, targets, controls) {
    target_weight <- as.numeric(names(mu) %in% targets) / length(targets)
    control_weight <- tabulate(match(controls, names(mu)), nbins = length(mu)) / length(controls)
    list(weight = setNames(target_weight - control_weight, names(mu)),
         baseline = mean(mu[targets]) - mean(mu[controls]),
         target_weight = target_weight, control_weight = control_weight)
}


make_legacy_reference <- function(expression, rows, targets, seed) {
    indicator <- Matrix::sparseMatrix(i = rows, j = rep(1L, length(rows)), x = 1 / length(rows),
                                      dims = c(ncol(expression), 1L))
    mu <- setNames(as.numeric(expression %*% indicator), rownames(expression))
    squared <- expression
    squared@x <- squared@x^2
    variance <- as.numeric(squared %*% indicator) - mu^2
    universe <- names(mu)[(mu > 0.1 & variance > 0) | names(mu) %in% targets]
    n_bins <- min(30L, max(1L, floor(length(universe) / 120L)))
    bins <- scalop::bin(x = sort(mu[universe]), breaks = n_bins)
    set.seed(seed)
    controls <- scalop::binmatch(targets, bins = bins, n = 100L, replace = FALSE)
    make_reference_weights(mu, targets, controls)
}

plot_bins <- c(-20L, 20L, 40L, 60L, 80L, 100L)

region_hex_map <- as.data.table(seurat_obj@misc$spatial_boundary$hex_region_summary)[HexID %chin% metadata$HexID & distance %in% plot_bins,
                                                                                     .(HexID, SampleName, RegionID, boundary_distance = distance)]
region_hex_map <- merge(region_hex_map, metadata[, .(HexID, SampleID, PatientID, Group)], by = "HexID", all.x = TRUE, sort = FALSE)
map_index <- match(region_hex_map$HexID, metadata$HexID)
sep_names <- c("N_SEP1", "N_SEP11", "N_SEP12", "N_SEP4", "N_SEP6", "N_SEP8")

sep_genes <- lapply(sep_gene_table[, ..sep_names], trimws)
keep_gene <- !grepl("^MT-|^RPL|^RPS", rownames(raw_counts))

all_targets <- unlist(sep_genes, use.names = FALSE)

target_indicator <- sparseMatrix(
    i = match(all_targets, rownames(raw_counts)),
    j = rep(seq_along(sep_genes), lengths(sep_genes)),
    x = 1, dims = c(nrow(raw_counts), length(sep_genes)),
    dimnames = list(rownames(raw_counts), sep_names)
)
sample_ids <- sort(unique(metadata$SampleID))
sample_columns <- split(seq_len(nrow(metadata)), metadata$SampleID)

first_counts <- downsample_cap(counts = raw_counts, replicate_id = 1L, sample_ids = sample_ids, sample_columns = sample_columns, target_umi = 500)
detected <- as.matrix(crossprod(first_counts > 0, target_indicator))
first_expression <- normalize_counts(counts = first_counts, keep_gene = keep_gene)

reference_rows <- lapply(sep_names, function(sep) which(detected[, sep] >= 5))

references <- setNames(lapply(seq_along(sep_names), function(i) {
    make_legacy_reference(first_expression, reference_rows[[i]], sep_genes[[i]], 42 + 100000L + i)
}), sep_names)

weights <- Matrix::Matrix(do.call(cbind, lapply(references, `[[`, "weight")), sparse = TRUE)
centers <- vapply(references, `[[`, numeric(1), "baseline")

print(data.table(SEP = sep_names, N_reference_hex = lengths(reference_rows), Center = unname(centers)))


needed <- sort(unique(map_index))
region_hex_map[, N_hex := .N, by = .(SampleID, RegionID, boundary_distance)]
region_hex_map[, `:=`(N_regions = uniqueN(RegionID), Bin = .GRP), by = .(SampleID, boundary_distance)]
bins <- unique(region_hex_map[, .(Bin, SampleID, PatientID, Group, Distance = boundary_distance)])

membership <- sparseMatrix(i = match(map_index, needed), j = region_hex_map$Bin,
                           x = 1 / (region_hex_map$N_hex * region_hex_map$N_regions), dims = c(length(needed), nrow(bins)))

bins[, Bin := NULL]

first_expression <- normalize_counts(counts = first_counts, keep_gene = keep_gene)
first_expression_needed <- first_expression[, needed, drop = FALSE]

score_draw <- function(rep_id) {
    expression <- if (rep_id == 1L) {
        first_expression_needed
    } else {
        capped <- downsample_cap(counts = raw_counts, replicate_id = rep_id,sample_ids = sample_ids, sample_columns = sample_columns, target_umi = 500L)
        normalize_counts(counts = capped[, needed, drop = FALSE], keep_gene = keep_gene)
    }
    hex_scores <- sweep(as.matrix(crossprod(expression, weights)), 2L, centers, "-")
    sample_values <- as.data.table(as.matrix(crossprod(membership, hex_scores)))
    
    melt(cbind(bins, sample_values), measure.vars = sep_names, variable.name = "SEP", value.name = "Score", variable.factor = FALSE)
}

replicates <- seq_len(20L)
n_workers <- 6L

draws <- if (.Platform$OS.type == "windows" || n_workers <= 1L) {
    lapply(replicates, score_draw)
} else {
    parallel::mclapply(replicates, score_draw, mc.cores = n_workers,
                       mc.preschedule = FALSE, mc.set.seed = FALSE)
}

sample_scores <- rbindlist(draws)[, .(Score = mean(Score)), by = .(SEP, Group, PatientID, SampleID, Distance)]
patient_scores <- sample_scores[, .(Score = mean(Score)), by = .(SEP, Group, PatientID, Distance)]
print(patient_scores[, .(N_patient_distance_scores = .N, N_patients = uniqueN(PatientID)), by = .(SEP, Group)])

patient_scores[, SEP := sub("^N_", "", SEP)]
patient_scores[, Z_score := (Score - mean(Score)) / sd(Score), by = SEP]

common_shift <- -min(patient_scores$Z_score)
patient_scores[, Shifted_Z := Z_score + common_shift]

group_scores <- patient_scores[, .(Mean_shifted_Z = mean(Shifted_Z), Mean_Z = mean(Z_score), 
                                   Mean_SigScore = mean(Score), N_patients = .N), by = .(SEP, Group, Distance)]

group_scores[, Common_shift := common_shift]
seps <- c("SEP4", "SEP6", "SEP8", "SEP1", "SEP11", "SEP12")
groups <- c("NR-like", "R-like")

plot_data <- copy(group_scores)
plot_data[, `:=`(SEP = factor(SEP, levels = seps), Group = factor(Group, levels = groups))]
plot_data[, Signature_set := factor(fifelse(SEP %in% c("SEP4", "SEP6", "SEP8"),
                                            "Endothelial", "Fibroblast"),
                                    levels = c("Endothelial", "Fibroblast"))]

sep_colors <- c(SEP1 = "#d7b0b0", SEP11 = "#0496ff", SEP12 = "#900b1d", SEP4 = "#EC6D3C", SEP6 = "#a8dadc", SEP8 = "#ff78b3")

plot_dataA = subset(plot_data,Signature_set == 'Fibroblast')
Fibro_plot <- ggplot(plot_dataA, aes(Distance, Mean_shifted_Z, color = SEP, linetype = Group, shape = Group, group = interaction(SEP, Group))) +
    geom_line(linewidth = 0.55) + 
    geom_point(size = 3,  fill = "white") +
    scale_x_continuous(breaks = plot_bins) + 
    scale_y_continuous(expand = expansion(mult = c(.06, .08))) +
    scale_color_manual(name = "SEP", values = sep_colors) +
    
    scale_linetype_manual(name = "Group", values = c("NR-like" = "solid", "R-like" = "dashed")) +
    scale_shape_manual(name = "Group", values = c("NR-like" = 16, "R-like" = 21)) +
    guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid", shape = 16)),
           linetype = guide_legend(order = 2), shape = guide_legend(order = 2)) +
    labs(title = "Fibroblast",
         x = "Distance to tumor boundary (µm)", 
         y = "Patient-equal mean shifted Z score") +
    theme_few() + 
    theme(legend.position = "right", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),  
          legend.key.width = grid::unit(1.2, "cm"),plot.caption = element_text(hjust = 0))

ggsave('/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC/Figure/ST/Boundary/Figure5G/Fibroblast.pdf',Fibro_plot,width = 6,height = 4.5)
library(ggthemes)
Endo_plot <- ggplot(plot_dataB, aes(Distance, Mean_shifted_Z, color = SEP, linetype = Group, shape = Group, group = interaction(SEP, Group))) +
    geom_line(linewidth = 0.55) + 
    geom_point(size = 3,  fill = "white") +
    scale_x_continuous(breaks = plot_bins) + 
    scale_y_continuous(expand = expansion(mult = c(.06, .08))) +
    scale_color_manual(name = "SEP", values = sep_colors) +
    
    scale_linetype_manual(name = "Group", values = c("NR-like" = "solid", "R-like" = "dashed")) +
    scale_shape_manual(name = "Group", values = c("NR-like" = 16, "R-like" = 21)) +
    guides(color = guide_legend(order = 1, override.aes = list(linetype = "solid", shape = 16)),
           linetype = guide_legend(order = 2), shape = guide_legend(order = 2)) +
    labs(title = "Endothelial",
         x = "Distance to tumor boundary (µm)", 
         y = "Patient-equal mean shifted Z score") +
    theme_few() + 
    theme(legend.position = "right", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),  
          legend.key.width = grid::unit(1.2, "cm"),plot.caption = element_text(hjust = 0))

ggsave('/Volumes/B/Projects/2025Project/SpatialAnalysis/Project_ICC/Figure/ST/Boundary/Figure5G/Endothelial.pdf',Endo_plot,width = 6,height = 4.5)

comparison_distances <- c(-20L, 20L, 40L, 60L, 80L, 100L)
paired_means <- group_scores[Distance %in% comparison_distances & Group %chin% groups,
                             .(SEP, Distance, Group = factor(Group, levels = groups), Mean_shifted_Z)]

curve_pairs <- dcast(paired_means, SEP + Distance ~ Group,
                     value.var = "Mean_shifted_Z", drop = c(TRUE, FALSE))

sep_curve_paired_t <- curve_pairs[, {
    complete <- is.finite(`R-like`) & is.finite(`NR-like`)
    result <- list(N_pairs = sum(complete), Mean_diff_R_minus_NR = NA_real_,
                   T_statistic = NA_real_, P_value = NA_real_)
    
    if (result$N_pairs == length(comparison_distances)) {
        result$Mean_diff_R_minus_NR <- mean(`R-like` - `NR-like`)
        fit <- tryCatch(stats::t.test(`R-like`, `NR-like`, paired = TRUE, alternative = "two.sided"),
                        error = identity)
        if (inherits(fit, "error")) {
            result$Status <- conditionMessage(fit)
        } else if (!is.finite(fit$p.value)) {
            result$Status <- "Undefined paired t-test"
        } else {
            result$T_statistic <- unname(fit$statistic)
            result$P_value <- fit$p.value
        }
    }
    result
}, by = SEP]
sep_curve_paired_t[, P_Holm := stats::p.adjust(P_value, method = "holm", n = length(seps))]
print(sep_curve_paired_t)

comparison_distances <- c(-20L, 20L, 40L, 60L, 80L)
paired_means <- group_scores[Distance %in% comparison_distances & Group %chin% groups,
                             .(SEP, Distance, Group = factor(Group, levels = groups), Mean_shifted_Z)]
stopifnot(setequal(paired_means$SEP, seps),
          !anyDuplicated(paired_means[, .(SEP, Distance, Group)]))
curve_pairs <- dcast(paired_means, SEP + Distance ~ Group,
                     value.var = "Mean_shifted_Z", drop = c(TRUE, FALSE))

sep_curve_paired_t_to_80 <- curve_pairs[, {
    complete <- is.finite(`R-like`) & is.finite(`NR-like`)
    result <- list(N_pairs = sum(complete), Mean_diff_R_minus_NR = NA_real_, T_statistic = NA_real_, P_value = NA_real_)
    if (result$N_pairs == length(comparison_distances)) {
        result$Mean_diff_R_minus_NR <- mean(`R-like` - `NR-like`)
        fit <- tryCatch(stats::t.test(`R-like`, `NR-like`, paired = TRUE, alternative = "two.sided"),
                        error = identity)
        if (inherits(fit, "error")) {
            result$Status <- conditionMessage(fit)
        } else if (!is.finite(fit$p.value)) {
            result$Status <- "Undefined paired t-test"
        } else {
            result$T_statistic <- unname(fit$statistic)
            result$P_value <- fit$p.value
        }
    }
    result
}, by = SEP]

sep_curve_paired_t_to_80[, P_Holm := stats::p.adjust(P_value, method = "holm", n = length(seps))]
setorder(sep_curve_paired_t_to_80, SEP)
print(sep_curve_paired_t_to_80)
