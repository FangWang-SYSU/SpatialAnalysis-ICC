#!/usr/bin/env Rscript

library(data.table)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: filter_tiles.R <target_file> <tiles_dir>\n")
  quit(status = 1)
}

target_file <- args[1]
tiles_dir <- args[2]

# Load barcodes
aa <- fread(target_file, header = FALSE)

# List .gz files
setwd(tiles_dir)
tile_files <- list.files(pattern = "\\.gz$")

# Match barcodes
tile_List <- mclapply(tile_files, function(x) {
  data <- tryCatch({
    fread(cmd = paste("zcat", shQuote(x)))
  }, error = function(e) data.table(cell_bc = character(0)))
  
  intersect(data$cell_bc, aa$V1)
}, mc.cores = 30)

names(tile_List) <- tile_files

# Print results with match > 2
for (i in seq_along(tile_List)) {
  matched <- tile_List[[i]]
  if (length(matched) > 5) {
    cat(sprintf("%-40s  %5d matched barcodes\n", names(tile_List)[i], length(matched)))
  }
}
