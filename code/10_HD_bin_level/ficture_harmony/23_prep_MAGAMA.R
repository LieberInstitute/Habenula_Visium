# MAGAMA on FICTURE cluster markers (k=17)
# prepare input data

# Find MeanRatio marker genes
library(SpatialExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(spatialLIBD)
library(duckplyr)

spe_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'spe','y_clean_spe.rds'
)
cluster_path = here(
  'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
  'bin_level_clusters_batch.parquet'
)
out_dir = here(
    'processed-data', '10_HD_bin_level', 'new_samples2','ficture_harmony', 'MAGMA', 'gene_sets'
)

mean_ratio_threshold = 1.05
max_num_genes = 200

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

cluster_df = read_parquet_duckdb(cluster_path) |>
  select(sample_id, barcode, FICTURE_k17) |>
  collect()

spe$ficture_cluster = tibble(sample_id = spe$sample_id, barcode = spe$barcode) |>
  left_join(cluster_df, by = c("sample_id", "barcode")) |>
  pull(FICTURE_k17)

dt <- as.data.table( fread(cluster_path) )

if (names(dt)[1] == "V1" && as.character(dt[1, 1]) == "gene") {
  new_names <- as.character(unlist(dt[1, ]))
  setnames(dt, new_names)
  dt <- dt[-1]
}

setnames(dt, 1, "gene_id")

# remove "DEPRECATED_" prefix if exists
dt[, gene_id := sub("^DEPRECATED_", "", gene_id)]

# factor 
factor_cols <- setdiff(names(dt), "gene_id")

# change to numeric
dt[, (factor_cols) := lapply(.SD, as.numeric), .SDcols = factor_cols]

# avoid division by zero
eps <- 1e-8

# calculate mean ratio for each factor and select top genes
marker_list <- lapply(factor_cols, function(fc) {
  other_cols <- setdiff(factor_cols, fc)
  
  tmp <- copy(dt)
  tmp[, target := get(fc)]
  tmp[, other_mean := rowMeans(.SD, na.rm = TRUE), .SDcols = other_cols]
  tmp[, MeanRatio := target / pmax(other_mean, eps)]
  
  tmp <- tmp[is.finite(MeanRatio) & !is.na(MeanRatio)]
  tmp <- tmp[MeanRatio > mean_ratio_threshold]
  tmp <- tmp[order(-MeanRatio, -target)]
  tmp <- head(tmp, max_num_genes)
  
  tmp[, set_id := paste0("Factor_", fc)]
  tmp[, .(set_id, gene_id, MeanRatio, target, other_mean)]
})

marker_stats <- rbindlist(marker_list, use.names = TRUE, fill = TRUE)

# factor marker counts
print(table(marker_stats$set_id))

# save results
fwrite(marker_stats, file.path(out_dir, "factor_mean_ratio_top200.tsv"), sep = "\t")