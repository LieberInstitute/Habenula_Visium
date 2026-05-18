library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)
library(duckplyr)
library(scater)
library(BiocSingular)
library(BiocParallel)

k = 10
spe_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'spe_raw.rds'
)
cell_bins_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'cellular_bins.csv.gz'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.parquet'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation', 'reduced_dims'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

# Load cellular bin IDs to exclude
cellular_bin_ids = read_csv_duckdb(cell_bins_path) |>
    select(bin_id) |>
    collect()

# Load FICTURE cluster assignments: collect first, then build key in R
ficture_df = read_parquet_duckdb(ficture_path) |>
    select(sample_id, barcode, ficture_cluster = paste0("FICTURE_k", k)) |>
    filter(!is.na(ficture_cluster) & ficture_cluster != "NA") |>
    collect() |>
    mutate(key = paste0(sample_id, "_", barcode)) |>
    select(key, ficture_cluster)

# Load SPE, subset to non-cellular bins, and merge ficture_cluster into colData
spe = readRDS(spe_path)
spe = spe[, !spe$key %in% cellular_bin_ids$bin_id]
spe$ficture_cluster = deframe(ficture_df)[spe$key]

spe_pb = registration_pseudobulk(
    spe, var_registration = "ficture_cluster", var_sample_id = "sample_id"
)

spe_pb = runPCA(
    spe_pb, ncomponents = 10, BSPARAM = IrlbaParam(),
    BPPARAM = MulticoreParam(num_cores)
)

p = plotReducedDim(spe_pb, dimred = "PCA", colour_by = "ficture_cluster")
pdf(file.path(plot_dir, sprintf('PCA_factor_k%d_no_cellular.pdf', k)))
print(p)
dev.off()

p = plotReducedDim(spe_pb, dimred = "PCA", colour_by = "sample_id")
pdf(file.path(plot_dir, sprintf('PCA_sample_id_k%d_no_cellular.pdf', k)))
print(p)
dev.off()

session_info()
