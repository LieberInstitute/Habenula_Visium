#   The extracellular registration stats are indeed highly noisy (there are
#   almost no significant markers for any cluster at any k). Check if the data
#   looks ok in low-dimensional space compared to cellular data

library(here)
library(SpatialExperiment)
library(sessioninfo)
library(tidyverse)
library(duckplyr)
library(scater)
library(BiocSingular)
library(BiocParallel)
library(qs2)

k = 10
dataset = c('extra', 'all')[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

if (dataset == 'extra') {
    spe_path = here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
        'spe_filtered.rds'
    )
    cluster_path = here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
        'ficture_plotting', 'extracellular.parquet'
    )
} else {
    spe_path = here(
        'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
        'spe_raw.rds'
    )
    cluster_path = here(
        'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
        'bin_level_clusters_batch.parquet'
    )
}
out_path = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation', 'reduced_dims', sprintf('spe_%s.qs2', dataset)
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation', 'reduced_dims'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

spe = runPCA(
    spe, ncomponents = 10, BSPARAM = IrlbaParam(),
    BPPARAM = MulticoreParam(num_cores)
)

#   PCA took a while to compute; save in case something goes wrong or I want to
#   edit a plot
assays(spe) = list()
qs_save(spe, out_path)

ficture_df = read_parquet_duckdb(cluster_path, prudence = 'stingy') |>
    dplyr::rename(factor_K1 = paste0('k', k)) |>
    select(bin_key, factor_K1)

spe$ficture_cluster = tibble(
        bin_key = paste(colnames(spe), spe$sample_id, sep = "_")
    ) |>
    left_join(ficture_df, by = 'bin_key') |>
    collect() |>
    pull(factor_K1)

p = plotReducedDim(spe, dimred = "PCA", colour_by = "factor_K1")
png(
    file.path(plot_dir, sprintf('PCA_factor_k%d_%s.png', k, dataset)),
    width = 400, height = 400
)
print(p)
dev.off()

p = plotReducedDim(spe, dimred = "PCA", colour_by = "sample_id")
png(
    file.path(plot_dir, sprintf('PCA_sample_id_k%d_%s.png', k, dataset)),
    width = 400, height = 400
)
print(p)
dev.off()

session_info()
