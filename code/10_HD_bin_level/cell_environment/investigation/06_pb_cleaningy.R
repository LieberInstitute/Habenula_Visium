#   The pseudobulked raw counts for the extracellular data are definitely
#   problematic. What about the pseudobulked batch-corrected counts?

#   The extracellular registration stats are indeed highly noisy (there are
#   almost no significant markers for any cluster at any k). Check if the data
#   looks ok in low-dimensional space compared to cellular data

library(here)
library(SpatialExperiment)
library(sessioninfo)
library(tidyverse)
library(scater)
library(BiocSingular)
library(BiocParallel)
library(duckplyr)
library(spatialLIBD)

k = 10
dataset = c('extra', 'all')[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

if (dataset == 'extra') {
    spe_path = here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
        'spe_filtered.rds'
    )
    chunks_path =  here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
        'cleaningY', 'temp_chunks', '%d.rds'
    )
    ficture_path = here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
        'ficture_plotting', 'extracellular.parquet'
    )
    num_chunks = 50
} else {
    spe_path = here(
        'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
        'spe', 'y_clean_spe.rds'
    )
    ficture_path = here(
        'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
        'bin_level_clusters_batch.parquet'
    )
}
svg_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out',
    'merged_SVGs.txt'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation', 'reduced_dims'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

spe = readRDS(spe_path)

if (dataset == 'extra') {
    #   Read in cleaningY chunks of counts and merge
    cleaned_counts = list()
    for (i in seq_len(num_chunks)) {
        cleaned_counts[[i]] = readRDS(sprintf(chunks_path, i))
    }
    cleaned_counts = do.call(rbind, cleaned_counts)

    #   Attach to SPE, pretending they're raw counts
    stopifnot(identical(dim(spe), dim(cleaned_counts)))
    assays(spe)$counts = cleaned_counts
    rm(cleaned_counts)
    gc()
}

#   Subset to SVGs for speed. We aren't guaranteed SVGs provide good signal in
#   the extracellular bins, but we need to reduce runtime somehow
svg = readLines(svg_path)
stopifnot(all(svg %in% rownames(spe)))
spe = spe[svg, ]

ficture_df = read_parquet_duckdb(ficture_path, prudence = 'stingy') |>
    dplyr::rename(ficture_cluster = paste0('k', k)) |>
    filter(!is.na(ficture_cluster)) |>
    select(bin_key, ficture_cluster)

spe$ficture_cluster = tibble(
        bin_key = paste(colnames(spe), spe$sample_id, sep = "_")
    ) |>
    left_join(ficture_df, by = 'bin_key') |>
    collect() |>
    pull(ficture_cluster)
spe = spe[, !is.na(spe$ficture_cluster)]

spe_pb = registration_pseudobulk(
    spe, var_registration = "ficture_cluster", var_sample_id = "sample_id"
)

spe_pb = runPCA(
    spe_pb, ncomponents = 10, BSPARAM = IrlbaParam(),
    BPPARAM = MulticoreParam(num_cores)
)

p = plotReducedDim(spe_pb, dimred = "PCA", colour_by = "ficture_cluster")
pdf(file.path(plot_dir, sprintf('cleaningY_PCA_factor_k%d_%s.pdf', k, dataset)))
print(p)
dev.off()

p = plotReducedDim(spe_pb, dimred = "PCA", colour_by = "sample_id")
pdf(
    file.path(
        plot_dir, sprintf('cleaningY_PCA_sample_id_k%d_%s.pdf', k, dataset)
    )
)
print(p)
dev.off()

session_info()
