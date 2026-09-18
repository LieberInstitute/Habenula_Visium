#   Register each extracellular FICTURE run against Banksy cell types (not 
#   clusters). We'll later (in another script) ask if FICTURE clusters nearby cells
#   agree transcriptionally with the Bansky cell types (cellular data)

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(duckplyr)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'spe_filtered.rds'
)
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
pseudo_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'pseudobulk_spe', sprintf('%s.rds', k)
)
model_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'modeling_results', sprintf('%s.rds', k)
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

ficture_df = read_parquet_duckdb(cluster_path, prudence = 'stingy') |>
    dplyr::rename(factor_K1 = paste0('k', k)) |>
    select(bin_key, factor_K1)

spe$ficture_cluster = tibble(
        bin_key = paste(colnames(spe), spe$sample_id, sep = "_"),
        idx = seq_along(colnames(spe))
    ) |>
    left_join(ficture_df, by = 'bin_key') |>
    #   This is critical, as duckplyr does not naturally preserve row order
    arrange(idx) |>
    pull(factor_K1)

message(
    sprintf(
        "Dropping %.1f%% of bins without a defined FICTURE cluster",
        mean(is.na(spe$ficture_cluster)) * 100
    )
)
spe = spe[, !is.na(spe$ficture_cluster)]

rowData(spe)$gene_id = rownames(spe)

#   Pseudobulk by donor + FICTURE cluster
model_results = registration_wrapper(
    spe,
    var_registration = 'ficture_cluster',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'symbol',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
