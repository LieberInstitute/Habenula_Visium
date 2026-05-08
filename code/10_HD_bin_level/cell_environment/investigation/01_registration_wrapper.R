#   Run registration_wrapper() by k on cellular data labeled by most-dominant
#   extracellular FICTURE cluster

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(duckplyr)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'crawdad', sprintf('input_cells_k%d.csv.gz', k)
)
pseudo_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'investigation', 'registration', sprintf('%s_pb.rds', k)
)
model_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'investigation', 'registration', sprintf('%s_modeling.rds', k)
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

ficture_df = read_csv_duckdb(ficture_path, prudence = 'stingy') |>
    select(cell_key, ficture_cluster)

spe$ficture_cluster = tibble(cell_key = spe$key) |>
    left_join(ficture_df, by = 'cell_key') |>
    pull(ficture_cluster)

message(
    sprintf(
        "Dropping %.1f%% of cells without a defined FICTURE cluster",
        mean(is.na(spe$ficture_cluster)) * 100
    )
)
spe = spe[, !is.na(spe$ficture_cluster)]

#   Pseudobulk by tissue section + FICTURE cluster
model_results = registration_wrapper(
    spe,
    var_registration = 'ficture_cluster',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
