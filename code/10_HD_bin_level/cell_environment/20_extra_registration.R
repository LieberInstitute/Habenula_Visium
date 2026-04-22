#   Register each extracellular FICTURE run against Banksy clusters (not cell
#   types). We'll later (in another script) ask if FICTURE clusters nearby cells
#   agree transcriptionally with the Bansky clusters (cellular data)

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(duckplyr)

k = c(4, 10, 20)[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'spe_filtered.rds'
)
cluster_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_outputs', 'cleaningy', sprintf('k_%d', k), 'analysis',
    sprintf('nF%d.d_12', k), 'cleaningy_joined_input.tsv.gz'
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

ficture_df = read_csv_duckdb(cluster_path, prudence = 'stingy') |>
    distinct(barcode, sample_id, factor_K1) |>
    filter(factor_K1 != 'NA')

spe$ficture_cluster = tibble(
        barcode = colnames(spe), sample_id = spe$sample_id
    ) |>
    left_join(ficture_df, by = c('barcode', 'sample_id')) |>
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
