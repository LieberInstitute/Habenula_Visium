#   Begin spatial registration process of Banksy clusters by running
#   registration_wrapper()

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.rds'
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'pb_spe.rds'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'modeling_results.rds'
)
tissue_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'tissue_key_map.csv.gz'
)

dir.create(dirname(pseudo_path), showWarnings = FALSE)

spe = readRDS(spe_path)

#   Add in info about tissue section. Pseudobulking is done by tissue section
#   as was done for the cellular data
tissue_df = read_csv(tissue_path, show_col_types = FALSE)
spe$tissue_section = tibble(key = spe$key) |>
    left_join(tissue_df, by = 'key') |>
    pull(tissue_section)
stopifnot(!any(is.na(spe$tissue_section)))

model_results = registration_wrapper(
    spe,
    var_registration = 'cell_type',
    var_sample_id = 'tissue_section',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
