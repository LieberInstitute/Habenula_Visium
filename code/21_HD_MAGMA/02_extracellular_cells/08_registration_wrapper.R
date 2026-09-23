#   Begin spatial registration process of Banksy clusters by running
#   registration_wrapper()

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(qs2)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.qs2'
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'pb_spe.rds'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'modeling_results.rds'
)

dir.create(dirname(pseudo_path), showWarnings = FALSE)

spe = qs_read(spe_path)

model_results = registration_wrapper(
    spe,
    var_registration = 'cell_type',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
