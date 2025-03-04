library(here)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)

#   Get Leiden resolution from array task ID
res = (seq_len(10) / 10)[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
res_neat = sub('\\.', '_', as.character(res))

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_8',
    sprintf('leiden_res0_8.csv', res_neat)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'pseudobulk_spe', sprintf('%s.rds', res_neat)
)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$banksy = cluster_df$banksy_lambda0_8[match(spe$key, cluster_df$key)]

#   Pseudobulk
spe_pseudo = registration_pseudobulk(
    spe,
    var_registration = 'banksy',
    var_sample_id = 'sample_id',
    pseudobulk_rds_file = out_path
)

session_info()
