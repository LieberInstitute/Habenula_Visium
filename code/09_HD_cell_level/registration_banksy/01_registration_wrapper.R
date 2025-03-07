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
    sprintf('leiden_res%s.csv', res_neat)
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'pseudobulk_spe', sprintf('%s.rds', res_neat)
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'registration_banksy',
    'modeling_results', sprintf('%s.rds', res_neat)
)

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE)

#   Load and bring counts into memory to speed up computations. Despite the huge
#   size of the data, the memory footprint is manageable due to the extreme
#   sparsity of the data
spe = loadHDF5SummarizedExperiment(spe_dir)
assays(spe)$counts = as(assays(spe)$counts, "dgCMatrix")

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$banksy = cluster_df$banksy_lambda0_8[match(spe$key, cluster_df$key)]

#   Pseudobulk
model_results = registration_wrapper(
    spe,
    var_registration = 'banksy',
    var_sample_id = 'sample_id',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
