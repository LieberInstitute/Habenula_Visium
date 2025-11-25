#   Begin spatial registration process of Banksy clusters by running
#   registration_wrapper()

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)

res = c(seq_len(20) / 10, 4, 8)[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
res_neat = sub('\\.', '_', as.character(res))
lambda_neat = 'lambda0_2'

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', lambda_neat,
    sprintf('leiden_res%s.csv', res_neat)
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'pseudobulk_spe', lambda_neat, sprintf('%s.rds', res_neat)
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'modeling_results', lambda_neat, sprintf('%s.rds', res_neat)
)

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$banksy = cluster_df[[paste0('banksy_', lambda_neat)]][
    match(spe$key, cluster_df$key)
]
spe$banksy = factor(spe$banksy, levels = sort(unique(spe$banksy)))

#   Pseudobulk
model_results = registration_wrapper(
    spe,
    var_registration = 'banksy',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
