#   Begin spatial registration process of Banksy clusters by running
#   registration_wrapper()

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)

res = c(seq_len(20) / 10, 4, 8)[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
res_neat = sub('\\.', '_', as.character(res))

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    sprintf('leiden_res%s.csv', res_neat)
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'pseudobulk_spe', sprintf('%s.rds', res_neat)
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', sprintf('%s.rds', res_neat)
)

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE, recursive = TRUE)

spe = readRDS(spe_path)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$banksy = cluster_df$banksy[match(spe$key, cluster_df$key)]
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
