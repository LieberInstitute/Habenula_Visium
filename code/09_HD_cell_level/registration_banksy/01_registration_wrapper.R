#   Begin spatial registration process of Banksy clusters by running
#   registration_wrapper()

library(here)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)

res = c(seq_len(20) / 10, 4, 8)[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
res_neat = sub('\\.', '_', as.character(res))
lambda_neat = 'lambda0_2'

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', lambda_neat,
    sprintf('leiden_res%s.csv', res_neat)
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'pseudobulk_spe', lambda_neat, sprintf('%s.rds', res_neat)
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'modeling_results', lambda_neat, sprintf('%s.rds', res_neat)
)
good_samples = c(
    "H1-W369TJK_D1_9090", "H1-MVPY9BW_A1_8433", "H1-MVPY9BW_D1_8667"
)

dir.create(dirname(pseudo_path), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(model_path), showWarnings = FALSE, recursive = TRUE)

#   Load and bring counts into memory to speed up computations. Despite the huge
#   size of the data, the memory footprint is manageable due to the extreme
#   sparsity of the data
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id %in% good_samples]
assays(spe)$counts = as(assays(spe)$counts, "dgCMatrix")

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
