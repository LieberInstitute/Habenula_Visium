# combine the cell type and then do the registration

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'spe_norm_filtered_split.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy', 'leiden_res1_8.csv'
)
cell_type = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy', 'cluster_annotation.csv'
)
pseudo_path_broad = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'pseudobulk_spe', 'broad.rds'
)
pseudo_path_fine = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'pseudobulk_spe', 'fine.rds'
)
model_path_broad = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', 'broad.rds'
)
model_path_fine = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', 'fine.rds'
)

spe = readRDS(spe_path)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
cell_type_df = read_csv(cell_type, show_col_types = FALSE)
cluster_df = cluster_df %>%
    left_join(
        cell_type_df,
        by = c("banksy" = "cluster")
    )
stopifnot(all(spe$key %in% cluster_df$key))
spe$broad_cell_type = cluster_df$broad_cell_type[match(spe$key, cluster_df$key)]
spe$fine_cell_type = cluster_df$fine_cell_type[match(spe$key, cluster_df$key)]
spe$broad_cell_type= factor(spe$broad_cell_type, levels = sort(unique(spe$broad_cell_type)))
spe$fine_cell_type= factor(spe$fine_cell_type, levels = sort(unique(spe$fine_cell_type)))

model_results_broad = registration_wrapper(
    spe,
    var_registration = 'broad_cell_type',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path_broad
)

model_results_fine = registration_wrapper(
    spe,
    var_registration = 'fine_cell_type',
    var_sample_id = 'sample_id',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path_fine
)

saveRDS(model_results_broad, model_path_broad)
saveRDS(model_results_fine, model_path_fine)

session_info()
