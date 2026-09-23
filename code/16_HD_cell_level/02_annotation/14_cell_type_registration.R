#   Like 01_registration but for cell types annotated for the optimal resolution
#   rather than the clusters themselves. In other words, merge clusters with the
#   same annotation and reregister

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
pseudo_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'pseudobulk_spe', '1_8_cell_types.rds'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8_cell_types.rds'
)

dir.create(dirname(pseudo_path), showWarnings = FALSE)
dir.create(dirname(model_path), showWarnings = FALSE)

spe = readRDS(spe_path)

#   Merge in cell types
anno_df = read_csv(anno_path, show_col_types = FALSE)
spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    ) |>
    pull(cell_type)
stopifnot(!any(is.na(spe$cell_type)))

spe = spe[, spe$cell_type != 'Drop']

#   Pseudobulk
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
