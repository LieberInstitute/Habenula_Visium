library(tidyverse)
library(here)
library(sessioninfo)

extra_model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'modeling_results.rds'
)
cell_model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8_cell_types.rds'
)

prep_model_stats = function(model_path, col_suffix) {
    readRDS(model_path)$enrichment |>
        as_tibble() |>
        select(ensembl, gene, matches('^(fdr|t_stat)_')) |>
        pivot_longer(
            cols = matches('^(fdr|t_stat)_'),
            names_to = c(".value", "cell_type"),
            names_pattern = "^(fdr|t_stat)_(.*)$"
        ) |>
        dplyr::rename(
            gene_id = ensembl, gene_name = gene
        ) |>
        dplyr::rename(!!sprintf('t_stat_%s', col_suffix) := t_stat) |>
        select(-fdr)  
}

cell_df = prep_model_stats(cell_model_path, 'cell')
extra_df = prep_model_stats(extra_model_path, 'extra')

stat_df = inner_join(
    cell_df, extra_df, by = c('gene_id', 'gene_name', 'cell_type')
)
