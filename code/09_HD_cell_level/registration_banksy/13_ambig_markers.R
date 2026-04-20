#   Find top markers in an easy-to-browse format for clusters whose cell types
#   were unclear from spatial registration

library(here)
library(tidyverse)
library(sessioninfo)

plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'annotation'
)
model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8.rds'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'ambig_markers.csv'
)
sig_cutoff = 0.05
num_top_markers = 50
ambig_clusters = c(2, 13, 14, 16)

readRDS(model_path)$enrichment |>
    as_tibble() |>
    select(ensembl, gene, matches('^(fdr|t_stat)_')) |>
    pivot_longer(
        cols = matches('^(fdr|t_stat)_'),
        names_to = c(".value", "cluster"),
        names_pattern = "^(fdr|t_stat)_X(\\d+)$"
    ) |>
    mutate(cluster = as.integer(cluster)) |>
    dplyr::rename(gene_id = ensembl, gene_name = gene) |>
    filter(cluster %in% ambig_clusters, fdr < sig_cutoff, t_stat > 0) |>
    select(cluster, gene_id, gene_name, fdr, t_stat) |>
    group_by(cluster) |>
    arrange(fdr) |>
    slice_head(n = num_top_markers) |>
    arrange(cluster, fdr) |>
    write_csv(out_path)

session_info()
