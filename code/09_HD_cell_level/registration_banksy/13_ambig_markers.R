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
    select(ensembl, gene, matches('^fdr_')) |>
    pivot_longer(
        cols = starts_with("fdr_"),
        names_to = "banksy",
        values_to = "fdr"
    ) |>
    mutate(cluster = as.integer(sub('fdr_X', '', banksy))) |>
    filter(cluster %in% ambig_clusters, fdr < sig_cutoff) |>
    dplyr::rename(gene_id = ensembl, gene_name = gene) |>
    select(cluster, gene_id, gene_name, fdr) |>
    arrange(cluster, fdr) |>
    group_by(cluster) |>
    slice_head(n = num_top_markers) |>
    write_csv(out_path)

session_info()
