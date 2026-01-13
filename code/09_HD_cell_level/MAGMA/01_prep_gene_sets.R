library(tidyverse)
library(here)
library(sessioninfo)

model_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'modeling_results', 'lambda0_2', '1_7.rds'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'MAGMA', 'gene_sets',
    'enrichment_markers.tsv'
)
sig_cutoff = 0.05

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

readRDS(model_path)$enrichment |>
    as_tibble() |>
    select(ensembl, matches('^fdr_X')) |>
    pivot_longer(
        cols = starts_with("fdr_X"),
        names_to = "banksy",
        values_to = "fdr"
    ) |>
    filter(fdr < sig_cutoff) |>
    mutate(set_id = as.integer(sub('fdr_X', '', banksy))) |>
    dplyr::rename(gene_id = ensembl) |>
    select(set_id, gene_id) |>
    write_tsv(out_path)

session_info()
