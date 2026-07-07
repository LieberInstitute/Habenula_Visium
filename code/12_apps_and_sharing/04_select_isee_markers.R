library(here)
library(tidyverse)
library(sessioninfo)

marker_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA', 'gene_sets',
    'fine.tsv'
)

message("Top markers by cell type:")
read_tsv(marker_path, show_col_types = FALSE) |>
    group_by(set_id) |>
    arrange(desc(MeanRatio)) |>
    slice_head(n = 1) |>
    ungroup() |>
    arrange(set_id) |>
    pull(gene_name) |>
    dput()

session_info()
