library(here)
library(tidyverse)
library(sessioninfo)

marker_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA', 'gene_sets',
    'fine.tsv'
)
color_path = here('code', 'hd_colors.R')

source(color_path)

message("Top markers by cell type:")
read_tsv(marker_path, show_col_types = FALSE) |>
    mutate(
        set_id = case_when(
                set_id == 'Endo.microglia' ~ 'Endo/microglia',
                set_id == 'Excit.Thal.Inhib_LHb_4.2' ~ 'Excit.Thal/Inhib_LHb_4.2',
                set_id == 'LHb.4.Inhib_LHb_4.2' ~ 'LHb.4/Inhib_LHb_4.2',
                TRUE ~ set_id
            ) |>
            factor(levels = names(cell_type_colors))
    ) |>
    group_by(set_id) |>
    arrange(desc(MeanRatio)) |>
    slice_head(n = 1) |>
    ungroup() |>
    arrange(set_id) |>
    pull(gene_name) |>
    dput()

session_info()
