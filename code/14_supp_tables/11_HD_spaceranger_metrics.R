library(tidyverse)
library(here)
library(sessioninfo)
library(SpatialExperiment)
library(qs2)

spe_path = here(
    'processed-data', '12_apps_and_sharing', '01_prep_objects', 'spe_shiny.qs2'
)
out_path = here(
    'processed-data', '14_supp_tables', 'HD_spaceranger_metrics.csv'
)

spe = qs_read(spe_path)
spe_df = colData(spe) |>
    as_tibble() |>
    count(donor) |>
    dplyr::rename(`Sample ID` = donor, num_cells_post_QC = n)

metrics_paths = list.files(
        here('processed-data', '01_spaceranger', 'five_samples_10_2025'),
        full.names = TRUE
    ) |>
    file.path('outs', 'metrics_summary.csv')

lapply(metrics_paths, function(x) read_csv(x, show_col_types = FALSE)) |>
    bind_rows() |>
    mutate(
        `Sample ID` = sprintf('Br%s', str_extract(`Sample ID`, '[0-9]{4}$'))
    ) |>
    left_join(spe_df, by = "Sample ID") |>
    write_csv(out_path)

session_info()
