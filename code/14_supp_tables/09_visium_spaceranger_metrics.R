library(tidyverse)
library(here)
library(sessioninfo)
library(SpatialExperiment)
library(qs2)

spe_path = here(
    'processed-data', '12_apps_and_sharing', '05_prep_visium_shiny',
    'spe_shiny.qs2'
)
metrics_path = here(
    'processed-data', '01_spaceranger', '%s', 'outs',
    'metrics_summary.csv'
)
out_path = here(
    'processed-data', '14_supp_tables', 'visium_spaceranger_metrics.csv'
)

spe = qs_read(spe_path)
spe_df = colData(spe) |>
    as_tibble() |>
    count(sample_id) |>
    dplyr::rename(num_spots_post_QC = n)

lapply(
        unique(spe_df$sample_id),
        function(x) read_csv(sprintf(metrics_path, x), show_col_types = FALSE)
    ) |>
    bind_rows() |>
    left_join(spe_df, by = c("Sample ID" = "sample_id")) |>
    write_csv(out_path)

session_info()
