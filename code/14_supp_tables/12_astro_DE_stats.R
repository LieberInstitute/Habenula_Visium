#   In a main figure we have logFC oriented such that up in MHb is positive
#   logFC. Make the supplementary stats table consistent with this orientation

library(tidyverse)
library(here)
library(sessioninfo)

de_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'aggregated_DE_stats.csv.gz'
)
out_path = here('processed-data', '14_supp_tables', 'astro_DE_stats.csv')

read_csv(de_path) |>
    mutate(logFC = -1 * logFC, t = -1 * t) |>
    write_csv(out_path)

session_info()
