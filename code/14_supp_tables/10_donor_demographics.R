library(here)
library(tidyverse)
library(sessioninfo)

in_path = here(
    'raw-data', 'sample_info', 'hb_atlas_donor_demographics_full.csv'
)
out_path = here(
    'processed-data', '14_supp_tables', 'donor_demographics.csv'
)

read_csv(in_path, show_col_types = FALSE) |>
    select(
        donor, in_visium_HE, in_visium_HD, in_multiome, age, ancestry, sex, PMI
    ) |>
    write_csv(out_path)

session_info()
