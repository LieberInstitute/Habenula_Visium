#   The existing FICTURE data required to generate a manual plot in R of
#   FICTURE clusters is extremely large, and scattered in separate files.
#   Generate one CSV with donor ID, k 14 cluster, and spatial coordinates
#   to make plotting less computationally expensive later

library(here)
library(tidyverse)
library(sessioninfo)

cluster_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
coord_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'ficture_inputs' , 'cleany', 'input.tsv.gz'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'k14_delete_this.csv.gz'
)

#   dplyr/readr equivalent to the duckplyr version
read_table(coord_path) |>
    dplyr::distinct(sample_id, barcode, .keep_all = TRUE) |>
    select(sample_id, barcode, X, Y) |>
    left_join(
        read_csv(cluster_path) |>
            select(sample_id, barcode, FICTURE_k14),
        by = c('sample_id', 'barcode')
    ) |>
    filter(!is.na(FICTURE_k14), FICTURE_k14 != "NA") |>
    select(sample_id, X, Y, FICTURE_k14) |>
    write_csv(out_path)

session_info()
