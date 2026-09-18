#   The existing FICTURE data required to generate a manual plot in R of
#   FICTURE clusters is extremely large, and scattered in separate files.
#   Generate one parquet file with donor ID, k = 8 cluster, and spatial
#   coordinates to make plotting less computationally expensive later

library(here)
library(tidyverse)
library(sessioninfo)
library(duckplyr)

cluster_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
coord_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'ficture_inputs', 'cleany', 'input.tsv.gz'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'k8_cluster_coords.parquet'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = TRUE)

#   Directly read cluster and spatial coordinate info into DuckDB-managed
#   memory without touching R's memory
temp = read_csv_duckdb(coord_path, prudence = 'stingy') |>
    dplyr::distinct(sample_id, barcode, .keep_all = TRUE) |>
    select(sample_id, barcode, X, Y) |>
    left_join(
        read_csv_duckdb(cluster_path, prudence = 'stingy') |>
            select(sample_id, barcode, FICTURE_k8),
        by = c('sample_id', 'barcode')
    ) |>
    filter(!is.na(FICTURE_k8), FICTURE_k8 != "NA") |>
    dplyr::rename(x = X, y = Y, ficture_cluster = FICTURE_k8) |>
    select(sample_id, x, y, ficture_cluster) |>
    compute_parquet(out_path)

session_info()
