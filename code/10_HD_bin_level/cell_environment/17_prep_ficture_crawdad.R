#   We plan to run a version of CRAWDAD incorporating clusters both from Banksy
#   and FICTURE. This script generates the inputs for CRAWDAD with both types
#   of clusters

library(tidyverse)
library(here)
library(duckplyr)

extra_bin_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_bins.csv.gz'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_outputs', 'cleaningy', 'k_10', 'analysis', 'nF10.d_12',
    'cleaningy_joined_input.tsv.gz'
)
ficture_out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_outputs', 'cleaningy', 'k_10', 'cleaned_clusters.parquet'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

ficture_df = read_csv_duckdb(ficture_path, prudence = 'stingy') |>
    distinct(barcode, sample_id, factor_K1) |>
    filter(factor_K1 != 'NA') |>
    compute_parquet(ficture_out_path) |>
    collect()
