library(here)
library(tidyverse)
library(SpatialExperiment)
library(duckplyr)

spe_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'spe_raw.rds'
)
extra_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
cell_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'cellular.parquet'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

#   Was going to investigate whether zeros in extracellular counts were biased
#   towards (potentially different sets of) samples. Had AI explore a bit but
#   moved towards prioritizing the next script
