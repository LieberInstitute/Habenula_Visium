#   AnnDatas are by capture area and don't have info about tissue section.
#   Export a CSV mapping key to tissue section which will be used to add the
#   tissue info to a dataset-wide AnnData

library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)

spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'tissue_key_map.csv.gz'
)

spe = readRDS(spe_path)

tibble(key = spe$key, tissue_section = spe$sample_id) |>
    write_csv(out_path)

session_info()
