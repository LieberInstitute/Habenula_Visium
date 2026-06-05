#   AnnDatas are by capture area and don't have info about tissue section.
#   Export a CSV mapping key to tissue section which will be used to add the
#   tissue info to a dataset-wide AnnData. Actually also include fixed
#   spatial coordinates (the split CSV has anatomically arranged samples) in
#   the CSV

library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'tissue_key_map.csv.gz'
)

spe = readRDS(spe_path)

tibble(
        key = spe$key, tissue_section = spe$sample_id,
        pxl_col_in_fullres = spatialCoords(spe)[, 'pxl_col_in_fullres'],
        pxl_row_in_fullres = spatialCoords(spe)[, 'pxl_row_in_fullres']
    ) |>
    write_csv(out_path)

session_info()
