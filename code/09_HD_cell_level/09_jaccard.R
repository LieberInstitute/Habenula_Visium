library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)

spe_good_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
spe_bad_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
cluster_good_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_8',
    'leiden_res1.csv'
)
cluster_bad_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_8',
    'leiden_res1.csv'
)

spe_good = loadHDF5SummarizedExperiment(spe_good_dir)
cluster_good = read_csv(cluster_good_path, show_col_types = FALSE) |>
    dplyr::rename(cluster_good = banksy_lambda0_8)

cluster_good_df = tibble(
        x = round(spatialCoords(spe_good)[, 'pxl_col_in_fullres'], 1),
        y = round(spatialCoords(spe_good)[, 'pxl_row_in_fullres'], 1),
        sample_id = spe_good$sample_id,
        key = spe_good$key,
        segmentation_type = spe_good$labels_joint_source
    ) |>
    filter(segmentation_type == 'primary') |>
    left_join(cluster_good, by = 'key') |>
    select(x, y, sample_id, cluster_good)
    