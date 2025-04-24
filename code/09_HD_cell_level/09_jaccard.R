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

prep_clustering_results = function(spe_dir, cluster_path, cluster_colname) {
    spe = loadHDF5SummarizedExperiment(spe_dir)
    
    cluster_df = tibble(
            x = round(spatialCoords(spe)[, 'pxl_col_in_fullres'], 1),
            y = round(spatialCoords(spe)[, 'pxl_row_in_fullres'], 1),
            sample_id = spe$sample_id,
            key = spe$key,
            segmentation_type = spe$labels_joint_source
        ) |>
        filter(segmentation_type == 'primary') |>
        left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
        dplyr::rename(!!cluster_colname := banksy_lambda0_8) |>
        select(x, y, sample_id, {{ cluster_colname }})

    stopifnot(!any(is.na(cluster_df[[cluster_colname]])))
    
    return(cluster_df)
}

cluster_good_df = prep_clustering_results(
    spe_good_dir, cluster_good_path, 'cluster_good'
)
cluster_bad_df = prep_clustering_results(
    spe_bad_dir, cluster_bad_path, 'cluster_bad'
)
