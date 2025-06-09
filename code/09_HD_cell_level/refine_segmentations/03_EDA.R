#   Explore the spatial distributions of cell sizes and number of extracellular
#   bins. This is to better understand if these can be used as metrics for QC in
#   the samples with poor-quality H&E images

library(here)
library(tidyverse)
library(data.table)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
extra_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'extracellular_bins.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_4.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'refine_segmentations', 'EDA'
)

dir.create(plot_dir, showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$labels_joint_source == 'primary']

extra_df = fread(extra_path) |>
    as_tibble() |>
    mutate(key = paste(cell_id, sample_id, sep = '_')) |>
    filter(key %in% spe$key) |>
    select(key) |>
    group_by(key) |>
    summarize(num_neighbors = n())
    