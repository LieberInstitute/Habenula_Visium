library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(DeconvoBuddies)
library(sessioninfo)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
k2_path = here('processed-data', '09_HD_cell_level', 'banksy', 'k2.csv')

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Add banksy k2 clusters to SPE
k2 = read_csv(k2_path, show_col_types = FALSE)
spe$banksy_k2 = k2$banksy_lambda0.2[match(as.numeric(colnames(spe)), k2$key)]

markers = findMarkers_1vAll(
    spe, assay_name = "logcounts", cellType_col = "banksy_k2", mod = NULL
)

session_info()
