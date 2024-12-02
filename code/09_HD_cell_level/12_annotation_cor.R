library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

t_stat_paths = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'markers',
    sprintf('k%s.rds', 2:28)
)

t_stats = lapply(
    t_stat_paths,
    function(path) {
        readRDS(path) |>
            column_to_rownames('gene')
    }
)

session_info()
