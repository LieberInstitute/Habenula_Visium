#   This script has two goals:
#       1. Break up existing capture areas into the 2 constituent tissue pieces
#       2. Arrange constituent pieces anatomically consistently (dorsal up,
#          medial left)
#   Then just write a new SPE with 10 sample IDs

library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)
library(spatialLIBD)
library(dbscan)

spe_in_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'spe_norm_filtered.rds'
)
spe_out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'spe_norm_filtered_split.rds'
)
db_scan_eps = 300 # found through trial and error
plot_dir = here('plots', '09_HD_cell_level', 'new_samples2', 'split_spe')

dir.create(plot_dir, showWarnings = FALSE)

spe = readRDS(spe_in_path)

for (sample_id in unique(spe$sample_id)) {
    spe_sub = spe[, spe$sample_id == sample_id]
    
    #   Use DBSCAN to separate the two tissue pieces
    spe_sub$tissue_piece = dbscan(
        spatialCoords(spe_sub), eps = db_scan_eps
    )$cluster
    stopifnot(length(unique(spe_sub$tissue_piece)) == 2)

    p = vis_clus(
        spe_sub, clustervar = 'tissue_piece', point_size = 20, spatial = FALSE,
        is_stitched = TRUE
    )
    png(
        file.path(plot_dir, sprintf('%s.png', sample_id)),
        width = 1000, height = 1000
    )
    print(p)
    dev.off()
}