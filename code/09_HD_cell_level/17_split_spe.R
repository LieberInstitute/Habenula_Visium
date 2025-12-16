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
plot_dir = here('plots', '09_HD_cell_level', 'new_samples2', 'split_spe')
db_scan_eps = 300 # found through trial and error

#   Will apply SpatialExperiment::rotateObject() and mirrorObject(axis = 'v'),
#   in that order, according to the following plan by tissue piece and donor.
#   Kelsey provided the anatomical orientations over Slack
rotate_plan = c(
    'Br9090_1' = 0,
    'Br9090_2' = 0,
    'Br3942_1' = 0,
    'Br3942_2' = 180,
    'Br9902_1' = 60,
    'Br9902_2' = 240,
    'Br8433_1' = 0,
    'Br8433_2' = 0,
    'Br8667_1' = 180,
    'Br8667_2' = 0
)
mirror_plan = c(
    'Br9090_1' = FALSE,
    'Br9090_2' = FALSE,
    'Br3942_1' = TRUE,
    'Br3942_2' = TRUE,
    'Br9902_1' = TRUE,
    'Br9902_2' = TRUE,
    'Br8433_1' = FALSE,
    'Br8433_2' = FALSE,
    'Br8667_1' = FALSE,
    'Br8667_2' = FALSE
)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(file.path(plot_dir, 'section_identity'), showWarnings = FALSE)

spe = readRDS(spe_in_path)

spe_list = list()
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
        file.path(plot_dir, 'section_identity', sprintf('%s.png', sample_id)),
        width = 1000, height = 1000
    )
    print(p)
    dev.off()

    for (piece in c(1, 2)) {
        spe_piece = spe_sub[, spe_sub$tissue_piece == piece]
        spe_piece$sample_id = sprintf(
            'Br%s_%d', str_extract(sample_id, '[0-9]{4}$'), piece
        )

        if (rotate_plan[[spe_piece$sample_id]] != 0) {
            spe_piece = rotateObject(
                spe_piece, degrees = rotate_plan[[spe_piece$sample_id]]
            )
        }
        if (mirror_plan[[spe_piece$sample_id]]) {
            spe_piece = mirrorObject(spe_piece, axis = 'v')
        }

        spe_list[[spe_piece$sample_id]] = spe_piece
    }
}
