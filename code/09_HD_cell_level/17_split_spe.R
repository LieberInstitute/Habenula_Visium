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
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
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
    'Br9902_1' = 50,
    'Br9902_2' = 235,
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

handle_Br9902 = function(spe_piece, this_sample_id) {
    #   Br9902 requires rotations that are not multiples of 90 degrees, which is
    #   not supported by rotateObject(). I'll just rotate the coordinates;
    #   handling the image in a meaningful way is possible but quite complex and
    #   difficult (we probably won't plot the background images for the
    #   manuscript)
    radians = rotate_plan[[this_sample_id]] * pi / 180

    #   Determine the matrix by which left-multiplication represents
    #   rotation. Then apply rotation about the origin
    rotation_mat <- matrix(
        c(
            cos(radians), sin(radians), -1 * sin(radians), cos(radians)
        ),
        nrow = 2
    )

    #   Get the dimensions of the midpoint of the "rectangle" containing the set
    #   of spatialCoords within the object
    dim_mid = dim(imgRaster(spe_piece)) / scaleFactors(spe_piece)[1] / 2

    #   Rotate about the center of the image
    trans_vec = rep(dim_mid, each = ncol(spe_piece))
    new_coords = t(rotation_mat %*% t(spatialCoords(spe_piece) - trans_vec)) + trans_vec |>
        round() |>
        as.integer()
    dimnames(new_coords) = dimnames(spatialCoords(spe_piece))
    spatialCoords(spe_piece) = new_coords

    return(spe_piece)
}

dir.create(plot_dir, showWarnings = FALSE)
dir.create(file.path(plot_dir, 'section_identity'), showWarnings = FALSE)
dir.create(file.path(plot_dir, 'final_orientation'), showWarnings = FALSE)

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
        this_sample_id = sprintf(
            'Br%s_%d', str_extract(sample_id, '[0-9]{4}$'), piece
        )
        spe_piece$sample_id = this_sample_id

        if (grepl('^Br9902', this_sample_id)) {
            spe_piece = handle_Br9902(spe_piece, this_sample_id)
        } else {
            if (rotate_plan[[this_sample_id]] != 0) {
                spe_piece = rotateObject(
                    spe_piece, degrees = rotate_plan[[this_sample_id]]
                )
            }
            if (mirror_plan[[this_sample_id]]) {
                spe_piece = mirrorObject(spe_piece, axis = 'v')
            }
        }

        spe_list[[this_sample_id]] = spe_piece
    }
}

spe_split = do.call(cbind, spe_list)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe_split$key %in% cluster_df$key))
spe_split$banksy = cluster_df$banksy_lambda0_2[
    match(spe_split$key, cluster_df$key)
]
spe_split$banksy = factor(spe_split$banksy, levels = sort(unique(spe_split$banksy)))

#   Plot the Banksy clusters over the tissue in the split SPE
for (sample_id in unique(spe_split$sample_id)) {
    p = vis_clus(
            spe_split, sampleid = sample_id, clustervar = 'banksy',
            point_size = 20, spatial = TRUE, is_stitched = TRUE
        ) +
        guides(fill = guide_legend(override.aes = list(size = 8)))
    png(
        file.path(plot_dir, 'final_orientation', sprintf('%s.png', sample_id)),
        width = 1000, height = 1000
    )
    print(p)
    dev.off()
}

saveRDS(spe_split, spe_out_path)

session_info()
