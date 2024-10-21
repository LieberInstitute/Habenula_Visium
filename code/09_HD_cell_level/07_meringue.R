library(MERINGUE)
library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(spatialLIBD)
library(tidyverse)

sample_id = 'H1-W369TJK_D1_9090'
spe_in_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
plot_dir = here('plots', '09_HD_cell_level', 'meringue')
out_path = here(
    'processed-data', '09_HD_cell_level',
    paste0(sample_id, '_meringue_clusters.csv')
)

dir.create(plot_dir, showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_in_dir)
rownames(spatialCoords(spe)) = colnames(spe)
spe$exclude_overlapping = FALSE

#   Grab a corner (1% of the total area) for testing
# a = spatialCoords(spe)[,1] > 0.1 * min(spatialCoords(spe)[, 1]) + 0.9 * max(spatialCoords(spe)[, 1])
# b = spatialCoords(spe)[,2] > 0.1 * min(spatialCoords(spe)[, 2]) + 0.9 * max(spatialCoords(spe)[, 2])
# spe = spe[, a & b]

#   Form network of spatial neighbors
message(Sys.time(), ' | Forming network of spatial neighbors')
W <- getSpatialNeighbors(spatialCoords(spe), filterDist = 500)
png(
    file.path(
        plot_dir, sprintf('neighbor_network_%s.png', sample_id)
    ),
    width = 1500, height = 1500
)
plotNetwork(spatialCoords(spe), W)
dev.off()

message(Sys.time(), ' | Performing spatially aware clustering')
spe$meringue_cluster = getSpatiallyInformedClusters(
    reducedDims(spe)$PCA, W = W, k = 50
)

png(
    file.path(
        plot_dir, sprintf('meringue_clusters_%s.png', sample_id)
    ),
    width = 1500, height = 1500
)
vis_clus(
        spe, clustervar = 'meringue_cluster', is_stitched = TRUE, point_size = 1
    ) +
    guides(fill = guide_legend(override.aes = list(size = 4)))
dev.off()

#   Write clustering results to CSV
write_csv(
    tibble(key = colnames(spe), meringue_cluster = spe$meringue_cluster),
    out_path
)

session_info()
