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
svg_path = here(
    'processed-data', '10_HD_bin_level', 'nnSVG_out', 'H1-W369TJK_D1_9090.csv'
)
out_path = here(
    'processed-data', '09_HD_cell_level',
    sprintf('meringue_patterns_%s.csv', sample_id)
)
top_n = 200

svg = read_csv(svg_path, show_col_types = FALSE) |>
    arrange(rank) |>
    slice_head(n = top_n) |>
    pull(gene_id)

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

#   'getSpatialPatterns' is prohibitively slow. Use SVGs from nnSVG instead
message(Sys.time(), ' | Computing spatial cross-correlation')
scc = spatialCrossCorMatrix(mat = assays(spe)$logcounts[svg,], weight = W)

#   Group spatial patterns and plot them
message(Sys.time(), ' | Grouping and plotting spatial patterns')
pdf(file.path(plot_dir, sprintf('spatial_patterns_%s.pdf', sample_id)))
ggroup <- groupSigSpatialPatterns(
    pos = spatialCoords(spe), 
    mat = as(assays(spe)$logcounts[svg,], "dgCMatrix"), 
    scc = scc, 
    power = 1, 
    hclustMethod = 'ward.D', 
    deepSplit = 2,
    zlim=c(-1.5,1.5)
)
dev.off()

#   Save genes paired with any patterns
write_csv(
    tibble(
        gene_id = names(ggroup$groups), pattern_id = factor(ggroup$groups)
    ),
    out_path
)

session_info()
