library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(Banksy)
library(cowplot)
library(scater)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
plot_dir = here('plots', '09_HD_cell_level', 'banksy')
out_path = here('processed-data', '09_HD_cell_level', 'banksy_clusters.csv')
random_seed = 0

#   Corresponding to "cell typing" and "domain segmentation"
lambda = c(0.2, 0.8)

dir.create(plot_dir, showWarnings = FALSE)
set.seed(random_seed)

spe = loadHDF5SummarizedExperiment(spe_dir)
spe$exclude_overlapping = FALSE

message(Sys.time(), ' | Running computeBanksy')
spe = computeBanksy(
    spe, assay_name = "logcounts", compute_agf = TRUE, seed = random_seed
)

message(Sys.time(), ' | Running PCA and UMAP on embedding')
spe = runBanksyPCA(spe, use_agf = TRUE, lambda = lambda, seed = random_seed)
spe = runBanksyUMAP(spe, use_agf = TRUE, lambda = lambda, seed = random_seed)

message(Sys.time(), ' | Performing clustering')
spe = clusterBanksy(spe, use_agf = TRUE, lambda = lambda, seed = random_seed)
spe = connectClusters(spe)

#   Get the names of the cluster columns and reducedDimNames() corresponding to
#   each lambda value
cluster_names = sapply(
        sprintf('^clust.*lam%s', lambda),
        function(pattern) {
            colnames(colData(spe))[grep(pattern, colnames(colData(spe)))]
        }
    ) |>
    unname()
rd_names = sapply(
        sprintf('^UMAP.*lam%s', lambda),
        function(pattern) {
            reducedDimNames(spe)[grep(pattern, reducedDimNames(spe))]
        }
    ) |>
    unname()

#   Plot clusters and UMAP
for (i in seq_len(length(lambda))) {
    #   First plot the clusters spatially
    p = vis_clus(
            spe, clustervar = cluster_names[i], is_stitched = TRUE,
            point_size = 1, spatial = FALSE
        ) +
        guides(fill = guide_legend(override.aes = list(size = 4)))
    png(
        file.path(plot_dir, sprintf('clusters_lambda%s.png', lambda[i])),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()

    #   First plot the clusters spatially
    p = plotReducedDim(
            spe, dimred = rd_names[i], point_size = 0.6,
            colour_by = cluster_names[i]
        ) +
        theme_bw(base_size = 20) +
        guides(color = guide_legend(override.aes = list(size = 4)))
    png(
        file.path(plot_dir, sprintf('UMAP_lambda%s.png', lambda[i])),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

#   Export clusters to CSV
cluster_df = colData(spe) |>
    as_tibble() |>
    mutate(key = colnames(spe)) |>
    select(key, any_of(cluster_names))
colnames(cluster_df) = c('key', sprintf('banksy_lambda%s', lambda))

write_csv(cluster_df, out_path)

session_info()
