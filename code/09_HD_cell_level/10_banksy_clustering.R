library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(Banksy)
library(cowplot)
library(scater)

k = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_banksy')
out_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', sprintf('k%s.csv', k)
)
plot_dir = here('plots', '09_HD_cell_level', 'banksy', sprintf('k%s', k))

dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Infer value(s) of lambda from the reducedDimNames
lambda = reducedDimNames(spe)[grep('^UMAP.*lam', reducedDimNames(spe))] |>
    str_extract('lam(.*)$', group = 1) |>
    as.numeric()

message(Sys.time(), ' | Performing clustering')
spe = clusterBanksy(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed, algo = "kmeans",
    kmeans.centers = k
)
spe = connectClusters(spe)

#   Get the names of the cluster and UMAP columns corresponding to each lambda
#   value
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

#   Plot clusters and colored UMAP for each lambda
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

    #   Then UMAP colored by cluster
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
