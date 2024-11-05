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
out_dir = here('processed-data', '09_HD_cell_level', 'spe_banksy')

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

#   Get the names of the cluster columns corresponding to each lambda value
rd_names = sapply(
        sprintf('^UMAP.*lam%s', lambda),
        function(pattern) {
            reducedDimNames(spe)[grep(pattern, reducedDimNames(spe))]
        }
    ) |>
    unname()

#   Plot UMAP for each lambda
for (i in seq_len(length(lambda))) {
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

message(Sys.time(), ' | Saving Banksy embedding (in full SPE object)')
saveHDF5SummarizedExperiment(
    spe, dir = out_dir, replace = TRUE, as.sparse = TRUE
)

session_info()
