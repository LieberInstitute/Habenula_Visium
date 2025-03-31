library(here)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(Banksy)
library(harmony)
library(cowplot)
library(scater)
library(tidyverse)
library(getopt)

lambda = 0.2

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'harmony_debugging',
    sprintf('spe_rewrite_lambda%s', sub('\\.', '_', as.character(lambda)))
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'harmony_debugging',
    'spe_rewrite2.rds'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'harmony_debugging',
    paste0('lambda', sub('\\.', '_', as.character(lambda)))
)

random_seed = 0

set.seed(random_seed)

spe = loadHDF5SummarizedExperiment(spe_dir)
sample_ids = unique(spe$sample_id)

#   Trim object to make lightweight
assays(spe) = list()
reducedDims(spe)$PCA = NULL
reducedDims(spe)$PCA_M1_lam0.2 = NULL
reducedDims(spe)$HARMONY = NULL
reducedDims(spe)$UMAP_HARMONY = NULL

################################################################################
#   Run Harmony on embedding, with UMAP after
################################################################################

message(Sys.time(), " | Running Harmony...")
pdf(file.path(plot_dir, "harmony_convergence_rewrite2.pdf"))
spe = RunHarmony(
    spe, group.by.vars = "sample_id", reduction.use = 'PCA_M1_lam0.2',
    plot_convergence = TRUE, kmeans_init_nstart = 20, kmeans_init_iter_max = 100
)
dev.off()

message(Sys.time(), ' | Running UMAP on Harmony-corrected embedding')
spe = runBanksyUMAP(
    spe,  dimred = "HARMONY", use_agf = TRUE, lambda = lambda,
    seed = random_seed
)

p = plot_grid(
    plotReducedDim(
            spe, sprintf("UMAP_M1_lam%s", lambda), point_size = 0.6,
            point_alpha = 0.5, color_by = "sample_id"
        ) +
        theme(legend.position = "none"),
    plotReducedDim(
            spe, "UMAP_HARMONY", point_size = 0.6, point_alpha = 0.5,
            color_by = "sample_id"
        ) +
       guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))),
    nrow = 1,
    rel_widths = c(1, 1.2)
)
png(file.path(plot_dir, 'harmony_umap_rewrite2.png'), width = 1500, height = 750)
print(p)
dev.off()

################################################################################
#   Save SpatialExperiment
################################################################################

message(Sys.time(), ' | Saving full SPE object')
saveRDS(spe, out_path)

session_info()
