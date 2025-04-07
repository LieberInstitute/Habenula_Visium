library(here)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(Banksy)
library(harmony)
library(cowplot)
library(scater)

lambda = 0.2
spe_dir = here(
    'processed-data', '09_HD_cell_level',
    sprintf('spe_banksy_lambda%s', sub('\\.', '_', as.character(lambda)))
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'harmony_debugging',
    sprintf('spe_temp_UMAP_lambda%s.rds', sub('\\.', '_', as.character(lambda)))
)
plot_path = here(
    'plots', '09_HD_cell_level', 'harmony_debugging', 'temp_banksy_umap.png'
)
random_seed = 0

dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(dirname(plot_path), showWarnings = FALSE)
set.seed(random_seed)

spe = loadHDF5SummarizedExperiment(spe_dir)
reducedDims(spe)$UMAP_M1_lam0.2 = NULL

message(Sys.time(), " | Running UMAP on Banksy embedding...")
spe <- runBanksyUMAP(spe, use_agf = TRUE, lambda = lambda)
message(Sys.time(), " | Running UMAP on HARMONY-corrected embedding...")
spe <- runBanksyUMAP(spe, dimred = "HARMONY")

#   Make object smaller on disk
message(Sys.time(), " | Saving SPE...")
assays(spe) = list()
saveRDS(spe, out_path)

message(Sys.time(), " | Trying to plot UMAPs...")
p = plot_grid(
    plotReducedDim(
            spe, "UMAP_M1_lam0.2", point_size = 0.6, point_alpha = 0.5,
            color_by = "sample_id"
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
png(plot_path, width = 1500, height = 750)
print(p)
dev.off()
