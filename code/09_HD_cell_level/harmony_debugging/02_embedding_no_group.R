library(here)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(Banksy)
library(harmony)
library(cowplot)
library(scater)
library(getopt)

#   Corresponding to "cell typing" and "domain segmentation"
lambda = c(0.2, 0.8)[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
out_dir = here(
    'processed-data', '09_HD_cell_level', 'harmony_debugging',
    sprintf('spe_no_group_lambda%s', sub('\\.', '_', as.character(lambda)))
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'nnSVG_out', 'merged_SVGs.txt'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'harmony_debugging',
    paste0('lambda', sub('\\.', '_', as.character(lambda)))
)

random_seed = 0

dir.create(plot_dir, showWarnings = FALSE)
set.seed(random_seed)

#   Load and subset to SVGs to avoid exceeding maximum number
#   of rows in a data table created internally. See
#   https://github.com/prabhakarlab/Banksy/issues/38#issuecomment-2310220881 and
#   https://github.com/prabhakarlab/Banksy_py/issues/12#issuecomment-2268114768
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[readLines(svg_path),]
spe$exclude_overlapping = FALSE

#   Split big SPE into a list (one object per sample). Then run 'computeBanksy'
#   and merge again
message(Sys.time(), ' | Running computeBanksy on each sample')
spe = do.call(
    cbind,
    lapply(
        unique(spe$sample_id),
        function(x) {
            computeBanksy(
                spe[, spe$sample_id == x], assay_name = "logcounts",
                compute_agf = TRUE, seed = random_seed
            )
        }
    )
)

message(Sys.time(), ' | Running (non-grouped) PCA on embedding')
spe = runBanksyPCA(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed
)

message(Sys.time(), ' | Running UMAP on embedding')
spe = runBanksyUMAP(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed
)

message(Sys.time(), " | Running Harmony...")
reducedDims(spe)$PCA = reducedDims(spe)[[sprintf('PCA_M1_lam%s', lambda)]]
reducedDims(spe)[[sprintf('PCA_M1_lam%s', lambda)]] = NULL
pdf(file.path(plot_dir, "harmony_convergence_no_group.pdf"))
spe = RunHarmony(
    spe, group.by.vars = "sample_id", plot_convergence = TRUE
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
png(
    file.path(plot_dir, 'harmony_umap_no_group.png'), width = 1500, height = 750
)
print(p)
dev.off()

message(Sys.time(), ' | Saving full SPE object')
saveHDF5SummarizedExperiment(
    spe, dir = out_dir, replace = TRUE, as.sparse = TRUE
)

session_info()
