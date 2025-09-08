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

#   Corresponding to "cell typing" and "domain segmentation"
lambda = c(0.2, 0.8)[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
out_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    sprintf('spe_banksy_lambda%s', sub('\\.', '_', as.character(lambda)))
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'nnSVG_out',
    'merged_SVGs.txt'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'banksy',
    paste0('lambda', sub('\\.', '_', as.character(lambda)))
)

random_seed = 0
buffer_prop = 0.5

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(random_seed)

#   Load and subset to SVGs to avoid exceeding maximum number
#   of rows in a data table created internally. See
#   https://github.com/prabhakarlab/Banksy/issues/38#issuecomment-2310220881 and
#   https://github.com/prabhakarlab/Banksy_py/issues/12#issuecomment-2268114768
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[readLines(svg_path),]
spe$exclude_overlapping = FALSE

sample_ids = unique(spe$sample_id)

################################################################################
#   Stagger spatial coordinates to fit each sample in a unique range
################################################################################

message(Sys.time(), ' | Staggering spatial coordinates')
coords = spatialCoords(spe) |>
    as_tibble() |>
    mutate(sample_id = factor(spe$sample_id)) |>
    rename(sdimx = pxl_col_in_fullres, sdimy = pxl_row_in_fullres)

#   Find a range of X values slightly larger than any particular sample
x_size = coords |>
    group_by(sample_id) |>
    summarize(x_diff = max(sdimx) - min(sdimx)) |>
    pull(x_diff) |>
    max()
x_size = (1 + buffer_prop) * x_size

#   Separate samples by placing each sample into the same Y range and adjacent
#   X ranges
spatialCoords(spe) = coords |>
    group_by(sample_id) |>
    mutate(
        sdimx = sdimx - min(sdimx) + x_size * (match(cur_group()$sample_id, sample_ids) - 1),
        sdimy = sdimy - min(sdimy)
    ) |>
    ungroup() |>
    select(sdimx, sdimy) |>
    as.matrix()

################################################################################
#   Compute Banksy embedding
################################################################################

message(Sys.time(), ' | Running computeBanksy on full dataset')
spe = computeBanksy(
    spe, assay_name = "logcounts", compute_agf = TRUE, seed = random_seed
)

message(Sys.time(), ' | Running PCA on embedding')
spe = runBanksyPCA(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed
)

################################################################################
#   Run Harmony on embedding, with UMAP before and after
################################################################################

message(Sys.time(), ' | Running UMAP on embedding')
spe = runBanksyUMAP(
    spe, use_agf = TRUE, lambda = lambda, seed = random_seed
)

message(Sys.time(), " | Running Harmony...")
reducedDims(spe)$PCA = reducedDims(spe)[[sprintf('PCA_M1_lam%s', lambda)]]
reducedDims(spe)[[sprintf('PCA_M1_lam%s', lambda)]] = NULL
pdf(file.path(plot_dir, "harmony_convergence.pdf"))
spe = RunHarmony(
    spe, group.by.vars = "sample_id", plot_convergence = TRUE
)
dev.off()

message(Sys.time(), ' | Running UMAP on Harmony-corrected embedding')
spe = runBanksyUMAP(
    spe,  dimred = "HARMONY", use_agf = TRUE, lambda = lambda,
    seed = random_seed
)

################################################################################
#   Explore effect of Harmony on UMAP
################################################################################

spe$lot = ifelse(
    spe$sample_id == 'H1-W369TJK_D1_9090', 'Lot 1', 'Lot 2'
)

#   All samples together, colored by sample ID and lot (separate plots)
for (color_var in c('sample_id', 'lot')) {
    p = plot_grid(
        plotReducedDim(
                spe, sprintf("UMAP_M1_lam%s", lambda), point_size = 0.6,
                point_alpha = 0.5, color_by = color_var
            ) +
            theme_classic(base_size = 18) +
            theme(legend.position = "none"),
        plotReducedDim(
                spe, "UMAP_HARMONY", point_size = 0.6, point_alpha = 0.5,
                color_by = color_var
            ) +
            theme_classic(base_size = 18) +
            guides(
                color = guide_legend(override.aes = list(size = 4, alpha = 1))
            ),
            nrow = 1,
            rel_widths = c(1, 1.2)
    )
    png(
        file.path(plot_dir, sprintf('harmony_umap_%s_together.png', color_var)),
        width = 1200, height = 600
    )
    print(p)
    dev.off()
}

#   Faceted by sample, colored by lot
p = plot_grid(
    plotReducedDim(
            spe, sprintf("UMAP_M1_lam%s", lambda), point_size = 0.6,
            point_alpha = 0.5, color_by = 'lot'
        ) +
        facet_wrap(~ spe$sample_id, nrow = 1) +
        theme_bw(base_size = 18) +
        guides(
            color = guide_legend(override.aes = list(size = 4, alpha = 1))
        ),
    plotReducedDim(
            spe, "UMAP_HARMONY", point_size = 0.6, point_alpha = 0.5,
            color_by = 'lot'
        ) +
        facet_wrap(~ spe$sample_id, nrow = 1) +
        theme_bw(base_size = 18) +
        guides(
            color = guide_legend(override.aes = list(size = 4, alpha = 1))
        ),
        nrow = 2
)
png(
    file.path(plot_dir, 'harmony_umap_apart.png'),
    width = 1500, height = 750
)
print(p)
dev.off()

message(Sys.time(), ' | Saving full SPE object')
saveHDF5SummarizedExperiment(
    spe, dir = out_dir, replace = TRUE, as.sparse = TRUE
)

session_info()
