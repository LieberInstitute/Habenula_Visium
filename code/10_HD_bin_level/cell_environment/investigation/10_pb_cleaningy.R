#   Is there signal between clusters after regressing out sample ID from the
#   psuedobulked extracellular data? These tests show there seems to only be
#   noise left over

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(jaffelab)
library(Matrix)
library(BiocSingular)
library(BiocParallel)
library(scater)

k = 10
spe_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'pseudobulk_spe', sprintf('%s.rds', k)
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out',
    'merged_SVGs.txt'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation', 'reduced_dims'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))

spe = readRDS(spe_path)

#   Regress out sample ID
mod = with(colData(spe), model.matrix(~ sample_id))
assays(spe)$logcounts = cleaningY(assays(spe)$logcounts, mod, P = 1)

#   Subset to SVGs in an attempt to clean up signal. Note results don't change
#   much when interactively testing without subsetting to SVGs
svg = readLines(svg_path)
stopifnot(all(svg %in% rownames(spe)))

spe = runPCA(
    spe, ncomponents = 10, subset_row = svg, BSPARAM = IrlbaParam(),
    BPPARAM = MulticoreParam(num_cores)
)
spe = runUMAP(
    spe, ncomponents = 2, subset_row = svg, BPPARAM = MulticoreParam(num_cores)
)

for (dim_red in c("PCA", "UMAP")) {
    p = plotReducedDim(spe, dimred = dim_red, colour_by = "ficture_cluster")
    pdf(
        file.path(
            plot_dir,
            sprintf('pb_cleaningY_%s_factor_k%d_extra.pdf', dim_red, k)
        )
    )
    print(p)
    dev.off()

    p = plotReducedDim(spe, dimred = dim_red, colour_by = "sample_id")
    pdf(
        file.path(
            plot_dir,
            sprintf('pb_cleaningY_%s_sample_id_k%d_extra.pdf', dim_red, k)
        )
    )
    print(p)
    dev.off()
}

