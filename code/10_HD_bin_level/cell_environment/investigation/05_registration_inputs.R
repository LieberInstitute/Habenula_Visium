#   The extracellular registration stats are indeed highly noisy (there are
#   almost no significant markers for any cluster at any k). Check if the data
#   looks ok in low-dimensional space compared to cellular data

library(here)
library(SpatialExperiment)
library(sessioninfo)
library(tidyverse)
library(scater)
library(BiocSingular)
library(BiocParallel)

k = 10
dataset = c('extra', 'all')[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

if (dataset == 'extra') {
    spe_path = here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
        'registration', 'pseudobulk_spe', sprintf('%s.rds', k)
    )
    cluster_var = 'ficture_cluster'
} else {
    spe_path = here(
        'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
        'registration', 'pseudobulk_spe', 'cleaning_y', sprintf('%d.rds', k)
    )
    cluster_var = 'ficture'
}
svg_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out',
    'merged_SVGs.txt'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment', 
    'investigation', 'reduced_dims'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
dir.create(plot_dir, showWarnings = FALSE)

spe = readRDS(spe_path)

#   Subset to SVGs for speed. We aren't guaranteed SVGs provide good signal in
#   the extracellular bins, but we need to reduce runtime somehow
svg = readLines(svg_path)
stopifnot(all(svg %in% rownames(spe)))

spe = runPCA(
    spe, ncomponents = 10, subset_row = svg, BSPARAM = IrlbaParam(),
    BPPARAM = MulticoreParam(num_cores)
)

p = plotReducedDim(spe, dimred = "PCA", colour_by = cluster_var)
pdf(file.path(plot_dir, sprintf('PCA_factor_k%d_%s.pdf', k, dataset)))
print(p)
dev.off()

p = plotReducedDim(spe, dimred = "PCA", colour_by = "sample_id")
pdf(file.path(plot_dir, sprintf('PCA_sample_id_k%d_%s.pdf', k, dataset)))
print(p)
dev.off()

session_info()
