library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(BiocParallel)
library(scran)
library(spatialLIBD)
library(sessioninfo)

spe_norm_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
spe_raw_dir = here('processed-data', '09_HD_cell_level', 'spe_raw')
plot_dir = here('plots', '09_HD_cell_level')

num_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

################################################################################
#   Drop empty genes and do cell-level QC
################################################################################

spe = loadHDF5SummarizedExperiment(spe_bin_dir)

#   Filter SPE: drop cells with 0 counts for all genes, and drop genes with 0
#   counts in every cell
message(Sys.time(), " - Filtering genes and spots")
spe <- spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]

#   There are no problematic-looking cells based on 'sum_umi' or 'sum_gene'
message("No outliers in 'sum_umi' or 'sum_gene':")
table(isOutlier(spe$sum_umi, type = "lower"))
table(isOutlier(spe$sum_gene, type = "lower"))

#   Explore cells with high mitochondrial ratios
spe$high_expr_chrM_ratio = isOutlier(spe$expr_chrM_ratio, type = "higher")
message(
    sprintf(
        "Dropping %s of %s cells (%s%%) with high mitochondrial ratios (with these stats):",
        sum(spe$high_expr_chrM_ratio),
        ncol(spe),
        round(100 * sum(spe$high_expr_chrM_ratio) / ncol(spe), 1)
    )
)
summary(spe$expr_chrM_ratio[spe$high_expr_chrM_ratio])

#   Plot bad cells based on mitochondrial ratio
p <- vis_clus(
    spe, clustervar = 'high_expr_chrM_ratio', is_stitched = TRUE, point_size = 1,
    spatial = FALSE
)
png(file.path(plot_dir, "high_expr_chrM_ratio.png"), width = 1500, height = 1500)
print(p)
dev.off()

#   Drop cells with high mitochondrial ratios
spe = spe[, !spe$high_expr_chrM_ratio]
spe$high_expr_chrM_ratio = NULL

################################################################################
#   Log normalization
################################################################################

message(Sys.time(), " - Running quickCluster()")
spe$scran_quick_cluster <- quickCluster(
    spe, BPPARAM = MulticoreParam(num_cores)
)

message("Quick cluster table:")
table(spe$scran_quick_cluster)

message("sizeFactors() summary:")
summary(sizeFactors(spe))

message(Sys.time(), " - Running logNormCounts()")
spe <- logNormCounts(spe)

#   Save normalized object
message(Sys.time(), " - Saving normalized SPE")
spe <- saveHDF5SummarizedExperiment(
    spe, dir = spe_norm_dir, replace = TRUE, as.sparse = TRUE
)

session_info()