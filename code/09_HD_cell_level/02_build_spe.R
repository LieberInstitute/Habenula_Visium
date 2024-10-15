library(Giotto)
library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(BiocParallel)
library(scran)
library(spatialLIBD)
library(sessioninfo)

sample_id = 'H1-W369TJK_D1_9090'
ad_in_path = here(
    'processed-data', '09_HD_cell_level', sprintf('%s.h5ad', sample_id)
)
spe_bin_dir = here('processed-data', '10_HD_bin_level', 'spe_raw')
spe_norm_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
spe_raw_dir = here('processed-data', '09_HD_cell_level', 'spe_raw')

num_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

################################################################################
#   Form basic SPE mostly from the bin2cell-output AnnData
################################################################################

message(Sys.time(), " - Forming basic SPE from the bin2cell-output AnnData...")

#   Convert to Giotto then SpatialExperiment
g = anndataToGiotto(anndata_path = ad_in_path)
spe = giottoToSpatialExperiment(g)[[1]]

#   Fix sample ID
spe$sample_id = sample_id

#   Use more standard names for spatialCoords and assays
colnames(spatialCoords(spe)) = c("pxl_col_in_fullres", "pxl_row_in_fullres")
names(assays(spe)) = "counts"

#   Fix negation of spatialCoords introduced during conversion here:
#   https://github.com/drieslab/GiottoClass/blob/73e36b94e7782499ee28dd0f8caf169958c02f26/inst/python/ad2g.py#L170
spatialCoords(spe)[, 'pxl_row_in_fullres'] = -1 * spatialCoords(spe)[, 'pxl_row_in_fullres']

#   Transfer over image-related data from the bin-level object (so far, it
#   doesn't seem there is a simpler method for retaining this data when
#   converting from AnnData). Also transfer rowData
spe_bin = loadHDF5SummarizedExperiment(spe_bin_dir)
imgData(spe) = imgData(spe_bin)

stopifnot(all(rownames(spe) %in% rownames(spe_bin)))
rowData(spe) = rowData(spe_bin[rownames(spe),])

message(Sys.time(), " - Adding spatialLIBD metrics...")

#   Add QC-related variables normally added by spatialLIBD. Taken from
#   https://github.com/LieberInstitute/spatialLIBD/blob/c82c789d8538fe52e90d33af852a92d23b2368c4/R/read10xVisiumWrapper.R#L124-L129
#   Note that seqnames(spe) is incorrectly defined, but no 'seqnames<-'
#   method exists in SpatialExperiment, hence the workaround when computing
#   'is_mito'
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)
rowData(spe)$gene_search <- paste0(
    rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
)
is_mito <- which(seqnames(spe_bin[rownames(spe),]) == 'chrM')
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi
spe$ManualAnnotation <- "NA"

rm(spe_bin)
gc()

#   Save now to allow assays to become HDF5-backed in hopes of driving memory
#   down during log-normalization
message(Sys.time(), " - Saving raw SPE to move assays to HDF5")
spe <- saveHDF5SummarizedExperiment(
    spe, dir = spe_raw_dir, replace = TRUE, as.sparse = TRUE
)

################################################################################
#   Filter, log normalize, and save
################################################################################

#   Filter SPE: drop cells with 0 counts for all genes, and drop genes with 0
#   counts in every cell
message(Sys.time(), " - Filtering genes and spots")
spe <- spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]

message(Sys.time(), " - Running quickCluster()")
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores)
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
