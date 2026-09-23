#   We already have the extracellular data by cell in the AnnData for LIANA+.
#   Just import into R as a SpatialExperiment

library(Giotto)
library(here)
library(tidyverse)
library(SpatialExperiment)
library(BiocParallel)
library(scran)
library(sessioninfo)
library(qs2)

ad_in_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'adata', 'liana_ready', 'extracellular.h5ad'
)
spe_bin_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
spe_out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.qs2'
)

dir.create(dirname(spe_out_path), recursive = TRUE, showWarnings = FALSE)

################################################################################
#   Main
################################################################################

spe_bin = readRDS(spe_bin_path)

#   Convert to Giotto then SpatialExperiment
g = anndataToGiotto(anndata_path = ad_in_path)
spe = giottoToSpatialExperiment(g)[[1]]

spe = spe[, spe$cell_type != 'Drop']

#   Use more standard names for spatialCoords and assays
colnames(spatialCoords(spe)) = c("pxl_col_in_fullres", "pxl_row_in_fullres")
names(assays(spe)) = "counts"

#   Fix negation of spatialCoords introduced during conversion here:
#   https://github.com/drieslab/GiottoClass/blob/73e36b94e7782499ee28dd0f8caf169958c02f26/inst/python/ad2g.py#L170
spatialCoords(spe)[, 'pxl_row_in_fullres'] = -1 * spatialCoords(spe)[, 'pxl_row_in_fullres']

#   Circumvent hard restrictions on changing sample_id
temp = spatialCoords(spe)
sce = as(spe, "SingleCellExperiment")
sce$sample_id = sce$tissue_section
spe = as(sce, "SpatialExperiment")
spatialCoords(spe) = temp

#   Transfer over image-related data from the bin-level object (so far, it
#   doesn't seem there is a simpler method for retaining this data when
#   converting from AnnData). Also transfer rowData. For unclear reasons,
#   some genes are dropped during creation of the bin-level object, and so
#   this step may drop genes in the cell-level object where rowData doesn't
#   exist
stopifnot(setequal(spe$sample_id, spe_bin$sample_id))
imgData(spe) = imgData(spe_bin)

num_dropped = length(setdiff(rownames(spe), rownames(spe_bin)))
if (num_dropped > 0) {
    warning(
        sprintf(
            "Dropping %s of %s genes not present in 'spe_bin'...",
            num_dropped, nrow(spe)
        )
    )
}
spe = spe[intersect(rownames(spe), rownames(spe_bin))]
rowData(spe) = rowData(spe_bin[rownames(spe),])

#   Add QC-related variables normally added by spatialLIBD. Taken from
#   https://github.com/LieberInstitute/spatialLIBD/blob/c82c789d8538fe52e90d33af852a92d23b2368c4/R/read10xVisiumWrapper.R#L124-L129
#   Note that seqnames(spe) is incorrectly defined, but no 'seqnames<-'
#   method exists in SpatialExperiment, hence the workaround when computing
#   'is_mito'
spe$sum_umi <- colSums(counts(spe))
spe$sum_gene <- colSums(counts(spe) > 0)
is_mito <- which(seqnames(spe_bin[rownames(spe),]) == 'chrM')
spe$expr_chrM <- colSums(counts(spe)[is_mito, , drop = FALSE])
spe$expr_chrM_ratio <- spe$expr_chrM / spe$sum_umi
spe$ManualAnnotation <- "NA"

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large)
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)

#   Save
message(Sys.time(), " - Saving raw SPE")
qs_save(spe, spe_out_path)

message("Memory usage:")
gc()
session_info()
