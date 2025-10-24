#   Take the sample-specific AnnDatas from bin2cell and produce a basic
#   dataset-wide SpatialExperiment (bring into R and merge)

library(Giotto)
library(here)
library(tidyverse)
library(SpatialExperiment)
library(BiocParallel)
library(scran)
library(sessioninfo)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = read_csv(sample_info_path)
sample_ids = sample_info$sample_id

ad_in_paths = here(
    'processed-data', '09_HD_cell_level', 'new_samples',
    sprintf('%s.h5ad', sample_ids)
)
spe_bin_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'spe_raw.rds'
)
spe_raw_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_raw.rds'
)

################################################################################
#   Functions
################################################################################

#   Given one sample ID, the path to a single-sample AnnData from bin2cell, and
#   a potentially multi-sample bin-level SpatialExperiment, return a
#   single-sample SpatialExperiment formed from the AnnData
anndata_to_spe = function(sample_id, ad_in_path, spe_bin) {
    message(
        Sys.time(),
        sprintf(
            " - Forming basic SPE from the bin2cell-output AnnData for sample %s...",
            sample_id
        )
    )

    #   Convert to Giotto then SpatialExperiment
    g = anndataToGiotto(anndata_path = ad_in_path)
    spe = giottoToSpatialExperiment(g)[[1]]

    #   Fix sample ID and key
    spe$sample_id = sample_id
    spe$key = paste(colnames(spe), sample_id, sep = "_")
    colnames(spe) = spe$key

    #   Use more standard names for spatialCoords and assays
    colnames(spatialCoords(spe)) = c("pxl_col_in_fullres", "pxl_row_in_fullres")
    names(assays(spe)) = "counts"

    #   Fix negation of spatialCoords introduced during conversion here:
    #   https://github.com/drieslab/GiottoClass/blob/73e36b94e7782499ee28dd0f8caf169958c02f26/inst/python/ad2g.py#L170
    spatialCoords(spe)[, 'pxl_row_in_fullres'] = -1 * spatialCoords(spe)[, 'pxl_row_in_fullres']

    #   Transfer over image-related data from the bin-level object (so far, it
    #   doesn't seem there is a simpler method for retaining this data when
    #   converting from AnnData). Also transfer rowData. For unclear reasons,
    #   some genes are dropped during creation of the bin-level object, and so
    #   this step may drop genes in the cell-level object where rowData doesn't
    #   exist
    imgData(spe) = imgData(spe_bin)

    stopifnot(all(rownames(spe_bin) %in% rownames(spe)))
    num_dropped = length(setdiff(rownames(spe), rownames(spe_bin)))
    if (num_dropped > 0) {
        warning(
            sprintf(
                "Dropping %s of %s genes for %s not present in 'spe_bin'...",
                num_dropped, nrow(spe), sample_id
            )
        )
    }
    spe = spe[intersect(rownames(spe), rownames(spe_bin))]
    rowData(spe) = rowData(spe_bin[rownames(spe),])

    message(
        Sys.time(),
        sprintf(" - Adding spatialLIBD metrics for sample %s...", sample_id)
    )

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

    gc()
    return(spe)
}

################################################################################
#   Main
################################################################################

spe_bin = readRDS(spe_bin_path)
stopifnot(setequal(sample_ids, unique(spe_bin$sample_id)))

#   Individually build single-sample SPEs from the individual AnnDatas, then
#   merge
spe_list = list()
for (i in seq_len(length(sample_ids))) {
    spe_list[[sample_ids[[i]]]] = anndata_to_spe(
        sample_ids[[i]], ad_in_paths[[i]],
        spe_bin[, spe_bin$sample_id == sample_ids[[i]]]
    )
}

message(Sys.time(), " - Merging SPEs across samples")
gene_sets = unname(lapply(spe_list, rownames))
stopifnot(length(unique(gene_sets)) == 1)
spe = do.call(cbind, spe_list)

#   Save
message(Sys.time(), " - Saving raw SPE")
saveRDS(spe, spe_raw_path)

message("Memory usage:")
gc()
session_info()
