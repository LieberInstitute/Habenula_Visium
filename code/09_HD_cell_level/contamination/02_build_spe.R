#   Take the sample-specific AnnDatas from bin2cell and produce a basic
#   dataset-wide SpatialExperiment (bring into R and merge)

library(Giotto)
library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)
library(qs2)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = read_csv(sample_info_path)
sample_ids = sample_info$sample_id

ad_in_paths = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'bin2cell', sprintf('%s.h5ad', sample_ids)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'spe', 'raw.qs2'
)

dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   Given one sample ID and the path to a single-sample AnnData from bin2cell,
#   return a single-sample SpatialExperiment formed from the AnnData
anndata_to_spe = function(sample_id, ad_in_path) {
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

    gc()
    return(spe)
}

################################################################################
#   Main
################################################################################

#   Individually build single-sample SPEs from the individual AnnDatas, then
#   merge
spe_list = list()
for (i in seq_len(length(sample_ids))) {
    spe_list[[sample_ids[[i]]]] = anndata_to_spe(
        sample_ids[[i]], ad_in_paths[[i]]
    )
}

message(Sys.time(), " - Merging SPEs across samples")
gene_sets = unname(lapply(spe_list, rownames))
stopifnot(length(unique(gene_sets)) == 1)
spe = do.call(cbind, spe_list)

#   Save
message(Sys.time(), " - Saving raw SPE")
qs_save(spe, out_path)

message("Memory usage:")
gc()
session_info()
