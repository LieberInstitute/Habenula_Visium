library(SpatialExperiment)
library(here)
library(tidyverse)
library(sessioninfo)
library(HDF5Array)
library(Matrix)
library(nnSVG)

sample_id = 'H1-W369TJK_D1_9090'
spe_norm_dir = here(
    'processed-data', '10_HD_bin_level', 'rasterized',
    sprintf('spe_%s_1x_standard_res', sample_id)
)
out_path <- here(
    "processed-data", '10_HD_bin_level', "nnSVG_out",
    paste0(sample_id, ".csv")
)

set.seed(0)
dir.create(dirname(out_path), showWarnings = FALSE)

message(Sys.time(), " | Loading SpatialExperiment")
spe <- loadHDF5SummarizedExperiment(spe_norm_dir)
names(assays(spe)) = 'logcounts'

#-------------------------------------------------------------------------------
#   Filter lowly expressed and mitochondrial genes
#-------------------------------------------------------------------------------

#   Rather than rasterizing the raw counts and using it to invoke filter_genes()
#   with similar parameters to Visium standard, we'll manually filter
#   mitochondrial genes and apply an expression cutoff the selects a similar
#   absolute number of genes (~3000) as previous Visium standard datasets,
#   using the rasterized logcounts

message(Sys.time(), " | Filtering genes")

#   Drop mitochondrial genes
drop_vec = grepl('(^MT-)|(^mt-)', rowData(spe)$gene_name)
message(sprintf("Dropping %s mitochondrial genes.", length(which(drop_vec))))
spe = spe[!drop_vec,]

#   Drop genes where 70% or less of spots have >0 logcounts
spe = spe[rowMeans(assays(spe)$logcounts > 0) > 0.7,]

message("Dimensions of spe after filtering:")
print(dim(spe))

#-------------------------------------------------------------------------------
#   Run nnSVG and export results
#-------------------------------------------------------------------------------

message(Sys.time(), " | Running nnSVG")
spe <- nnSVG(spe)

message(Sys.time(), " | Exporting results")
write_csv(as_tibble(rowData(spe)), out_path)

session_info()
