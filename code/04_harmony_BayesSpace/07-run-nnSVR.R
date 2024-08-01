library(spatialLIBD)
library(SpatialExperiment)
library(here)
library(tidyverse)
library(jaffelab)
library(sessioninfo)
library(spatialNAcUtils)
library(HDF5Array)
library(nnSVG)
library(scran)
library(scuttle)

spe_dir = here(
  'processed-data', '05_harmony_BayesSpace', 'spe_filtered_hdf5'
)

message(Sys.time(), ' | Loading SpatialExperiment')
spe = loadHDF5SummarizedExperiment(spe_dir)
sample_id = unique(spe$sample_id)[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
out_path = here(
  'processed-data', '05_harmony_BayesSpace', 'nnSVG_out',
  paste0(sample_id, '.csv')
)

set.seed(0)
dir.create(dirname(out_path), showWarnings = FALSE)

#-------------------------------------------------------------------------------
#   Subset to this sample, filter lowly expressed and mitochondrial genes, and
#   take spots with at least some nonzero counts. Also drop overlapping spots
#-------------------------------------------------------------------------------

message(Sys.time(), ' | Filtering genes, spots, and to ', sample_id)
spe = spe[, (spe$sample_id == sample_id) & !spe$exclude_overlapping]

#   Bring fully into memory to speed up computations later
assays(spe) = list(counts = assays(spe)$counts)
spe = realize(spe)

spe = filter_genes(
  spe, 
  filter_genes_ncounts = 3, 
  filter_genes_pcspots = 0.5, 
  filter_mito = TRUE
)
spe = spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]
message("Dimensions of spe after filtering:")
print(dim(spe))

#-------------------------------------------------------------------------------
#   Recompute logcounts (library-size normalization as recommended in
#   https://bioconductor.org/packages/release/bioc/vignettes/nnSVG/inst/doc/nnSVG.html)
#-------------------------------------------------------------------------------

message(Sys.time(), ' | Re-computing logcounts')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)

#-------------------------------------------------------------------------------
#   Run nnSVG and export results
#-------------------------------------------------------------------------------

message(Sys.time(), ' | Running nnSVG')
spe = nnSVG(spe)

message(Sys.time(), ' | Exporting results')
write.csv(rowData(spe), out_path, row.names = FALSE, quote = FALSE)

session_info()