library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocParallel)

spe_in_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
hvg_out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'HVGs.txt'
)

num_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
set.seed(0)

spe <- loadHDF5SummarizedExperiment(spe_in_dir)

#   Get top 10% highly variable genes
## From
## https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#3_Variance_modelling
message(Sys.time(), " | Computing HVGs...")
dec <- modelGeneVar(
    spe,  block = spe$sample_id, BPPARAM = MulticoreParam(num_cores)
)
top_hvgs <- getTopHVGs(dec, prop = 0.1)

#   Save HVGs
message(Sys.time(), " | Savings HVGs...")
writeLines(top_hvgs, con = hvg_out_path)

session_info()
