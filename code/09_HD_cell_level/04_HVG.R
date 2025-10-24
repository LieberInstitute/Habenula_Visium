library(sessioninfo)
library(here)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocParallel)

spe_in_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_norm_filtered.rds'
)
hvg_out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'HVGs.txt'
)

num_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
set.seed(0)

spe <- readRDS(spe_in_path)

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
