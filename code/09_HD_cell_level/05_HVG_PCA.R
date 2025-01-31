library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocSingular)
library(BiocParallel)

spe_in_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
hvg_out_path = here('processed-data', '09_HD_cell_level', 'HVGs.txt')
num_pcs = 50

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

#   Perform PCA (subsetting by HVGs). Use IrlbaParam() for speed and memory,
#   inspired by https://pachterlab.github.io/voyager/articles/vig6_merfish.html#pca-for-larger-datasets
message(Sys.time(), " | Running PCA...")
spe = runPCA(
    spe, subset_row = top_hvgs, ncomponents = num_pcs, BSPARAM = IrlbaParam()
)

#   Save PCs (in place) and HVGs
message(Sys.time(), " | Savings PCs and HVGs...")
writeLines(top_hvgs, con = hvg_out_path)
quickResaveHDF5SummarizedExperiment(spe)

session_info()
