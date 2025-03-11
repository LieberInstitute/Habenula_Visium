library(sessioninfo)
library(here)
library(HDF5Array)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocSingular)
library(BiocParallel)

spe_in_dir = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'spe_raw'
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'nnSVG_out', 'merged_SVGs.txt'
)
num_pcs = 50

set.seed(0)

spe <- loadHDF5SummarizedExperiment(spe_in_dir)

#   Perform PCA (subsetting by SVGs). Use IrlbaParam() for speed and memory,
#   inspired by https://pachterlab.github.io/voyager/articles/vig6_merfish.html#pca-for-larger-datasets
message(Sys.time(), " | Running PCA...")
spe = runPCA(
    spe, subset_row = readLines(svg_path), ncomponents = num_pcs,
    BSPARAM = IrlbaParam()
)

#   Save PCs (in place) 
message(Sys.time(), " | Savings PCs...")
quickResaveHDF5SummarizedExperiment(spe)

session_info()
