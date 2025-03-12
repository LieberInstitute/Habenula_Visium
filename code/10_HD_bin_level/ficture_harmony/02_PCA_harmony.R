library(sessioninfo)
library(here)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocSingular)
library(BiocParallel)
library(harmony)
library(tidyverse)

spe_in_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'spe_raw.rds'
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'nnSVG_out', 'merged_SVGs.txt'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'harmony_embedding.csv.gz'
)
num_pcs = 50

num_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
set.seed(0)

spe <- readRDS(spe_in_path)

#   Filter raw SPE: drop bins with 0 counts for all genes, and drop genes with
#   0 counts in every bin
message(Sys.time(), ' | Filtering bins and genes...')
spe <- spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large)
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)

#   Perform PCA (subsetting by SVGs). Use IrlbaParam() for speed and memory,
#   inspired by https://pachterlab.github.io/voyager/articles/vig6_merfish.html#pca-for-larger-datasets
message(Sys.time(), " | Running PCA...")
spe = runPCA(
    spe, subset_row = readLines(svg_path), ncomponents = num_pcs,
    BPPARAM = MulticoreParam(num_cores), qBSPARAM = IrlbaParam()
)
message(Sys.time(), " | Running Harmony...")
spe = RunHarmony(
    spe, group.by.vars = "sample_id", dims.use = "PCA", ncores = num_cores
)

#   Write just the Harmony embedding
reducedDims(spe)$HARMONY |>
    rownames_to_column('key') |>
    mutate(X = spatialCoords(spe)[, 1], Y = spatialCoords(spe)[, 2]) |>
    as_tibble() |>
    write_csv(out_path)

session_info()
