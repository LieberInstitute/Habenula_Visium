library(sessioninfo)
library(here)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocSingular)
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

set.seed(0)

spe <- readRDS(spe_in_path)

#   Perform PCA (subsetting by SVGs). Use IrlbaParam() for speed and memory,
#   inspired by https://pachterlab.github.io/voyager/articles/vig6_merfish.html#pca-for-larger-datasets
message(Sys.time(), " | Running PCA...")
spe = runPCA(
    spe, subset_row = readLines(svg_path), ncomponents = num_pcs,
    BSPARAM = IrlbaParam()
)
message(Sys.time(), " | Running Harmony...")
spe = RunHarmony(spe, group.by.vars = "sample_id", dims.use = "PCA")

#   Write just the Harmony embedding
reducedDims(spe)$HARMONY |>
    rownames_to_column('key') |>
    mutate(X = spatialCoords(spe)[, 1], Y = spatialCoords(spe)[, 2]) |>
    as_tibble() |>
    write_csv(out_path)

session_info()
