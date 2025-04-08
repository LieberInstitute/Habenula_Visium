library(sessioninfo)
library(here)
library(SpatialExperiment)
library(scran)
library(scater)
library(BiocSingular)
library(BiocParallel)
library(crescendo)
library(tidyverse)

spe_in_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'spe_harmony.rds'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'harmony_embedding.csv.gz'
)
random_seed = 0

num_cores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
set.seed(random_seed)

#   Load SPE and drop logcounts to save memory
message(Sys.time(), ' | Loading SPE...')
spe <- readRDS(spe_in_path)
assays(spe)$logcounts = NULL
gc()

#   Make colnames unique
spe$key = paste(spe$sample_id, colnames(spe), sep = '_')
colnames(spe) = spe$key

message(Sys.time(), ' | Computing total UMI per cell...')
spe$nUMI = colSums(assays(spe)$counts)

message(Sys.time(), ' | Running crescendo...')
result = crescendo(
    Ycounts = assays(spe)$counts,
    meta_data = colData(spe) |> as.data.frame(),
    R = t(reducedDims(spe)$HARMONY),
    batch_var = "sample_id",
    prop = 0.01,
    seed = random_seed,
    mc.cores = num_cores,
    verbose = TRUE
)

session_info()
