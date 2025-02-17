library(here)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(Banksy)
library(getopt)

#   Corresponding to "cell typing" and "domain segmentation"
lambda = c(0.2, 0.8)[as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
out_dir = here(
    'processed-data', '09_HD_cell_level',
    sprintf('spe_banksy_%s', sub('\\.', '_', as.character(lambda)))
)
svg_path = here('TODO')

random_seed = 0

set.seed(random_seed)

#   Load and subset to SVGs to avoid exceeding maximum number
#   of rows in a data table created internally. See
#   https://github.com/prabhakarlab/Banksy/issues/38#issuecomment-2310220881 and
#   https://github.com/prabhakarlab/Banksy_py/issues/12#issuecomment-2268114768
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[readLines(svg_path),]
spe$exclude_overlapping = FALSE

#   Split big SPE into a list (one object per sample). Then run 'computeBanksy'
#   and merge again
message(Sys.time(), ' | Running computeBanksy on each sample')
spe = do.call(
    cbind,
    lapply(
        unique(spe$sample_id),
        function(x) {
            computeBanksy(
                spe[, spe$sample_id == x], assay_name = "logcounts",
                compute_agf = TRUE, seed = random_seed
            )
        }
    )
)

message(Sys.time(), ' | Running PCA on embedding')
spe = runBanksyPCA(
    spe, use_agf = TRUE, lambda = lambda, group = 'sample_id',
    seed = random_seed
)

message(Sys.time(), ' | Running UMAP on embedding')
spe = runBanksyUMAP(
    spe, use_agf = TRUE, lambda = lambda, group = 'sample_id',
    seed = random_seed
)

message(Sys.time(), ' | Saving Banksy embedding (in full SPE object)')
saveHDF5SummarizedExperiment(
    spe, dir = out_dir, replace = TRUE, as.sparse = TRUE
)

session_info()
