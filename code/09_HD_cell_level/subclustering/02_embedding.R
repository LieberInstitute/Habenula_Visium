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
    'processed-data', '09_HD_cell_level', 'subclustering',
    sprintf('spe_banksy_lambda%s', sub('\\.', '_', as.character(lambda)))
)
svg_path = here(
    'processed-data', '10_HD_bin_level', 'nnSVG_out', 'merged_SVGs.txt'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_2',
    'leiden_res0_9.csv'
)

random_seed = 0

set.seed(random_seed)

spe = loadHDF5SummarizedExperiment(spe_dir)
spe$exclude_overlapping = FALSE

#   Read in clustering results to determine which cells are habenula
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
is_habenula = cluster_df$banksy_lambda0_2[match(spe$key, cluster_df$key)] == 1

#   Subset to SVGs to avoid exceeding maximum number
#   of rows in a data table created internally. See
#   https://github.com/prabhakarlab/Banksy/issues/38#issuecomment-2310220881 and
#   https://github.com/prabhakarlab/Banksy_py/issues/12#issuecomment-2268114768.
#   Also subset to habenula cells
spe = spe[readLines(svg_path), is_habenula]


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
