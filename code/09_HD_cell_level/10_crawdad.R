library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(crawdad)
library(rjson)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_norm_filtered.rds'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'banksy', 'lambda0_2',
    'leiden_res1_3_subset.csv'
)
cor_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'registration_banksy',
    'lambda0_2', 'cor_vs_snRNAseq_fine_subset.rds'
)
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
plot_path = here(
    'plots', '09_HD_cell_level', 'new_samples', 'crawdad', 'dot_plot_%s.pdf'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'crawdad',
    '%s_results.csv'
)
scales = c(100, 200, 500, 1000, 5000)
random_seed = 0
cor_index = 13

num_cores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
dir.create(dirname(plot_path), showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

sample_info = read_csv(sample_info_path)
task_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
sample_id = sample_info$sample_id[task_id]
spaceranger_dir = sample_info$spaceranger_dir[task_id]

spe = readRDS(spe_path)
spe = spe[, spe$sample_id == sample_id]

#   Ultimately, we'll be converting spatial coordinates to units of microns,
#   which is more interpretable than pixels
scalefactors_path = here(
    spaceranger_dir, 'outs', 'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
micron_per_px = fromJSON(file = scalefactors_path)[['microns_per_pixel']]

#   Compute a reference table matching clusters to fine cell types
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        layer_label = ifelse(
            layer_confidence == 'good', layer_label, paste0('C',cluster)
        )
    )

#   Gather spatial coordinates and Banksy clusters
cell_df = tibble(
        key = spe$key,
        x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * micron_per_px,
        y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * micron_per_px
    ) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = factor(
            anno_df$layer_label[match(banksy_lambda0_2, anno_df$cluster)]
        )
    ) |>
    select(x, y, cell_type) |>
    as.data.frame()
stopifnot(!any(is.na(cell_df$cell_type)))

pos_df = toSF(pos = select(cell_df, c(x, y)), cellTypes = cell_df$cell_type)

#   Shuffle cell-type assignments to create null background
shuffle_list = makeShuffledCells(
    pos_df, scales = scales, seed = random_seed, ncores = num_cores,
    verbose = TRUE
)

#   Main Z-score calculation for each reference-neighbor pair
results = findTrends(
    pos_df, shuffleList = shuffle_list, returnMeans = FALSE, ncores = num_cores,
    verbose = TRUE
)

#   Reformat and compute multiple-testing-corrected Z-score
results = meltResultsList(results, withPerms = TRUE)
z_sig = correctZBonferroni(results)

#   Main dot plot figure
pdf(sprintf(plot_path, sample_id))
vizColocDotplot(
    results, zSigThresh = z_sig, zScoreLimit = 2 * z_sig, dotSizes = c(2, 10)
)
dev.off()

#   Export results
write_csv(results, sprintf(out_path, sample_id))

message('Memory usage:')
gc()

session_info()
