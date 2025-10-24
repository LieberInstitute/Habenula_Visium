#   Try running CRAWDAD as if we only had one sample ("the BayesSpace trick"),
#   by putting cells from all samples in the same space in adjacent sections,
#   separated by the largest spatial scale used in the analysis. The idea is
#   described more at https://github.com/JEFworks-Lab/CRAWDAD/issues/34 (the
#   authors hadn't responded at the time of writing this script)

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
plot_dir = here('plots', '09_HD_cell_level', 'new_samples', 'crawdad')
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'crawdad',
    'combined_results.csv'
)
scales = c(100, 200, 500, 1000, 5000)
random_seed = 0
cor_index = 13

num_cores = as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

custom_dotplot = function(result_df, z_sig, filename, color_var = 'Z') {
    p = ggplot(
            result_df, aes(x = reference, y = neighbor, color = !!sym(color_var), size = scale)
        ) +
        geom_point() +
        scale_color_gradientn(
            colors = c('blue', 'white', 'white', 'red'),
            values = rescale(
                c(min(result_df$Z), -1 * z_sig, z_sig, max(result_df$Z))
            )
        ) +
        scale_radius(
            trans = 'reverse',
            breaks = seq(
                min(result_df$scale), max(result_df$scale), length.out = 3
            ),
            range = c(2, 15)
        ) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    pdf(file.path(plot_dir, filename), width = 9)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

sample_info = read_csv(sample_info_path)
sample_ids = sample_info$sample_id[1:3]

spe = readRDS(spe_path)

#   Compute a reference table matching clusters to fine cell types
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        layer_label = ifelse(
            layer_confidence == 'good', layer_label, paste0('C',cluster)
        )
    )

#   Read in cell coordinates (in microns) for each sample
cell_df_list = list()
for (sample_id in sample_ids) {
    #   Ultimately, we'll be converting spatial coordinates to units of microns,
    #   which is more interpretable than pixels
    spaceranger_dir = sample_info$spaceranger_dir[sample_info$sample_id == sample_id]
    scalefactors_path = here(
        spaceranger_dir, 'outs', 'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
    )
    micron_per_px = fromJSON(file = scalefactors_path)[['microns_per_pixel']]

    small_spe = spe[, spe$sample_id == sample_id]

    offset = max(scales) * (match(sample_id, sample_ids) - 1)
    cell_df_list[[sample_id]] = tibble(
            key = small_spe$key,
            x = spatialCoords(small_spe)[, 'pxl_col_in_fullres'] *
                micron_per_px,
            y = spatialCoords(small_spe)[, 'pxl_row_in_fullres'] *
                micron_per_px
        ) |>
        #   Place each sample in adjacent (but spaced-apart) regions
        mutate(x = x - min(x) + offset, y = y - min(y))
}

#   Gather spatial coordinates and Banksy clusters
cell_df = do.call(rbind, cell_df_list) |>
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

#   Export results
write_csv(results, out_path)

#   Main dot plot figure
pdf(file.path(plot_dir, 'dot_plot_combined_artificial.pdf'))
vizColocDotplot(
    results, zSigThresh = z_sig, zScoreLimit = 2 * z_sig, dotSizes = c(2, 10)
)
dev.off()

#   Also plot a version visually similar to the one in 11_crawdad_plot.R
results = results |>
    #   First average Z-scores across permutations
    group_by(neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= z_sig) |>
    group_by(neighbor, reference) |>
    filter(scale == min(scale)) |>
    ungroup() |>
    #   Cap Z-score at twice the magnitude of the significance threshold
    mutate(Z = sign(Z) * pmin(abs(Z), z_sig * 2))

custom_dotplot(
    results, z_sig, filename = 'dot_plot_combined_artificial_custom.pdf'
)

message('Memory usage:')
gc()

session_info()
