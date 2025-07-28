library(here)
library(tidyverse)
library(crawdad)
library(spatialLIBD)
library(HDF5Array)
library(scales)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
result_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad',
    '%s_results.csv'
)
spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_3_subset.csv'
)
cor_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'lambda0_2', 'cor_vs_snRNAseq_fine_subset.rds'
)
cell_type_colors = c(
    Microglia = "#2F97FF", MHb.2 = "#FFA239", Other = "#DFE1DD"
)
plot_dir = here('plots', '09_HD_cell_level', 'probe_fix', 'crawdad')
cor_index = 13

dir.create(file.path(plot_dir, 'spatial_plots'), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

custom_dotplot = function(result_df, z_sig, filename) {
    p = ggplot(
            result_df, aes(x = reference, y = neighbor, color = Z, size = scale)
        ) +
        geom_point() +
        scale_color_gradientn(
            colors = c('blue', '#CECECE', '#CECECE', 'red'),
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
#   CRAWDAD-specific plots
################################################################################

sample_ids = readLines(sample_id_path)[1:3]

result_list = list()
for (sample_id in sample_ids) {
    result_list[[sample_id]] = sprintf(result_paths, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        mutate(sample_id = sample_id)

    result_list[[sample_id]]$Z_sig = correctZBonferroni(
        result_list[[sample_id]]
    )
}

result_df = do.call(rbind, result_list) |>
    #   First average Z-scores across permutations
    group_by(sample_id, neighbor, scale, reference, Z_sig) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= Z_sig) |>
    group_by(sample_id, neighbor, reference) |>
    filter(scale == min(scale)) |>
    ungroup()

#   Significance thresholds should only be determined by number of clusters,
#   which should be equal in all samples. Grab the single cutoff
z_sig = unname(unlist(lapply(result_list, function(x) x$Z_sig)))
stopifnot(length(unique(z_sig)) == 1)
z_sig = unique(z_sig)

result_df = result_df |>
    #   Retain pairs where all samples are significant, and sign of Z scores
    #   agree across samples
    group_by(reference, neighbor) |>
    filter(all(Z > 0) | all(Z < 0)) |>
    filter(n() == length(unique(sample_ids))) |>
    #   Take the mean Z-score and scale across samples
    group_by(neighbor, reference) |>
    summarize(scale = mean(scale), Z = mean(Z)) |>
    ungroup() |>
    #   Cap Z-score at twice the magnitude of the significance threshold
    mutate(Z = sign(Z) * pmin(abs(Z), z_sig * 2))

custom_dotplot(result_df, z_sig, 'dot_plot_combined.pdf')

#   Plot Z-scores vs scale for a particularly interesting cell-type pair
p = do.call(rbind, result_list) |>
    #   Average Z-scores across permutations
    group_by(sample_id, neighbor, scale, reference, Z_sig) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Improve plot appearance
    mutate(
        sample_id = paste0('Br', str_extract(sample_id, '[0-9]{4}$')),
        facet_anno = sprintf("Ref: %s\nNeighbor: %s", reference, neighbor)
    ) |>
    #   Focus on a particular pair (and its reverse)
    filter(
        ((neighbor == 'MHb.2') & (reference == 'Microglia')) |
        ((neighbor == 'Microglia') & (reference == 'MHb.2'))
    ) |>
    ggplot(aes(x = scale, y = Z, color = sample_id, group = sample_id)) +
        geom_line() +
        geom_point() +
        facet_wrap(~ facet_anno, nrow = 1) +
        theme_bw(base_size = 20) +
        labs(x = 'Scale (Microns)', color = 'Sample ID')

pdf(file.path(plot_dir, 'z_scores_MHb_microglia.pdf'), width = 10, height = 5)
print(p)
dev.off()

################################################################################
#   Spatial distribution of particular cell types
################################################################################

spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id %in% sample_ids]

#   Compute a reference table matching clusters to fine cell types
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        layer_label = ifelse(
            layer_confidence == 'good', layer_label, paste0('C',cluster)
        )
    )

#   Annotate Banksy clusters with cell type
anno_join_df = tibble(key = spe$key) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$layer_label[
            match(banksy_lambda0_2, anno_df$cluster)
        ],
        cell_type_pair = ifelse(
            cell_type %in% c("Microglia", "MHb.2"), cell_type, 'Other'
        )
    )
stopifnot(!any(is.na(anno_join_df$cell_type)))
spe$cell_type = anno_join_df$cell_type
spe$cell_type_pair = anno_join_df$cell_type_pair

for (sample_id in sample_ids) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        #   Plot the cell-type pair spatially in each sample
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = 'cell_type_pair',
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = cell_type_colors
            ) +
                guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(
        file.path(plot_dir, 'spatial_plots', sprintf('%s_pair.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()

    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        #   Plot all clusters spatially
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = 'cell_type',
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = cell_type_colors
            ) +
                guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(
        file.path(plot_dir, 'spatial_plots', sprintf('%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

session_info()
