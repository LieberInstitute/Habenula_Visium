#   Explore how clustering results line up before and after filtering bad
#   probes. One problematic cluster (composed of many bad probes) was very
#   spatially scattered-- is this true after filtering bad probes?

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)

spe_good_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
spe_bad_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
cluster_good_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_8',
    'leiden_res1.csv'
)
cluster_bad_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_8',
    'leiden_res1.csv'
)
svg_good_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'nnSVG_out',
    'merged_SVGs.txt'
)
svg_bad_path = here(
    'processed-data', '10_HD_bin_level', 'nnSVG_out', 'merged_SVGs.txt'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_8',
    'leiden_res1'
)

prep_clustering_results = function(spe_dir, cluster_path, cluster_colname) {
    spe = loadHDF5SummarizedExperiment(spe_dir)
    
    cluster_df = tibble(
            x = round(spatialCoords(spe)[, 'pxl_col_in_fullres'], 1),
            y = round(spatialCoords(spe)[, 'pxl_row_in_fullres'], 1),
            sample_id = spe$sample_id,
            key = spe$key,
            segmentation_type = spe$labels_joint_source
        ) |>
        filter(segmentation_type == 'primary') |>
        left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
        dplyr::rename(!!cluster_colname := banksy_lambda0_8) |>
        select(x, y, sample_id, {{ cluster_colname }})

    stopifnot(!any(is.na(cluster_df[[cluster_colname]])))
    
    return(cluster_df)
}

cluster_good_df = prep_clustering_results(
    spe_good_dir, cluster_good_path, 'cluster_good'
)
cluster_bad_df = prep_clustering_results(
    spe_bad_dir, cluster_bad_path, 'cluster_bad'
)
cluster_df = inner_join(
    cluster_good_df, cluster_bad_df, by = c('x', 'y', 'sample_id')
)

#   Compute the Jaccard index for each combination of "good" and "bad" clusters
jaccard_df_list = list()
for (good_val in unique(cluster_df$cluster_good)) {
    for (bad_val in unique(cluster_df$cluster_bad)) {
        intersect_size = cluster_df |>
            filter((cluster_good == good_val) & (cluster_bad == bad_val)) |>
            nrow()
        union_size = cluster_df |>
            filter((cluster_good == good_val) | (cluster_bad == bad_val)) |>
            nrow()
        jaccard_df_list[[length(jaccard_df_list) + 1]] = tibble(
            cluster_good = good_val,
            cluster_bad = bad_val,
            jaccard_index = intersect_size / union_size
        )
    }
}

#   Plot a heatmap of Jaccard indices
p = do.call(rbind, jaccard_df_list) |>
    ggplot(
        aes(
            x = cluster_good, y = cluster_bad,
            fill = jaccard_index
        )
    ) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_size = 25) +
    labs(
        x = 'Cluster: Filtered Probe Set', y = 'Cluster: Full Probe Set',
        fill = 'Jaccard\nIndex'
    )
pdf(file.path(plot_dir, 'jaccard_index.pdf'), width = 9)
print(p)
dev.off()

#   Banksy results are extremely similar. Is it because the input genes are
#   mostly the same?
svg_good = readLines(svg_good_path)
svg_bad = readLines(svg_bad_path)
stopifnot(length(svg_good) == length(svg_bad))
message(sprintf('SVG intersect proportion: %.2f', mean(svg_good %in% svg_bad)))

session_info()
