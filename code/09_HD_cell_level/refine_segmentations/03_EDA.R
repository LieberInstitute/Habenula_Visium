#   Explore the spatial distributions of cell sizes and number of extracellular
#   bins. This is to better understand if these can be used as metrics for QC in
#   the samples with poor-quality H&E images

library(here)
library(tidyverse)
library(data.table)
library(spatialLIBD)
library(HDF5Array)
library(cowplot)
library(paletteer)
library(sessioninfo)

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
extra_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'extracellular_bins.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_4.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'refine_segmentations', 'EDA'
)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Load and prep data
################################################################################

#   We'll only deal with primary cells, since the set of secondary cells changes
#   between the full-data SPE and the extracellular analysis. Also, only the
#   primary cells are affected by bad H&E images
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$labels_joint_source == 'primary']

#   For each cell, compute number of extracellular bins
extra_df = fread(extra_path) |>
    as_tibble() |>
    mutate(key = paste(cell_id, sample_id, sep = '_')) |>
    filter(key %in% spe$key) |>
    select(key) |>
    group_by(key) |>
    summarize(num_neighbors = n())

#   Join this metric with the SPE
spe$num_neighbors = colData(spe) |>
    as_tibble() |>
    left_join(extra_df, by = 'key') |>
    pull(num_neighbors)
spe$num_neighbors[is.na(spe$num_neighbors)] = 0

################################################################################
#   Check spatial distribution of number of cells and number of extracellular
#   bins
################################################################################

#   Check distribution of cells' number of constituent bins by sample
p = colData(spe)[, c('sample_id', 'bin_count')] |>
    as_tibble() |>
    filter(bin_count < 80) |>
    ggplot(aes(x = bin_count, fill = sample_id)) +
        geom_density(alpha = 0.25) +
        theme_bw(base_size = 20) +
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        labs(x = 'Cell Size (Bins)', y = 'Density', fill = 'Sample ID')
pdf(file.path(plot_dir, 'cell_size_density.pdf'), width = 10, height = 6)
print(p)
dev.off()

spe$any_neighbors = spe$num_neighbors > 0

plot_list_neighbor = list()
plot_list_size = list()
for (sample_id in unique(spe$sample_id)) {
    plot_list_neighbor[[sample_id]] = vis_clus(
        spe, sampleid = sample_id, clustervar = 'any_neighbors',
        is_stitched = TRUE, point_size = 10
    )
    plot_list_size[[sample_id]] = vis_gene(
        spe, sampleid = sample_id, geneid = 'bin_count', cap_percentile = 0.99,
        is_stitched = TRUE, point_size = 10
    )
}

png(file.path(plot_dir, 'any_neighbors_spatial.png'), width = 4000, height = 800)
plot_grid(plotlist = plot_list_neighbor, nrow = 1)
dev.off()

png(file.path(plot_dir, 'cell_size_spatial.png'), width = 4000, height = 800)
plot_grid(plotlist = plot_list_size, nrow = 1)
dev.off()

################################################################################
#   Are certain clusters specific to samples with bad H&E images?
################################################################################

banksy_df = read_csv(banksy_path, show_col_types = FALSE)
cluster_colors = paletteer_d(
    "Polychrome::palette36", length(unique(banksy_df$banksy_lambda0_2))
)

p = colData(spe)[, c('sample_id', 'key')] |>
    as_tibble() |>
    left_join(banksy_df, by = 'key') |>
    ggplot(aes(x = sample_id, fill = factor(banksy_lambda0_2))) +
        geom_bar(position = 'fill') +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = 'Sample ID', y = 'Proportion of Cells', fill = 'Banksy Cluster') +
        scale_fill_manual(values = cluster_colors)
pdf(file.path(plot_dir, 'banksy_by_sample.pdf'), height = 10)
print(p)
dev.off()

################################################################################
#   Are cells with no extracellular bins specific to certain Banksy clusters?
################################################################################

#   From the spatial plots, we know that cells with no extracellular bins have
#   a spatial pattern, including being slightly less likely to occur in habenula.
#   Let's confirm this statistically with a chi-squared test.

cluster_df = colData(spe)[, 'key', drop = FALSE] |>
    as_tibble() |>
    left_join(extra_df, by = 'key') |>
    left_join(banksy_df, by = 'key') |>
    mutate(
        banksy_cluster = factor(banksy_lambda0_2),
        num_neighbors = ifelse(is.na(num_neighbors), 0, num_neighbors)
    ) |>
    select(-banksy_lambda0_2)

no_df = cluster_df |>
    filter(num_neighbors == 0) |>
    group_by(banksy_cluster) |>
    summarize(num_cells_no = n())

cluster_df = cluster_df |>
    group_by(banksy_cluster) |>
    summarize(num_cells_total = n()) |>
    left_join(no_df, by = 'banksy_cluster')
    
chi_result = chisq.test(
    cluster_df$num_cells_no,
    p = cluster_df$num_cells_total / sum(cluster_df$num_cells_total)
)

if (chi_result$p.value < 0.05) {
    message(
        'Some Banksy clusters are enriched in cells with no extracellular bins.'
    )
}

session_info()
