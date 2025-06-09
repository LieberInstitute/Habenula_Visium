#   Explore the spatial distributions of cell sizes and number of extracellular
#   bins. This is to better understand if these can be used as metrics for QC in
#   the samples with poor-quality H&E images

library(here)
library(tidyverse)
library(data.table)
library(spatialLIBD)
library(HDF5Array)
library(cowplot)
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

plot_list = list()
for (sample_id in unique(spe$sample_id)) {
    plot_list[[sample_id]] = vis_clus(
        spe, sampleid = sample_id, clustervar = 'any_neighbors',
        is_stitched = TRUE, point_size = 10
    )
}

png(file.path(plot_dir, 'any_neighbors_spatial.png'), width = 4000, height = 800)
plot_grid(plotlist = plot_list, nrow = 1)
dev.off()
