library(here)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_8',
    'leiden_res1.csv'
)
plot_path = here(
    'plots', '09_HD_cell_level', 'registration_banksy', 'habenula_spot_plot.png'
)
sample_id = 'H1-MVPY9BW_A1_8433'
habenula_clusters = c(4, 6, 11)

#   Load ans subset to a good-looking sample
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id == sample_id]

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$banksy = cluster_df$banksy_lambda0_8[match(spe$key, cluster_df$key)]

spe$banksy = ifelse(spe$banksy %in% habenula_clusters, spe$banksy, 'Other')

p = vis_clus(
        spe, sampleid = sample_id, clustervar = 'banksy',
        is_stitched = TRUE, point_size = 20, spatial = FALSE,
        colors = colors
    ) +
    guides(fill = guide_legend(override.aes = list(size = 10)))
png(plot_path, width = 1500, height = 1500)
print(p)
dev.off()

session_info()