library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_norm_filtered.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'banksy', 'lambda0_2',
    'leiden_res1_4.csv'
)
plot_path = here(
    'plots', '09_HD_cell_level', 'new_samples', 'registration_banksy',
    'habenula_spot_plot.png'
)
sample_id = 'H1-W369TJK_D1_9090'
habenula_clusters = c(2, 11)
cluster_colors = c('2' = '#0150B8', '11' = '#C23853', 'Other' = '#ACB3B6')

#   Load and subset to a good-looking sample
spe = readRDS(spe_path)
spe = spe[, spe$sample_id == sample_id]

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$banksy = cluster_df$banksy_lambda0_2[match(spe$key, cluster_df$key)]

spe$banksy = factor(
    ifelse(
        spe$banksy %in% habenula_clusters, as.character(spe$banksy), 'Other'
    ),
    levels = c(as.character(habenula_clusters), 'Other')
)

p = vis_clus(
        spe, sampleid = sample_id, clustervar = 'banksy',
        is_stitched = TRUE, point_size = 20, spatial = FALSE,
        colors = cluster_colors
    ) +
    guides(fill = guide_legend(override.aes = list(size = 15)))
png(plot_path, width = 1500, height = 1500)
print(p)
dev.off()

session_info()
