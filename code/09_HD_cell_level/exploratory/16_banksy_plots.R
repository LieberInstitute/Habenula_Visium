library(here)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(paletteer)
library(tidyverse)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_8',
    'leiden_res0_8.csv'
)
plot_path = here(
    'plots', '09_HD_cell_level', 'banksy', 'lambda0_8', 'leiden_res0_8',
    'manual_clusters_H1-MVPY9BW_A1_8433.png'
)
colors = c(
    paletteer_d("Polychrome::palette36", 14)[c(1, 3:10, 12:14)], '#2222AAFF', '#222200FF'
)
colors[8] = '#008020FF'

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(all(spe$key %in% cluster_df$key))
spe$banksy = cluster_df$banksy_lambda0_8[match(spe$key, cluster_df$key)]

#   Replot with 
p = vis_clus(
        spe, sampleid = 'H1-MVPY9BW_A1_8433', clustervar = 'banksy',
        is_stitched = TRUE, point_size = 20, spatial = FALSE,
        colors = colors
    ) +
    guides(fill = guide_legend(override.aes = list(size = 10)))
png(plot_path, width = 1500, height = 1500)
print(p)
dev.off()

session_info()
