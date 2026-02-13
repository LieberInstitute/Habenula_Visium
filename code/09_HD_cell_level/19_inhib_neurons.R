#   At the Banksy clustering resolution = 4, we get a decent Inhib.Thal
#   cluster that doesn't exist at the optimal resolution. Quickly investigate
#   where this cluster is located. Note resolution 8 has an even more conclusive
#   Inhib.Thal cluster, but it's about ~400 cells and not even visible in plots

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

cor_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cor_vs_snRNAseq_fine.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res4.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here('plots', '09_HD_cell_level', 'no_secondary', 'inhib_neurons')
cor_index = 21
cluster_colors = c('Inhib.Thal' = '#00FF51', 'Other' = '#A7ABA8')

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vis_banksy = function(spe, clustervar, sample_id, plot_path) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = clustervar,
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = cluster_colors
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(plot_path, width = 1500, height = 1500)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)

inhib_cluster = annotate_registered_clusters(
        readRDS(cor_path)[[cor_index]], cutoff_merge_ratio = 0.1
    ) |>
    filter(layer_label == 'Inhib.Thal') |>
    pull(cluster) |>
    as.numeric()

spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    mutate(cell_type = ifelse(banksy == inhib_cluster, 'Inhib.Thal', 'Other')) |>
    pull(cell_type)
stopifnot(!any(is.na(spe$cell_type)))

for (sample_id in unique(spe$sample_id)) {
    vis_banksy(
        spe, clustervar = 'cell_type', sample_id = sample_id,
        plot_path = file.path(plot_dir, paste0(sample_id, '.png'))
    )
}

session_info()
