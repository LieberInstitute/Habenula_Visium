#   Check whether manual annotation of habenula/thalamus from the Shiny app
#   looks reasonable. Then annotate the optimal Banksy resolution with cell
#   types using registration results against the Yalcinbas fine snRNA-seq data

library(here)
library(tidyverse)
library(spatialLIBD)
library(Polychrome)

cor_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'lambda0_2', 'cor_vs_snRNAseq_fine.rds'
)
cor_index = 17
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
hb_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'hb_thal_manual_anno.csv.gz'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'annotation'
)
manual_anno = c(
    '18' = 'Astrocyte', # From mapping to multiome 'C.20.Astrocyte'
    '24' = 'Oligo',     # From mapping to multiome 'C.02.Oligo'
    '12' = 'TODO'       # Need to check manual hb/thal anno in HD (referencing multiome would be circular!)
)
region_colors = c(
    habenula = "#B1092D", thalamus = "#0B52C4", other = "#DFE1DD"
)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(file.path(plot_dir, 'region'), showWarnings = FALSE)
dir.create(file.path(plot_dir, 'banksy_as_is'), showWarnings = FALSE)

################################################################################
#   Add Banksy clusters and region annotation to the SPE. Plot both
################################################################################

spe = readRDS(spe_path)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(setequal(spe$key, cluster_df$key))
spe$banksy = factor(
    as.character(cluster_df$banksy_lambda0_2[match(spe$key, cluster_df$key)]),
    levels = as.character(sort(unique(cluster_df$banksy_lambda0_2)))
)
cluster_colors = palette36.colors(length(levels(spe$banksy)))
names(cluster_colors) = levels(spe$banksy)

#   Add manual annotation of habenula/thalamus/other
hb_anno_df = read_csv(hb_anno_path, show_col_types = FALSE) |>
    dplyr::rename(key = spot_name, region_anno = ManualAnnotation) |>
    select(key, region_anno)
spe$region_anno = tibble(key = spe$key) |>
    left_join(hb_anno_df, by = "key") |>
    pull(region_anno) |>
    replace_na('other')

#   Plot the region annotation on each sample to make sure it worked
for (sample_id in unique(spe$sample_id)) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p_region = vis_clus(
                spe, sampleid = sample_id, clustervar = 'region_anno',
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = region_colors
            ) +
                guides(fill = guide_legend(override.aes = list(size = 8)))
        p_clus = vis_clus(
                spe, sampleid = sample_id, clustervar = 'banksy',
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = cluster_colors
            ) +
                guides(fill = guide_legend(override.aes = list(size = 8)))
        
    }
    png(
        file.path(plot_dir, 'region', sprintf('%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p_region)
    dev.off()
    png(
        file.path(plot_dir, 'banksy_as_is', sprintf('%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p_clus)
    dev.off()
}

################################################################################
#   Annotate clusters with cell types
################################################################################

anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble()
