#   Check whether manual annotation of habenula/thalamus from the Shiny app
#   looks reasonable. Then annotate the optimal Banksy resolution with cell
#   types using registration results against the Yalcinbas fine snRNA-seq data

library(here)
library(tidyverse)
library(spatialLIBD)
library(Polychrome)
library(sessioninfo)

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
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'cluster_annotation.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'annotation'
)
manual_anno = c(
    '18' = 'Astrocyte', # From mapping to multiome 'C.20.Astrocyte'
    '24' = 'Oligo',     # From mapping to multiome 'C.02.Oligo'
    # From manual hb/thal anno in HD (referencing multiome would be circular!)
    '12' = 'Excit.Thal',
    '16' = 'LHb.2.7',   # Just rewriting HD results ('LHb.7/LHb.2')
    '22' = 'LHb.1.3.4'  # Same here
)
region_colors = c(
    habenula = "#B1092D", thalamus = "#0B52C4", other = "#DFE1DD"
)

dir.create(plot_dir, showWarnings = FALSE)
for (subdir in c('region', 'banksy_as_is', 'banksy_C4_C27', 'anno_broad', 'anno_fine')) {
    dir.create(file.path(plot_dir, subdir), showWarnings = FALSE)
}


################################################################################
#   Functions
################################################################################

#   vis_clus with HD settings and saving to file
vis_clus_HD = function(spe, sampleid, clustervar, plot_dir, ...) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sampleid, clustervar = clustervar,
                is_stitched = TRUE, point_size = 20, spatial = FALSE, ...
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    
    png(
        file.path(plot_dir, sprintf('%s.png', sampleid)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

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
    vis_clus_HD(
        spe, sample_id, 'region_anno', file.path(plot_dir, 'region'),
        colors = region_colors
    )

    #   The reason for plotting Banksy results again is to have consistent
    #   colors across samples; some samples lack certain clusters and the
    #   default original plots improperly assigned colors to clusters
    vis_clus_HD(
        spe, sample_id, 'banksy', file.path(plot_dir, 'banksy_as_is'),
        colors = cluster_colors
    )
}

################################################################################
#   Annotate clusters with cell types
################################################################################

message('Cluster 12 locates in thalamus based on manual annotation:')
table(spe$region_anno[spe$banksy == '12'])

#   These clusters registered with poor layer confidence and look ambiguous.
#   Where are they located?
cluster_colors = region_colors
names(cluster_colors) = c('4', '27', 'other')
spe$temp = case_when(
    spe$banksy %in% c('4', '27') ~ spe$banksy,
    TRUE ~ 'other'
)

#   Check locations of clusters 4 and 27
for (sample_id in unique(spe$sample_id)) {
    vis_clus_HD(
        spe, sample_id, 'temp', file.path(plot_dir, 'banksy_C4_C27'),
        colors = cluster_colors
    )
}

#   Create mapping of cluster to cell type
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        fine_cell_type = case_when(
            cluster %in% names(manual_anno) ~ manual_anno[cluster],
            layer_confidence == 'poor' ~ 'Ambig',
            cluster == '12' ~ 'Excit.Thal',
            TRUE ~ layer_label
        ),
        broad_cell_type = str_replace(
            fine_cell_type, '(\\..*$|Excit\\.|Inhib\\.)', ''
        )
    ) |>
    select(cluster, broad_cell_type, fine_cell_type) |>
    arrange(as.integer(cluster))

write_csv(anno_df, out_path)

#   Now plot the broad and fine annotations on each sample
spe$anno_broad = factor(
    anno_df$broad_cell_type[match(spe$banksy, anno_df$cluster)],
    levels = unique(anno_df$broad_cell_type)
)
spe$anno_fine = factor(
    anno_df$fine_cell_type[match(spe$banksy, anno_df$cluster)],
    levels = unique(anno_df$fine_cell_type)
)

for (sample_id in unique(spe$sample_id)) {
    for (anno_level in c('anno_broad', 'anno_fine')) {
        vis_clus_HD(spe, sample_id, anno_level, file.path(plot_dir, anno_level))
    }
}

session_info()
