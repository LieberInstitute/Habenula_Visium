#   Check whether manual annotation of habenula/thalamus from the Shiny app
#   looks reasonable. Then annotate the optimal Banksy resolution with cell
#   types using (mostly) registration results against the mid-res multiome data

library(here)
library(tidyverse)
library(spatialLIBD)
library(Polychrome)
library(sessioninfo)

cor_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'all_genes', 'cor_vs_multiome_mid.rds'
)
cor_index = 18
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
hb_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'hb_thal_manual_anno.csv.gz'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'annotation'
)
demo_samples = c("Br9090_1", "Br8433_1")
region_colors = c(
    habenula = "#B1092D", thalamus = "#0B52C4", other = "#DFE1DD"
)
fine_colors = c(
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    MHb.1 = '#FF00FF',
    MHb.1.2 = '#f183e6',
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.1.3.4 = '#004F2D',
    LHb.4 = '#84DCC6',
    Inhib_LHb_4.2 = '#1E19AD',
    Excit.Thal = '#9e4ad1'
)
broad_colors = c(
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    MHb = '#FF00FF',
    LHb = '#004F2D',
    Excit.Thal = '#9e4ad1'
)
ambig_1_colors = c(
    '8' = '#1E19AD',
    '11' = '#758E4F',
    '15' = '#F6AE2D',
    'Astrocyte' = '#AF1D1D',
    'Other' = '#939393'
)
ambig_2_colors = c(
    '2' = '#758E4F',
    '13' = '#F6AE2D',
    'Astrocyte' = '#AF1D1D',
    'Endo' = '#1E19AD',
    'Other' = '#939393'
)
ambig_3_colors = c(
    '5' = '#F6AE2D',
    'LHb.2.7' = '#758E4F',
    'LHb.1.3.4' = '#490437',
    'LHb.4' = '#AF1D1D',
    'Inhib_LHb_4.2' = '#1E19AD',
    'Other' = '#939393'
)

#   Just rewriting HD results (stuff like 'LHb.7/LHb.2' => 'LHb.2.7') 
manual_anno = c('19' = 'LHb.2.7', '23' = 'LHb.1.3.4')

dir.create(plot_dir, showWarnings = FALSE)
for (subdir in c('region', 'anno_broad', 'anno_fine', 'ambig_1', 'ambig_2', 'ambig_3')) {
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
                is_stitched = TRUE, point_size = 25, spatial = FALSE, ...
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
    as.character(cluster_df$banksy[match(spe$key, cluster_df$key)]),
    levels = as.character(sort(unique(cluster_df$banksy)))
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
}

################################################################################
#   Investigate ambiguous clusters
################################################################################

spe$ambig_1 = case_when(
        spe$banksy %in% c(8, 11, 15) ~ as.character(spe$banksy),
        spe$banksy %in% c(6, 16, 12, 24, 14, 18) ~ 'Astrocyte',
        TRUE ~ 'Other'
    ) |>
    factor(levels = names(ambig_1_colors))
spe$ambig_2 = case_when(
        spe$banksy %in% c(2, 13) ~ as.character(spe$banksy),
        spe$banksy %in% c(6, 16, 12, 24, 14, 18) ~ 'Astrocyte',
        spe$banksy %in% c(22, 26, 21, 17) ~ 'Endo',
        TRUE ~ 'Other'
    ) |>
    factor(levels = names(ambig_2_colors))
spe$ambig_3 = case_when(
        spe$banksy == 19 ~ 'LHb.2.7',
        spe$banksy == 23 ~ 'LHb.1.3.4',
        spe$banksy %in% c(7, 25) ~ 'LHb.4',
        spe$banksy == 9 ~ 'Inhib_LHb_4.2',
        spe$banksy == 5 ~ '5',
        TRUE ~ 'Other'
    ) |>
    factor(levels = names(ambig_3_colors))

for (this_sample_id in demo_samples) {
    vis_clus_HD(
        spe, this_sample_id, 'ambig_1', file.path(plot_dir, 'ambig_1'),
        colors = ambig_1_colors
    )
    vis_clus_HD(
        spe, this_sample_id, 'ambig_2', file.path(plot_dir, 'ambig_2'),
        colors = ambig_2_colors
    )
    vis_clus_HD(
        spe, this_sample_id, 'ambig_3', file.path(plot_dir, 'ambig_3'),
        colors = ambig_3_colors
    )
}

################################################################################
#   Annotate clusters with cell types
################################################################################

#   Create mapping of cluster to cell type
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        fine_cell_type = case_when(
            cluster %in% names(manual_anno) ~ manual_anno[cluster],
            TRUE ~ layer_label
        ),
        broad_cell_type = str_replace(fine_cell_type, '\\.[0-9]+.*$', '')
    ) |>
    select(cluster, broad_cell_type, fine_cell_type) |>
    arrange(as.integer(cluster))

write_csv(anno_df, out_path)

#   Now plot the broad and fine annotations on each sample
spe$anno_broad = factor(
    anno_df$broad_cell_type[match(spe$banksy, anno_df$cluster)],
    levels = names(broad_colors)
)
spe$anno_fine = factor(
    anno_df$fine_cell_type[match(spe$banksy, anno_df$cluster)],
    levels = names(fine_colors)
)

for (sample_id in unique(spe$sample_id)) {
    vis_clus_HD(
        spe, sample_id, 'anno_broad', file.path(plot_dir, 'anno_broad'),
        colors = broad_colors
    )
    vis_clus_HD(
        spe, sample_id, 'anno_fine', file.path(plot_dir, 'anno_fine'),
        colors = fine_colors
    )
}

session_info()
