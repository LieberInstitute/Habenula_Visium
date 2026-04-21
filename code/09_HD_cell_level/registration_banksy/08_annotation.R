#   Check whether manual annotation of habenula/thalamus from the Shiny app
#   looks reasonable. Then annotate the optimal Banksy resolution with cell
#   types using (mostly) registration results against the mid-res multiome data

library(here)
library(tidyverse)
library(spatialLIBD)
library(Polychrome)
library(sessioninfo)
library(cowplot)

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
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'annotation'
)
demo_samples = c("Br9090_1", "Br8433_1", "Br9902_1")
fine_colors = c(
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    'Endo/microglia' = '#a4511b',
    MHb.1 = '#FF00FF',
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.4 = '#84DCC6',
    "LHb.4/Inhib_LHb_4.2" = '#004F2D',
    "Excit.Thal/Inhib_LHb_4.2" = '#9faefb',
    Excit_LHb = '#6aff00',
    Excit.Thal = '#9e4ad1',
    Subependymal = '#4c00ff',
    Ependymal = '#0e005c'
)
broad_colors = c(
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    'Endo/microglia' = '#a4511b',
    MHb = '#FF00FF',
    LHb = '#004F2D',
    Excit.Thal = '#9e4ad1',
    Subependymal = '#4c00ff',
    Ependymal = '#0e005c'
)
ambig_colors = c(
    'LHb' = '#758E4F',
    'Excit.Thal' = '#AF1D1D',
    '9' = '#1E19AD',
    'Other' = '#aeaeae'
)
ambig_2_colors = c(
    'LHb' = '#758E4F',
    'Excit.Thal' = '#AF1D1D',
    '25' = '#1E19AD',
    'Other' = '#aeaeae'
)

#   We looked at spatial plots, registration against the fine multiome data, and
#   top markers to manually resolve some ambiguous or hard-to-label clusters
manual_anno = c(
    '15' = 'MHb.2',
    '8' = 'MHb.1',
    '23' = 'Excit_LHb',
    '5' = 'LHb.4/Inhib_LHb_4.2',
    '13' = 'Subependymal',
    '14' = 'Subependymal',
    '2' = 'Endo/microglia',
    '9' = 'Excit.Thal/Inhib_LHb_4.2',
    '16' = 'Drop'
)

dir.create(plot_dir, showWarnings = FALSE)
for (subdir in c('anno_broad', 'anno_fine', 'ambig', 'ambig_2', 'faceted')) {
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
#   Add Banksy clusters to the SPE and plot
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

#   Plot each cluster individually in all samples
all_samples = unique(spe$sample_id)[grepl('_1$', unique(spe$sample_id))]
for (this_cluster in levels(spe$banksy)) {
    p_list = list()
    for (this_sample_id in all_samples) {
        spe$temp = ifelse(spe$banksy == this_cluster, this_cluster, 'Other') |>
            factor(levels = c(this_cluster, 'Other'))

        #   Run twice to overcome a bug with different behavior on the first
        #   plot
        for (i in seq_len(2)) {
            p_list[[this_sample_id]] = vis_clus(
                    spe, sampleid = this_sample_id, clustervar = 'temp',
                    is_stitched = TRUE, point_size = 10, spatial = FALSE,
                    colors = c('#AF1D1D', '#aeaeae')
                ) +
                guides(fill = guide_legend(override.aes = list(size = 5)))
        }
    }
    p = plot_grid(plotlist = p_list, nrow = 1)
    png(
        file.path(plot_dir, 'faceted', sprintf('%s.png', this_cluster)),
        width = 2000, height = 400
    )
    print(p)
    dev.off()
}

################################################################################
#   Investigate ambiguous clusters
################################################################################

spe$ambig = case_when(
        spe$banksy %in% c(19, 23, 7, 25, 5) ~ 'LHb',
        spe$banksy == 20 ~ 'Excit.Thal',
        spe$banksy == 9 ~ '9',
        TRUE ~ 'Other'
    ) |>
    factor(levels = names(ambig_colors))
spe$ambig_2 = case_when(
        spe$banksy %in% c(19, 23, 7, 5) ~ 'LHb',
        spe$banksy == 20 ~ 'Excit.Thal',
        spe$banksy == 25 ~ '25',
        TRUE ~ 'Other'
    ) |>
    factor(levels = names(ambig_2_colors))

for (this_sample_id in demo_samples) {
    vis_clus_HD(
        spe, this_sample_id, 'ambig', file.path(plot_dir, 'ambig'),
        colors = ambig_colors
    )
    vis_clus_HD(
        spe, this_sample_id, 'ambig_2', file.path(plot_dir, 'ambig_2'),
        colors = ambig_2_colors
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
        broad_cell_type = case_when(
            grepl('^MHb', fine_cell_type) ~ 'MHb',
            grepl('LHb', fine_cell_type) & !grepl('^Excit\\.Thal', fine_cell_type) ~ 'LHb',
            TRUE ~ fine_cell_type
        )
    ) |>
    select(cluster, broad_cell_type, fine_cell_type) |>
    arrange(as.integer(cluster))

write_csv(anno_df, out_path)

#   Highly ambiguous, sample-specific, and small cluster
spe = spe[, spe$banksy != 16]

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
