#   Check whether manual annotation of habenula/thalamus from the Shiny app
#   looks reasonable. Then annotate the optimal Banksy resolution with cell
#   types using registration results against the Yalcinbas fine snRNA-seq data

library(here)
library(tidyverse)
library(spatialLIBD)
library(Polychrome)
library(sessioninfo)

cor_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cor_vs_snRNAseq_fine.rds'
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
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.1.3.4 = '#004F2D',
    LHb.4 = '#84DCC6',
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

#   Just rewriting HD results (stuff like 'LHb.7/LHb.2' => 'LHb.2.7') 
manual_anno = c('19' = 'LHb.2.7', '23' = 'LHb.1.3.4')

dir.create(plot_dir, showWarnings = FALSE)
for (subdir in c('region', 'anno_broad', 'anno_fine')) {
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
