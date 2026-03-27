library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

sample_id = 'Br9090_1'
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
hb_thal_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'hb_thal_manual_anno.csv.gz'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'misc_paper_figs'
)
cell_type_colors = c(
    Excit.Thal = '#9e4ad1',
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    MHb.1 = '#FF00FF',
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.1.3.4 = '#004F2D',
    LHb.4 = '#84DCC6'
)
cell_type_colors_2 = c(
    MHb.1 = '#FF00FF',
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.1.3.4 = '#004F2D',
    LHb.4 = '#84DCC6',
    other_hb = '#635A69',
    non_hb = '#C4C4C4'
)
region_colors_1 = c(
    LHb = '#93151d',
    MHb = '#93151d',
    other = '#C4C4C4'
)
region_colors_2 = c(
    MHb = '#2546D8',
    LHb = '#93151d',
    other_hb = '#635A69',
    non_hb = '#C4C4C4'
)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vis_clus_hd = function(
        spe, clustervar, sample_id, plot_path, colors = region_colors_1
    ) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = clustervar,
                is_stitched = TRUE, point_size = 30, spatial = FALSE,
                colors = colors
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

#   Cell type- cluster map and manual annotation of brain region (hb + thal)
anno_df = read_csv(ct_anno_path, show_col_types = FALSE)
region_df = read_csv(hb_thal_anno_path, show_col_types = FALSE) |>
    dplyr::rename(key = spot_name)

cell_type_df = tibble(key = spe$key, tissue_id = spe$sample_id) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    left_join(region_df, by = 'key') |>
    mutate(
        cell_type = factor(
            anno_df$fine_cell_type[
                match(as.character(banksy), anno_df$cluster)
            ],
            levels = names(cell_type_colors)
        ),
        cell_type_lhb = ifelse(grepl('^LHb', cell_type), 'LHb', 'other'),
        cell_type_mhb = ifelse(grepl('^MHb', cell_type), 'MHb', 'other'),
        donor = paste0('Br', str_extract(key, '[0-9]{4}$')),
        region_anno = replace_na(ManualAnnotation, 'other'),
        cell_type_hb = factor(
            case_when(
                (region_anno == 'habenula') & (cell_type_lhb == 'LHb') ~ 'LHb',
                (region_anno == 'habenula') & (cell_type_mhb == 'MHb') ~ 'MHb',
                region_anno == 'habenula' ~ 'other_hb',
                TRUE ~ 'non_hb'
            ),
            levels = names(region_colors_2)
        ),
        cell_type_custom = factor(
            case_when(
                grepl('^[ML]Hb', cell_type) ~ cell_type,
                region_anno == 'habenula' ~ 'other_hb',
                TRUE ~ 'non_hb'
            ),
            levels = names(cell_type_colors_2)
        )
    )
stopifnot(!any(is.na(cell_type_df$cell_type)))

spe$cell_type_lhb = cell_type_df$cell_type_lhb
spe$cell_type_mhb = cell_type_df$cell_type_mhb
spe$cell_type_hb = cell_type_df$cell_type_hb
spe$cell_type_custom = cell_type_df$cell_type_custom

#   Order donors by proportion of habenula, but keep tissue sections ordered by 1, 2
donor_order = cell_type_df |>
    group_by(donor) |>
    summarize(prop_hb = mean(grepl('^[ML]Hb', cell_type))) |>
    arrange(prop_hb) |>
    pull(donor)
cell_type_df$tissue_id = factor(
    cell_type_df$tissue_id,
    levels = paste0(rep(donor_order, each = 2), c('_1', '_2'))
)

#   Composition barplot of cell types by tissue section
p = ggplot(cell_type_df, aes(x = tissue_id, fill = cell_type)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = cell_type_colors) +
    labs(x = "Tissue Section", y = "Proportion", fill = "Cell Type") +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(file.path(plot_dir, 'cell_type_tissue_section_barplot.pdf'))
print(p)
dev.off()

#   Simple plot to show how we subset to habenula
vis_clus_hd(
    spe = spe, clustervar = 'cell_type_hb', sample_id = sample_id,
    colors = region_colors_2, plot_path = file.path(plot_dir, 'Hb_anno.png')
)

#   Plot for our internal interest: how much of data-driven habenula clusters
#   lie outside manually annotated habenula region?
for (tissue_id in paste0(unique(cell_type_df$donor), '_1')) {
    vis_clus_hd(
        spe = spe, clustervar = 'cell_type_custom', sample_id = tissue_id,
        colors = cell_type_colors_2,
        plot_path = file.path(
            plot_dir, sprintf('data_driven_vs_manual_hb_%s.png', tissue_id)
        )
    )
}

#   Plot of where medial and lateral habenula are spatially in the sample
spe_hb = spe[, cell_type_df$region_anno == 'habenula']
vis_clus_hd(
    spe = spe_hb, clustervar = 'cell_type_mhb', sample_id = sample_id,
    plot_path = file.path(plot_dir, 'MHb_anno.png')
)
vis_clus_hd(
    spe = spe_hb, clustervar = 'cell_type_lhb', sample_id = sample_id,
    plot_path = file.path(plot_dir, 'LHb_anno.png')
)

session_info()
