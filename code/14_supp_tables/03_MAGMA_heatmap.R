library(tidyverse)
library(here)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(sessioninfo)

out_path = here('processed-data', '14_supp_tables', 'magma_sig_sets.csv')
hd_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'heatmap_data.csv'
)
hd_extra_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'heatmap_data.csv'
)
k8_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'MAGMA', 'heatmap_data.csv'
)
k17_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', 'heatmap_data_k17.csv'
)
multiome_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/heatmap_data.csv'
multiome_map_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/cell_type_map.csv'
hd_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
gwas_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/gwas_info.csv'
plot_dir = here('plots', '14_supp_tables')
cluster_levels = c(
    "Hb", "MHb", "LHb", "MHb_A", "MHb_B", "MHb_C", "MHb_D", "Excit_LHb",
    "LHb_A", "LHb_B", "LHb_C", "GABA_LHb_C.1", "GABA_LHb_C.2",
    "Excit.Thal/GABA_LHb_C.2", "Excit.Thal", "Inhib.Thal", "Ependymal",
    "Subependymal", "Astrocyte", "Endo", "Microglia", "Oligo", "OPC",
    sprintf("Factor_%d", 0:16)
)
gwas_category_levels = c('P', 'Psychiatric', 'Substance Use')
gwas_levels = read_csv(gwas_path, show_col_types = FALSE)$manuscript_name

################################################################################
#   Functions
################################################################################

fix_map_df = function(map_path) {
    map_df = rbind(
        read_csv(map_path, show_col_types = FALSE) |>
            select(old_cell_type, new_cell_type),
        tibble(
            old_cell_type = c('Hb', 'MHb', 'LHb'),
            new_cell_type = c('Hb', 'MHb', 'LHb')
        )
    )
    return(map_df)
}

read_cell_df = function(cell_path, map_df) {
    cell_df = read_csv(cell_path, show_col_types = FALSE) |>
        left_join(map_df, by = c('cell_type' = 'old_cell_type')) |>
        filter(
            (cell_type_group == 'fine') | grepl('^[ML]?Hb$', cell_type)
        ) |>
        select(new_cell_type, gwas_group, neg_log_p, p_label) |>
        dplyr::rename(cell_type = new_cell_type)
    return(cell_df)
}

read_ficture_df = function(ficture_path) {
    ficture_df = read_csv(ficture_path, show_col_types = FALSE) |>
        select(cell_type, gwas_group, neg_log_p, p_label) |>
        mutate(cell_type = str_replace(cell_type, '^X', 'Factor_'))
    return(ficture_df)
}

add_magma_labels = function(magma_df) {
    magma_df |>
        mutate(
            broad_cell_type = case_when(
                grepl('[ML]Hb', cell_type) ~ 'habenula',
                grepl('Thal', cell_type) ~ 'thalamus',
                grepl('^Factor_', cell_type) ~ 'unassigned',
                TRUE ~ 'glia'
            ),
            gwas_category = factor(
                case_when(
                    grepl('^[ACOS]UD', gwas_group) ~ 'Substance Use',
                    gwas_group == 'p_factor_Grotzinger' ~ 'P',
                    TRUE ~ 'Psychiatric'
                ),
                levels = gwas_category_levels
            ),
            gwas_group = factor(gwas_group, levels = gwas_levels)
        )
}

magma_heatmap = function(magma_df) {
    marker_df = magma_df |>
        group_by(cell_type, dataset) |>
        mutate(
            neg_log_fdr = -log10(p.adjust(10^(-1 * neg_log_p), method = 'fdr')),
            fdr_label = case_when(neg_log_fdr > -log10(0.05) ~ '*', TRUE ~ '')
        ) |>
        ungroup() |>
        add_magma_labels() |>
        mutate(cell_type = factor(cell_type, levels = cluster_levels))

    p = marker_df |>
        ggplot(
                aes(
                    x = gwas_group, y = cell_type, fill = neg_log_fdr,
                    label = fdr_label
                )
            ) +
            geom_tile() +
            geom_text(size = 6) +
            scale_fill_viridis_c() +
            facet_grid(dataset ~ gwas_category, scales = "free", space = "free") +
            theme_bw(base_size = 15) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
            labs(x = "GWAS Trait", y = "Cell Type", fill = "-log10(FDR)")
    return(p)
}

magma_heatmap_complex = function(magma_df) {
    marker_df = magma_df |>
        group_by(cell_type, dataset) |>
        mutate(
            neg_log_fdr = -log10(p.adjust(10^(-1 * neg_log_p), method = 'fdr')),
            fdr_label = case_when(neg_log_fdr > -log10(0.05) ~ '*', TRUE ~ '')
        ) |>
        ungroup() |>
        add_magma_labels() |>
        mutate(cell_type = factor(cell_type, levels = cluster_levels))

    row_df = marker_df |>
        distinct(dataset, cell_type, broad_cell_type) |>
        mutate(row_id = paste(dataset, cell_type, sep = '::')) |>
        arrange(dataset, cell_type) |>
        select(row_id, dataset, cell_type, broad_cell_type)

    col_df = marker_df |>
        distinct(gwas_group, gwas_category) |>
        arrange(gwas_category, gwas_group) |>
        filter(!is.na(gwas_group))

    heatmap_df = marker_df |>
        mutate(row_id = paste(dataset, cell_type, sep = '::')) |>
        select(row_id, gwas_group, neg_log_fdr) |>
        pivot_wider(names_from = gwas_group, values_from = neg_log_fdr)

    label_df = marker_df |>
        mutate(row_id = paste(dataset, cell_type, sep = '::')) |>
        select(row_id, gwas_group, fdr_label) |>
        pivot_wider(names_from = gwas_group, values_from = fdr_label)

    heatmap_mat = heatmap_df |>
        column_to_rownames('row_id') |>
        as.matrix()
    heatmap_mat = heatmap_mat[row_df$row_id, col_df$gwas_group, drop = FALSE]

    label_mat = label_df |>
        column_to_rownames('row_id') |>
        as.matrix()
    label_mat = label_mat[row_df$row_id, col_df$gwas_group, drop = FALSE]
    label_mat[is.na(label_mat)] = ''

    row_annotation_colors = c(
        habenula = '#0072B2',
        thalamus = '#CC79A7',
        glia = '#E69F00',
        unassigned = '#999999'
    )

    row_ha = rowAnnotation(
        `Broad cell type` = row_df$broad_cell_type,
        col = list(`Broad cell type` = row_annotation_colors),
        show_annotation_name = TRUE,
        annotation_name_gp = gpar(fontsize = 12)
    )

    col_fun = colorRamp2(
        seq(
            min(heatmap_mat, na.rm = TRUE),
            max(heatmap_mat, na.rm = TRUE),
            length.out = 256
        ),
        viridis::viridis(256)
    )

    Heatmap(
        heatmap_mat,
        name = '-log10(FDR)',
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = row_df$dataset,
        column_split = col_df$gwas_category,
        left_annotation = row_ha,
        row_labels = row_df$cell_type,
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        column_names_rot = 90,
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (label_mat[i, j] != '') {
                grid.text(label_mat[i, j], x, y, gp = gpar(fontsize = 12))
            }
        },
        heatmap_legend_param = list(title = '-log10(FDR)')
    )
}

################################################################################
#   Main
################################################################################

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

hd_map_df = fix_map_df(hd_map_path)
multiome_map_df = fix_map_df(multiome_map_path)

hd_cell_df = read_cell_df(hd_cell_path, hd_map_df) |>
    #   Not enough genes for MAGMA to give reliable results
    filter(!(cell_type %in% c('Endo/microglia', 'LHb_C/GABA_LHb_C.2'))) |>
    mutate(dataset = 'HD Cellular')

hd_extra_df = read_cell_df(hd_extra_path, hd_map_df) |>
    #   Not enough genes for MAGMA to give reliable results
    filter(!(cell_type %in% c('Endo/microglia', 'LHb_C/GABA_LHb_C.2'))) |>
    mutate(dataset = 'HD Extracellular')

multiome_df = read_cell_df(multiome_path, multiome_map_df) |>
    mutate(dataset = 'Multiome Cellular')

k17_df = read_ficture_df(k17_path) |>
    mutate(dataset = 'Extracellular k = 17')

p = bind_rows(hd_cell_df, multiome_df) |>
    filter(!grepl('^[ML]?Hb$', cell_type)) |>
    magma_heatmap_complex()
pdf(file.path(plot_dir, 'MAGMA_heatmap_main_fine_cellular.pdf'), width = 8, height = 10)
draw(p)
dev.off()

p = bind_rows(hd_extra_df, k17_df) |>
    filter(!grepl('^[ML]?Hb$', cell_type)) |>
    magma_heatmap_complex()
pdf(
    file.path(plot_dir, 'MAGMA_heatmap_main_fine_extracellular.pdf'),
    width = 8, height = 10
)
draw(p)
dev.off()

p = bind_rows(hd_cell_df, hd_extra_df, multiome_df) |>
    filter(grepl('^[ML]?Hb$', cell_type)) |>
    mutate(
        dataset = case_when(
            dataset == 'HD Cellular' ~ 'HD Cell.',
            dataset == 'HD Extracellular' ~ 'HD Extra.',
            dataset == 'Multiome Cellular' ~ 'Multiome',
            TRUE ~ dataset
        )
    ) |>
    magma_heatmap()
pdf(file.path(plot_dir, 'MAGMA_heatmap_main_broad.pdf'), width = 8, height = 6)
print(p)
dev.off()

k8_df = read_ficture_df(k8_path) |>
    mutate(dataset = 'All-Bin k = 8')

p = bind_rows(k8_df, k17_df) |>
    magma_heatmap()
pdf(file.path(plot_dir, 'MAGMA_heatmap_supp.pdf'), width = 8, height = 10)
print(p)
dev.off()

#   Export significant sets, which makes filtering of gene-level stats easier
#   in a different script
bind_rows(hd_cell_df, hd_extra_df, multiome_df, k8_df, k17_df) |>
    group_by(cell_type, dataset) |>
    filter(p.adjust(10^(-1 * neg_log_p), method = 'fdr') < 0.05) |>
    ungroup() |>
    select(cell_type, gwas_group, dataset) |>
    mutate(
        dataset = case_when(
            dataset == 'HD Cellular' ~ 'HD_cell_types',
            dataset == 'HD Extracellular' ~ 'HD_extracellular',
            dataset == 'Multiome Cellular' ~ 'Multiome_cell_types',
            dataset == 'All-Bin k = 8' ~ 'HD_FICTURE_all_bin_k8',
            dataset == 'Extracellular k = 17' ~ 'HD_FICTURE_extracellular_k17',
            TRUE ~ dataset
        )
    ) |>
    write_csv(out_path)

session_info()
