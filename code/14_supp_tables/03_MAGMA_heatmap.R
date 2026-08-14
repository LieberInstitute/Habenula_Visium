library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

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

magma_heatmap = function(magma_df) {
    marker_df = magma_df |>
        group_by(cell_type, dataset) |>
        mutate(
            neg_log_fdr = -log10(p.adjust(10^(-1 * neg_log_p), method = 'fdr'))
        ) |>
        ungroup() |>
        mutate(
            cell_type = factor(cell_type, levels = cluster_levels),
            gwas_category = case_when(
                grepl('^[ACOS]UD', gwas_group) ~ 'Substance Use',
                gwas_group == 'p_factor_Grotzinger' ~ 'P-Factor',
                TRUE ~ 'Psychiatric'
            )
        )

    p = marker_df |>
        ggplot(
                aes(
                    x = gwas_group, y = cell_type, fill = neg_log_p, label = p_label
                )
            ) +
            geom_tile() +
            geom_text(size = 6) +
            scale_fill_viridis_c() +
            facet_grid(dataset ~ gwas_category, scales = "free", space = "free") +
            theme_bw(base_size = 15) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
            labs(x = "GWAS Trait", y = "Cell Type", fill = "-log10(p)")
    return(p)
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

p = bind_rows(hd_cell_df, hd_extra_df, multiome_df) |>
    magma_heatmap()
pdf(file.path(plot_dir, 'MAGMA_heatmap_main.pdf'), width = 8, height = 15)
print(p)
dev.off()

p = bind_rows(hd_cell_df, hd_extra_df, multiome_df) |>
    mutate(
        broad_cell_type = case_when(
            cell_type == 'Excit.Thal/GABA_LHb_C.2' ~ 'Other',
            grepl('[ML]?Hb', cell_type) ~ 'Habenula',
            grepl('Thal', cell_type) ~ 'Thalamus',
            TRUE ~ 'Glia'
        )
    ) |>
    filter(broad_cell_type != 'Other') |>
    group_by(broad_cell_type, gwas_group, dataset) |>
    summarise(prop_sig = mean(neg_log_p > -log10(0.05), na.rm = TRUE)) |>
    ungroup() |>
    mutate(
        gwas_category = case_when(
            grepl('^[ACOS]UD', gwas_group) ~ 'Substance Use',
            gwas_group == 'p_factor_Grotzinger' ~ 'P-Factor',
            TRUE ~ 'Psychiatric'
        )
    ) |>
    ggplot(aes(x = gwas_group, y = prop_sig, fill = broad_cell_type)) +
        geom_col(position = "fill") +
        facet_grid(dataset ~ gwas_category, scales = "free", space = "free") +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "GWAS Trait", y = "Prop. Significant", fill = "Broad Cell Type")
pdf(file.path(plot_dir, 'MAGMA_stacked_barplot.pdf'), height = 8)
print(p)
dev.off()

k8_df = read_ficture_df(k8_path) |>
    mutate(dataset = 'All-Bin k = 8')

k17_df = read_ficture_df(k17_path) |>
    mutate(dataset = 'Extracellular k = 17')

p = bind_rows(k8_df, k17_df) |>
    magma_heatmap()
pdf(file.path(plot_dir, 'MAGMA_heatmap_supp.pdf'), width = 8, height = 10)
print(p)
dev.off()

session_info()
