library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

hd_cell_path = here(
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
plot_dir = here('plots', '14_supp_tables', 'plots')
cluster_levels = c(
    "MHb_A", "MHb_B", "MHb_C", "MHb_D", "Excit_LHb", "LHb_A", "LHb_B", "LHb_C",
    "GABA_LHb_C.1", "GABA_LHb_C.2", "Excit.Thal/GABA_LHb_C.2", "Excit.Thal",
    "Inhib.Thal", "Ependymal", "Subependymal", "Astrocyte", "Endo", "Microglia",
    "Oligo", "OPC", sprintf("Factor_%d", 0:16)
)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

hd_cell_df = read_csv(hd_cell_path, show_col_types = FALSE) |>
    filter(
        cell_type_group == 'fine',
        #   Not enough genes for MAGMA to give reliable results
        !(cell_type %in% c('Endo/microglia', 'LHb_C/GABA_LHb_C.2'))) |>
    left_join(
        read_csv(hd_map_path, show_col_types = FALSE),
        by = c('cell_type' = 'old_cell_type')
    ) |>
    select(new_cell_type, gwas_group, neg_log_p, p_label) |>
    dplyr::rename(cell_type = new_cell_type) |>
    mutate(dataset = 'HD_cell_types')

k8_df = read_csv(k8_path, show_col_types = FALSE) |>
    select(cell_type, gwas_group, neg_log_p, p_label) |>
    mutate(
        cell_type = str_replace(cell_type, '^X', 'Factor_'),
        dataset = 'HD_FICTURE_all_bin_k8'
    )

k17_df = read_csv(k17_path, show_col_types = FALSE) |>
    select(cell_type, gwas_group, neg_log_p, p_label) |>
    mutate(
        cell_type = str_replace(cell_type, '^X', 'Factor_'),
        dataset = 'HD_FICTURE_extracellular_k17'
    )

multiome_df = read_csv(multiome_path, show_col_types = FALSE) |>
    filter(cell_type_group == 'fine') |>
    left_join(
        read_csv(multiome_map_path, show_col_types = FALSE),
        by = c('cell_type' = 'old_cell_type')
    ) |>
    select(new_cell_type, gwas_group, neg_log_p, p_label) |>
    dplyr::rename(cell_type = new_cell_type) |>
    mutate(dataset = 'Multiome_cell_types')

marker_df = bind_rows(hd_cell_df, k8_df, k17_df, multiome_df) |>
    group_by(cell_type, gwas_group) |>
    mutate(
        cell_type = ifelse(
            !grepl('^Factor_', cell_type) & (n() > 1),
            sprintf(
                '%s (%s)', cell_type, str_extract(dataset, '^(HD|Multiome)')
            ),
            cell_type
        )
    )

cluster_levels = c(
    t(outer(cluster_levels, c(' (HD)', ' (Multiome)', ''), FUN = paste0))
)
cluster_levels = cluster_levels[cluster_levels %in% marker_df$cell_type]

p = marker_df |>
    mutate(
        cell_type = factor(cell_type, levels = cluster_levels),
        gwas_category = case_when(
            grepl('^[ACOS]UD', gwas_group) ~ 'Substance Use',
            gwas_group == 'p_factor_Grotzinger' ~ 'P-Factor',
            TRUE ~ 'Psychiatric'
        ),
        cell_type_category = case_when(
                grepl('^Factor_', cell_type) & (dataset == 'HD_FICTURE_all_bin_k8') ~ 'All-Data Factors',
                grepl('^Factor_', cell_type) & (dataset == 'HD_FICTURE_extracellular_k17') ~ 'Extracellular Factors',
                TRUE ~ 'Cell Types'
            ) |>
            factor(
                levels = c(
                    'Cell Types', 'All-Data Factors', 'Extracellular Factors'
                )
            )
    ) |>
    ggplot(
            aes(
                x = gwas_group, y = cell_type, fill = neg_log_p, label = p_label
            )
        ) +
        geom_tile() +
        geom_text(size = 6) +
        scale_fill_viridis_c() +
        facet_grid(
            cell_type_category ~ gwas_category,
            scales = "free", space = "free"
        ) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "GWAS Trait", y = "Cell Type", fill = "-log10(p)")
pdf(file.path(plot_dir, 'MAGMA_heatmap.pdf'), width = 10, height = 25)
print(p)
dev.off()

session_info()
