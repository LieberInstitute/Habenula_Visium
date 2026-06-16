library(tidyverse)
library(here)
library(sessioninfo)

source('/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/05_03_annotation_adjustments/celltype_colors.R')

extra_model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spatial_registration', 'modeling_results.rds'
)
cell_model_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'modeling_results', '1_8_cell_types.rds'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    't_stat_comparison'
)
cell_type_colors = c(
    Excit.Thal = "#4d55b7",
    Excit.Thal.Inhib_LHb_4.2 = "#2d1e6a",
    LHb.4 = "#00607A",
    LHb.4.Inhib_LHb_4.2 = "#003d4e",
    Inhib.Thal = "#9a9fe7",
    Astrocyte = "#532222",
    OPC = "#829454",
    Oligo =  "#384a08",
    Microglia = "#141b02",
    LHb.2.7 = "#6C9FA9",
    Endo = "#d95f02",
    Endo.microglia = "#8d3e01",
    Excit_LHb = "#A8B8BC",
    MHb.1 = "#A86A9A",
    MHb.2 = "#BCA6B6",
    Ependymal = "#f5a105ff",
    Subependymal = "#976d1d"
) 

dir.create(plot_dir, showWarnings = FALSE)

prep_model_stats = function(model_path, col_suffix) {
    readRDS(model_path)$enrichment |>
        as_tibble() |>
        select(ensembl, gene, matches('^(fdr|t_stat)_')) |>
        pivot_longer(
            cols = matches('^(fdr|t_stat)_'),
            names_to = c(".value", "cell_type"),
            names_pattern = "^(fdr|t_stat)_(.*)$"
        ) |>
        dplyr::rename(
            gene_id = ensembl, gene_name = gene
        ) |>
        dplyr::rename(!!sprintf('t_stat_%s', col_suffix) := t_stat) |>
        select(-fdr)  
}

cell_df = prep_model_stats(cell_model_path, 'cell')
extra_df = prep_model_stats(extra_model_path, 'extra')

stat_df = inner_join(
    cell_df, extra_df, by = c('gene_id', 'gene_name', 'cell_type')
)

facet_labels = stat_df |>
    group_by(cell_type) |>
    summarise(
        slope = coef(lm(t_stat_extra ~ t_stat_cell))[['t_stat_cell']],
        intercept = coef(lm(t_stat_extra ~ t_stat_cell))[['(Intercept)']],
        cor = cor(t_stat_cell, t_stat_extra, use = 'complete.obs')
    ) |>
    mutate(
        label_eq = sprintf(
            'y = %.2fx %s %.2f', slope, ifelse(intercept >= 0, '+', '-'),
            abs(intercept)
        ),
        label_cor = sprintf('cor = %.2f', cor)
    )

p = ggplot(stat_df, aes(x = t_stat_cell, y = t_stat_extra, color = cell_type)) +
    geom_point(alpha = 0.2) +
    geom_smooth(
        method = 'lm', se = FALSE, linetype = 'dotted', color = 'black',
        linewidth = 0.7
    ) +
    geom_text(
        data = facet_labels, aes(label = label_eq),
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.8,
        color = 'black', size = 3, inherit.aes = FALSE
    ) +
    geom_text(
        data = facet_labels, aes(label = label_cor),
        x = -Inf, y = Inf, hjust = -0.1, vjust = 3.4,
        color = 'black', size = 3, inherit.aes = FALSE
    ) +
    facet_wrap(~cell_type, nrow = 4) +
    scale_color_manual(values = cell_type_colors) +
    coord_fixed() +
    theme_bw(base_size = 13) +
    guides(color = 'none') +
    labs(
        x = 'Cellular t-stat',
        y = 'Extracellular t-stat',
    )
png(file.path(plot_dir, 't_stat_comparison.png'), width = 7, height = 7, units = 'in', res = 300)
print(p)
dev.off()

session_info()
