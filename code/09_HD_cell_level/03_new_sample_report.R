#   Explore how the new 2 samples compare to others in terms of notable metrics

library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)

spe_cell_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_norm_filtered.rds'
)
spe_bin_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'spe_norm_filtered.rds'
)
sample_colors = c(
    'H1-W369TJK_D1_9090' = '#FF4A1C',
    'H1-MVPY9BW_A1_8433' = '#03045E',
    'H1-MVPY9BW_D1_8667' = '#0077B6',
    'H1-6FX4YN3_A1_3942' = '#98DFAF',
    'H1-6FX4YN3_D1_9902' = '#5FB49C'
)
metric_names = c('sum_umi', 'sum_gene', 'expr_chrM_ratio')
plot_dir = here('plots', '09_HD_cell_level', 'new_samples', 'QC', 'after')

spe_cell = readRDS(spe_cell_path)
spe_bin = readRDS(spe_bin_path)

stopifnot(setequal(spe_cell$sample_id, names(sample_colors)))
names(sample_colors) = paste0(
    'Br', str_extract(names(sample_colors), '[0-9]{4}$')
)

################################################################################
#   Boxplots of QC metrics across cells, by sample
################################################################################

#   Boxplots of each metric by sample
p = colData(spe_cell)[, c('sample_id', metric_names)] |>
    as_tibble() |>
    pivot_longer(
        cols = all_of(metric_names),
        names_to = 'metric',
        values_to = 'value'
    ) |>
    mutate(
        sample_id = factor(
            paste0('Br', str_extract(sample_id, '[0-9]{4}$')),
            levels = c('Br9090', 'Br8433', 'Br8667', 'Br3942', 'Br9902')
        )
    ) |>
    #   While this throws off the boxplot components (e.g. the median), there
    #   appears to be no way to limit extreme values in a data-driven way in a
    #   faceted plot using ggplot2!
    group_by(metric) |>
    filter(value < quantile(value, 0.98)) |>
    ungroup() |>
    ggplot(aes(x = sample_id, y = value, fill = sample_id)) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = sample_colors) +
        facet_wrap(~ metric, scales = 'free_y') +
        theme_bw(base_size = 20) +
        theme(
            legend.position = 'none',
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        ) +
        labs(x = 'Sample ID', y = 'Value')
pdf(file.path(plot_dir, 'metrics_by_sample.pdf'), width = 9, height = 5)
print(p)
dev.off()

################################################################################
#   Barplots of number of cells by sample
################################################################################

metric_df = colData(spe_cell)[, c('sample_id', 'labels_joint_source')] |>
    as_tibble() |>
    group_by(sample_id, labels_joint_source) |>
    summarize(n_cells = n())

metric_df = rbind(
        metric_df,
        metric_df |>
            group_by(sample_id) |>
            summarize(n_cells = sum(n_cells)) |>
            mutate(labels_joint_source = 'total')
    ) |>
    ungroup() |>
    dplyr::rename(seg_type = labels_joint_source) |>
    left_join(
        tibble(sample_id = spe_bin$sample_id) |>
            group_by(sample_id) |>
            summarize(n_bins = n()),
        by = 'sample_id'
    ) |>
    mutate(
        n_cells_scaled = n_cells / n_bins,
        sample_id = factor(
            paste0('Br', str_extract(sample_id, '[0-9]{4}$')),
            levels = c('Br9090', 'Br8433', 'Br8667', 'Br3942', 'Br9902')
        )
    )

for (metric_name in c('n_cells', 'n_cells_scaled')) {
    p = metric_df |>
        ggplot(aes(x = sample_id, y = .data[[metric_name]], fill = sample_id)) +
            geom_bar(stat = "identity") +
            facet_wrap(~ seg_type, scales = 'free_y') +
            scale_fill_manual(values = sample_colors) +
            theme_bw(base_size = 20) +
            theme(
                legend.position = 'none',
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            ) +
            labs(x = 'Sample ID', y = metric_name, fill = 'Segmentation Type')
    pdf(
        file.path(plot_dir, paste0(metric_name, '_by_sample.pdf')),
        width = 9,
        height = 5
    )
    print(p)
    dev.off()
}

session_info()
