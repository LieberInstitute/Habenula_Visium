#   Explore how the new 2 samples compare to others in terms of notable metrics

library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)
library(spatialLIBD)

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

p = colData(spe_cell)[, c('sample_id', metric_names)] |>
    as_tibble() |>
    pivot_longer(
        cols = all_of(metric_names),
        names_to = 'metric',
        values_to = 'value'
    ) |>
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
pdf(file.path(plot_dir, 'metrics_by_sample.pdf'))
print(p)
dev.off()
