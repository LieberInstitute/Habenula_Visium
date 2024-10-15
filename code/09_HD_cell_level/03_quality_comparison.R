#   Compare the DLPFC, HPC, and habenula Visium HD samples on various quality
#   metrics

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)

dlpfc_hpc_sample_info_path = "/dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/raw-data/sample_info/sample_info.csv"
spe_dlpfc_hpc_bin_dir = "/dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/processed-data/06_bin_level/spe_norm"
spe_dlpfc_hpc_cell_dir = "/dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/processed-data/07_cell_level/spe_norm"
spe_habenula_bin_dir = here('processed-data', '10_HD_bin_level', 'spe_norm')
spe_habenula_cell_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
bin_plot_dir = here('plots', '10_HD_bin_level')
cell_plot_dir = here('plots', '09_HD_cell_level')

dir.create(bin_plot_dir, showWarnings = FALSE)

dlpfc_hpc_sample_info = read_csv(
    dlpfc_hpc_sample_info_path, show_col_types = FALSE
)

################################################################################
#   Bin-level metrics
################################################################################

message(Sys.time(), " - Gathering bin-level metrics")

spe_dlpfc_hpc_bin = loadHDF5SummarizedExperiment(spe_dlpfc_hpc_bin_dir)
spe_habenula_bin = loadHDF5SummarizedExperiment(spe_habenula_bin_dir)

#   Gather DLPFC and HPC bin-level metrics
bin_metrics = colData(spe_dlpfc_hpc_bin) |>
    as_tibble() |>
    mutate(
        region = dlpfc_hpc_sample_info$region[
            match(sample_id, dlpfc_hpc_sample_info$sample_id)
        ]
    ) |>
    select(region, sum_umi, sum_gene, expr_chrM_ratio)

#   Add habenula metrics
bin_metrics = bin_metrics |>
    rbind(
        colData(spe_habenula_bin) |>
            as_tibble() |>
            mutate(region = "habenula") |>
            select(region, sum_umi, sum_gene, expr_chrM_ratio)
    )

#   Plot all 3 metric as boxplots in a single row
plot_list = list()
for (metric in c('sum_umi', 'sum_gene', 'expr_chrM_ratio')) {
    plot_list[[metric]] = ggplot(
            bin_metrics, aes(x = region, y = !!sym(metric), color = region)
        ) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw(base_size = 15) +
        guides(color = "none") +
        labs(title = metric) +
        coord_cartesian(
            ylim = c(0, boxplot.stats(bin_metrics[[metric]])$stats[5] * 1.1)
        )
}
pdf(file.path(bin_plot_dir, 'quality_comparison.pdf'), width = 10, height = 5)
plot_grid(plotlist = plot_list, nrow = 1)
dev.off()

