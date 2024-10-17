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
spe_dlpfc_hpc_dir = "/dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/processed-data/07_cell_level/spe_raw"
spe_habenula_dir = here('processed-data', '09_HD_cell_level', 'spe_raw')
plot_dir = here('plots', '09_HD_cell_level')
wm_genes <- c("MBP", "GFAP", "PLP1", "AQP4")

################################################################################
#   Compare quality metrics among DLPFC, HPC, and habenula
################################################################################

dlpfc_hpc_sample_info = read_csv(
    dlpfc_hpc_sample_info_path, show_col_types = FALSE
)

spe_dlpfc_hpc = loadHDF5SummarizedExperiment(spe_dlpfc_hpc_dir)
spe_habenula = loadHDF5SummarizedExperiment(spe_habenula_dir)
spe_dlpfc_hpc = spe_dlpfc_hpc[, spe_dlpfc_hpc$in_tissue]
spe_habenula = spe_habenula[, spe_habenula$in_tissue]

#   Gather DLPFC and HPC cell-level metrics
metrics = colData(spe_dlpfc_hpc) |>
    as_tibble() |>
    mutate(
        region = dlpfc_hpc_sample_info$region[
            match(sample_id, dlpfc_hpc_sample_info$sample_id)
        ]
    ) |>
    select(region, sum_umi, sum_gene, expr_chrM_ratio)

#   Add habenula metrics
metrics = metrics |>
    rbind(
        colData(spe_habenula) |>
            as_tibble() |>
            mutate(region = "habenula") |>
            select(region, sum_umi, sum_gene, expr_chrM_ratio)
    )

#   Plot all 3 metric as boxplots in a single row
plot_list = list()
for (metric in c('sum_umi', 'sum_gene', 'expr_chrM_ratio')) {
    #   Get the top of the highest top whisker
    y_max = metrics |>
        group_by(region) |>
        summarize(top = boxplot.stats(!!sym(metric))$stats[5]) |>
        summarize(top = max(top)) |>
        pull(top)

    plot_list[[metric]] = ggplot(
            metrics, aes(x = region, y = !!sym(metric), color = region)
        ) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw(base_size = 15) +
        guides(color = "none") +
        labs(title = metric) +
        coord_cartesian(ylim = c(0, y_max))
}
pdf(file.path(plot_dir, 'quality_comparison.pdf'), width = 10, height = 5)
plot_grid(plotlist = plot_list, nrow = 1)
dev.off()

################################################################################
#   Plot quality metrics and white matter spatially for habenula
################################################################################

spe_habenula$exclude_overlapping = FALSE

#   Individually plot several quality metrics
for (metric in c('sum_umi', 'sum_gene', 'expr_chrM_ratio')) {
    p <- vis_gene(
        spe_habenula, geneid = metric, point_size = 1, is_stitched = TRUE,
        spatial = FALSE
    )

    png(
        file.path(plot_dir, paste0(metric, ".png")),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

#   Plot a combination of white-matter genes
wm_genes <- rownames(spe_habenula)[
    match(wm_genes, rowData(spe_habenula)$gene_name)
]
p <- vis_gene(
    spe_habenula, geneid = wm_genes, multi_gene_method = "pca", is_stitched = TRUE,
    point_size = 1, spatial = FALSE
)

png(file.path(plot_dir, "white_matter.png"), width = 1500, height = 1500)
print(p)
dev.off()

session_info()
