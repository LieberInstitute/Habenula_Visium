#   Compare the DLPFC, HPC, and habenula Visium HD samples on various quality
#   metrics

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)
library(scran)
library(viridis)

dlpfc_hpc_sample_info_path = "/dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/raw-data/sample_info/sample_info.csv"
spe_dlpfc_hpc_dir = "/dcs05/lieber/lcolladotor/Visium_HD_DLPFC_pilot_LIBD4100/Visium_HD_DLPFC_pilot/processed-data/07_cell_level/spe_raw"
spe_habenula_dir = here('processed-data', '09_HD_cell_level', 'spe_raw')
spe_habenula_norm_dir = here('processed-data', '09_HD_cell_level', 'spe_norm')
plot_dir = here('plots', '09_HD_cell_level')
marker_genes = list(
    white_matter = c("MBP", "GFAP", "PLP1", "AQP4"),
    habenula = c("POU4F1", "GPR151", "CHRNB4", "HTR2C"),
    thalamus = c("LYPD6B", "ADARB2", "RORB")
)

################################################################################
#   Compare quality metrics among DLPFC, HPC, and habenula
################################################################################

dlpfc_hpc_sample_info = read_csv(
    dlpfc_hpc_sample_info_path, show_col_types = FALSE
)

spe_dlpfc_hpc = loadHDF5SummarizedExperiment(spe_dlpfc_hpc_dir)
spe_habenula = loadHDF5SummarizedExperiment(spe_habenula_dir)

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

plot_list_comparison = list()
plot_list_violin = list()
dir.create(file.path(plot_dir, 'QC'), showWarnings = FALSE)
for (metric in c('sum_umi', 'sum_gene', 'expr_chrM_ratio')) {
    #   Get the top of the highest top whisker
    y_max = metrics |>
        group_by(region) |>
        summarize(top = boxplot.stats(!!sym(metric))$stats[5]) |>
        summarize(top = max(top)) |>
        pull(top)

    #   Boxplot comparing each region
    plot_list_comparison[[metric]] = ggplot(
            metrics, aes(x = region, y = !!sym(metric), color = region)
        ) +
        geom_boxplot(outlier.shape = NA) +
        theme_bw(base_size = 15) +
        guides(color = "none") +
        labs(title = metric) +
        coord_cartesian(ylim = c(0, y_max))
    
    plot_list_violin[[metric]] = metrics |>
        filter(region == "habenula") |>
        ggplot(aes(x = region, y = !!sym(metric))) +
            geom_violin() +
            labs(title = metric) +
            theme_bw(base_size = 15)
}

#   Plot all 3 metric as boxplots (with all regions) in a single row
pdf(
    file.path(plot_dir, 'QC', 'quality_comparison.pdf'),
    width = 10, height = 5
)
plot_grid(plotlist = plot_list_comparison, nrow = 1)
dev.off()

#   Plot all 3 metric as violin plots for just habenula in a single row
pdf(file.path(plot_dir, 'QC', 'habenula_violin.pdf'))
plot_grid(plotlist = plot_list_violin, nrow = 1)
dev.off()

################################################################################
#   Spatial plots
################################################################################

spe_habenula = loadHDF5SummarizedExperiment(spe_habenula_norm_dir)
spe_habenula$exclude_overlapping = FALSE

#-------------------------------------------------------------------------------
#   Quality metrics
#-------------------------------------------------------------------------------

#   Individually plot several quality metrics
for (metric in c('sum_umi', 'sum_gene', 'expr_chrM_ratio')) {
    p <- vis_gene(
        spe_habenula, geneid = metric, point_size = 1, is_stitched = TRUE,
        spatial = FALSE
    )

    png(
        file.path(plot_dir, 'QC', paste0(metric, ".png")),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()

    #   Do a version where high outliers are removed
    if (metric != "expr_chrM_ratio") {
        is_outlier = isOutlier(spe_habenula[[metric]], type = "higher")
        p <- vis_gene(
            spe_habenula[, !is_outlier], geneid = metric, point_size = 1,
            is_stitched = TRUE, spatial = FALSE
        )

        png(
            file.path(plot_dir, 'QC', paste0(metric, "_no_outliers.png")),
            width = 1500, height = 1500
        )
        print(p)
        dev.off()
    }
}

#-------------------------------------------------------------------------------
#   Marker genes
#-------------------------------------------------------------------------------

#   All markers should be measured
stopifnot(all(unlist(marker_genes) %in% rowData(spe_habenula)$gene_name))

#   Convert from gene symbol to Ensembl ID
marker_genes = lapply(
    marker_genes,
    function(x) {
        rownames(spe_habenula)[
            match(x, rowData(spe_habenula)$gene_name)
        ]
    }
)

#   Plot a combination of white-matter genes
dir.create(file.path(plot_dir, 'marker_genes'), showWarnings = FALSE)
for (marker_name in names(marker_genes)) {
    p <- vis_gene(
        spe_habenula, geneid = marker_genes[[marker_name]],
        multi_gene_method = "pca", is_stitched = TRUE,
        point_size = 1, spatial = FALSE
    )

    png(
        file.path(plot_dir, 'marker_genes', sprintf("%s.png", marker_name)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

session_info()
