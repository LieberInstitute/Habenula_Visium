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
sr_metrics_path = '/dcs04/lieber/lcolladotor/with10x_LIBD001/ranger_metrics/processed-data/spaceranger/03_merged_metrics/merged_metrics.csv'
plot_dir = here('plots', '09_HD_cell_level')
marker_genes = list(
    white_matter = c("MBP", "GFAP", "PLP1", "AQP4"),
    habenula = c("POU4F1", "GPR151", "CHRNB4", "HTR2C"),
    thalamus = c("LYPD6B", "ADARB2", "RORB")
)

################################################################################
#   Compare quality metrics among DLPFC, HPC, and habenula
################################################################################

#-------------------------------------------------------------------------------
#   Perform in-silico dissection of the habenula using markers
#-------------------------------------------------------------------------------

spe_habenula = loadHDF5SummarizedExperiment(spe_habenula_norm_dir)
spe_habenula$exclude_overlapping = FALSE

habenula_markers = rownames(spe_habenula)[
    match(marker_genes[['habenula']], rowData(spe_habenula)$gene_name)
]

#   Use the Z-score method (https://github.com/LieberInstitute/spatialLIBD/blob/05e78863bbfae5addee2c80aca5a39c091d39145/R/multi_gene_z_score.R#L16)
#   for combining several habenula genes into a single expression metric. Then
#   use it to infer which cells belong to the habenula
x = assays(spe_habenula[habenula_markers,])$logcounts
spe_habenula$habenula_exp = unname(colMeans((x - rowMeans(x)) / rowSds(x)))
spe_habenula$is_habenula = spe_habenula$habenula_exp > 0

dir.create(file.path(plot_dir, 'QC'), showWarnings = FALSE)
p <- vis_clus(
        spe_habenula, clustervar = 'is_habenula', point_size = 1,
        is_stitched = TRUE, spatial = FALSE
    ) +
    guides(fill = guide_legend(override.aes = list(size = 4)))

png(file.path(plot_dir, 'QC', "is_habenula.png"), width = 1500, height = 1500)
print(p)
dev.off()

#   Form into a tibble for joining with the raw SPE later
habenula_df = colData(spe_habenula) |>
    as_tibble() |>
    mutate(key = colnames(spe_habenula)) |>
    select(key, habenula_exp, is_habenula)

#-------------------------------------------------------------------------------
#   Gather metrics from habenula, DLPFC, and HPC into a single tibble
#-------------------------------------------------------------------------------

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

#   Gather habenula metrics. Include some duplicate rows: include the full
#   habenula sample ("habenula") and the subset we think is actually habenula
#   ("habenula_subset") as separate rows
habenula_cd = colData(spe_habenula) |>
    as_tibble() |>
    mutate(key = colnames(spe_habenula)) |>
    left_join(habenula_df, by = "key")

habenula_cd = rbind(
        habenula_cd |> 
            mutate(region = "habenula"),
        habenula_cd |>
            filter(is_habenula) |>
            mutate(region = "habenula_subset")
    ) |>
    select(region, sum_umi, sum_gene, expr_chrM_ratio)

metrics = rbind(metrics, habenula_cd)

#-------------------------------------------------------------------------------
#   Explore metrics in habenula and other regions using boxplots and violin
#   plots
#-------------------------------------------------------------------------------

plot_list_comparison = list()
plot_list_violin = list()
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
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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

#-------------------------------------------------------------------------------
#   Plot similar spaceranger metrics from standard Visium datasets of the same
#   regions as a reference 
#-------------------------------------------------------------------------------

sr_metrics = read_csv(sr_metrics_path, show_col_types = FALSE) |>
    filter(
        study_name %in% c('spatial_hpc/spaceranger_2022-04-12_SPag033122', 'spatial_hpc/spaceranger_novaseq', 'Habenula_Visium', 'spatialDLPFC')
    ) |>
    mutate(
        study_name = ifelse(
            grepl('^spatial_hpc', study_name), 'spatial_hpc', study_name
        )
    )

plot_list = list()
for (metric in c('Median.UMI.Counts.per.Spot', 'Median.Genes.per.Spot')) {
    plot_list[[metric]] = ggplot(
            sr_metrics,
            aes(x = study_name, y = !!sym(metric), color = study_name)
        ) +
        geom_boxplot() +
        theme_bw(base_size = 15) +
        guides(color = "none") +
        labs(title = metric)
}

pdf(
    file.path(plot_dir, 'QC', 'spaceranger_standard_visium.pdf'),
    width = 10, height = 5
)
plot_grid(plotlist = plot_list, nrow = 1)
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
