#   Drop problematic bins from the 8um bin-level data to form a final QC'd
#   object

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)
library(scran)

plot_dir = here('plots', '10_HD_bin_level', 'new_samples', 'QC')
spe_in_dir = here('processed-data', '10_HD_bin_level', 'new_samples', 'spe_norm')
spe_out_dir = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'spe_norm_filtered'
)

spe = loadHDF5SummarizedExperiment(spe_in_dir)
spe$exclude_overlapping = FALSE

dir.create(
    file.path(plot_dir, 'before'), recursive = TRUE, showWarnings = FALSE
)
dir.create(file.path(plot_dir, 'after'), showWarnings = FALSE)
dir.create(file.path(plot_dir, 'artifacts'), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   Find the average capped sum UMI value within a range of 'array_coord' values
#   and return a tibble with this value
umi_within_window = function(lower_threshold, cd_small, window = 10) {
    mean_umi = cd_small |>
        filter(
            array_coord >= lower_threshold,
            array_coord < lower_threshold + window
        ) |>
        pull(sum_umi_capped) |>
        mean()

    temp = tibble(lower_threshold = lower_threshold, mean_umi = mean_umi)
    return(temp)
}

#   With a sliding window, scan average UMI within each window across a range of
#   starting values. The idea is that the artifact produces a sudden increase in
#   average UMI counts at the boundary and nowhere else
scan_window = function(cd_small, start_range, window = 10) {
    umi_df = do.call(
        rbind,
        lapply(
            start_range,
            umi_within_window, cd_small = cd_small, window = window
        )
    )
    return(umi_df)
}

spatial_qc_plots = function(spe, plot_dir) {
    #   Plot the spatial distribution of UMI and gene counts, colored to best
    #   see low values of each. Also plot mitochondrial ratio normally
    for (sample_id in unique(spe$sample_id)) {
        for (metric in c('sum_umi_capped', 'sum_gene_capped')) {
            p = vis_gene(
                    spe, sampleid = sample_id, geneid = metric,
                    is_stitched = TRUE, point_size = 1, spatial = TRUE
                ) +
                scale_color_viridis_c(direction = -1) +
                scale_fill_viridis_c(direction = -1)
            
            png(
                file.path(plot_dir, sprintf('%s_%s.png', metric, sample_id)),
                width = 1000, height = 1000
            )
            print(p)
            dev.off()
        }

        p = vis_gene(
            spe, sampleid = sample_id, geneid = 'expr_chrM_ratio',
            is_stitched = TRUE, point_size = 1, spatial = TRUE
        )
        
        png(
            file.path(plot_dir, sprintf('expr_chrM_ratio_%s.png', sample_id)),
            width = 1000, height = 1000
        )
        print(p)
        dev.off()
    }
}

################################################################################
#   Explore QC metrics spatially
################################################################################

#   The severe sparsity of the counts combined with high-valued outliers makes
#   the MAD method implemented by scran not meaningful for 'sum_umi' and
#   'sum_gene'. In particular, lower outliers must be negative, which is
#   obviously not possible. Fine tuning the numbers of MADs representing "low
#   quality" is also difficult to do in a data-driven way
stopifnot(median(spe$sum_umi) - 3 * sd(spe$sum_umi) < 0)
stopifnot(median(spe$sum_gene) - 3 * sd(spe$sum_gene) < 0)

#   Calculate a version of 'sum_umi' and 'sum_gene' that's capped at the
#   median value (by sample) as a maximum. This will enable more easily
#   visualizing low values of these metrics in the presence of high-valued
#   outliers
temp = colData(spe) |>
    as_tibble() |>
    select(sample_id, sum_umi, sum_gene) |>
    group_by(sample_id) |>
    mutate(
        sum_umi_capped = ifelse(
            sum_umi < median(sum_umi), sum_umi, median(sum_umi)
        ),
        sum_gene_capped = ifelse(
            sum_gene < median(sum_gene), sum_gene, median(sum_gene)
        )
    )
spe$sum_umi_capped = temp$sum_umi_capped
spe$sum_gene_capped = temp$sum_gene_capped

#   Save QC plots of several metrics
spatial_qc_plots(spe, file.path(plot_dir, 'before'))

################################################################################
#   Remove artifact in H1-MVPY9BW_A1_8433
################################################################################

#   Visually show where the UMI boundary is. In the plot, the boundary occurs at
#   783, which indicates bad coordinates start at 793 since a sliding window of
#   10 is used
p = colData(spe) |>
    as_tibble() |>
    filter(sample_id == 'H1-MVPY9BW_A1_8433') |>
    select(array_col, sum_umi_capped) |>
    dplyr::rename(array_coord = array_col) |>
    scan_window(760:810) |>
    ggplot(aes(x = lower_threshold, y = mean_umi)) +
        geom_line() +
        theme_bw(base_size = 25) +
        geom_vline(xintercept = 783) +
        labs(x = "Window Start", y = "Mean UMI in Window")

pdf(file.path(plot_dir, 'artifacts', 'H1-MVPY9BW_A1_8433_array_QC.pdf'))
print(p)
dev.off()

#   Array col and pixel col move in the same direction
small_spe = spe[, spe$sample_id == 'H1-MVPY9BW_A1_8433']
stopifnot(
    abs(1 - cor(small_spe$array_col, spatialCoords(small_spe)[, 'pxl_col_in_fullres']))
    < 1e-3
)

#   The cell-level SPE will also have problematic cells in the same region as
#   observed in the bin-level, but it doesn't have array coordinates. Instead,
#   find the pixel col of the boundary (since there's a direct correspondence
#   with array col) and print in the log for use in the cell-level script
spe_boundary = small_spe[, small_spe$array_col == 793]
message(
    sprintf(
        "The artifact in 'H1-MVPY9BW_A1_8433' occurs at values of 'pxl_col_in_fullres' above %s",
        round(median(spatialCoords(spe_boundary)[, 'pxl_col_in_fullres']))
    )
)

#   Filter out the artifact
spe = spe[, (spe$sample_id != 'H1-MVPY9BW_A1_8433') | (spe$array_col <= 793)]

################################################################################
#   Remove artifacts in H1-XQQD7C7_A1_8518
################################################################################

small_spe = spe[, spe$sample_id == 'H1-XQQD7C7_A1_8518']

#-------------------------------------------------------------------------------
#   Remove the array-row-related artifact (horizontal)
#-------------------------------------------------------------------------------

p = colData(small_spe) |>
    as_tibble() |>
    select(array_row, sum_umi_capped) |>
    dplyr::rename(array_coord = array_row) |>
    scan_window(740:780, window = 1) |>
    ggplot(aes(x = lower_threshold, y = mean_umi)) +
        geom_line() +
        theme_bw(base_size = 25) +
        geom_vline(xintercept = 764) +
        labs(x = "Window Start", y = "Mean UMI in Window")

pdf(
    file.path(
        plot_dir, 'artifacts', 'H1-XQQD7C7_A1_8518_array_QC_horizontal.pdf'
    )
)
print(p)
dev.off()

#   Array row and pixel row move in the same direction
stopifnot(
    abs(1 - cor(small_spe$array_row, spatialCoords(small_spe)[, 'pxl_row_in_fullres']))
    < 1e-3
)

#   Find the pixel row of the boundary (since there's a direct correspondence
#   with array row) and print in the log for use in the cell-level script
spe_boundary = small_spe[, small_spe$array_row == 764]
message(
    sprintf(
        "The artifact in 'H1-XQQD7C7_A1_8518' occurs at values of 'pxl_row_in_fullres' above %s",
        round(median(spatialCoords(spe_boundary)[, 'pxl_row_in_fullres']))
    )
)

#   Filter out the artifact
spe = spe[, (spe$sample_id != 'H1-XQQD7C7_A1_8518') | (spe$array_row <= 764)]

#-------------------------------------------------------------------------------
#   Remove the artifact on the right and an edge-related low-UMI strip
#-------------------------------------------------------------------------------

#   Considering UMI counts below the median, check the distribution. It's
#   clearly bimodal, and a cutoff of 5 UMI separates out the problematic bins
lower_umi = small_spe$sum_umi_capped[
    small_spe$sum_umi_capped < max(small_spe$sum_umi_capped)
]
p = ggplot(tibble(lower_umi = lower_umi), aes(x = lower_umi)) +
    geom_density() +
    geom_vline(xintercept = 5) +
    theme_bw(base_size = 15) +
    labs(x = "Total UMI Among Lower Half")
pdf(file.path(plot_dir, 'artifacts', 'H1-XQQD7C7_A1_8518_umi_cutoff.pdf'))
print(p)
dev.off()

#   Plot which bins are filtered out by the cutoff
spe$low_umi = spe$sum_umi < 5
p = vis_clus(
    spe, sampleid = 'H1-XQQD7C7_A1_8518', clustervar = 'low_umi',
    is_stitched = TRUE, point_size = 1, spatial = TRUE
)

png(
    file.path(plot_dir, 'artifacts', 'H1-XQQD7C7_A1_8518_low_umi.png'),
    width = 1000, height = 1000
)
print(p)
dev.off()

#   Filter out the artifacts
spe = spe[, (spe$sample_id != 'H1-XQQD7C7_A1_8518') | !spe$low_umi]
spe$low_umi = NULL

#-------------------------------------------------------------------------------
#   Remove the array-row-related artifact (vertical)
#-------------------------------------------------------------------------------

#   This artifact is trickier, because it actually can be defined by higher
#   average UMI, particularly at the top of the tissue (though the artifact
#   is present for the full range)

small_spe = spe[, spe$sample_id == 'H1-XQQD7C7_A1_8518']

p = colData(small_spe) |>
    as_tibble() |>
    #   Grab the very top of the tissue, where the artifact is most prominent
    filter(array_row < 200) |>
    select(array_col, sum_umi_capped) |>
    dplyr::rename(array_coord = array_col) |>
    scan_window(125:165, window = 5) |>
    ggplot(aes(x = lower_threshold, y = mean_umi)) +
        geom_line() +
        theme_bw(base_size = 25) +
        geom_vline(xintercept = 137) +
        labs(x = "Window Start", y = "Mean UMI in Window")

pdf(
    file.path(plot_dir, 'artifacts', 'H1-XQQD7C7_A1_8518_array_QC_vertical.pdf')
)
print(p)
dev.off()

#   Array col and pixel col move in the opposite direction
stopifnot(
    abs(-1 - cor(small_spe$array_col, spatialCoords(small_spe)[, 'pxl_col_in_fullres']))
    < 1e-3
)

#   Find the pixel col of the boundary (since there's a direct correspondence
#   with array col) and print in the log for use in the cell-level script
spe_boundary = small_spe[, small_spe$array_col == 142]
message(
    sprintf(
        "The artifact in 'H1-XQQD7C7_A1_8518' occurs at values of 'pxl_col_in_fullres' above %s",
        round(median(spatialCoords(spe_boundary)[, 'pxl_col_in_fullres']))
    )
)

#   Filter out the artifact
spe = spe[, (spe$sample_id != 'H1-XQQD7C7_A1_8518') | (spe$array_col > 142)]

################################################################################
#   Save object with problematic bins removed
################################################################################

#   Save QC plots of several metrics after filtering
spatial_qc_plots(spe, file.path(plot_dir, 'after'))

saveHDF5SummarizedExperiment(
    spe, dir = spe_out_dir, replace = TRUE, as.sparse = TRUE
)

session_info()
