#   Using QC info collected on the 8um bin-level data, filter out problematic
#   cells and log normalize counts. Save a final SpatialExperiment

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(scran)
library(sessioninfo)
library(spatialLIBD)

spe_raw_path = here('processed-data', '09_HD_cell_level', 'new_samples', 'spe_raw.rds')
spe_norm_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples', 'spe_norm_filtered.rds'
)
plot_dir = here('plots', '09_HD_cell_level', 'new_samples', 'QC')
min_umi_cutoff = 10
H1_MVPY9BW_A1_8433_artifact = 34223
H1_XQQD7C7_A1_8518_horizontal_artifact = 28588
H1_XQQD7C7_A1_8518_vertical_artifact = 32766

dir.create(
    file.path(plot_dir, 'before'), recursive = TRUE, showWarnings = FALSE
)
dir.create(file.path(plot_dir, 'after'), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

spatial_qc_plots = function(spe, plot_dir) {
    #   Plot the spatial distribution of UMI and gene counts, colored to best
    #   see low values of each. Also plot mitochondrial ratio normally
    for (sample_id in unique(spe$sample_id)) {
        for (metric in c('sum_umi_capped', 'sum_gene_capped')) {
            p = vis_gene(
                    spe, sampleid = sample_id, geneid = metric,
                    is_stitched = TRUE, point_size = 1, spatial = TRUE,
                    assayname = 'counts'
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
            is_stitched = TRUE, point_size = 1, spatial = TRUE,
            assayname = 'counts'
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
#   Load and calculate metrics for plotting
################################################################################

message(Sys.time(), " | Loading and computing QC metrics")
spe = readRDS(spe_raw_path)
spe$exclude_overlapping = FALSE

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

#   Print the median 'sum_umi' value by sample, and what percentage of cells
#   will be dropped by the minimum cutoff
message(
    sprintf(
        "The following table shows median UMI, and percentage of cells dropped using a minimum cutoff of %s UMI:",
        min_umi_cutoff
    )
)
colData(spe) |>
    as_tibble() |>
    select(sample_id, sum_umi) |>
    group_by(sample_id) |>
    summarize(
        med_umi = median(sum_umi),
        perc_dropped = 100 * mean(sum_umi < min_umi_cutoff)
    ) |>
    print()

spatial_qc_plots(spe, file.path(plot_dir, 'before'))

################################################################################
#   Log normalize and perform basic filtering
################################################################################

#   Filter SPE: drop cells with low counts for all genes, and drop genes with 0
#   counts in every cell
message(Sys.time(), " | Filtering genes and spots")
spe <- spe[rowSums(assays(spe)$counts) > 0, spe$sum_umi >= min_umi_cutoff]

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large)
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)

################################################################################
#   Filter out problematic cells
################################################################################

#   Filter out the artifact in H1-MVPY9BW_A1_8433 using info gained in bin-level
#   QC
spe = spe[
    ,
    (spe$sample_id != 'H1-MVPY9BW_A1_8433') |
    (spatialCoords(spe)[, 'pxl_col_in_fullres'] <= H1_MVPY9BW_A1_8433_artifact)
]

#   Filter out the artifacts in H1-XQQD7C7_A1_8518 using info gained in
#   bin-level QC
spe = spe[
    ,
    (spe$sample_id != 'H1-XQQD7C7_A1_8518') |
    (
        (spatialCoords(spe)[, 'pxl_row_in_fullres'] <= H1_XQQD7C7_A1_8518_horizontal_artifact) &
        (spatialCoords(spe)[, 'pxl_col_in_fullres'] <= H1_XQQD7C7_A1_8518_vertical_artifact)
    )
]

spatial_qc_plots(spe, file.path(plot_dir, 'after'))

#   Save normalized object
message(Sys.time(), " | Saving normalized and QCd SPE")
saveRDS(spe, spe_norm_path)

session_info()
