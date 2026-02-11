#   Filter out problematic cells, perform log normalization, and import
#   annotated Banksy clusters

library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)
library(spatialLIBD)

spe_raw_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_raw.rds'
)
spe_norm_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.rds'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'MAGMA', 'extracellular',
    'spe_checks'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
min_umi_cutoff = 10
H1_MVPY9BW_A1_8433_artifact = 34326

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

#   vis_clus with HD settings and saving to file
vis_clus_HD = function(spe, sampleid, clustervar, plot_dir, ...) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sampleid, clustervar = clustervar,
                is_stitched = TRUE, point_size = 25, spatial = FALSE, ...
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    
    png(
        file.path(plot_dir, sprintf('%s.png', sampleid)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_raw_path)
spe$exclude_overlapping = FALSE

#   Filter out the artifact in H1-MVPY9BW_A1_8433 using info gained in bin-level
#   QC
spe = spe[
    ,
    (spe$sample_id != 'H1-MVPY9BW_A1_8433') |
    (spatialCoords(spe)[, 'pxl_col_in_fullres'] <= H1_MVPY9BW_A1_8433_artifact)
]

#   Filter SPE: drop cells with low counts for all genes, and drop genes with 0
#   counts in every cell
message(Sys.time(), " | Filtering genes and spots")
spe <- spe[rowSums(assays(spe)$counts) > 0, spe$sum_umi >= min_umi_cutoff]

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large)
message(Sys.time(), ' | Performing log normalization...')
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)

#   Add in annotated cell types
message(Sys.time(), " | Adding annotated cell types...")
anno_df = read_csv(ct_anno_path, show_col_types = FALSE)
spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    ) |>
    pull(cell_type)
message(
    sprintf(
        "Dropping %.2f%% of cells without a Banksy cluster",
        100 * mean(is.na(spe$cell_type))
    )
)
spe = spe[, !is.na(spe$cell_type)]

#   Now plot cell types on the off chance there was an issue matching keys
#   between cellular and extracellular objects
for (sample_id in unique(spe$sample_id)) {
    vis_clus_HD(spe, sample_id, 'cell_type', plot_dir)
}

saveRDS(spe, spe_norm_path)

session_info()
