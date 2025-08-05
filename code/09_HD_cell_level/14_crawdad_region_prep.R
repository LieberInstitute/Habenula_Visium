#   We're interested in running CRAWDAD on only the habenula or only thalamus.
#   This script reads in annotations of the regions from Shiny, and otherwise
#   prepares input data for CRAWDAD

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(crawdad)
library(rjson)
library(sessioninfo)

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region', 'region_anno'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_3_subset.csv'
)
cor_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'lambda0_2', 'cor_vs_snRNAseq_fine_subset.rds'
)
scalefactors_path = here(
    'processed-data', '01_spaceranger', 'probe_fix', '%s', 'outs',
    'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region',
    'input_cells.csv.gz'
)
habenula_anno_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region',
    'habenula_shiny_annotations.csv.gz'
)
thalamus_anno_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region',
    'thalamus_shiny_annotations.csv.gz'
)
cor_index = 13
region_colors = c(
    habenula = "#B1092D", thalamus = "#0B52C4", other = "#DFE1DD"
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

sample_ids = readLines(sample_id_path)[1:3]

spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id %in% sample_ids]

#   Ultimately, we'll be converting spatial coordinates to units of microns,
#   which is more interpretable than pixels
micron_per_px = c()
for (sample_id in sample_ids) {
    micron_per_px = c(
        micron_per_px,
        fromJSON(file = sprintf(scalefactors_path, sample_id))[['microns_per_pixel']]
    )
}
scale_df = tibble(
    sample_id = sample_ids,
    micron_per_px = micron_per_px
)

#   Compute a reference table matching clusters to fine cell types
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        layer_label = ifelse(
            layer_confidence == 'good', layer_label, paste0('C',cluster)
        )
    )

#   Read in Shiny annotations of habenula and thalamus regions and attach to
#   the SpatialExperiment object
region_df = rbind(
        read_csv(habenula_anno_path, show_col_types = FALSE),
        read_csv(thalamus_anno_path, show_col_types = FALSE)
    ) |>
    dplyr::rename(key = spot_name) |>
    #   It's possible for a cell to be annotated as habenula and thalamus. In
    #   this case, call it habenula
    group_by(key) |>
    arrange(ManualAnnotation) |>
    slice_head(n = 1) |>
    ungroup()

spe$region_anno = tibble(key = spe$key) |>
    left_join(region_df, by = 'key') |>
    pull(ManualAnnotation) |>
    replace_na('other')

#   Plot the region annotation on each sample to make sure it worked
for (sample_id in sample_ids) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = 'region_anno',
                is_stitched = TRUE, point_size = 20, spatial = FALSE,
                colors = region_colors
            ) +
                guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(
        file.path(plot_dir, sprintf('%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

#   Gather spatial coordinates, Banksy clusters, and region annotations into a
#   single CSV for input to CRAWDAD
cell_df = tibble(
        key = spe$key,
        sample_id = spe$sample_id,
        region_anno = spe$region_anno,
        x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * scale_df$micron_per_px[
            match(sample_id, sample_ids)
        ],
        y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * scale_df$micron_per_px[
            match(sample_id, sample_ids)
        ]
    ) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = factor(
            anno_df$layer_label[match(banksy_lambda0_2, anno_df$cluster)]
        )
    )

#   Signal to drop combinations of cell type and region that consitute less than
#   1% of the region's cells
cell_counts_df = cell_df |>
    group_by(region_anno, cell_type) |>
    summarize(n = n()) |>
    group_by(region_anno) |>
    mutate(drop = n < sum(n) * 0.01) |>
    ungroup()

#   Clean up and export
cell_df |>
    left_join(cell_counts_df, by = c('region_anno', 'cell_type')) |>
    select(key, x, y, sample_id, region_anno, cell_type, drop) |>
    write_csv(out_path)

message('Memory usage:')
gc()

session_info()
