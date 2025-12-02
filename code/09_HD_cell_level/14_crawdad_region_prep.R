#   We're interested in running CRAWDAD on only the habenula or only thalamus.
#   This script reads in annotations of the regions from Shiny, and otherwise
#   prepares input data for CRAWDAD

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(crawdad)
library(rjson)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
cor_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'lambda0_2', 'cor_vs_snRNAseq_fine.rds'
)
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'crawdad', 'region',
    'input_cells.csv.gz'
)
hb_thal_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'hb_thal_manual_anno.csv.gz'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'cluster_annotation.csv'
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
spe = readRDS(spe_path)

#   Ultimately, we'll be converting spatial coordinates to units of microns,
#   which is more interpretable than pixels
micron_per_px = c()
for (sample_id in sample_info$sample_id) {
    spaceranger_dir = sample_info$spaceranger_dir[
        sample_info$sample_id == sample_id
    ]
    scalefactors_path = here(
        spaceranger_dir, 'outs', 'binned_outputs', 'square_002um', 'spatial',
        'scalefactors_json.json'
    )
    micron_per_px = c(
        micron_per_px,
        fromJSON(file = scalefactors_path)[['microns_per_pixel']]
    )
}
scale_df = tibble(
    sample_id = sample_info$sample_id,
    micron_per_px = micron_per_px
)

#   Mapping from clusters to cell types
anno_df = read_csv(ct_anno_path, show_col_types = FALSE)

#   Read in Shiny annotations of habenula and thalamus regions and attach to
#   the SpatialExperiment object
region_df = read_csv(hb_thal_anno_path, show_col_types = FALSE) |>
    dplyr::rename(key = spot_name)

spe$region_anno = tibble(key = spe$key) |>
    left_join(region_df, by = 'key') |>
    pull(ManualAnnotation) |>
    replace_na('other')

#   Gather spatial coordinates, Banksy clusters, and region annotations into a
#   single CSV for input to CRAWDAD
cell_df = tibble(
        key = spe$key,
        sample_id = spe$sample_id,
        region_anno = spe$region_anno,
        x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * scale_df$micron_per_px[
            match(sample_id, sample_info$sample_id)
        ],
        y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * scale_df$micron_per_px[
            match(sample_id, sample_info$sample_id)
        ]
    ) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy_lambda0_2), anno_df$cluster)
        ]
    ) |>
    filter(cell_type != 'Ambig')

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
