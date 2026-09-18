#   Prepare input data for CRAWDAD

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(crawdad)
library(rjson)
library(sessioninfo)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'crawdad',
    'input_cells.csv.gz'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)

dir.create(dirname(out_path), showWarnings = FALSE)

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
spe = readRDS(spe_path)

#   Ultimately, we'll be converting spatial coordinates to units of microns,
#   which is more interpretable than pixels
micron_per_px = c()
for (tissue_id in sample_info$tissue_id) {
    spaceranger_dir = sample_info$spaceranger_dir[
        sample_info$tissue_id == tissue_id
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
    sample_id = sample_info$tissue_id,
    micron_per_px = micron_per_px
)

#   Mapping from clusters to cell types
anno_df = read_csv(ct_anno_path, show_col_types = FALSE)

#   Gather spatial coordinates and Banksy clusters into a single CSV for input
#   to CRAWDAD
tibble(
        key = spe$key,
        sample_id = spe$sample_id,
        x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * scale_df$micron_per_px[
            match(sample_id, sample_info$tissue_id)
        ],
        y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * scale_df$micron_per_px[
            match(sample_id, sample_info$tissue_id)
        ]
    ) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    ) |>
    filter(cell_type != 'Drop') |>
    select(key, x, y, sample_id, cell_type) |>
    write_csv(out_path)

message('Memory usage:')
gc()

session_info()
