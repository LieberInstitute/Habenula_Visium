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
plot_dir = here('plots', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region')
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
    'habenula_shiny_annotations.csv'
)
thalamus_anno_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region',
    'thalamus_shiny_annotations.csv'
)
cor_index = 13

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

sample_ids = readLines(sample_id_path)

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

#   Gather spatial coordinates, Banksy clusters, and region annotations into a
#   single CSV for input to CRAWDAD
tibble(
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
    ) |>
    select(x, y, sample_id, region_anno, cell_type) |>
    write_csv(out_path)

message('Memory usage:')
gc()

session_info()
