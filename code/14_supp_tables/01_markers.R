#   Gather all mean-ratio markers throughout the habenula atlas project into
#   one CSV

library(here)
library(tidyverse)
library(sessioninfo)

hd_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA', 'gene_sets',
    'fine.tsv'
)
hd_extra_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'gene_sets', 'fine.tsv'
)
k8_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2','ficture_harmony',
    'MAGMA', 'gene_sets', 'ficturek8.tsv'
)
k17_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', 'gene_sets', 'k17.tsv'
)
multiome_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/gene_sets/fine.tsv'
multiome_map_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/cell_type_map.csv'
hd_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
out_path = here('processed-data', '14_supp_tables', 'markers.csv')

dir.create(dirname(out_path), showWarnings = FALSE)

hd_map_df = read_csv(hd_map_path, show_col_types = FALSE) |>
    mutate(old_cell_type = str_replace(old_cell_type, '/', '.'))

hd_cell_df = read_tsv(hd_cell_path, show_col_types = FALSE) |>
    left_join(hd_map_df, by = c('set_id' = 'old_cell_type')) |>
    select(new_cell_type, gene_id, gene_name, MeanRatio) |>
    dplyr::rename(set_id = new_cell_type) |>
    mutate(dataset = 'HD_cell_types')

hd_extra_df = read_tsv(hd_extra_path, show_col_types = FALSE) |>
    left_join(hd_map_df, by = c('set_id' = 'old_cell_type')) |>
    select(new_cell_type, gene_id, gene_name, MeanRatio) |>
    dplyr::rename(set_id = new_cell_type) |>
    mutate(dataset = 'HD_extracellular')

k8_df = read_tsv(k8_path, show_col_types = FALSE) |>
    mutate(
        set_id = str_replace(set_id, '^X', 'Factor_'),
        dataset = 'HD_FICTURE_all_bin_k8'
    )

k17_df = read_tsv(k17_path, show_col_types = FALSE) |>
    mutate(
        set_id = str_replace(set_id, '^X', 'Factor_'),
        dataset = 'HD_FICTURE_extracellular_k17'
    )

multiome_df = read_tsv(multiome_path, show_col_types = FALSE) |>
    left_join(read_csv(multiome_map_path, show_col_types = FALSE), by = c('set_id' = 'old_cell_type')) |>
    select(new_cell_type, gene_id, gene_name, mean_ratio) |>
    dplyr::rename(set_id = new_cell_type, MeanRatio = mean_ratio) |>
    mutate(dataset = 'Multiome_cell_types')

marker_df = bind_rows(hd_cell_df, hd_extra_df, k8_df, k17_df, multiome_df) |>
    select(dataset, set_id, gene_id, gene_name, MeanRatio) |>
    dplyr::rename(mean_ratio = MeanRatio) |>
    arrange(dataset, set_id, desc(mean_ratio))

stopifnot(!any(is.na(marker_df)))

write_csv(marker_df, out_path)

session_info()
