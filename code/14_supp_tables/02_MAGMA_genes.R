#   Gather all mean-ratio markers throughout the habenula atlas project into
#   one CSV

library(here)
library(tidyverse)
library(sessioninfo)

hd_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'top_genes.csv'
)
hd_extra_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'top_genes.csv'
)
k8_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'MAGMA', 'top_genes.csv'
)
k17_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', 'top_genes_k17.csv'
)
multiome_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/top_genes.csv'
multiome_map_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/cell_type_map.csv'
hd_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')
out_path = here('processed-data', '14_supp_tables', 'magma_top_genes.csv')

dir.create(dirname(out_path), showWarnings = FALSE)

hd_cell_df = read_csv(hd_cell_path, show_col_types = FALSE) |>
    filter(cell_type_res == 'fine') |>
    left_join(
        read_csv(hd_map_path, show_col_types = FALSE),
        by = c('cell_type' = 'old_cell_type')
    ) |>
    select(new_cell_type, gwas, gene_id, gene_name, p) |>
    dplyr::rename(cell_type = new_cell_type) |>
    mutate(dataset = 'HD_cell_types')

hd_extra_df = read_csv(hd_extra_path, show_col_types = FALSE) |>
    filter(cell_type_res == 'fine') |>
    left_join(
        read_csv(hd_map_path, show_col_types = FALSE),
        by = c('cell_type' = 'old_cell_type')
    ) |>
    select(new_cell_type, gwas, gene_id, gene_name, p) |>
    dplyr::rename(cell_type = new_cell_type) |>
    mutate(dataset = 'HD_extracellular')

k8_df = read_csv(k8_path, show_col_types = FALSE) |>
    mutate(
        cell_type = str_replace(cell_type, '^X', 'Factor_'),
        dataset = 'HD_FICTURE_all_bin_k8'
    )

k17_df = read_csv(k17_path, show_col_types = FALSE) |>
    mutate(
        cell_type = str_replace(cell_type, '^X', 'Factor_'),
        dataset = 'HD_FICTURE_extracellular_k17'
    )

multiome_df = read_csv(multiome_path, show_col_types = FALSE) |>
    filter(cell_type_res == 'fine') |>
    left_join(
        read_csv(multiome_map_path, show_col_types = FALSE),
        by = c('cell_type' = 'old_cell_type')
    ) |>
    select(new_cell_type, gwas, gene_id, gene_name, p) |>
    dplyr::rename(cell_type = new_cell_type) |>
    mutate(dataset = 'Multiome_cell_types')

marker_df = bind_rows(hd_cell_df, hd_extra_df, k8_df, k17_df, multiome_df) |>
    select(dataset, cell_type, gwas, gene_id, gene_name, p) |>
    dplyr::rename(set_id = cell_type) |>
    arrange(dataset, set_id, p)

stopifnot(!any(is.na(marker_df)))

write_csv(marker_df, out_path)

session_info()
