# Data preprocessing for LIANA
# 1. Create a cell dataframe with x, y coordinates and cell types
# 2. Save the dataframe as a CSV file for LIANA input

library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(rjson)
library(sessioninfo)
library(zellkonverter)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_5.csv'
)
cor_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'lambda0_2', 'cor_vs_snRNAseq_fine.rds'
)
scalefactors_path = here(
    'processed-data', '01_spaceranger', "five_samples_10_2025",'%s', 'outs',
    'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'liana',"adata"
)

ct_anno_path = here(
  "processed-data", "09_HD_cell_level", "new_samples2", "registration_banksy",
  "cluster_annotation.csv"
)
anno_df = read_csv(ct_anno_path, show_col_types = FALSE) |>
  mutate(cluster = as.character(cluster)) 

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
sample_ids = sample_info$sample_id

cell_df_list <- vector("list", length(sample_ids))

for(i in seq_along(sample_ids)){
  sample_id = sample_ids[i]

  spe = readRDS(spe_path)
  spe = spe[, spe$sample_id == sample_id]

  micron_per_px = fromJSON(file = sprintf(scalefactors_path, sample_id))[["microns_per_pixel"]]

  banksy_df = read_csv(banksy_path, show_col_types = FALSE) |>
    mutate(banksy_lambda0_2 = as.character(banksy_lambda0_2)) 

  cell_df = tibble(
    key = spe$key,
    sample_id = spe$sample_id,
    x = spatialCoords(spe)[, "pxl_col_in_fullres"] * micron_per_px,
    y = spatialCoords(spe)[, "pxl_row_in_fullres"] * micron_per_px
  ) |>
    left_join(banksy_df, by = "key") |>
    mutate(
      cell_type = anno_df$fine_cell_type[
        match(banksy_lambda0_2, anno_df$cluster)
      ]
    ) |>
    filter(!is.na(cell_type), cell_type != "Ambig") |>
    select(key, sample_id, x, y, cell_type) |>
    as.data.frame()

  ##filter cell_type NA
  stopifnot(!any(is.na(cell_df$cell_type)))

  cell_df_list[[i]] <- cell_df
}

cell_df_merged <- do.call(rbind, cell_df_list)

## add array_row/array_col back
spe = readRDS(spe_path)
cell_df_merged <- left_join(
  cell_df_merged,
  as.data.frame(colData(spe)[, c("key", "array_row", "array_col")]),
  by = "key"
)

write.csv(
  cell_df_merged,
  file = file.path(out_path, "cell_df.csv"),
  row.names = FALSE
)

# ct_anno_path = here(
#     'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
#     'cluster_annotation.csv'
# )

# scales = c(100, 200, 500, 1000, 5000)
# random_seed = 0
# cor_index = 15

# #   Mapping from clusters to cell types
# anno_df = read_csv(ct_anno_path, show_col_types = FALSE)

# cell_df = tibble(
#         key = spe$key,
#         sample_id = spe$sample_id,
#         region_anno = spe$region_anno,
#         x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * scale_df$micron_per_px[
#             match(sample_id, sample_info$sample_id)
#         ],
#         y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * scale_df$micron_per_px[
#             match(sample_id, sample_info$sample_id)
#         ]
#     ) |>
#     left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
#     mutate(
#         cell_type = anno_df$fine_cell_type[
#             match(as.character(banksy_lambda0_2), anno_df$cluster)
#         ]
#     ) |>
#     filter(cell_type != 'Ambig')


# #   Compute a reference table matching clusters to fine cell types
# anno_df = readRDS(cor_path)[[cor_index]] |>
#     annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
#     as_tibble() |>
#     mutate(
#         layer_label = ifelse(
#             layer_confidence == 'good', layer_label, paste0('C',cluster)
#         )
#     )

# cell_df_list <- vector("list", 5)

# sample_info = read_csv(sample_info_path, show_col_types = FALSE)
# sample_ids = sample_info$sample_id

# for(i in 1:5){
# sample_id = sample_ids[i]
# spe = readRDS(spe_path)
# spe = spe[, spe$sample_id == sample_id]
# micron_per_px = fromJSON(file = sprintf(scalefactors_path, sample_id))[['microns_per_pixel']]

# #   Gather spatial coordinates and Banksy clusters
# cell_df = tibble(
#         key = spe$key,
#         x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * micron_per_px,
#         y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * micron_per_px
#     ) |>
#     left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
#     mutate(
#         cell_type = factor(
#             anno_df$layer_label[match(banksy_lambda0_2, anno_df$cluster)]
#         )
#     ) |>
#     select(key,x,y,cell_type) |>
#     as.data.frame()
# stopifnot(!any(is.na(cell_df$cell_type)))

# cell_df$sample_id = sample_id
# cell_df_list[[i]] <- cell_df
# }
# cell_df_merged <- do.call(rbind, cell_df_list)

# spe = readRDS(spe_path)
# cell_df_merged <- left_join(cell_df_merged,
#     as.data.frame(colData(spe)[,c("key","array_row","array_col")]),
#     by = 'key'
# )

# write.csv(
#     cell_df_merged, file = file.path(out_path, 'cell_df.csv'),
#     row.names = FALSE
# )

# *************************************************************************
# spe_list <- vector("list", 3)

# for(i in 1:3){
# sample_id = readLines(sample_id_path)[i]
# spe = loadHDF5SummarizedExperiment(spe_dir)
# spe = spe[, spe$sample_id == sample_id]
# micron_per_px = fromJSON(file = sprintf(scalefactors_path, sample_id))[['microns_per_pixel']]

# #   Gather spatial coordinates and Banksy clusters
# cell_df = tibble(
#         key = spe$key,
#         x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * micron_per_px,
#         y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * micron_per_px
#     ) |>
#     left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
#     mutate(
#         cell_type = factor(
#             anno_df$layer_label[match(banksy_lambda0_2, anno_df$cluster)]
#         )
#     ) |>
#     select(key, cell_type) |>
#     as.data.frame()
# stopifnot(!any(is.na(cell_df$cell_type)))

# colData(spe)$cell_type=left_join(
#     as.data.frame(colData(spe)),
#     cell_df,
#     by = 'key'
# )$cell_type
# colData(spe)$cluster=left_join(
#     as.data.frame(colData(spe)),
#     cell_df,
#     by = 'key'
# )$banksy_lambda0_2
# stopifnot(!any(is.na(colData(spe)$cell_type)))

# spe_list[[i]] <- spe
# }
# spe_merged <- do.call(cbind, spe_list)

# writeH5AD(spe_merged, "out.h5ad")