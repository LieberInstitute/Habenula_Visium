library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(rjson)
library(sessioninfo)
library(zellkonverter)

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', 'lambda0_2',
    'leiden_res1_5_subset.csv'
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
    'processed-data', '10_HD_bin_level', 'LIANA'
)

scales = c(100, 200, 500, 1000, 5000)
random_seed = 0
cor_index = 15

#   Compute a reference table matching clusters to fine cell types
anno_df = readRDS(cor_path)[[cor_index]] |>
    annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
    as_tibble() |>
    mutate(
        layer_label = ifelse(
            layer_confidence == 'good', layer_label, paste0('C',cluster)
        )
    )

cell_df_list <- vector("list", 3)

for(i in 1:3){
sample_id = readLines(sample_id_path)[i]
spe = loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id == sample_id]
micron_per_px = fromJSON(file = sprintf(scalefactors_path, sample_id))[['microns_per_pixel']]

#   Gather spatial coordinates and Banksy clusters
cell_df = tibble(
        key = spe$key,
        x = spatialCoords(spe)[, 'pxl_col_in_fullres'] * micron_per_px,
        y = spatialCoords(spe)[, 'pxl_row_in_fullres'] * micron_per_px
    ) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = factor(
            anno_df$layer_label[match(banksy_lambda0_2, anno_df$cluster)]
        )
    ) |>
    select(key,x,y,cell_type) |>
    as.data.frame()
stopifnot(!any(is.na(cell_df$cell_type)))

cell_df$sample_id = sample_id
cell_df_list[[i]] <- cell_df
}
cell_df_merged <- do.call(rbind, cell_df_list)

spe = loadHDF5SummarizedExperiment(spe_dir)
cell_df_merged <- left_join(cell_df_merged,
    as.data.frame(colData(spe)[,c("key","array_row","array_col")]),
    by = 'key'
)

write.csv(
    cell_df_merged, file = file.path(out_path, 'cell_df.csv'),
    row.names = FALSE
)

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