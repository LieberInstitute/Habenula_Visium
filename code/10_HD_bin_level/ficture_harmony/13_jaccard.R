#   Explore how clustering results line up before and after filtering normalized
#   probes. One problematic cluster (composed of many normalized probes) was very
#   spatially scattered-- is this true after filtering normalized probes?

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(data.table)

spe_cleany_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/probe_fix/ficture_harmony/spe/y_clean_spe.rds"
# spe_normalized_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/probe_fix/spe_norm"

cluster_cleany_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
cluster_normalized_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'bin_level_clusters_normalized.csv.gz'
)

plot_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/plots/10_HD_bin_level/probe_fix/ficture_harmony"

get_coord_df <- function(spe) {
  tibble(
    x = round(spatialCoords(spe)[, "pxl_col_in_fullres"], 1),
    y = round(spatialCoords(spe)[, "pxl_row_in_fullres"], 1),
    sample_id = spe$sample_id,
    barcode = spe$barcode
  )
}

spe1 = readRDS(spe_cleany_dir)
spe1$barcode = colnames(spe1)

coord1 <- get_coord_df(spe1)

a1 <- fread(cluster_cleany_path)
a2 <- fread(cluster_normalized_path)

df1 <- left_join(coord1, a1, by = c("sample_id", "barcode"))
df2 <- left_join(coord1, a2, by = c("sample_id", "barcode"))

cluster_cols_1 <- setdiff(names(a1), c("sample_id", "barcode"))
cluster_cols_2 <- setdiff(names(a2), c("sample_id", "barcode"))

common_cluster_cols <- intersect(cluster_cols_1, cluster_cols_2)

out_dir <- file.path(plot_dir, "jaccard_cross_spe")
dir.create(out_dir, showWarnings = FALSE)

# Pairwise combination: clustering columns of df1 vs clustering columns of df2
for (col in common_cluster_cols) {
  
  # Align by coordinates + sample_id
  merged_df <- inner_join(
    df1 %>% select(x, y, sample_id, !!col := all_of(col)),
    df2 %>% select(x, y, sample_id, !!col := all_of(col)),
    by = c("x", "y", "sample_id")
  ) %>% drop_na()

  # calculate Jaccard Index
  jaccard_list <- list()
  for (val1 in unique(merged_df[[paste0(col, ".x")]])) {
  for (val2 in unique(merged_df[[paste0(col, ".y")]])) {
    intersect_n <- sum(merged_df[[paste0(col, ".x")]] == val1 & merged_df[[paste0(col, ".y")]] == val2)
    union_n <- sum(merged_df[[paste0(col, ".x")]] == val1 | merged_df[[paste0(col, ".y")]] == val2)
      jaccard_list[[length(jaccard_list) + 1]] <- tibble(
        cluster_1 = val1,
        cluster_2 = val2,
        jaccard_index = intersect_n / union_n
      )
    }
  }

  jaccard_df <- bind_rows(jaccard_list)

  # plot
  p <- ggplot(jaccard_df, aes(factor(cluster_1), factor(cluster_2), fill = jaccard_index)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_bw(base_size = 12) +
    labs(
      x = paste0("Cluster from cleany: ", col),
      y = paste0("Cluster from lib-size normalized: ", col),
      fill = "Jaccard Index"
    )

  # save PDF
  ggsave(
    filename = file.path(out_dir, paste0("jaccard_", col, "_common.pdf")),
    plot = p,
    width = 8, height = 6
  )
}

# prep_clustering_results = function(spe_dir, cluster_path, cluster_colname) {
#     spe = readRDS(spe_dir)  
#     cluster_df = tibble(
#             x = round(spatialCoords(spe)[, 'pxl_col_in_fullres'], 1),
#             y = round(spatialCoords(spe)[, 'pxl_row_in_fullres'], 1),
#             sample_id = spe$sample_id
#         ) |>
#         left_join(fread(cluster_path), by = c("sample_id","barcode")) |>
#     return(cluster_df)
# }

# cluster_cleany_df = prep_clustering_results(
#     spe_cleany_dir, cluster_cleany_path, 'cluster_cleany'
# )

# prep_clustering_results = function(spe_dir, cluster_path, cluster_colname) {
#     spe = loadHDF5SummarizedExperiment(spe_dir)
    
#     cluster_df = tibble(
#             x = round(spatialCoords(spe)[, 'pxl_col_in_fullres'], 1),
#             y = round(spatialCoords(spe)[, 'pxl_row_in_fullres'], 1),
#             sample_id = spe$sample_id,
#             key = spe$key,
#             segmentation_type = spe$labels_joint_source
#         ) |>
#         filter(segmentation_type == 'primary') |>
#         left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
#         dplyr::rename(!!cluster_colname := banksy_lambda0_8) |>
#         select(x, y, sample_id, {{ cluster_colname }})

#     stopifnot(!any(is.na(cluster_df[[cluster_colname]])))
    
#     return(cluster_df)
# }

# cluster_normalized_df = prep_clustering_results(
#     spe_normalized_dir, cluster_normalized_path, 'cluster_normalized'
# )

# cluster_df = inner_join(
#     cluster_cleany_df, cluster_normalized_df, by = c('x', 'y', 'sample_id')
# )

# #   Compute the Jaccard index for each combination of "cleany" and "normalized" clusters
# jaccard_df_list = list()
# for (cleany_val in unique(cluster_df$cluster_cleany)) {
#     for (normalized_val in unique(cluster_df$cluster_normalized)) {
#         intersect_size = cluster_df |>
#             filter((cluster_cleany == cleany_val) & (cluster_normalized == normalized_val)) |>
#             nrow()
#         union_size = cluster_df |>
#             filter((cluster_cleany == cleany_val) | (cluster_normalized == normalized_val)) |>
#             nrow()
#         jaccard_df_list[[length(jaccard_df_list) + 1]] = tibble(
#             cluster_cleany = cleany_val,
#             cluster_normalized = normalized_val,
#             jaccard_index = intersect_size / union_size
#         )
#     }
# }

# #   Plot a heatmap of Jaccard indices
# p = do.call(rbind, jaccard_df_list) |>
#     ggplot(
#         aes(
#             x = cluster_cleany, y = cluster_normalized,
#             fill = jaccard_index
#         )
#     ) +
#     geom_tile() +
#     scale_fill_viridis_c() +
#     coord_cartesian(expand = FALSE) +
#     theme_bw(base_size = 25) +
#     labs(
#         x = 'Cluster: cleany clusters', y = 'Cluster: normalized clusters',
#         fill = 'Jaccard\nIndex'
#     )
# pdf(file.path(plot_dir, 'jaccard_index.pdf'), width = 9)
# print(p)
# dev.off()

session_info()
