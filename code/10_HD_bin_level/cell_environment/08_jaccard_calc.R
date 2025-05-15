#   For a specific value of k (determined by the array task), compare FICTURE
#   results when using just the extracellular region vs the full tissue, as
#   spatially they look similar. There may not necessarily be a benefit to
#   subsetting to the extracellular region. Just compute the plot for this k;
#   in the next script, we'll plot a multi-page PDF with all results

library(here)
library(tidyverse)
library(data.table)
library(SpatialExperiment)
library(sessioninfo)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

ficture_extra_path = here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
        'ficture_outputs', 'normalized', 'k_%d', 'analysis', 'nF%d.d_12',
        'normalized_joined_input.tsv.gz'
    ) |>
    sprintf(k, k)
ficture_all_path = here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
        'ficture_outputs', 'normalized', 'k_%d', 'analysis', 'nF%d.d_12',
        'normalized_joined_input.tsv.gz'
    ) |>
    sprintf(k, k)
out_path = here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
        'temp_jaccard', 'k_%d.rds'
    ) |>
    sprintf(k)
ficture_colnames = c('sample_id', 'barcode', 'factor_K1')

dir.create(dirname(out_path), showWarnings = FALSE)

message(Sys.time(), ' | Reading in extracellular results...')
ficture_extra = fread(ficture_extra_path, select = ficture_colnames) |>
    as_tibble() |>
    group_by(sample_id, barcode) |>
    slice_head(n = 1) |>
    ungroup() |>
    filter(!is.na(factor_K1)) |>
    dplyr::rename(cluster_extra = factor_K1)

message(Sys.time(), ' | Reading in full-tissue results...')
ficture_all = fread(ficture_all_path, select = ficture_colnames) |>
    as_tibble() |>
    group_by(sample_id, barcode) |>
    slice_head(n = 1) |>
    ungroup() |>
    filter(!is.na(factor_K1)) |>
    dplyr::rename(cluster_all = factor_K1)

message(Sys.time(), ' | Joining...')
ficture = inner_join(ficture_extra, ficture_all, by = c('sample_id', 'barcode'))

#   Compute the Jaccard index for each combination of extracellular and
#   full-tissue clusters
jaccard_df_list = list()
for (extra_val in unique(ficture$cluster_extra)) {
    for (all_val in unique(ficture$cluster_all)) {
        intersect_size = ficture |>
            filter((cluster_extra == extra_val) & (cluster_all == all_val)) |>
            nrow()
        union_size = ficture |>
            filter((cluster_extra == extra_val) | (cluster_all == all_val)) |>
            nrow()
        jaccard_df_list[[length(jaccard_df_list) + 1]] = tibble(
            cluster_extra = extra_val,
            cluster_all = all_val,
            jaccard_index = intersect_size / union_size
        )
    }
}

#   Heatmap of Jaccard indices
p = do.call(rbind, jaccard_df_list) |>
    ggplot(
        aes(
            x = cluster_extra, y = cluster_all,
            fill = jaccard_index
        )
    ) +
    geom_tile() +
    scale_fill_viridis_c() +
    coord_cartesian(expand = FALSE) +
    theme_bw(base_size = 25) +
    labs(
        x = 'Extracellular Cluster', y = 'Full-Tissue Cluster',
        fill = 'Jaccard\nIndex', title = paste('k =', k)
    )
saveRDS(p, out_path)

session_info()
