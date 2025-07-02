#   Do clusters become more sample-specific as k increases?

library(here)
library(data.table)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

ficture_cluster_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'bin_level_clusters_normalized.csv.gz'
)
ficture_cleany_cluster_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
banksy_cluster_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', '%s',
    'leiden_res%s.csv'
)
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'ficture_harmony')

all_banksy_res = c(seq_len(20) / 10, 4, 8)
all_banksy_lambda = c(0.2, 0.8)

sample_cutoff = 0.5

################################################################################
#   Read in and clean clustering results
################################################################################

ficture_df = fread(ficture_cluster_path) |>
    as_tibble() |>
    mutate(sample_id = factor(sample_id)) |>
    pivot_longer(
        cols = matches('^FICTURE_k'),
        names_to = 'k', values_to = 'cluster'
    ) |>
    filter(!is.na(cluster))

#   Compute scalars for each sample, such that all add to 1 and each
#   represents the relative number of bins comprising each sample
size_df = ficture_df |>
    group_by(sample_id, k) |>
    summarize(sample_size_scalar = n()) |>
    group_by(k) |>
    mutate(
        sample_size_scalar = length(unique(sample_id)) * sample_size_scalar /
            sum(sample_size_scalar)
    ) |>
    ungroup()

ficture_df = ficture_df |>
    left_join(size_df, by = c('sample_id', 'k')) |>
    group_by(sample_id, k, cluster) |>
    summarize(num_bins = n() / sample_size_scalar[1]) |>
    group_by(k, cluster) |>
    summarize(max_prop = max(num_bins) / sum(num_bins)) |>
    group_by(k) |>
    summarize(num_balanced = sum(max_prop <= sample_cutoff)) |>
    ungroup() |>
    mutate(
        k = as.integer(sub('^FICTURE_k', '', k)),
        method = 'ficture_normalized'
    )

ficture_cleany_df = fread(ficture_cleany_cluster_path) |>
    as_tibble() |>
    mutate(sample_id = factor(sample_id)) |>
    pivot_longer(
        cols = matches('^FICTURE_k'),
        names_to = 'k', values_to = 'cluster'
    ) |>
    filter(!is.na(cluster)) |>
    group_by(sample_id, k, cluster) |>
    summarize(num_bins = n()) |>
    group_by(k, cluster) |>
    summarize(max_prop = max(num_bins) / sum(num_bins)) |>
    group_by(k) |>
    summarize(num_balanced = sum(max_prop <= sample_cutoff)) |>
    ungroup() |>
    mutate(
        k = as.integer(sub('^FICTURE_k', '', k)),
        method = 'ficture_cleany'
    )

banksy_df_list = list()
for (lambda in all_banksy_lambda) {
    lambda_neat = paste0('lambda', sub('\\.', '_', as.character(lambda)))
    for (res in all_banksy_res) {
        res_neat = sub('\\.', '_', as.character(res))
        banksy_df_list[[length(banksy_df_list) + 1]] = sprintf(
                banksy_cluster_paths, lambda_neat, res_neat
            ) |>
            read_csv(show_col_types = FALSE) |>
            dplyr::rename(cluster = paste0('banksy_', lambda_neat)) |>
            mutate(
                sample_id = factor(sub('^[0-9]+_', '', key)),
                k = length(unique(cluster)),
                method = paste0('banksy_', lambda_neat)
            )
    }
}

banksy_df = do.call(rbind, banksy_df_list) |>
    group_by(method, sample_id, k, cluster) |>
    summarize(num_bins = n()) |>
    group_by(method, k, cluster) |>
    summarize(max_prop = max(num_bins) / sum(num_bins)) |>
    group_by(method, k) |>
    summarize(num_balanced = sum(max_prop <= sample_cutoff)) |>
    ungroup() |>
    select(k, num_balanced, method)

p = ggplot(
        rbind(ficture_df, ficture_cleany_df, banksy_df),
        aes(x = k, y = num_balanced, color = method)
    ) +
    geom_line() +
    labs(x = 'k', y = 'Number of Balanced Clusters', color = 'Method') +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'sample_specificity_Jun27.pdf'))
print(p)
dev.off()

session_info()
