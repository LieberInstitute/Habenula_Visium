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
#   Functions
################################################################################

calculate_balanced = function(cluster_df, method_name) {
    #   Compute scalars for each sample, such that all add to 1 and each
    #   represents the relative number of bins comprising each sample
    size_df = cluster_df |>
        group_by(sample_id, k) |>
        summarize(sample_size_scalar = n()) |>
        group_by(k) |>
        mutate(
            sample_size_scalar = length(unique(sample_id)) *
                sample_size_scalar / sum(sample_size_scalar)
        ) |>
        ungroup()

    cluster_df = cluster_df |>
        left_join(size_df, by = c('sample_id', 'k')) |>
        group_by(sample_id, k, cluster) |>
        summarize(num_bins = n() / sample_size_scalar[1]) |>
        group_by(k, cluster) |>
        summarize(max_prop = max(num_bins) / sum(num_bins)) |>
        group_by(k) |>
        summarize(num_balanced = sum(max_prop <= sample_cutoff)) |>
        ungroup() |>
        mutate(method = method_name) |>
        select(k, num_balanced, method)
    
    return(cluster_df)
}

################################################################################
#   Read in and clean clustering results
################################################################################

message(Sys.time(), ' - Cleaning normalized FICTURE results')
ficture_norm_df = fread(ficture_cluster_path) |>
    as_tibble() |>
    mutate(sample_id = factor(sample_id)) |>
    pivot_longer(
        cols = matches('^FICTURE_k'),
        names_to = 'k', values_to = 'cluster'
    ) |>
    filter(!is.na(cluster)) |>
    mutate(k = as.integer(sub('^FICTURE_k', '', k))) |>
    calculate_balanced(method_name = 'FICTURE_normalized')

message(Sys.time(), ' - Cleaning cleaningY FICTURE results')
ficture_clean_df = fread(ficture_cleany_cluster_path) |>
    as_tibble() |>
    mutate(sample_id = factor(sample_id)) |>
    pivot_longer(
        cols = matches('^FICTURE_k'),
        names_to = 'k', values_to = 'cluster'
    ) |>
    filter(!is.na(cluster)) |>
    mutate(k = as.integer(sub('^FICTURE_k', '', k))) |>
    calculate_balanced(method_name = 'FICTURE_cleaningY')

message(Sys.time(), ' - Cleaning Banksy results')
banksy_df_outer_list = list()
for (lambda in all_banksy_lambda) {
    lambda_neat = paste0('lambda', sub('\\.', '_', as.character(lambda)))
    banksy_df_inner_list = list()
    for (res in all_banksy_res) {
        res_neat = sub('\\.', '_', as.character(res))
        banksy_df_inner_list[[length(banksy_df_inner_list) + 1]] = sprintf(
                banksy_cluster_paths, lambda_neat, res_neat
            ) |>
            read_csv(show_col_types = FALSE) |>
            dplyr::rename(cluster = paste0('banksy_', lambda_neat)) |>
            mutate(
                sample_id = factor(sub('^[0-9]+_', '', key)),
                k = length(unique(cluster))
            )            
    }
    banksy_df_outer_list[[length(banksy_df_outer_list) + 1]] = do.call(
            rbind, banksy_df_inner_list
        ) |>
        calculate_balanced(method_name = paste0('banksy_', lambda_neat))
}
banksy_df = do.call(rbind, banksy_df_outer_list)

message(Sys.time(), ' - Plotting')

#   All clustering methods (for internal use)
p = ggplot(
        rbind(ficture_norm_df, ficture_clean_df, banksy_df),
        aes(x = k, y = num_balanced, color = method)
    ) +
    geom_abline(slope = 1, linetype = 'dashed') +
    geom_line() +
    labs(x = 'k', y = 'Number of Balanced Clusters', color = 'Method') +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'sample_specificity.pdf'), height = 5)
print(p)
dev.off()

#   FICTURE only (for a supplementary figure)
p = rbind(ficture_norm_df, ficture_clean_df) |>
    mutate(
        method = ifelse(
            method == 'FICTURE_normalized', 'Normalized Only', 'Batch Corrected'
        )
    ) |>
    ggplot(aes(x = k, y = num_balanced, color = method)) +
        geom_abline(slope = 1, linetype = 'dashed') +
        geom_line() +
        labs(x = 'k', y = 'Number of Balanced Clusters', color = 'Method') +
        theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'sample_specificity_ficture.pdf'), height = 5)
print(p)
dev.off()

session_info()
