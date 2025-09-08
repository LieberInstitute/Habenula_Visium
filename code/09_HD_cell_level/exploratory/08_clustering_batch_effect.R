#   Check how sample-specific Banksy clusters are as a function of k. This
#   script is no longer used; instead
#   code/10_HD_bin_level/ficture_harmony/12_sample_specificity.R
#   is used.

library(tidyverse)
library(here)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
banksy_paths = here(
    'processed-data', '09_HD_cell_level', 'banksy', 'lambda0_8',
    'leiden_res%s.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'clustering_batch_effect'
)
banksy_res = seq_len(10) / 10

dir.create(plot_dir, showWarnings = FALSE)

sample_ids = readLines(sample_id_path)

#   Read in all Banksy clustering results for each Leiden resolution
banksy_df_list = list()
for (this_res in banksy_res) {
    banksy_df_list[[length(banksy_df_list) + 1]] = sprintf(
            banksy_paths, sub('\\.', '_', as.character(this_res))
        ) |>
        read_csv(show_col_types = FALSE) |>
        rename(cluster = banksy_lambda0_8) |>
        mutate(
            method = "Banksy",
            k = length(unique(cluster)),
            sample_id = sub('^[0-9]+_', '', key)
        ) |>
        select(sample_id, method, k, cluster)
}

banksy_df = do.call(rbind, banksy_df_list) |>
    #   For each unique cluster value, compute the proportion of cells belonging
    #   to each sample
    group_by(sample_id, method, k, cluster) |>
    summarize(n = n()) |>
    group_by(method, k, cluster) |>
    mutate(prop = n / sum(n)) |>
    #   Take the maximum such proportion across samples
    summarize(max_prop = max(prop)) |>
    #   Now average to get a global metric of how sample-specific clusters tend
    #   to be for each method and k
    group_by(method, k) |>
    summarize(mean_max_prop = mean(max_prop))

p = ggplot(
        banksy_df, aes(x = k, y = mean_max_prop, color = method, group = method)
    ) +
    geom_line() +
    geom_point() +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 20) + 
    labs(
        x = 'Number of Clusters', y = 'Sample Polarization',
        color = 'Clustering\nMethod'
    )

pdf(file.path(plot_dir, 'mean_max_prop.pdf'), width = 10)
print(p)
dev.off()

session_info()
