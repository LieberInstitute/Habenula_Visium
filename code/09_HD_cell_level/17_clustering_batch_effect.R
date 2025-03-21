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
    group_by(sample_id, method, k, cluster) |>
    summarize(n = n()) |>
    group_by(method, k, cluster) |>
    mutate(prop = n / sum(n)) |>
    summarize(max_prop = max(prop)) |>
    group_by(method, k) |>
    summarize(mean_max_prop = mean(max_prop))

p = ggplot(
        banksy_df, aes(x = k, y = mean_max_prop, color = method, group = method)
    ) +
    geom_line() +
    theme_bw(base_size = 20)

pdf(file.path(plot_dir, 'mean_max_prop.pdf'))
print(p)
dev.off()

session_info()
