library(here)
library(tidyverse)
library(sessioninfo)

banksy_paths = here(
    'processed-data', '09_HD_cell_level', 'banksy', '%s', 'leiden_%s.csv'
)
ficture_path = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate',
    'ficture_merged.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'ficture_aggregate')
all_lambda = c(0.2, 0.8)
all_res = seq_len(10) / 10

ficture_df = read_csv(ficture_path, show_col_types = FALSE) |>
    select(key, FICTURE_cluster) |>
    filter(!is.na(FICTURE_cluster))

for (lambda in all_lambda) {
    lambda_neat = paste0('lambda', sub('\\.', '_', as.character(lambda)))
    plot_list = list()
    for (res in all_res) {
        res_neat = paste0('res', sub('\\.', '_', as.character(res)))

        #   Form a tibble with key, FICTURE, and Banksy clusters as columns
        cluster_df = sprintf(banksy_paths, lambda_neat, res_neat) |>
            read_csv(show_col_types = FALSE) |>
            rename(banksy_cluster = paste0('banksy_', lambda_neat)) |>
            inner_join(ficture_df, by = "key")

        #   Compute the Jaccard index for each pair of FICTURE and Banksy clusters
        jaccard_df_list = list()
        for (ficture_val in unique(cluster_df$FICTURE_cluster)) {
            for (banksy_val in unique(cluster_df$banksy_cluster)) {
                intersect_size = cluster_df |>
                    filter((FICTURE_cluster == ficture_val) & (banksy_cluster == banksy_val)) |>
                    nrow()
                union_size = cluster_df |>
                    filter((FICTURE_cluster == ficture_val) | (banksy_cluster == banksy_val)) |>
                    nrow()
                jaccard_df_list[[length(jaccard_df_list) + 1]] = tibble(
                    FICTURE_cluster = ficture_val,
                    banksy_cluster = banksy_val,
                    jaccard_index = intersect_size / union_size
                )
            }
        }

        #   Plot a heatmap of Jaccard indices for each pair of FICTURE and
        #   Banksy clusters
        plot_list[[length(plot_list) + 1]] = do.call(rbind, jaccard_df_list) |>
            ggplot(
                aes(
                    x = FICTURE_cluster, y = banksy_cluster,
                    fill = jaccard_index
                )
            ) +
            geom_tile() +
            scale_fill_viridis_c() +
            coord_cartesian(expand = FALSE) +
            theme_bw(base_size = 25) +
            labs(
                x = 'FICTURE Cluster', y = 'Banksy Cluster',
                fill = 'Jaccard\nIndex',
                title = sprintf('Banksy: lambda = %s, res = %s', lambda, res)
            )
    }

    pdf(file.path(plot_dir, sprintf('ficture_banksy_lambda_%s.pdf', lambda_neat)))
    print(plot_list)
    dev.off()
}

session_info()
