library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)
library(cowplot)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
banksy_paths = here(
    'processed-data', '09_HD_cell_level', 'banksy', '%s', 'leiden_%s.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'subclustering')
markers = c("POU4F1", "GPR151", "CHRNB4", "HTR2C")
all_lambda = c(0.2, 0.8)
all_res = seq_len(10) / 10

dir.create(plot_dir, showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Compute mean of Z-scored expression for habenula markers
markers = rownames(spe)[match(markers, rowData(spe)$gene_name)]
marker_mat = assays(spe)$logcounts[markers,]
spe$habenula_markers = colMeans(
    (marker_mat - rowMeans(marker_mat)) / rowSds(marker_mat)
)

for (lambda in all_lambda) {
    lambda_neat = paste0('lambda', sub('\\.', '_', as.character(lambda)))
    plot_list = list()
    for (res in all_res) {
        res_neat = paste0('res', sub('\\.', '_', as.character(res)))

        #   Read in Banksy clusters
        cluster_df = sprintf(banksy_paths, lambda_neat, res_neat) |>
            read_csv(show_col_types = FALSE) |>
            rename(banksy_cluster = paste0('banksy_', lambda_neat)) |>
            left_join(
                colData(spe) |>
                    as_tibble() |>
                    select(key, habenula_markers),
                by = "key"
            )

        cor_df_list = list()
        for (cluster_val in unique(cluster_df$banksy_cluster)) {
            cor_df_list[[length(cor_df_list) + 1]] = tibble(
                cluster_val = cluster_val,
                cor_val = cor(
                    cluster_df$banksy_cluster == cluster_val,
                    cluster_df$habenula_markers
                )
            )
        }

        #   Plot correlation of habenula markers with presence of each Banksy
        #   cluster
        plot_list[[length(plot_list) + 1]] = do.call(rbind, cor_df_list) |>
            ggplot(aes(x = cluster_val, y = res, fill = cor_val)) +
            geom_tile() +
            scale_fill_viridis_c() +
            theme_bw(base_size = 20) +
            labs(
                x = 'Banksy Cluster', y = '',
                fill = '',
                title = sprintf('Banksy: lambda = %s, res = %s', lambda, res)
            )
    }

    pdf(
        file.path(
            plot_dir, sprintf('banksy_habenula_lambda_%s.pdf', lambda_neat)
        )
    )
    print(plot_grid(plotlist = plot_list, ncol = 1))
    dev.off()
}

session_info()
