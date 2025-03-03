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
    cor_df_list = list()
    for (res in all_res) {
        res_neat = paste0('res', sub('\\.', '_', as.character(res)))

        #   Read in Banksy clusters
        cluster_df = sprintf(banksy_paths, lambda_neat, res_neat) |>
            read_csv(show_col_types = FALSE) |>
            dplyr::rename(banksy_cluster = paste0('banksy_', lambda_neat)) |>
            left_join(
                colData(spe) |>
                    as_tibble() |>
                    select(key, habenula_markers),
                by = "key"
            )

        for (cluster_val in unique(cluster_df$banksy_cluster)) {
            cor_df_list[[length(cor_df_list) + 1]] = tibble(
                cluster_val = cluster_val,
                cor_val = cor(
                    cluster_df$banksy_cluster == cluster_val,
                    cluster_df$habenula_markers
                ),
                res = res
            )
        }
    }

    #   Plot correlation of habenula markers with presence of each Banksy
    #   cluster
    p = do.call(rbind, cor_df_list) |>
        mutate(
            cluster_val = factor(
                cluster_val, levels = sort(unique(cluster_df$banksy_cluster))
            )
        ) |>
        ggplot(aes(x = cluster_val, y = 1, fill = cor_val)) +
        geom_tile() +
        facet_wrap(~paste('Resolution =', res), scales = 'free_x', ncol = 1) +
        scale_fill_viridis_c() +
        theme_bw(base_size = 10) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        labs(x = 'Banksy Cluster', y = '', fill = 'Cor.: Habenula Markers')
    pdf(
        file.path(plot_dir, sprintf('banksy_habenula_%s.pdf', lambda_neat)),
        width = 10
    )
    print(p)
    dev.off()

    message("Top correlations with habenula markers:")
    do.call(rbind, cor_df_list) |>
        arrange(desc(cor_val)) |>
        head() |>
        print()
}

session_info()
