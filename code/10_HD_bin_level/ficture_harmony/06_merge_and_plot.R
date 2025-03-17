#   Gather cell-level FICTURE results across samples and export. Plot the top
#   FICTURE cluster for each sample, then create composition plots to determine
#   the extent to which batch effects affect cluster assignments by FICTURE

library(here)
library(tidyverse)
library(spatialLIBD)
library(HDF5Array)
library(paletteer)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
ficture_paths = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'bin2cell_out',
    sprintf('%s.csv', readLines(sample_id_path))
)
ficture_path_out = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'ficture_merged.csv'
)
plot_dir = here('plots', '10_HD_bin_level', 'ficture_harmony')

dir.create(plot_dir, showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)

#   Read in and concatenate cell-level FICTURE results for all samples
ficture_df_list = list()
for (ficture_path in ficture_paths) {
    ficture_df_list[[ficture_path]] = read_csv(
        ficture_path, show_col_types = FALSE
    )
}
ficture_df = do.call(rbind, ficture_df_list)

#   Compute the top FICTURE cluster
ficture_df$FICTURE_cluster = ficture_df |>
    select(matches('^FICTURE_[0-9]+')) |>
    as.matrix() |>
    apply(1, function(x) ifelse(all(x == 0), NA, which.max(x) - 1))

write_csv(ficture_df, ficture_path_out)

#   Add top FICTURE cluster to colData
stopifnot(all(spe$key %in% ficture_df$key))
col_data = colData(spe) |>
    as_tibble() |>
    left_join(ficture_df |> select(key, FICTURE_cluster), by = "key") |>
    DataFrame()
rownames(col_data) = colnames(spe)
colData(spe) = col_data

#   Plot top FICTURE cluster for each sample
for (sample_id in unique(spe$sample_id)) {
    p = vis_clus(
            spe, sampleid = sample_id, clustervar = "FICTURE_cluster",
            is_stitched = TRUE, point_size = 20
        ) +
        guides(fill = guide_legend(override.aes = list(size = 8)))
    png(
        file.path(plot_dir, sprintf('%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

ficture_df = colData(spe) |>
    as_tibble() |>
    select(FICTURE_cluster, sample_id) |>
    filter(!is.na(FICTURE_cluster))

#   Distribution of samples by cluster
p = ggplot(ficture_df, aes(x = FICTURE_cluster, fill = sample_id)) +
    geom_bar(position = "fill") +
    theme_bw(base_size = 20) +
    labs(x = "FICTURE Cluster", y = "Proportion of Cells", fill = "Sample ID")
pdf(file.path(plot_dir, 'batch_effect_cluster.pdf'), width = 10)
print(p)
dev.off()

#   Distribution of clusters by sample
p = ggplot(ficture_df, aes(x = sample_id, fill = factor(FICTURE_cluster))) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = paletteer_d("Polychrome::palette36", 12)) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Sample ID", y = "Proportion of Cells", fill = "FICTURE\nCluster")
pdf(file.path(plot_dir, 'batch_effect_sample.pdf'), height = 10)
print(p)
dev.off()

session_info()
