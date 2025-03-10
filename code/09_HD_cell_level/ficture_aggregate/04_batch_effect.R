#   Use composition bar plots to assess the extent to which batch effects affect
#   cluster assignments by FICTURE at the cell level

library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
library(sessioninfo)

spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
ficture_path = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate',
    'ficture_merged.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'ficture_aggregate')

spe = loadHDF5SummarizedExperiment(spe_dir)
ficture_df = read_csv(ficture_path, show_col_types = FALSE) |>
    select(key, FICTURE_cluster) |>
    filter(!is.na(FICTURE_cluster)) |>
    inner_join(
        colData(spe) |>
            as_tibble() |>
            select(key, sample_id),
        by = "key"
    )

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
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "Sample ID", y = "Proportion of Cells", fill = "FICTURE\nCluster")
pdf(file.path(plot_dir, 'batch_effect_sample.pdf'), height = 10)
print(p)
dev.off()

session_info()
