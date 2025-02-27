#   Compare FICTURE output plots to the cell-level FICTURE clusters determined
#   through the bin2cell-based method to confirm that method works

library(here)
library(tidyverse)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
spe_dir = here('processed-data', '09_HD_cell_level', 'spe_norm_filtered')
ficture_paths = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate', 'bin2cell_out', 
    sprintf('%s.csv', readLines(sample_id_path))
)
plot_dir = here('plots', '09_HD_cell_level', 'ficture_aggregate')

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
    apply(1, function(x) which.max(x) - 1)

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

session_info()
