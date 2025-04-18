library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

ficture_cor_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'registration', 'cor_vs_snRNAseq_fine.rds'
)
ficture_cluster_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony',
    'bin_level_clusters.csv.gz'
)
marker_cor_val = 0.3
non_marker_cor_val = 0.15

process_cor_df = function(cor_df) {
    #   Tidy up and convert to long format
    cor_df = cor_df |>
        as.data.frame() |>
        rownames_to_column('cell_type') |>
        as_tibble() |>
        mutate(is_habenula = grepl('^[ML]Hb', cell_type)) |>
        pivot_longer(
            cols = -c(cell_type, is_habenula),
            names_to = 'cluster', values_to = 'cor_val'
        )
    
    #   Number of clusters registering only to habenula cell types
    num_pure_hb_clusters = cor_df |>
        group_by(cluster) |>
        summarize(
            only_hb = all(
                ifelse(cor_val > marker_cor_val, is_habenula, TRUE) &
                ifelse(!is_habenula, cor_val < non_marker_cor_val, TRUE)
            )
        ) |>
        filter(only_hb) |>
        nrow()
    
    #   Number of cell types having at least one cluster uniquely registering
    #   to them
    num_non_hb_cell_types = cor_df |>
        group_by(cluster) |>
        arrange(desc(cor_val)) |>
        filter(
            (cor_val[1] > marker_cor_val) &
            !is_habenula[1] &
            (cor_val[2] < non_marker_cor_val)
        ) |>
        filter(cor_val > marker_cor_val) |>
        pull(cell_type) |>
        unique() |>
        length()
    
    cor_df |>
        group_by(cluster) |>
        filter(
            all(
                ifelse(cor_val > marker_cor_val, is_habenula, TRUE) &
                ifelse(!is_habenula, cor_val < non_marker_cor_val, TRUE)
            )
        ) |>
}

ficture_cor = readRDS(ficture_cor_path)
