library(tidyverse)
library(here)
library(sessioninfo)
library(crawdad)

liana_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana2', 'table',
    'significant_interactions_across_donors_%s.0.csv'
)
sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
crawdad_paths = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'crawdad',
    'output', '%s_results.csv'
)
plot_dir = here('plots', '10_HD_bin_level', 'no_secondary', 'liana2')
bandwidths = c(1000, 2500, 3000, 3500, 5000, 6000, 7000)
min_num_signif = 5

sample_info = read_csv(sample_info_path, show_col_types = FALSE)
sample_ids = sample_info$tissue_id

crawdad_list = list()
for (sample_id in sample_ids) {
    crawdad_list[[sample_id]] = sprintf(crawdad_paths, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        mutate(sample_id = sample_id)
}

#   Get the Z-score significance threshold (same in all samples)
z_sig = do.call(rbind, crawdad_list) |>
    filter(sample_id == sample_ids[1]) |>
    correctZBonferroni()

crawdad_df = do.call(rbind, crawdad_list) |>
    filter(reference != 'Excit.Thal', neighbor != 'Excit.Thal') |>
    #   First average Z-scores across permutations
    group_by(sample_id, neighbor, scale, reference) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= z_sig) |>
    group_by(sample_id, neighbor, reference) |>
    filter(scale == min(scale)) |>
    #   Take the mean Z-score and scale across samples
    group_by(neighbor, reference) |>
    summarize(scale = mean(scale), Z = mean(Z)) |>
    ungroup() |>
    dplyr::rename(source = reference, target = neighbor) |>
    select(source, target, Z)

liana_df_list = list()
for (this_bandwidth in bandwidths) {
    liana_df_list[[as.character(this_bandwidth)]] = sprintf(
            liana_path, this_bandwidth
        ) |>
        read_csv(show_col_types = FALSE) |>
        group_by(source, target) |>
        #   We're technically counting the same pairs multiple times since there
        #   are multiple donors, but this still works out favorably for
        #   interpretation
        summarize(num_int = n()) |>
        ungroup() |>
        mutate(bandwidth = this_bandwidth)
}
liana_df = do.call(rbind, liana_df_list) |>
    left_join(crawdad_df, by = c('source', 'target')) |>
    mutate(Z = ifelse(is.na(Z), 0, Z))

p = liana_df |>
    group_by(bandwidth) |>
    #   What proportion of interactions are unexpected: namely non-neurons
    #   communicating to neurons (among cases where those cell types are
    #   generally dispersed according to CRAWDAD)?
    summarize(
        prop_bad_int = sum(
            num_int[(Z < 0) & !grepl('Hb', source) & grepl('Hb', target)]
        ) / sum(num_int)
    ) |>
    ggplot(aes(x = bandwidth, y = prop_bad_int)) +
        geom_point() +
        geom_line() +
        theme_bw(base_size = 15) +
        labs(x = 'LIANA Bandwidth', y = 'Proportion of unexpected interactions')
pdf(file.path(plot_dir, 'crawdad_bandwidth_check.pdf'))
print(p)
dev.off()

session_info()
