#   Since there are a huge number of spatial registration results to comb
#   through, this script intends to automate selecting interesting results.
#   In particular, one of our main goals is to find several distinct clusters
#   that partition the habenula. It's also nice to see clean matches against
#   non-habenula cell types.
#
#   Another intention of this script is to see if clusters become more
#   sample-specific as k increases

library(here)
library(data.table)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

ficture_cor_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration', 'cor_vs_snRNAseq_fine.rds'
)
ficture_cluster_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'bin_level_clusters.csv.gz'
)
banksy_cor_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    '%s', 'cor_vs_snRNAseq_fine.rds'
)
banksy_cluster_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', '%s',
    'leiden_res%s.csv'
)
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'ficture_harmony')

all_banksy_res = seq_len(20) / 10
all_banksy_lambda = c(0.2, 0.8)

marker_cor_val = 0.3
non_marker_cor_val = 0.15

sample_cutoff = 0.5

################################################################################
#   Functions
################################################################################

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
    num_hb_clusters = cor_df |>
        group_by(cluster) |>
        summarize(
            only_hb = all(
                    ifelse(cor_val > marker_cor_val, is_habenula, TRUE) &
                    ifelse(!is_habenula, cor_val < non_marker_cor_val, TRUE)
                ) &
                any(cor_val > marker_cor_val)
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
            (cor_val[2] < non_marker_cor_val),
            cor_val > marker_cor_val
        ) |>
        pull(cell_type) |>
        unique() |>
        length()
    
    #   Number of habenula cell types having at least one cluster registering
    #   to them
    num_hb_cell_types = cor_df |>
        group_by(cluster) |>
        filter(
            all(
                ifelse(cor_val > marker_cor_val, is_habenula, TRUE) &
                ifelse(!is_habenula, cor_val < non_marker_cor_val, TRUE)
            ),
            cor_val > marker_cor_val
        ) |>
        pull(cell_type) |>
        unique() |>
        length()
    
    summary_df = tibble(
        num_hb_clus = num_hb_clusters,
        num_non_hb_CT = num_non_hb_cell_types,
        num_hb_CT = num_hb_cell_types,
        k = length(unique(cor_df$cluster)),
    )

    return(summary_df)
}

################################################################################
#   Gather metrics across clustering results
################################################################################

#   Collect metrics for all FICTURE spatial registration results
ficture_cor = readRDS(ficture_cor_path)
ficture_df = do.call(rbind, lapply(ficture_cor, process_cor_df)) |>
    mutate(method = 'ficture', res = NA, lambda = NA)

#   Collect metrics for all Banksy spatial registration results
banksy_df_list = list()
for (lambda in all_banksy_lambda) {
    lambda_neat = sub('\\.', '_', as.character(lambda))
    banksy_cor = sprintf(banksy_cor_paths, sprintf('lambda%s', lambda_neat)) |>
        readRDS()

    for (i in seq_len(length(all_banksy_res))) {
        banksy_df_list[[length(banksy_df_list) + 1]] = process_cor_df(
                banksy_cor[[i]]
            ) |>
            mutate(
                method = 'banksy',
                res = all_banksy_res[i],
                lambda = lambda
            )
    }
}
banksy_df = do.call(rbind, banksy_df_list)

################################################################################
#   Explore top-ranking results
################################################################################

#   We're first prioritizing the ability of clustering to split the habenula.
#   Next, we consider how many habenula cell types are represented in habenula
#   clusters, and also how many non-habenula cell types
message('Top 5 clustering settings (non-Hb cell types first):')
rbind(ficture_df, banksy_df) |>
    arrange(
        desc(num_hb_clus),
        desc(num_non_hb_CT),
        desc(num_hb_CT)
    ) |>
    print(n = 5)

message('Top 5 clustering settings (Hb cell types first):')
rbind(ficture_df, banksy_df) |>
    arrange(
        desc(num_hb_clus),
        desc(num_hb_CT),
        desc(num_non_hb_CT)
    ) |>
    print(n = 5)

message('Top 3 clustering settings by method (non-Hb cell types first):')
rbind(ficture_df, banksy_df) |>
    group_by(method) |>
    arrange(
        desc(num_hb_clus),
        desc(num_non_hb_CT),
        desc(num_hb_CT)
    ) |>
    slice_head(n = 3) |>
    ungroup() |>
    print(n = 6)

message('Top 3 clustering settings by method (Hb cell types first):')
rbind(ficture_df, banksy_df) |>
    group_by(method) |>
    arrange(
        desc(num_hb_clus),
        desc(num_hb_CT),
        desc(num_non_hb_CT)
    ) |>
    slice_head(n = 3) |>
    ungroup() |>
    print(n = 6)

################################################################################
#   Read in and clean clustering results
################################################################################

ficture_df = fread(ficture_cluster_path) |>
    as_tibble() |>
    mutate(sample_id = factor(sample_id)) |>
    pivot_longer(
        cols = matches('^FICTURE_k'),
        names_to = 'k', values_to = 'cluster'
    ) |>
    filter(!is.na(cluster)) |>
    group_by(sample_id, k, cluster) |>
    summarize(num_bins = n()) |>
    group_by(k, cluster) |>
    summarize(max_prop = max(num_bins) / sum(num_bins)) |>
    group_by(k) |>
    summarize(num_balanced = sum(max_prop <= sample_cutoff)) |>
    ungroup() |>
    mutate(
        k = as.integer(sub('^FICTURE_k', '', k)),
        method = 'ficture'
    )

banksy_df_list = list()
for (lambda in all_banksy_lambda) {
    lambda_neat = paste0('lambda', sub('\\.', '_', as.character(lambda)))
    for (res in all_banksy_res) {
        res_neat = sub('\\.', '_', as.character(res))
        banksy_df_list[[length(banksy_df_list) + 1]] = sprintf(
                banksy_cluster_paths, lambda_neat, res_neat
            ) |>
            read_csv(show_col_types = FALSE) |>
            dplyr::rename(cluster = paste0('banksy_', lambda_neat)) |>
            mutate(
                sample_id = factor(sub('^[0-9]+_', '', key)),
                k = length(unique(cluster)),
                method = paste0('banksy_', lambda_neat)
            )
    }
}

banksy_df = do.call(rbind, banksy_df_list) |>
    group_by(method, sample_id, k, cluster) |>
    summarize(num_bins = n()) |>
    group_by(method, k, cluster) |>
    summarize(max_prop = max(num_bins) / sum(num_bins)) |>
    group_by(method, k) |>
    summarize(num_balanced = sum(max_prop <= sample_cutoff)) |>
    ungroup() |>
    select(k, num_balanced, method)

p = ggplot(
        rbind(ficture_df, banksy_df),
        aes(x = k, y = num_balanced, color = method)
    ) +
    geom_line() +
    labs(x = 'k', y = 'Number of Balanced Clusters', color = 'Method') +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'sample_specificity.pdf'))
print(p)
dev.off()

session_info()
