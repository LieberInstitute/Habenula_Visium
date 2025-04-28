#   Since there are a huge number of spatial registration results to comb
#   through, this script intends to automate selecting interesting results.
#   In particular, one of our main goals is to find several distinct clusters
#   that partition the habenula. It's also nice to see clean matches against
#   non-habenula cell types.

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)

cor_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    '%s', 'cor_vs_snRNAseq_fine.rds'
)
cluster_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'banksy', '%s',
    'leiden_res%s.csv'
)

all_res = c(seq_len(20) / 10, 4, 8)
all_lambda = c(0.2, 0.8)

cutoff_merge_ratio = 0.1

################################################################################
#   Functions
################################################################################

process_cor_df = function(cor_df) {
    #   Annotate and tidy up
    cor_df = cor_df |>
        annotate_registered_clusters(cutoff_merge_ratio = cutoff_merge_ratio) |>
        as_tibble() |>
        filter(layer_confidence == 'good')
    
    #   Number of clusters registering only to habenula cell types
    num_hb_clusters = cor_df |>
        filter(
            sapply(
                layer_label,
                function(x) all(grepl('^[ML]Hb', str_split(x, '/')[[1]]))
            )
        ) |>
        nrow()
    
    #   Number of cell types having at least one cluster uniquely registering
    #   to them
    num_non_hb_cell_types = cor_df |>
        filter(!grepl('/', layer_label), !grepl('^[ML]Hb', layer_label)) |>
        pull(layer_label) |>
        unique() |>
        length()
    
    #   Number of habenula cell types having at least one cluster registering
    #   to them
    temp = strsplit(paste(cor_df$layer_label, collapse = '/'), '/')[[1]]
    num_hb_cell_types = length(unique(temp[grepl('^[ML]Hb', temp)]))

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

#   Collect metrics for all Banksy spatial registration results
summary_df_list = list()
for (lambda in all_lambda) {
    lambda_neat = sub('\\.', '_', as.character(lambda))
    cor_df = sprintf(cor_paths, sprintf('lambda%s', lambda_neat)) |>
        readRDS()

    for (i in seq_len(length(all_res))) {
        summary_df_list[[length(summary_df_list) + 1]] = process_cor_df(
                cor_df[[i]]
            ) |>
            mutate(
                method = 'banksy',
                res = all_res[i],
                lambda = lambda
            )
    }
}
summary_df = do.call(rbind, summary_df_list)

################################################################################
#   Explore top-ranking results
################################################################################

#   We're first prioritizing the ability of clustering to split the habenula.
#   Next, we consider how many habenula cell types are represented in habenula
#   clusters, and also how many non-habenula cell types
message('Top 10 clustering settings:')
summary_df |>
    arrange(
        desc(num_hb_clus),
        desc(num_hb_CT),
        desc(num_non_hb_CT)
    ) |>
    print(n = 10)

session_info()
