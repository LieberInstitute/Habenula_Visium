#   Use sctype to try to resolve cluster 4, an ambiguous but large cluster we'd
#   rather not discard. Also check cluster 27

library(tidyverse)
library(here)
library(SpatialExperiment)
library(HGNChelper)
library(sessioninfo)

#   Source sc-type functions
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'spe_norm_filtered.rds'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
score_out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'registration_banksy', 'sctype_ambig_scores.csv'
)
ambig_clusters = c('4', '27')
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

spe = readRDS(spe_path)

#   Add in cluster assignments to 'spe'. Subset to only the ambiguous clusters
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(setequal(spe$key, cluster_df$key))
spe$banksy = factor(
    as.character(cluster_df$banksy_lambda0_2[match(spe$key, cluster_df$key)]),
    levels = as.character(sort(unique(cluster_df$banksy_lambda0_2)))
)
spe = spe[, spe$banksy %in% ambig_clusters]

gs_list = gene_sets_prepare(db_, "Brain")

#   Just make sure we have the necessary genes in our data
message("Proportion of marker genes present in our data:")
tibble(
        cell_type = rep(
            names(gs_list$gs_positive),
            times = sapply(gs_list$gs_positive, length)
        ),
        gene_name = unlist(gs_list$gs_positive)
    ) |>
    group_by(cell_type) |>
    summarize(representation = mean(gene_name %in% rowData(spe)$gene_name)) |>
    print(n = length(gs_list$gs_positive))

#   For memory, subset to marker genes (as gene symbols)
spe = spe[rowData(spe)$gene_name %in% unlist(gs_list$gs_positive), ]
rownames(spe) = rowData(spe)$gene_name

#   Score cells with brain cell types
score_mat = sctype_score(
    scRNAseqData = as.matrix(logcounts(spe)), scaled = TRUE,
    gs = gs_list$gs_positive
)

#   Compile scores for each cluster-and-cell-type combo as a tibble
score_df_list = list()
for (cluster in ambig_clusters) {
    cluster_index = spe$banksy == cluster
    sorted_scores = sort(rowMeans(score_mat[,cluster_index]), decreasing = TRUE)
    score_df_list[[cluster]] = tibble(
        cluster = cluster,
        cell_type = names(sorted_scores),
        score = sorted_scores
    )
}
score_df = bind_rows(score_df_list)

write_csv(score_df, score_out_path)

score_df |>
    group_by(cluster) |>
    slice_max(order_by = score, n = 1) |>
    ungroup()

#   Now check if scores are >0.25

session_info()
