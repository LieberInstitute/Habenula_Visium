#   Use sctype to try to resolve cluster 4, an ambiguous but large cluster we'd
#   rather not discard

library(tidyverse)
library(here)
library(SpatialExperiment)
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
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

spe = readRDS(spe_path)

#   Add in cluster assignments to 'spe'
cluster_df = read_csv(cluster_path, show_col_types = FALSE)
stopifnot(setequal(spe$key, cluster_df$key))
spe$banksy = factor(
    as.character(cluster_df$banksy_lambda0_2[match(spe$key, cluster_df$key)]),
    levels = as.character(sort(unique(cluster_df$banksy_lambda0_2)))
)

gs_list <- gene_sets_prepare(db_, "Brain")

es.max <- sctype_score(scRNAseqData = as.matrix(logcounts(sce)), 
                       scaled = TRUE, 
                       gs = gs_list$gs_positive)

## compile scores
cL_results <- purrr::map_dfr(clusters, function(cluster){
    cluster_index = sce[[cl_col]] == cluster
    
    es.max.cl = sort(rowSums(es.max[,cluster_index]), decreasing = !0)
    cL_resutls <- head(tibble(cluster = cluster, 
                              type = names(es.max.cl), 
                              scores = es.max.cl, 
                              ncells = sum(cluster_index)
    ),10)
    return(cL_resutls)
})

sctype_scores <-  cL_results |> 
    group_by(cluster) |> 
    top_n(n = 1, wt = scores)  |>
    mutate(confident = scores >= (ncells/4))