# # library(slurmjobs)
# # slurmjobs::job_single('01_compute_cor', create_shell = TRUE, memory = '20G', command = "01_compute_cor.R")
# 
# # To submit the job use: sbatch 01_compute_cor.sh

library("here")
library("purrr")
library("spatialLIBD")
library("sessioninfo")


## Input dir
dir_input <- here("processed-data", "05_layer_differential_expression", "modeling_results_BS")

## Create output directories
dir_rdata <- here("processed-data", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
## Plot dir 
dir_plot <- here("plots", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

## specify the number of BayesSpace k to use 
k <- seq(2,28)

## Load Registration Results 
bayesSpace_registration_fn <-
  map(k, ~ here(
    dir_input,
    paste0(
      "modeling_results_BayesSpace_k",
      sprintf("%02d", .x),
      ".Rdata"
    )
  ))
## Load the 3 model results (anova, enrichment, pairwise) for each domain in the BS
bayesSpace_registration <-
  lapply(bayesSpace_registration_fn, function(x) {
    get(load(x))
  })
stopifnot(is.list(bayesSpace_registration))
# names(bayesSpace_registration[[1]])
# [1] "anova"      "enrichment" "pairwise" 

## Select t-stats from the registration enrichment data (spatial domains)
registration_t_stats <-
  map(bayesSpace_registration, function(data) {
    x <- data$enrichment
    t_stats <- x[, grep("^t_stat_", colnames(x))]
    colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))
    return(t_stats)
  })
stopifnot(is.list(registration_t_stats))
# head(registration_t_stats[[1]])
# Sp02D01    Sp02D02
# ENSG00000237491 -1.6855500  1.6855500
# ENSG00000228794 -1.6513162  1.6513162
# ENSG00000223764 -0.8991561  0.8991561

registration_vars <-
    c("final_Annotations", "final_Annotations_broad")

message(' Processing Spatial Registration for ', length(k), ' BayesSpace k')

compute_cor <- function(current_var) {
    # Load input snRNA-seq data
    # testing: current_var = "final_Annotations"
    results_enrichment <-
        readRDS(here(
            "processed-data",
            "05_snRNA-seq_model_stats",
            paste0("enrichment_", current_var, ".rds")
        ))
    colnames(results_enrichment)
    modeling_res_enrichment <- list("enrichment" = results_enrichment)
    ## check out table
    # results_enrichment[1:5, 1:5]
    #                 t_stat_Sp02D01 t_stat_Sp02D02 p_value_Sp02D01 p_value_Sp02D02 fdr_Sp02D01
    # ENSG00000237491     -1.6855500      1.6855500       0.1263326       0.1263326   0.5432782
    # ENSG00000228794     -1.6513162      1.6513162       0.1332383       0.1332383   0.5432782
    # ENSG00000223764     -0.8991561      0.8991561       0.3921055       0.3921055   0.6688699
    # ENSG00000187634      0.5440598     -0.5440598       0.5996823       0.5996823   0.8063719
    # ENSG00000188976     -1.2522290      1.2522290       0.2422019       0.2422019   0.5724077

    lapply(
      registration_t_stats,
      layer_stat_cor, #results_enrichment
      modeling_results = modeling_res_enrichment,
      top_n = 100
    )
}

cor_fine <- compute_cor("final_Annotations")
# head(cor_fine[[1]])
#         Astrocyte       Endo Excit.Thal Inhib.Thal       LHb.1       LHb.2       LHb.3       LHb.4       LHb.5
# Sp02D01 -0.2038447  0.1249809 -0.1806295  -0.228601 -0.04530781  0.03818429  0.06463135  0.04138618  0.03397194
# Sp02D02  0.2038447 -0.1249809  0.1806295   0.228601  0.04530781 -0.03818429 -0.06463135 -0.04138618 -0.03397194

cor_broad <- compute_cor("final_Annotations_broad")
# head(cor_broad[[1]])
#         Astrocyte      Endo Excit.Thal Inhib.Thal        LHb          MHb    Microglia      Oligo         OPC
# Sp02D01 -0.2094872  0.123983 -0.1993256 -0.3184479  0.1457124 -0.009116221  0.003706224  0.5103645 -0.05427858
# Sp02D02  0.2094872 -0.123983  0.1993256  0.3184479 -0.1457124  0.009116221 -0.003706224 -0.5103645  0.05427858

## Annotate clusters / classify by layer confidence classes (good/poor)
# annotated_clusters_fine <-
#     lapply(cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
# annotated_clusters_fine[[5]]
# cluster layer_confidence     layer_label
# 1 Sp06D04             good     MHb.2/MHb.1
# 2 Sp06D02             good            Endo
# 3 Sp06D06             good           MHb.1
# 4 Sp06D01             good           Oligo
# 5 Sp06D03             poor Astrocyte/Endo*
# 6 Sp06D05             good       Astrocyte

# ## relaxed merging threshold 0.1 
# annotated_clusters_broad <-
#     lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
# # head(annotated_clusters_broad)

annotated_clusters_broad <-
  lapply(cor_broad, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)

## With default confidence and cutoff_merge_ratio 
annotated_clusters_fine <-
  lapply(cor_fine, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)

# ## Use annotation labels on the correlation matrices
# cor_fine <- mapply(function(cor, label_data) {
#     rownames(cor) <-
#         paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
#     return(cor)
# }, cor_fine, annotated_clusters_fine)
# 
# cor_broad <- mapply(function(cor, label_data) {
#     rownames(cor) <-
#         paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
#     return(cor)
# }, cor_broad, annotated_clusters_broad)

stopifnot(is.list(cor_fine))
head(cor_fine[[1]])
stopifnot(is.list(cor_broad))
head(cor_broad[[1]])

# ## Confidence marks "x" need to be re-loaded ?
# annotated_clusters_broad <-
#   lapply(cor_broad, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)
# 
# ## With default confidence and cutoff_merge_ratio 
# annotated_clusters_fine <-
#   lapply(cor_fine, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)

# data.frame(
#   "broad" = sort(rownames(cor_broad[[8]])),
#   "fine" = sort(rownames(cor_fine[[8]]))
# )

save(cor_fine,
    cor_broad,
    #file = file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq.Rdata")
    file = file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq_top100.Rdata")
)

##   Make heatmaps broad res

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_broadRes.pdf"))
for (i in seq_len(length(cor_broad))) {
  print(
    layer_stat_cor_plot(
      cor_broad[[i]], annotation = annotated_clusters_broad[[i]],
      heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
    )
  )
}
dev.off()

##   Make heatmaps fine res

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_fineRes.pdf"))
for (i in seq_len(length(cor_fine))) {
  print(
    layer_stat_cor_plot(
      cor_fine[[i]], annotation = annotated_clusters_fine[[i]],
      heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
    )
  )
}
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

## =============================================================================
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-07-02 14:09:24 EDT"
# > proc.time()
# user   system  elapsed 
# 164.899    4.289 8690.528 
# > options(width = 120)
# > session_info()
# (R 4.3.2)
# benchmarkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel             1.36.0      2023-10-24 [2] Bioconductor
# BiocSingular             1.18.0      2023-10-24 [2] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
# Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
# cachem                   1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
# cellranger               1.1.0       2016-07-27 [2] CRAN (R 4.3.2)
# cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# codetools                0.2-19      2023-02-01 [3] CRAN (R 4.3.2)
# colorspace               2.1-0       2023-01-23 [2] CRAN (R 4.3.2)
# config                   0.3.2       2023-08-30 [2] CRAN (R 4.3.2)
# cowplot                  1.1.3       2024-01-22 [2] CRAN (R 4.3.2)
# crayon                   1.5.2       2022-09-29 [2] CRAN (R 4.3.2)
# curl                     5.2.0       2023-12-08 [2] CRAN (R 4.3.2)
# data.table               1.15.0      2024-01-30 [2] CRAN (R 4.3.2)
# DBI                      1.2.1       2024-01-12 [2] CRAN (R 4.3.2)
# dbplyr                   2.4.0       2023-10-26 [2] CRAN (R 4.3.2)
# DelayedArray           * 0.28.0      2023-10-24 [2] Bioconductor
# DelayedMatrixStats       1.24.0      2023-10-24 [2] Bioconductor
# digest                   0.6.34      2024-01-11 [2] CRAN (R 4.3.2)
# doParallel               1.0.17      2022-02-07 [2] CRAN (R 4.3.2)
# dotCall64                1.1-1       2023-11-28 [2] CRAN (R 4.3.2)
# dplyr                  * 1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
# DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
# edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
# fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
# farver                   2.1.1       2022-07-06 [2] CRAN (R 4.3.2)
# fastmap                  1.1.1       2023-02-24 [2] CRAN (R 4.3.2)
# fields                   15.2        2023-08-17 [2] CRAN (R 4.3.2)
# filelock                 1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
# forcats                * 1.0.0       2023-01-29 [2] CRAN (R 4.3.2)
# foreach                  1.5.2       2022-02-02 [2] CRAN (R 4.3.2)
# generics                 0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb           * 1.38.5      2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData         1.2.11      2024-02-09 [2] Bioconductor
# GenomicAlignments        1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
# GenomicRanges          * 1.54.1      2023-10-29 [2] Bioconductor
# ggbeeswarm               0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2                * 3.4.4       2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel                  0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
# glue                     1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
# golem                    0.4.1       2023-06-05 [2] CRAN (R 4.3.2)
# gridExtra                2.3         2017-09-09 [2] CRAN (R 4.3.2)
# gtable                   0.3.4       2023-08-21 [2] CRAN (R 4.3.2)
# HDF5Array              * 1.30.0      2023-10-24 [2] Bioconductor
# here                   * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
# hms                      1.1.3       2023-03-21 [2] CRAN (R 4.3.2)
# htmltools                0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets              1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                   1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
# httr                     1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
# interactiveDisplayBase   1.40.0      2023-10-24 [2] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [2] Bioconductor
# irlba                    2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
# iterators                1.0.14      2022-02-05 [2] CRAN (R 4.3.2)
# jquerylib                0.1.4       2021-04-26 [2] CRAN (R 4.3.2)
# jsonlite                 1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
# KEGGREST                 1.42.0      2023-10-24 [2] Bioconductor
# labeling                 0.4.3       2023-08-29 [2] CRAN (R 4.3.2)
# later                    1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
# lattice                  0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval                 0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
# lifecycle                1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
# limma                    3.58.1      2023-10-31 [2] Bioconductor
# lobstr                 * 1.1.2       2022-06-22 [2] CRAN (R 4.3.2)
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# lubridate              * 1.9.3       2023-09-27 [2] CRAN (R 4.3.2)
# magick                   2.8.2       2023-12-20 [2] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
# Matrix                 * 1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
# matrixStats            * 1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# mime                     0.12        2021-09-28 [2] CRAN (R 4.3.2)
# munsell                  0.5.0       2018-06-12 [2] CRAN (R 4.3.2)
# paletteer                1.6.0       2024-01-21 [2] CRAN (R 4.3.2)
# pillar                   1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig                2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
# plotly                   4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
# png                      0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
# promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
# purrr                  * 1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
# R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
# rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
# RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                     1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
# RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
# readr                  * 2.1.5       2024-01-10 [2] CRAN (R 4.3.2)
# readxl                 * 1.4.3       2023-07-06 [2] CRAN (R 4.3.2)
# rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
# restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
# rhdf5                  * 2.46.1      2023-11-29 [2] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters             1.14.1      2023-11-06 [2] Bioconductor
# Rhdf5lib                 1.24.2      2024-02-07 [2] Bioconductor 3.18 (R 4.3.2)
# rjson                    0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
# rlang                    1.1.3       2024-01-10 [2] CRAN (R 4.3.2)
# rprojroot                2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools                2.18.0      2023-10-24 [2] Bioconductor
# RSQLite                  2.3.5       2024-01-21 [2] CRAN (R 4.3.2)
# rstudioapi               0.15.0      2023-07-07 [2] CRAN (R 4.3.2)
# rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer              1.62.0      2023-10-24 [2] Bioconductor
# S4Arrays               * 1.2.0       2023-10-24 [2] Bioconductor
# S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
# scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scater                   1.30.1      2023-11-16 [2] Bioconductor
# scuttle                  1.12.0      2023-10-24 [2] Bioconductor
# sessioninfo            * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
# shiny                    1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
# shinyWidgets             0.8.1       2024-01-10 [2] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0      2023-10-24 [2] Bioconductor
# spam                     2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray            * 1.2.3       2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats        1.14.0      2023-10-24 [2] Bioconductor
# SpatialExperiment      * 1.12.0      2023-10-24 [2] Bioconductor
# spatialLIBD            * 1.15.4      2024-05-01 [1] Github (LieberInstitute/spatialLIBD@77a5303)
# statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
# stringi                  1.8.3       2023-12-11 [2] CRAN (R 4.3.2)
# stringr                * 1.5.1       2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
# tibble                 * 3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                  * 1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect               1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse              * 2.0.0       2023-02-22 [2] CRAN (R 4.3.2)
# timechange               0.3.0       2024-01-18 [2] CRAN (R 4.3.2)
# tzdb                     0.4.0       2023-05-12 [2] CRAN (R 4.3.2)
# utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
# vctrs                    0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
# vipor                    0.4.7       2023-12-18 [2] CRAN (R 4.3.2)
# viridis                  0.6.5       2024-01-29 [2] CRAN (R 4.3.2)
# viridisLite              0.4.2       2023-05-02 [2] CRAN (R 4.3.2)
# withr                    3.0.0       2024-01-16 [2] CRAN (R 4.3.2)
# XML                      3.99-0.16.1 2024-01-22 [2] CRAN (R 4.3.2)
# xtable                   1.8-4       2019-04-21 [2] CRAN (R 4.3.2)
# XVector                  0.42.0      2023-10-24 [2] Bioconductor
# yaml                     2.3.8       2023-12-11 [2] CRAN (R 4.3.2)
# zlibbioc                 1.48.0      2023-10-24 [2] Bioconductor
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
