# # library(slurmjobs)
# # slurmjobs::job_single('01_compute_cor', create_shell = TRUE, memory = '20G', command = "01_compute_cor.R")
# 
# # To submit the job use: sbatch 01_compute_cor.sh

library("here")
library("purrr")
library("spatialLIBD")
library("ComplexHeatmap")
library("sessioninfo")
library("tidyverse")


## Input dir
rds_input <- here("processed-data", "05_snRNA-seq_model_stats", "enrichment_snRNA-multiome_v2.rds")
## Create output directories
dir_rdata <- here("processed-data", "07_spatial_registration_vs_multiome_snRNA-seq")
dir_plot <- here("plots", "07_spatial_registration_vs_multiome_snRNA-seq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

# ## specify the number of BayesSpace k to use 
# k <- seq(2,28)

#################### Fine clusters snRNAseq vs Multiome snRNAseq ##################################

## load multiome snRNAseq t-stats enrichment data

sn_multiome_data <- readRDS(rds_input)
head(sn_multiome_data$enrichment[5:10])
#                 t_stat_C.05.DD_LHb t_stat_C.06 t_stat_C.07.DD_MHb t_stat_C.08
# ENSG00000238009          1.2234449    1.759594         -0.5585619  2.26415183
# ENSG00000241860          0.2151696    1.234082         -0.2209549  0.07914740
# ENSG00000237491          0.4194972    1.935467         -0.7272817  1.36389851

results_enrichment_multiome <- readRDS(rds_input)$enrichment |>
  filter(!duplicated(ensembl))
rownames(results_enrichment_multiome)
class(results_enrichment_multiome) # [1] "data.frame"

## load snRNAseq t-stats enrichment data (fine resolution)

results_enrichment <-
  readRDS(here(
    "processed-data",
    "05_snRNA-seq_model_stats",
    paste0("enrichment_final_Annotations.rds")
  ))
colnames(results_enrichment)

## filter only enrichment t-stats

modeling_res_enrichment <- results_enrichment[, grep("^t_stat_", colnames(results_enrichment))]
colnames(modeling_res_enrichment) <- gsub("^t_stat_", "", colnames(modeling_res_enrichment))
# modeling_res_enrichment[1:3,]
# Astrocyte        Endo Excit.Thal Inhib.Thal      LHb.1
# ENSG00000238009  0.6873346  1.23174976  1.2820563   1.526893 -0.6502180
# ENSG00000241860  0.7626024 -2.42952913  1.4923036   2.026956  1.2102714
# ENSG00000237491 -0.2347861 -0.06842523  0.8264148   1.981078 -0.1328745
# modeling_res_enrichment <- list(enrichment = modeling_res_enrichment)
# class(modeling_res_enrichment)

cor_fine <- layer_stat_cor(
  stats = results_enrichment_multiome, # data.frame / y-axis
  modeling_results = list(enrichment = modeling_res_enrichment), # list / x-axis
  model_type = "enrichment",
  top_n = 100
)

table(rownames(results_enrichment_multiome) %in% rownames(modeling_res_enrichment))


## Old version to correlate multiome vs snRNAseq 

cor_fine <- layer_stat_cor(
  stats = modeling_res_enrichment, # data.frame
  modeling_results = sn_multiome_data,
  model_type = "enrichment",
  top_n = 100
)
head(modeling_res_enrichment)
## With default confidence and cutoff_merge_ratio 

annotated_clusters_fine <- annotate_registered_clusters(cor_fine, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)

# ## Use annotation labels on the correlation matrices
# rownames(cor_fine) <- paste0(rownames(cor_fine), " ~ ", 
#                              annotated_clusters_fine$layer_label[match(rownames(cor_fine), annotated_clusters_fine$cluster)])
# 
# ## With default confidence and cutoff_merge_ratio 
# annotated_clusters_fine <- annotate_registered_clusters(cor_fine, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)

head(cor_fine)
#           C.01        C.02       C.03        C.04  C.05.DD_LHb        C.06
# LHb.2 -0.1959626 -0.06687099 -0.1797187 -0.04842418  0.636548725  0.01036223
# LHb.7 -0.1068896 -0.01307742 -0.1562915 -0.02322458  0.735296169 -0.01500998
# LHb.6 -0.1717414  0.01985625 -0.1635756 -0.22551095  0.147718838 -0.22087233

save(cor_fine,
     file = file.path(dir_rdata, "cor_multiome_vs_snRNA-seq_top100.Rdata")
)


##   Make heatmaps fine clusters snRNAseq vs Multiome snRNAseq
                    
plt_name <- paste0("cor_top100_registration_snMultiome_snRNAseq_v2.pdf")
pdf(here(dir_plot, plt_name))

layer_stat_cor_plot(
  cor_fine, annotation = annotated_clusters_fine,
  heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10)
)

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
