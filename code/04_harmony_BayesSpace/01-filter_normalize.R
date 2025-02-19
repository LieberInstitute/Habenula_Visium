library("spatialLIBD")
library("here")
library("tidyverse")
library("scran")
library("purrr")
library("ggplot2")
library("BiocParallel")
library("scater")
library("scry")
library("BiocSingular")
library("sessioninfo")
library("HDF5Array")

## set path directories
dir_rdata <- here("processed-data", "04_harmony_BayesSpace")
# filtered_in_path <- here("processed-data", "02_build_spe", "spe_qc_low_lib_edge.rds")
spe_in_path <- here("processed-data", "02_build_spe", "spe_scran_spotsweeper.rds")
filtered_hdf5_dir <- file.path(dir_rdata, "spe_filtered_hdf5")
dir_plots <- here("plots", "04_harmony_BayesSpace")

num_cores <- 2
set.seed(20240223)

## Create output directories
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)


## load a filtered spe object
spe <- readRDS(spe_in_path)
spe
cat("Number of spots after removed any remaining empty spots and/or genes with zero counts:", dim(spe)[2], "\n")


## Check initial outliers detected by both scran and SpotSweeper

colnames(colData(spe))
# From scran: discard = (low_lib_size | low_n_features) | high_subsets_Mito_percent)         
message("scran_low_lib_size\t\t", table(spe$scran_low_lib_size)[[1]])
message("scran_low_n_features\t\t", table(spe$scran_low_n_features)[[1]])
message("scran_high_subsets_Mito_percent\t\t", table(spe$scran_high_subsets_Mito_percent)[[1]])
message("scran_discard\t\t", table(spe$scran_discard)[[1]])
# From SS: SpotSweeper local_outliers = sum_umi_outliers | sum_gene_outliers | expr_chrM_ratio_outliers
message("SpotSweeper_sum_umi_outliers\t\t", table(spe$sum_umi_outliers)[[2]])
message("SpotSweeper_sum_gene_outliers\t\t", table(spe$sum_gene_outliers)[[2]])
message("SpotSweeper_expr_chrM_ratio_outliers\t\t", table(spe$expr_chrM_ratio_outliers)[[2]])
message("SpotSweeper_\t\t", table(spe$local_outliers)[[2]])


###############################################################################
#   Automatic selection of spots detected by scran and SpotSweeper 
################################################################################


message("Merging scran and spot sweeper outliers to remove")
v_scran_discard_keys <- unlist(spe$key[spe$scran_discard==T])
length(v_scran_discard_keys)
v_spotS_outlier_keys <- unlist(spe$key[spe$local_outliers==T])
length(v_spotS_outlier_keys)
v_keys_bad_spots <- append(as.vector(v_scran_discard_keys), as.vector(v_spotS_outlier_keys))
message("Total ouliers merged: ", length(v_keys_bad_spots))
v_keys_bad_spots <- unique(v_keys_bad_spots)
df_bad_spots <- as.data.frame(v_keys_bad_spots)

message("Total unique outliers: ", nrow(df_bad_spots))

## Remove all the outliers detected automatically

head(df_bad_spots, 3)
# match de unique IDs and get the index row from spe
m <- match(previous_work$key, spe$key)
m <- match(df_bad_spots$v_keys_bad_spots, spe$key)
m
length(m) # for example: 284 manual annotations

spe <- spe[, !spe$key %in% spe$key[c(m)]] 

# for (bad_spot in v_keys_bad_spots) { 
#   #spe <- spe[, !spe$key == outlier] 
#   spe <- spe[, spe$key != bad_spot] 
#   }


message("Outliers removed!")
message(" - Current spots: ", length(spe$key))


## Double check any remaining empty spots and/or genes with zero counts
length(spe$key[spe$in_tissue])
spe <- spe[
  rowSums(assays(spe)$counts) > 0,
  (colSums(assays(spe)$counts) > 0)
]
spe[, colSums(counts(spe)) > 0] #33409
spe[rowSums(assays(spe)$counts) > 0] #33409

message("Number of spots after removed any remaining empty spots and/or genes with zero counts:", dim(spe)[2], "\n")

## Save new spe object with spots manually annotated drop
saveRDS(spe, file.path(dir_rdata, "spe_qcED_spatialLIBD.rds"))



################################################################################
#   Compute log-normalized counts
################################################################################



#   Filter SPE: take only spots in tissue, drop spots with 0 counts for all
#   genes, and drop genes with 0 counts in every spot
message(Sys.time(), " - Running quickCluster()")

# pre-clustering step where cells in each cluster are normalized separately and the size factors are rescaled to be comparable across clusters

Sys.time()
spe$scran_quick_cluster <- quickCluster(
    spe,
    BPPARAM = MulticoreParam(num_cores),
    block = spe$sample_id,
    block.BPPARAM = MulticoreParam(num_cores),
    #use.ranks=TRUE
    #min.mean = 0.1
)
Sys.time()

## Test to avoid warning in computeSumFactors() which generate error on logNormCounts() final step:
# (1) I applied `use.ranks=TRUE`, which removes low-abundance genes with many tied ranks, especially due to zeros, which may reduce the precision of the clustering
# (2) I applied `min.mean = 0.1` for UMI data - the function will automatically try to determine this from the data if min.mean=NULL.


print("Quick cluster table:")
table(spe$scran_quick_cluster)

## deconvolution size factors normalization

message(Sys.time(), " - Running computeSumFactors()")
Sys.time()
spe <- computeSumFactors(spe,
    clusters = spe$scran_quick_cluster,
    BPPARAM = MulticoreParam(num_cores)
)
Sys.time()

# Warning message: -> I used use.ranks=TRUE to perform quickCluster()
#   In .rescale_clusters(clust.profile, ref.col = ref.clust, min.mean = min.mean) :
#   inter-cluster rescaling factor for cluster 7 is not strictly positive,
# reverting to the ratio of average library sizes

message(Sys.time(), " - Running checking sizeFactors()")
summary(sizeFactors(spe))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1491  0.4805  1.0000  1.2669 27.5203 

## plot deconvolution size factor for each cell compared to the equivalent size factor derived from the library size

lib.sf <- librarySizeFactors(spe)

pdf(here(dir_plots, "Histogram_log10_size_factor.pdf"))
hist(log10(lib.sf), xlab="Log10[Size factor]", col='grey80')
dev.off()

pdf(here(dir_plots, "Deconvolution_size_factor.pdf"))
plot(lib.sf, sizeFactors(spe), xlab="Library size factor",
     ylab="Deconvolution size factor", pch=16, # log='xy',
     col=as.integer(factor(spe$sizeFactor)))
abline(a=0, b=1, col="red")
dev.off()

## run log normalization

message(Sys.time(), " - Running logNormCounts()")

spe <- logNormCounts(spe) # Error in .local(x, ...) : size factors should be positive
assayNames(spe)

# Save spe QCed with log-counts 

saveRDS(spe, file.path(dir_rdata, "spe_qcED_spatialLIBD_log_not_QCed.rds"))

# spe <- saveHDF5SummarizedExperiment(
#     spe,
#     dir = paste0(filtered_hdf5_dir, "_temp"), replace = TRUE
# )
# gc()



# ################################################################################


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-06-13 12:09:37 EDT"
# > proc.time()
# user   system  elapsed 
# 713.525   11.483 2943.889 
# > options(width = 120)
# > session_info()
# .8       2022-06-12 [2] CRAN (R 4.3.2)
# benchmarkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
# Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
# BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
# BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
# BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
# BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           * 1.36.0      2023-10-24 [2] Bioconductor
# BiocSingular           * 1.18.0      2023-10-24 [2] Bioconductor
# BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
# Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
# bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
# bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
# bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
# blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
# bluster                  1.12.0      2023-10-24 [2] Bioconductor
# bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
# cachem                   1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
# cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
# cluster                  2.1.6       2023-12-01 [3] CRAN (R 4.3.2)
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
# dqrng                    0.3.2       2023-11-29 [2] CRAN (R 4.3.2)
# DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
# edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
# ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
# fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
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
# igraph                   2.0.1.9008  2024-02-09 [2] Github (igraph/rigraph@39158c6)
# interactiveDisplayBase   1.40.0      2023-10-24 [2] Bioconductor
# IRanges                * 2.36.0      2023-10-24 [2] Bioconductor
# irlba                    2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
# iterators                1.0.14      2022-02-05 [2] CRAN (R 4.3.2)
# jquerylib                0.1.4       2021-04-26 [2] CRAN (R 4.3.2)
# jsonlite                 1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
# KEGGREST                 1.42.0      2023-10-24 [2] Bioconductor
# later                    1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
# lattice                  0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval                 0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
# lifecycle                1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
# limma                    3.58.1      2023-10-31 [2] Bioconductor
# locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
# lubridate              * 1.9.3       2023-09-27 [2] CRAN (R 4.3.2)
# magick                   2.8.2       2023-12-20 [2] CRAN (R 4.3.2)
# magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
# maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
# Matrix                 * 1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
# matrixStats            * 1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
# memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
# metapod                  1.10.1      2023-12-24 [2] Bioconductor 3.18 (R 4.3.2)
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
# rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer              1.62.0      2023-10-24 [2] Bioconductor
# S4Arrays               * 1.2.0       2023-10-24 [2] Bioconductor
# S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
# scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
# scater                 * 1.30.1      2023-11-16 [2] Bioconductor
# scran                  * 1.30.2      2024-01-22 [2] Bioconductor 3.18 (R 4.3.2)
# scry                   * 1.14.0      2023-10-24 [2] Bioconductor
# scuttle                * 1.12.0      2023-10-24 [2] Bioconductor
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
