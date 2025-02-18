library("ggplot2")
library("SpotSweeper")
library("spatialLIBD")
library("purrr")
library("here")
library("sessioninfo")
library("ggpubr")

## code adapted from https://github.com/LieberInstitute/LFF_spatial_ERC/blob/ad7546fddf48d4047bbc95ab7c980961ac0a3549/code/02_build_spe/04_SpotSweeper.R

plot_dir <- here("plots", "02_build_spe")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

data_dir <- here("processed-data", "02_build_spe")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)


set.seed(20241105)
# spe <- readRDS(here("processed-data", "02_build_spe", "spe_raw.rds"))
spe <- readRDS(here("processed-data", "02_build_spe", "spe_qc_low_lib_edge.rds"))
spe
colnames(colData(spe))
dim(spe)
# [1] 27028 44681

cat("Spots in tissue:", dim(spe)[2], "\n")
# Spots in tissue: 44681

lobstr::obj_size(spe)
# 5.96 GB

## SpotSweeper variables sum_umi, sum_gene and expr_chrM_ratio are created
message(Sys.time(), "Local Outliers - low sum_umi")
spe <- localOutliers(spe,
                     metric = "sum_umi",
                     direction = "lower",
                     log = TRUE)
table(spe$sum_umi_outliers)

message(Sys.time(), "Local Outliers - low sum_gene")
spe <- localOutliers(spe,
                     metric = "sum_gene",
                     direction = "lower",
                     log = TRUE)
table(spe$sum_gene_outliers)

message(Sys.time(), "Local Outliers - High expr_chrM_ratio")
spe <- localOutliers(spe,
                     metric = "expr_chrM_ratio",
                     direction = "higher",
                     log = FALSE)

## combine all outliers into "local_outliers" column
spe$local_outliers <- as.logical(spe$sum_umi_outliers) |
  as.logical(spe$sum_gene_outliers) |
  as.logical(spe$expr_chrM_ratio_outliers)

message("Local Outliers")
table(spe$local_outliers)

message("Local Outliers on edge")
spe$edge_spot
table(spe$local_outliers, spe$edge_spot)

#### find artifacts using SpotSweeper ####
## Only works one sample at a time

unique(spe$sample_id)
# [1] "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1" "V13B23-285_D1"
# [5] "V13B23-281_A1" "V13B23-281_B1" "V13B23-281_C1" "V13B23-281_D1"
# [9] "V14F07-340_A1" "V14F07-340_B1" "V14F07-340_C1" "V14F07-340_D1"
# [13] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1"
# colnames(colData(spe))

## Double check no NAs or zeros on the data due SpotSweeper request none zeros to compute findArtifacts()
if (any(colSums(counts(spe)) == 0)) {
    message("Removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}

message(Sys.time(), " - findArtifact ring 5")
## default 5 neighbors
artifact_df <- purrr::map_dfr(unique(spe$sample_id), function(samp){
  message(Sys.time(), " - ", samp)
  spe_temp <- findArtifacts(spe[,spe$sample_id == samp],
                            mito_percent = "expr_chrM_ratio",
                            mito_sum = "expr_chrM",
                            n_rings = 5,
                            name = "artifact"
  )

  return(as.data.frame(colData(spe_temp)[,c("sample_id", "key", "artifact")]))
})
## Add artifact to spe
identical(spe$key, artifact_df$key) #TRUE
spe$artifact <- artifact_df$artifact

# ## testing with 20 neighbors
# message(Sys.time(), " - findArtifact ring 20")
# artifact_df2 <- purrr::map_dfr(unique(spe$sample_id), function(samp){
#   message(Sys.time(), " - ", samp)
#   spe_temp <- findArtifacts(spe[,spe$sample_id == samp],
#                             mito_percent = "expr_chrM_ratio",
#                             mito_sum = "expr_chrM",
#                             n_rings = 20,
#                             name = "artifact_r20"
#   )
#
#   return(as.data.frame(colData(spe_temp)[,c("sample_id", "key", "artifact_r20")]))
# })
#
#
# ## Add artifact to spe
# identical(spe$key, artifact_df2$key) #TRUE
# spe$artifact <- artifact_df2$artifact_r20



#### save spot sweeper data ####

spotsweeper_data <- as.data.frame(colData(spe)[,c("sample_id", "key", "array_row", "array_col", "sum_umi_outliers", "sum_gene_outliers", "expr_chrM_ratio_outliers", "local_outliers", "artifact")])
head(spotsweeper_data)
write.csv(spotsweeper_data, file = here(data_dir, "SpotSweeper_outliers_detected.csv"))

point_size = 1.1


#### plotting ####

plot_all_spot_sweep <- function(spe, sample = unique(spe$sample_id)[1]){

  spe <- spe[,spe$sample_id == sample]

  # library size
  p1 <- plotQC(spe, metric = "sum_umi_log", outliers = "SpotSweeper_sum_umi_outliers", point_size = point_size) +
    ggtitle(paste(sample, "Sum UMI"))

  # unique genes
  p2 <- plotQC(spe, metric = "sum_gene_log", outliers = "SpotSweeper_sum_gene_outliers", point_size = point_size) +
    ggtitle("Sum Genes")

  # mitochondrial percent
  p3 <- plotQC(spe, metric = "expr_chrM_ratio", outliers = "SpotSweeper_expr_chrM_ratio_outliers", point_size = point_size) +
    ggtitle("ChrM Ratio")

  # all local outliers
  p4 <- plotQC(spe, metric = "sum_umi_log", outliers = "SpotSweeper_local_outliers", point_size = point_size, stroke = 0.75) +
    ggtitle("All Local Outliers")

  ## artifact
  p5 <- plotQC(spe, metric = "expr_chrM_ratio", outliers = "SpotSweeper_artifact", point_size = point_size, stroke = 0.75) + # sum_umi_log
    ggtitle("Artifact")

  plot_list <- list(p1, p2, p3, p4, p5)
  ggarrange(
    plotlist = plot_list,
    ncol = 3, nrow = 2,
    common.legend = FALSE)
}

message(Sys.time(), " - Plot SpotSweeper QC for all samples")
pdf(here(plot_dir, "SpotSweeper_ALL_QC_sample.pdf"), width = 12, height = 10)
purrr::map(sort(unique(spe$sample_id)), ~plot_all_spot_sweep(spe = spe, sample = .x))
dev.off()


## save a list of plots by specific metric

plt_plots <- function(plt_list, plot_name) {
  pdf(file = here(plot_dir, plot_name), width = 12, height = 10)
  # title_main <- "local_outliers - sum_umi"
  plt_main <- ggarrange(
    plotlist = plt_list,
    ncol = 4, nrow = 4,
    common.legend = FALSE,
    font.label=list(color="black",size=8)
  )
  print(plt_main)
  dev.off()
}


## Only plot log sum-umi local outliers

plot_all_local_ouliers_sum_umi <- function(spe, sample = unique(spe$sample_id)[1]){
  spe <- spe[,spe$sample_id == sample]
  plt1 <- plotQC(spe, metric = "sum_umi_log", outliers = "local_outliers",
                 point_size = point_size, stroke = 0.75) +  ggtitle(sample)
}

local_ouliers_plts <- purrr::map(sort(unique(spe$sample_id)), ~
                                   plot_all_local_ouliers_sum_umi(spe = spe, sample = .x))
plt_plots(local_ouliers_plts, "SpotSweeper_local_umi_ouliers.pdf")


## Only plot high Mito local outliers

plot_all_local_ouliers_mito <- function(spe, sample = unique(spe$sample_id)[1]){
  spe <- spe[,spe$sample_id == sample]
  plt1 <- plotQC(spe, metric = "expr_chrM_ratio", outliers = "local_outliers",
                 point_size = point_size, stroke = 0.75) +  ggtitle(sample)
}

local_ouliers_plts <- purrr::map(sort(unique(spe$sample_id)), ~
                                   plot_all_local_ouliers_mito(spe = spe, sample = .x))
plt_plots(local_ouliers_plts, "SpotSweeper_local_mito_ouliers.pdf")


## Only plot log sum-gene local outliers

plot_all_local_ouliers_gene <- function(spe, sample = unique(spe$sample_id)[1]){
  spe <- spe[,spe$sample_id == sample]
  plt1 <- plotQC(spe, metric = "sum_gene_log", outliers = "local_outliers",
                 point_size = point_size, stroke = 0.75) +  ggtitle(sample)
}

local_ouliers_plts <- purrr::map(sort(unique(spe$sample_id)), ~
                                   plot_all_local_ouliers_gene(spe = spe, sample = .x))
plt_plots(local_ouliers_plts, "SpotSweeper_local_gene_ouliers.pdf")


## Only plot artifacts detected using mito dev and ring 5 and 20

# colnames(colData(spe))
plot_all_local_artifacts <- function(spe, sample = unique(spe$sample_id)[1]){
  spe <- spe[,spe$sample_id == sample]
  plt1 <- plotQC(spe, metric = "expr_chrM_ratio", outliers = "artifact",
                 point_size = 0.3, stroke = 0.75) + ggtitle("Artifact")
}
local_ouliers_plts <- purrr::map(sort(unique(spe$sample_id)), ~
                                   plot_all_local_artifacts(spe = spe, sample = .x))
plt_plots(local_ouliers_plts, "SpotSweeper_mito_artifacts_r5.pdf")

# plot_all_local_artifacts_r20 <- function(spe, sample = unique(spe$sample_id)[1]){
#   spe <- spe[,spe$sample_id == sample]
#   plt1 <- plotQC(spe, metric = "expr_chrM_ratio", outliers = "artifact_r20",
#                  point_size = 0.3, stroke = 0.75) + ggtitle("Artifact_r20")
# }
# local_ouliers_plts <- purrr::map(sort(unique(spe$sample_id)), ~
#                                    plot_all_local_artifacts_r20(spe = spe, sample = .x))
# plt_plots(local_ouliers_plts, "SpotSweeper_mito_artifacts_r20.pdf")

colnames(colData(spe))[grep("outliers", colnames(colData(spe)))]
# [1] "sum_umi_outliers"         "sum_gene_outliers"
# [3] "expr_chrM_ratio_outliers" "local_outliers"


## Define order of samples for the grid plots

slide_order <- sort(unique(spe$sample_id))
sample_order <- unlist(sapply(slide_order, function(i) {
    sort(unique(spe$sample_id)[grepl(i, unique(spe$sample_id))])
}))
## Re-order samples
new_order <- unlist(lapply(sample_order, function(i) {
  which(spe$sample_id == i)
}))
stopifnot(all(seq_len(ncol(spe)) %in% new_order))
spe <- spe[, new_order]
unique(spe$sample_id)
# [1] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1"
# [5] "V13B23-285_A1" "V13B23-285_B1" "V13B23-285_C1" "V13B23-285_D1"
# [9] "V14F07-340_A1" "V14F07-340_B1" "V14F07-340_C1" "V14F07-340_D1"

## Note qc_low_lib_edge from scran have been removed from this object
data_dir <- here("processed-data", "02_build_spe")
saveRDS(spe, file.path(data_dir, "spe_scran_spotsweeper.rds"))


message("Process completed!!!")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



# > ## Reproducibility information
#   > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2025-02-05 13:31:49 EST"
# > proc.time()
# user   system  elapsed
# 689.115    4.828 2503.095
# > options(width = 120)
# > session_info()
# (R 4.4.2)
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.4.0)
# benchmarkme            1.0.8     2022-06-12 [2] CRAN (R 4.4.0)
# benchmarkmeData        1.0.4     2020-04-23 [2] CRAN (R 4.4.0)
# Biobase              * 2.66.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocFileCache          2.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocGenerics         * 0.52.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocIO                 1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocManager            1.30.25   2024-08-28 [2] CRAN (R 4.4.1)
# BiocNeighbors          2.0.1     2024-11-28 [2] Bioconductor 3.20 (R 4.4.2)
# BiocParallel           1.40.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocSingular           1.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocVersion            3.20.0    2024-05-01 [2] Bioconductor 3.20 (R 4.4.0)
# Biostrings             2.74.1    2024-12-16 [2] Bioconductor 3.20 (R 4.4.2)
# bit                    4.5.0.1   2024-12-03 [2] CRAN (R 4.4.2)
# bit64                  4.6.0-1   2025-01-16 [2] CRAN (R 4.4.2)
# bitops                 1.0-9     2024-10-03 [2] CRAN (R 4.4.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.4.0)
# broom                  1.0.7     2024-09-26 [2] CRAN (R 4.4.1)
# bslib                  0.9.0     2025-01-30 [2] CRAN (R 4.4.2)
# cachem                 1.1.0     2024-05-16 [2] CRAN (R 4.4.0)
# car                    3.1-3     2024-09-27 [2] CRAN (R 4.4.1)
# carData                3.0-5     2022-01-06 [2] CRAN (R 4.4.0)
# circlize               0.4.16    2024-02-20 [2] CRAN (R 4.4.0)
# cli                    3.6.3     2024-06-21 [2] CRAN (R 4.4.1)
# clue                   0.3-66    2024-11-13 [2] CRAN (R 4.4.2)
# cluster                2.1.8     2024-12-11 [3] CRAN (R 4.4.2)
# codetools              0.2-20    2024-03-31 [3] CRAN (R 4.4.2)
# colorspace             2.1-1     2024-07-26 [2] CRAN (R 4.4.1)
# ComplexHeatmap         2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# config                 0.3.2     2023-08-30 [2] CRAN (R 4.4.0)
# cowplot                1.1.3     2024-01-22 [2] CRAN (R 4.4.0)
# crayon                 1.5.3     2024-06-20 [2] CRAN (R 4.4.1)
# curl                   6.2.0     2025-01-23 [2] CRAN (R 4.4.2)
# data.table             1.16.4    2024-12-06 [2] CRAN (R 4.4.2)
# DBI                    1.2.3     2024-06-02 [2] CRAN (R 4.4.0)
# dbplyr                 2.5.0     2024-03-19 [2] CRAN (R 4.4.0)
# DelayedArray           0.32.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# digest                 0.6.37    2024-08-19 [2] CRAN (R 4.4.1)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.4.0)
# dplyr                  1.1.4     2023-11-17 [2] CRAN (R 4.4.0)
# DT                     0.33      2024-04-04 [2] CRAN (R 4.4.0)
# edgeR                  4.4.2     2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
# escheR                 1.6.0     2024-10-29 [1] Bioconductor 3.20 (R 4.4.2)
# ExperimentHub          2.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# farver                 2.1.2     2024-05-13 [2] CRAN (R 4.4.0)
# fastmap                1.2.0     2024-05-15 [2] CRAN (R 4.4.0)
# filelock               1.0.3     2023-12-11 [2] CRAN (R 4.4.0)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.4.0)
# Formula                1.2-5     2023-02-24 [2] CRAN (R 4.4.0)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.42.3    2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
# GenomeInfoDbData       1.2.13    2024-10-01 [2] Bioconductor
# GenomicAlignments      1.42.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# GenomicRanges        * 1.58.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.4.0)
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.4.0)
# ggplot2              * 3.5.1     2024-04-23 [2] CRAN (R 4.4.0)
# ggpubr               * 0.6.0     2023-02-10 [2] CRAN (R 4.4.0)
# ggrepel                0.9.6     2024-09-07 [2] CRAN (R 4.4.1)
# ggsignif               0.6.4     2022-10-13 [2] CRAN (R 4.4.0)
# GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.4.0)
# glue                   1.8.0     2024-09-30 [2] CRAN (R 4.4.1)
# golem                  0.5.1     2024-08-27 [2] CRAN (R 4.4.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.4.0)
# gtable                 0.3.6     2024-10-25 [2] CRAN (R 4.4.2)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.4.0)
# htmltools              0.5.8.1   2024-04-04 [2] CRAN (R 4.4.0)
# htmlwidgets            1.6.4     2023-12-06 [2] CRAN (R 4.4.0)
# httpuv                 1.6.15    2024-03-26 [2] CRAN (R 4.4.0)
# httr                   1.4.7     2023-08-15 [2] CRAN (R 4.4.0)
# IRanges              * 2.40.1    2024-12-05 [2] Bioconductor 3.20 (R 4.4.2)
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.4.0)
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.4.0)
# jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.4.0)
# jsonlite               1.8.9     2024-09-20 [2] CRAN (R 4.4.1)
# KEGGREST               1.46.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.4.0)
# later                  1.4.1     2024-11-27 [2] CRAN (R 4.4.2)
# lattice                0.22-6    2024-03-20 [3] CRAN (R 4.4.2)
# lazyeval               0.2.2     2019-03-15 [2] CRAN (R 4.4.0)
# lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.4.0)
# limma                  3.62.2    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# locfit                 1.5-9.11  2025-02-03 [2] CRAN (R 4.4.2)
# magick                 2.8.5     2024-09-20 [2] CRAN (R 4.4.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.4.0)
# MASS                   7.3-64    2025-01-04 [3] CRAN (R 4.4.2)
# Matrix                 1.7-2     2025-01-23 [3] CRAN (R 4.4.2)
# MatrixGenerics       * 1.18.1    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# matrixStats          * 1.5.0     2025-01-07 [2] CRAN (R 4.4.2)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.4.0)
# mime                   0.12      2021-09-28 [2] CRAN (R 4.4.0)
# munsell                0.5.1     2024-04-01 [2] CRAN (R 4.4.0)
# paletteer              1.6.0     2024-01-21 [2] CRAN (R 4.4.0)
# pillar                 1.10.1    2025-01-07 [2] CRAN (R 4.4.2)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.4.0)
# plotly                 4.10.4    2024-01-13 [2] CRAN (R 4.4.0)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.4.0)
# promises               1.3.2     2024-11-28 [2] CRAN (R 4.4.2)
# purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.4.0)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.4.0)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.4.0)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.4.0)
# Rcpp                   1.0.14    2025-01-12 [2] CRAN (R 4.4.2)
# RCurl                  1.98-1.16 2024-07-11 [2] CRAN (R 4.4.1)
# rematch2               2.1.2     2020-05-01 [2] CRAN (R 4.4.0)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.4.0)
# rjson                  0.2.23    2024-09-16 [2] CRAN (R 4.4.1)
# rlang                  1.1.5     2025-01-17 [2] CRAN (R 4.4.2)
# rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.4.0)
# Rsamtools              2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# RSQLite                2.3.9     2024-12-03 [2] CRAN (R 4.4.2)
# rstatix                0.7.2     2023-02-01 [2] CRAN (R 4.4.0)
# rstudioapi             0.17.1    2024-10-22 [2] CRAN (R 4.4.2)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.4.0)
# rtracklayer            1.66.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# S4Arrays               1.6.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# S4Vectors            * 0.44.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# sass                   0.4.9     2024-03-15 [2] CRAN (R 4.4.0)
# ScaledMatrix           1.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# scales                 1.3.0     2023-11-28 [2] CRAN (R 4.4.0)
# scater                 1.34.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# scuttle                1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.4.0)
# shape                  1.4.6.1   2024-02-23 [2] CRAN (R 4.4.0)
# shiny                  1.10.0    2024-12-14 [2] CRAN (R 4.4.2)
# shinyWidgets           0.8.7     2024-09-23 [2] CRAN (R 4.4.1)
# SingleCellExperiment * 1.28.1    2024-11-10 [2] Bioconductor 3.20 (R 4.4.2)
# SparseArray            1.6.1     2025-01-19 [2] Bioconductor 3.20 (R 4.4.2)
# spatialEco             2.0-2     2023-11-17 [1] CRAN (R 4.4.2)
# SpatialExperiment    * 1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# spatialLIBD          * 1.19.6    2025-02-04 [1] Github (LieberInstitute/spatialLIBD@27db42d)
# SpotSweeper          * 1.2.0     2024-11-13 [1] Bioconductor 3.20 (R 4.4.2)
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.4.0)
# SummarizedExperiment * 1.36.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# terra                  1.8-15    2025-01-24 [2] CRAN (R 4.4.2)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.4.0)
# tidyr                  1.3.1     2024-01-24 [2] CRAN (R 4.4.0)
# tidyselect             1.2.1     2024-03-11 [2] CRAN (R 4.4.0)
# UCSC.utils             1.2.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# vctrs                  0.6.5     2023-12-01 [2] CRAN (R 4.4.0)
# vipor                  0.4.7     2023-12-18 [2] CRAN (R 4.4.0)
# viridis                0.6.5     2024-01-29 [2] CRAN (R 4.4.0)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.4.0)
