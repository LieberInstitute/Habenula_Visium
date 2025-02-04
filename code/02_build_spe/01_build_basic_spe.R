library("spatialLIBD")
library("tidyverse")
library("here")
library("lobstr")
library("sessioninfo")



## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
if (!dir.exists(dir_rdata)) {
    dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
}

## Define the donor info using information from
## https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/raw-data/Visium_SPG_AD_ITG_MasterExcelSummarySheet.xlsx

## Sample ID
Sid <- paste0("S", rep(1:16), "_v") 
## Slide number and slide numver
array_id <- strsplit(c("A1 B1 C1 D1 A1 B1 C1 D1 A1 B1 C1 D1 A1 B1 C1 D1"), "\\s+")[[1]]
sample_slide_id <- strsplit(c("V13B23-285 V13B23-285 V13B23-285 V13B23-285 V13B23-281 V13B23-281 V13B23-281 V13B23-281 V14F07-340 V14F07-340 V14F07-340 V14F07-340 V13B23-280 V13B23-280 V13B23-280 V13B23-280"), "\\s+")[[1]]
sample_slide_id <- paste0(sample_slide_id, "_", array_id)
## Brain ID
brain_id <- c(rep("Br8518", 4), rep("Br6522", 4), rep("Br9090",4), rep("Br9037",4)) 
brain_area = strsplit(c("PR5 AR6 AR6 AR6 PR6 PR6 PR6 PR6 PL6 PL6 PL6 PL6 AL5 AL5 AL5 AL5"), "\\s+")[[1]]

## Define some info for the samples
sample_info <- data.frame(
  sample_id_short = c(Sid), # S10_v
  sample_id = c(sample_slide_id), #V14F07-340_D1
  brain_id = c(brain_id),
  brain_area = c(brain_area),
  age = c(rep(41.3, 8), rep(57.5, 4), rep(51.4, 4)), 
  sex = c(rep("F", 4), rep("M", 12)),  
  ethnicity  = c(rep("EA/CAUC", 12), rep("AA", 4)), 
  pmi = c(rep(10.5, 4), rep(30.5, 4), rep(26.5, 4), rep(18.5, 4)), 
  diagnosis = c(rep("Control", 16)), 
  rin = c(rep(6.3, 4), rep(7.4, 4), rep(6.1, 4), rep(7.8, 4))
)
sample_info$sample_id_short <- sprintf("S%02d_Hb_V", parse_number(sample_info$sample_id_short))

sample_info$sample_path <-
    file.path(
        here::here("processed-data", "01_spaceranger"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Combine sample info with the donor info
sample_info[c(colnames(sample_info))] #"sample_id", "subject", "age", "sex", "race", "pmi", "diagnosis", "rin"

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    #sample_info$sample_id_short,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    #reference_gtf = NULL
    reference_gtf = "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
)
Sys.time()
# 2024-02-22 14:43:49.823746 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2024-02-22 14:43:56.531453 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2024-02-22 14:43:56.687378 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2024-02-22 14:43:56.86425 rtracklayer::import: reading the reference GTF file
# 2024-02-22 14:44:28.499147 adding gene information to the SPE object
# 2024-02-22 14:44:28.523581 adding information used by spatialLIBD
# [1] "2024-02-22 14:44:28 EST"

# class: SpatialExperiment
# dim: 36601 24960

# spe@int_colData$reducedDims
# colnames(spe)
# rownames(spe)

colnames(colData(spe))
head(spe$sample_id)
tail(spe$sample_id)
head(assays(spe)$counts[, 1:10]) # colnames=samples; rownames=genes

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

# head(colData(spe))

# ## Read in cell counts and segmentation results
# segmentations_list <-
#     lapply(sample_info$sample_id, function(sampleid) {
#         file <-
#             here(
#                 "processed-data",
#                 "spaceranger",
#                 sampleid,
#                 "outs",
#                 "spatial",
#                 "tissue_spot_counts.csv"
#             )
#         if (!file.exists(file)) {
#             return(NULL)
#         }
#         x <- read.csv(file)
#         x$key <- paste0(x$barcode, "_", sampleid)
#         return(x)
#     })
# ## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <-
#     Reduce(function(...) {
#         merge(..., all = TRUE)
#     }, segmentations_list[lengths(segmentations_list) > 0])
#
# ## Add the information
# segmentation_match <- match(spe$key, segmentations$key)
# segmentation_info <-
#     segmentations[segmentation_match, -which(
#         colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
#     )]
# colData(spe) <- cbind(colData(spe), segmentation_info)


message("Initial number of spots:", dim(spe)[2], "\n")

## Remove genes with no data
expr <- which(rowSums(counts(spe)) > 0)
message("Number of genes with counts:", length(expr))

no_expr <- which(rowSums(counts(spe)) == 0)
message("Number of genes with no counts:", length(no_expr))
# Number of genes with no counts: 13162
# Merged samples: Number of genes with no counts: 10315

message("% genes with counts:", (length(expr) / nrow(spe)) * 100)
message("% genes with no counts:", (length(no_expr) / nrow(spe)) * 100)


spe <- spe[-no_expr, ]

## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw <- spe

## Size in Gb
lobstr::obj_size(spe_raw)
# 1.81 GB

saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))

message("Saved spe_raw.rds")

## Now drop the spots outside the tissue
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
message("Spots in tissue:", dim(spe)[2], "\n")

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}
dim(spe)

lobstr::obj_size(spe)
# Merged samples: 1.16 MB

saveRDS(spe, file.path(dir_rdata, "spe.rds"))

message("Saved spe.rds")

message("Completed!!!")

# library("slurmjobs")
# job_single(
#   name = "01_build_basic_spe", memory = "50G", cores = 2, create_shell = TRUE
# )

# ==============================================================================

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2025-02-04 10:08:11 EST"
# > proc.time()
# user   system  elapsed 
# 162.239    3.853 1682.249 
# > options(width = 120)
# > session_info()
# CRAN (R 4.4.0)
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
# bslib                  0.9.0     2025-01-30 [2] CRAN (R 4.4.2)
# cachem                 1.1.0     2024-05-16 [2] CRAN (R 4.4.0)
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
# DelayedMatrixStats     1.28.1    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# digest                 0.6.37    2024-08-19 [2] CRAN (R 4.4.1)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.4.0)
# dplyr                * 1.1.4     2023-11-17 [2] CRAN (R 4.4.0)
# dqrng                  0.4.1     2024-05-28 [2] CRAN (R 4.4.0)
# DropletUtils           1.26.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# DT                     0.33      2024-04-04 [2] CRAN (R 4.4.0)
# edgeR                  4.4.2     2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
# ExperimentHub          2.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# fastmap                1.2.0     2024-05-15 [2] CRAN (R 4.4.0)
# filelock               1.0.3     2023-12-11 [2] CRAN (R 4.4.0)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.4.0)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.4.0)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.42.3    2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
# GenomeInfoDbData       1.2.13    2024-10-01 [2] Bioconductor
# GenomicAlignments      1.42.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# GenomicRanges        * 1.58.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.4.0)
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.4.0)
# ggplot2              * 3.5.1     2024-04-23 [2] CRAN (R 4.4.0)
# ggrepel                0.9.6     2024-09-07 [2] CRAN (R 4.4.1)
# GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.4.0)
# glue                   1.8.0     2024-09-30 [2] CRAN (R 4.4.1)
# golem                  0.5.1     2024-08-27 [2] CRAN (R 4.4.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.4.0)
# gtable                 0.3.6     2024-10-25 [2] CRAN (R 4.4.2)
# HDF5Array              1.34.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.4.0)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.4.0)
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
# later                  1.4.1     2024-11-27 [2] CRAN (R 4.4.2)
# lattice                0.22-6    2024-03-20 [3] CRAN (R 4.4.2)
# lazyeval               0.2.2     2019-03-15 [2] CRAN (R 4.4.0)
# lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.4.0)
# limma                  3.62.2    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# lobstr               * 1.1.2     2022-06-22 [2] CRAN (R 4.4.0)
# locfit                 1.5-9.11  2025-02-03 [2] CRAN (R 4.4.2)
# lubridate            * 1.9.4     2024-12-08 [2] CRAN (R 4.4.2)
# magick                 2.8.5     2024-09-20 [2] CRAN (R 4.4.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.4.0)
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
# prettyunits            1.2.0     2023-09-24 [2] CRAN (R 4.4.0)
# promises               1.3.2     2024-11-28 [2] CRAN (R 4.4.2)
# purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.4.0)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.4.0)
# R.oo                   1.27.0    2024-11-01 [2] CRAN (R 4.4.2)
# R.utils                2.12.3    2023-11-18 [2] CRAN (R 4.4.0)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.4.0)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.4.0)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.4.0)
# Rcpp                   1.0.14    2025-01-12 [2] CRAN (R 4.4.2)
# RCurl                  1.98-1.16 2024-07-11 [2] CRAN (R 4.4.1)
# readr                * 2.1.5     2024-01-10 [2] CRAN (R 4.4.0)
# rematch2               2.1.2     2020-05-01 [2] CRAN (R 4.4.0)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.4.0)
# rhdf5                  2.50.2    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# rhdf5filters           1.18.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# Rhdf5lib               1.28.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# rjson                  0.2.23    2024-09-16 [2] CRAN (R 4.4.1)
# rlang                  1.1.5     2025-01-17 [2] CRAN (R 4.4.2)
# rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.4.0)
# Rsamtools              2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# RSQLite                2.3.9     2024-12-03 [2] CRAN (R 4.4.2)
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
# sparseMatrixStats      1.18.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# SpatialExperiment    * 1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# spatialLIBD          * 1.19.6    2025-02-04 [1] Github (LieberInstitute/spatialLIBD@27db42d)
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.4.0)
# stringi                1.8.4     2024-05-06 [2] CRAN (R 4.4.0)
# stringr              * 1.5.1     2023-11-14 [2] CRAN (R 4.4.0)
# SummarizedExperiment * 1.36.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.4.0)
# tidyr                * 1.3.1     2024-01-24 [2] CRAN (R 4.4.0)
# tidyselect             1.2.1     2024-03-11 [2] CRAN (R 4.4.0)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.4.0)
# timechange             0.3.0     2024-01-18 [2] CRAN (R 4.4.0)
# tzdb                   0.4.0     2023-05-12 [2] CRAN (R 4.4.0)
# UCSC.utils             1.2.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# vctrs                  0.6.5     2023-12-01 [2] CRAN (R 4.4.0)
# vipor                  0.4.7     2023-12-18 [2] CRAN (R 4.4.0)
# viridis                0.6.5     2024-01-29 [2] CRAN (R 4.4.0)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.4.0)
# withr                  3.0.2     2024-10-28 [2] CRAN (R 4.4.2)
# XML                    3.99-0.18 2025-01-01 [2] CRAN (R 4.4.2)
# xtable                 1.8-4     2019-04-21 [2] CRAN (R 4.4.0)
# XVector                0.46.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# yaml                   2.3.10    2024-07-26 [2] CRAN (R 4.4.1)
# zlibbioc               1.52.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# 
# [1] /users/csoto/R/4.4.x
# [2] /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/library