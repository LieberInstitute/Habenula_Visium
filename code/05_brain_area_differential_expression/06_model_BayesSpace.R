########################################################################
## Create Bayes-Space model 
## Authors
# copied from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/07_layer_differential_expression/03_model_BayesSpace.R
## Adapted: CSC
## Data: XX
## For 60 to 80k spots: $srun --pty --mem=30GB --x11 bash
########################################################################


k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## read input arguments
# args = commandArgs(trailingOnly = TRUE)
# k <- args[2]
# 13 = 1 Hb domain
# 21 = 2 Hb domains
# 26 = 3 Hb domains

## If empty or testing
if (is.na(k)) {
  k <- 13
}

k_nice <- sprintf("%02d", k)

library("here")
library("sessioninfo")
library("spatialLIBD")

# > packageVersion("spatialLIBD")
# [1] ‘1.21.5’

## output directory
dir_rdata <- here("processed-data","05_brain_area_differential_expression")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here("plots","05_brain_area_differential_expression")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load spe_pseudo data
sce_pseudo <-
  readRDS(
    file.path(
      dir_rdata,
      paste0("sce_pseudo_PCA_brain_area_k", k_nice, ".rds")
    )
  )
## verification
sce_pseudo
colnames(colData(sce_pseudo))
# [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"    
# [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
# [9] "nspots"          "pmi"             "rin"             "sample_id"      
# [13] "sex"             "sum_umi"   
table(sce_pseudo[["brain_id"]])
table(sce_pseudo[["brain_area2"]])
# G0 G1 G2 G3 G4 
# 2  9  9  9  6 
levels(sce_pseudo$BayesSpace)

# ## rename annotated BayesSpace names/levels to fix error when computing 'registration_stats_pairwise' - we need syntactically valid names
# old_bs_names <- sce_pseudo$BayesSpace
# new_bs_names <- gsub("/", "_",  # substitute cell-type '/' separator
#                      gsub(" ~ ", ".", sce_pseudo$BayesSpace)) # substitute ' ~ '
# new_bs_names <- gsub("\\*$", "",  new_bs_names) # substitute '*' added to significant cell-types
# unique(new_bs_names)
# # [1] "Sp13D01.Oligo"      "Sp13D02.Endo"       "Sp13D03.Endo"
# # [4] "Sp13D04.Astrocyte"  "Sp13D06.Astrocyte"  "Sp13D07.Astrocyte"
# # [7] "Sp13D08.OPC_Astroc" "Sp13D09.Astrocyte_" "Sp13D10.Oligo"
# # [10] "Sp13D11.Habenula"   "Sp13D12.Oligo"      "Sp13D13.Oligo"
# 
# level_map <- setNames(new_bs_names, old_bs_names)
# # Match original BayesSpace labels in sce to new names
# sce_pseudo$BayesSpace <- level_map[as.character(sce_pseudo$BayesSpace)]
# levels(sce_pseudo$BayesSpace)
# [1] "Sp11D01.Oligo"      "Sp11D02.Microglia"  "Sp11D03.Endo"      
# [4] "Sp11D04.Astrocyte"  "Sp11D05.Excit.Thal" "Sp11D06.Astrocyte" 
# [7] "Sp11D07.Astrocyte"  "Sp11D08.Astrocyte"  "Sp11D09.Oligo"     
# [10] "Sp11D10.Habenula"   "Sp11D11.Oligo"

## To avoid having to change parameters later on
sce_pseudo$registration_variable <- sce_pseudo$BayesSpace
sce_pseudo$registration_sample_id <- sce_pseudo$sample_id

## Drop unused levels
## - avoid error when some levels may not actually be present in the data, making the model non-identifiable in registration_model
table(sce_pseudo$registration_variable)
level_counts <- table(sce_pseudo$registration_variable)
if (any(level_counts == 0)) {
    sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)
    cat("After dropping unused levels:\n")
    print(table(sce_pseudo$registration_variable))
}

#covars <- c("sample_id")   # add sex, age when we have more than 2 classes
#covars <- NULL
covars <- c("brain_id")
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "all"
#colData(sce_pseudo)

## Taken from spatialLIBD::registration_wrapper()
## https://github.com/LieberInstitute/spatialLIBD/blob/master/R/registration_wrapper.R
registration_mod <-
  registration_model(sce_pseudo, covars = covars)

# head(registration_mod)
# colnames(registration_mod)
# rownames(registration_mod)
# colData(sce_pseudo)

message("Registration model done!")

block_cor <-
  registration_block_cor(sce_pseudo, registration_model = registration_mod)

results_enrichment <-
  registration_stats_enrichment(
    sce_pseudo,
    block_cor = block_cor,
    covars = covars,
    gene_ensembl = gene_ensembl,
    gene_name = gene_name
  )
results_pairwise <-
  registration_stats_pairwise(
    sce_pseudo,
    registration_model = registration_mod,
    block_cor = block_cor,
    gene_ensembl = gene_ensembl,
    gene_name = gene_name
  )
if (k >= 3) {
  results_anova <-
    registration_stats_anova(
      sce_pseudo,
      block_cor = block_cor,
      covars = covars,
      gene_ensembl = gene_ensembl,
      gene_name = gene_name,
      suffix = suffix
    )
} else {
  results_anova <- NULL
}

modeling_results <- list(
  "anova" = results_anova,
  "enrichment" = results_enrichment,
  "pairwise" = results_pairwise
)

## Save the final results
dir_rdata <- here("processed-data","05_brain_area_differential_expression", "modeling_results_BS")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
save(
  modeling_results,
  file = file.path(
    dir_rdata,
    paste0("modeling_results_BayesSpace_k", k_nice, ".Rdata")
  )
)

message(' Modeling pseudo-bulk bayes space results completed!')

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

## Reproducibility information
# > Sys.time()
# [1] "2025-05-26 11:58:10 EDT"
# > proc.time()
# user   system  elapsed 
# 21.375    1.103 1069.358 
# > options(width = 120)
# > session_info()
# 20 (R 4.4.2)
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.4.0)
# benchmarkme            1.0.8     2022-06-12 [2] CRAN (R 4.4.0)
# benchmarkmeData        1.0.4     2020-04-23 [2] CRAN (R 4.4.0)
# Biobase              * 2.66.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocFileCache          2.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocGenerics         * 0.52.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocIO                 1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocManager            1.30.25   2024-08-28 [2] CRAN (R 4.4.1)
# BiocNeighbors          2.0.1     2024-11-28 [2] Bioconductor 3.20 (R 4.4.2)
# BiocParallel           1.40.2    2025-04-10 [2] Bioconductor
# BiocSingular           1.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BiocVersion            3.20.0    2024-05-01 [2] Bioconductor 3.20 (R 4.4.0)
# Biostrings             2.74.1    2024-12-16 [2] Bioconductor 3.20 (R 4.4.2)
# bit                    4.6.0     2025-03-06 [2] CRAN (R 4.4.3)
# bit64                  4.6.0-1   2025-01-16 [2] CRAN (R 4.4.2)
# bitops                 1.0-9     2024-10-03 [2] CRAN (R 4.4.1)
# blob                   1.2.4     2023-03-17 [2] CRAN (R 4.4.0)
# bslib                  0.9.0     2025-01-30 [2] CRAN (R 4.4.2)
# cachem                 1.1.0     2024-05-16 [2] CRAN (R 4.4.0)
# circlize               0.4.16    2024-02-20 [2] CRAN (R 4.4.0)
# cli                    3.6.5     2025-04-23 [2] CRAN (R 4.4.3)
# clue                   0.3-66    2024-11-13 [2] CRAN (R 4.4.2)
# cluster                2.1.8     2024-12-11 [3] CRAN (R 4.4.3)
# codetools              0.2-20    2024-03-31 [3] CRAN (R 4.4.3)
# colorspace             2.1-1     2024-07-26 [2] CRAN (R 4.4.1)
# ComplexHeatmap         2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# config                 0.3.2     2023-08-30 [2] CRAN (R 4.4.0)
# cowplot                1.1.3     2024-01-22 [2] CRAN (R 4.4.0)
# crayon                 1.5.3     2024-06-20 [2] CRAN (R 4.4.1)
# curl                   6.2.2     2025-03-24 [2] CRAN (R 4.4.3)
# data.table             1.17.2    2025-05-12 [2] CRAN (R 4.4.3)
# DBI                    1.2.3     2024-06-02 [2] CRAN (R 4.4.0)
# dbplyr                 2.5.0     2024-03-19 [2] CRAN (R 4.4.0)
# DelayedArray           0.32.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# dichromat              2.0-0.1   2022-05-02 [2] CRAN (R 4.4.0)
# digest                 0.6.37    2024-08-19 [2] CRAN (R 4.4.1)
# doParallel             1.0.17    2022-02-07 [2] CRAN (R 4.4.0)
# dplyr                  1.1.4     2023-11-17 [2] CRAN (R 4.4.0)
# DT                     0.33      2024-04-04 [2] CRAN (R 4.4.0)
# edgeR                  4.4.2     2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
# ExperimentHub          2.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# farver                 2.1.2     2024-05-13 [2] CRAN (R 4.4.0)
# fastmap                1.2.0     2024-05-15 [2] CRAN (R 4.4.0)
# filelock               1.0.3     2023-12-11 [2] CRAN (R 4.4.0)
# foreach                1.5.2     2022-02-02 [2] CRAN (R 4.4.0)
# generics               0.1.4     2025-05-09 [2] CRAN (R 4.4.3)
# GenomeInfoDb         * 1.42.3    2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
# GenomeInfoDbData       1.2.13    2024-10-01 [2] Bioconductor
# GenomicAlignments      1.42.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# GenomicRanges        * 1.58.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# GetoptLong             1.0.5     2020-12-15 [2] CRAN (R 4.4.0)
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.4.0)
# ggplot2                3.5.2     2025-04-09 [2] CRAN (R 4.4.3)
# ggrepel                0.9.6     2024-09-07 [2] CRAN (R 4.4.1)
# GlobalOptions          0.1.2     2020-06-10 [2] CRAN (R 4.4.0)
# glue                   1.8.0     2024-09-30 [2] CRAN (R 4.4.1)
# golem                  0.5.1     2024-08-27 [2] CRAN (R 4.4.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.4.0)
# gtable                 0.3.6     2024-10-25 [2] CRAN (R 4.4.2)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.4.0)
# htmltools              0.5.8.1   2024-04-04 [2] CRAN (R 4.4.0)
# htmlwidgets            1.6.4     2023-12-06 [2] CRAN (R 4.4.0)
# httpuv                 1.6.16    2025-04-16 [2] CRAN (R 4.4.3)
# httr                   1.4.7     2023-08-15 [2] CRAN (R 4.4.0)
# IRanges              * 2.40.1    2024-12-05 [2] Bioconductor 3.20 (R 4.4.2)
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.4.0)
# iterators              1.0.14    2022-02-05 [2] CRAN (R 4.4.0)
# jquerylib              0.1.4     2021-04-26 [2] CRAN (R 4.4.0)
# jsonlite               2.0.0     2025-03-27 [2] CRAN (R 4.4.3)
# KEGGREST               1.46.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# later                  1.4.2     2025-04-08 [2] CRAN (R 4.4.3)
# lattice                0.22-6    2024-03-20 [3] CRAN (R 4.4.3)
# lazyeval               0.2.2     2019-03-15 [2] CRAN (R 4.4.0)
# lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.4.0)
# limma                  3.62.2    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# locfit                 1.5-9.12  2025-03-05 [2] CRAN (R 4.4.3)
# magick                 2.8.6     2025-03-23 [2] CRAN (R 4.4.3)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.4.0)
# Matrix                 1.7-2     2025-01-23 [3] CRAN (R 4.4.3)
# MatrixGenerics       * 1.18.1    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# matrixStats          * 1.5.0     2025-01-07 [2] CRAN (R 4.4.2)
# memoise                2.0.1     2021-11-26 [2] CRAN (R 4.4.0)
# mime                   0.13      2025-03-17 [2] CRAN (R 4.4.3)
# paletteer              1.6.0     2024-01-21 [2] CRAN (R 4.4.0)
# pillar                 1.10.2    2025-04-05 [2] CRAN (R 4.4.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.4.0)
# plotly                 4.10.4    2024-01-13 [2] CRAN (R 4.4.0)
# png                    0.1-8     2022-11-29 [2] CRAN (R 4.4.0)
# promises               1.3.2     2024-11-28 [2] CRAN (R 4.4.2)
# purrr                  1.0.4     2025-02-05 [2] CRAN (R 4.4.2)
# R6                     2.6.1     2025-02-15 [2] CRAN (R 4.4.2)
# rappdirs               0.3.3     2021-01-31 [2] CRAN (R 4.4.0)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.4.0)
# Rcpp                   1.0.14    2025-01-12 [2] CRAN (R 4.4.2)
# RCurl                  1.98-1.17 2025-03-22 [2] CRAN (R 4.4.3)
# rematch2               2.1.2     2020-05-01 [2] CRAN (R 4.4.0)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.4.0)
# rjson                  0.2.23    2024-09-16 [2] CRAN (R 4.4.1)
# rlang                  1.1.6     2025-04-11 [2] CRAN (R 4.4.3)
# rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.4.0)
# Rsamtools              2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# RSQLite                2.3.11    2025-05-04 [2] CRAN (R 4.4.3)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.4.0)
# rtracklayer            1.66.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# S4Arrays               1.6.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# S4Vectors            * 0.44.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# sass                   0.4.10    2025-04-11 [2] CRAN (R 4.4.3)
# ScaledMatrix           1.14.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# scales                 1.4.0     2025-04-24 [2] CRAN (R 4.4.3)
# scater                 1.34.1    2025-03-03 [2] Bioconductor 3.20 (R 4.4.3)
# scuttle                1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# sessioninfo          * 1.2.3     2025-02-05 [2] CRAN (R 4.4.2)
# shape                  1.4.6.1   2024-02-23 [2] CRAN (R 4.4.0)
# shiny                  1.10.0    2024-12-14 [2] CRAN (R 4.4.2)
# shinyWidgets           0.9.0     2025-02-21 [2] CRAN (R 4.4.3)
# SingleCellExperiment * 1.28.1    2024-11-10 [2] Bioconductor 3.20 (R 4.4.2)
# SparseArray            1.6.2     2025-02-20 [2] Bioconductor 3.20 (R 4.4.3)
# SpatialExperiment    * 1.16.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# spatialLIBD          * 1.21.5    2025-05-16 [1] Github (LieberInstitute/spatialLIBD@aff00db)
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.4.0)
# SummarizedExperiment * 1.36.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.4.0)
# tidyr                  1.3.1     2024-01-24 [2] CRAN (R 4.4.0)
# tidyselect             1.2.1     2024-03-11 [2] CRAN (R 4.4.0)
# UCSC.utils             1.2.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# vctrs                  0.6.5     2023-12-01 [2] CRAN (R 4.4.0)
# vipor                  0.4.7     2023-12-18 [2] CRAN (R 4.4.0)
# viridis                0.6.5     2024-01-29 [2] CRAN (R 4.4.0)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.4.0)
# XML                    3.99-0.18 2025-01-01 [2] CRAN (R 4.4.2)
# xtable                 1.8-4     2019-04-21 [2] CRAN (R 4.4.0)
# XVector                0.46.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# yaml                   2.3.10    2024-07-26 [2] CRAN (R 4.4.1)
# zlibbioc               1.52.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# 
# [1] /users/csoto/R/4.4.x
# [2] /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/library

