# copied from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/07_layer_differential_expression/03_model_BayesSpace.R

# library(slurmjobs)
# slurmjobs::job_single('03_model_BayesSpace', create_shell = TRUE, memory = '30G', command = "03_model_BayesSpace.R")

# To submit the job use: sbatch 03_model_BayesSpace.sh

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## For testing
if (is.na(k)) {
  k <- 9
}


library("here")
library("sessioninfo")
library("spatialLIBD")

## output directory
dir_rdata <- here("processed-data","05_layer_differential_expression")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_rdata)) ## Check that it was created successfully
dir_plots <- here("plots","05_layer_differential_expression")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(dir_plots))

## load spe_pseudo data
sce_pseudo <-
  readRDS(
    file.path(
      dir_rdata,
      paste0("sce_pseudo_BayesSpace_k", sprintf("%02d", k), ".rds")
    )
  )
sce_pseudo
colnames(colData(sce_pseudo))

## To avoid having to change parameters later on
sce_pseudo$registration_variable <- sce_pseudo$BayesSpace
sce_pseudo$registration_sample_id <- sce_pseudo$sample_id

## Set arguments used in spatialLIBD::registration_wrapper()
covars <- c("brain_id")   # add sex, age when we have more than 2 classes
# I set covars = Null to avoid error when running "registration_block_cor", it gets next warming and other error when running "registration_stats_pairwise"
# Warning message:
#     In limma::duplicateCorrelation(logcounts(sce_pseudo), registration_model,  :
#     Block factor already encoded in the design matrix: setting intrablock correlation to zero.

#covars <- NULL
gene_ensembl <- "gene_id"
gene_name <- "gene_name"
suffix <- "all"

#colData(sce_pseudo)

## Taken from spatialLIBD::registration_wrapper()
## https://github.com/LieberInstitute/spatialLIBD/blob/master/R/registration_wrapper.R
table(sce_pseudo[["brain_id"]])

# defines the statistical model that will be used for computing the block correlation as well as pairwise statistics
# The output of model.matrix() which you can inspect to verify that your sample-level covariates are being properly modeled
registration_mod <-
  registration_model(sce_pseudo, covars = covars)

head(registration_mod)
colnames(registration_mod)
rownames(registration_mod)
colData(sce_pseudo)

# computes the block correlation at sample ID level
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
dir_rdata <- here("processed-data","05_layer_differential_expression", "modeling_results_BS")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
save(
  modeling_results,
  file = file.path(
    dir_rdata,
    paste0("modeling_results_BayesSpace_k", sprintf("%02d", k), ".Rdata")
  )
)

message(' Modeling pseudo-bulk bayes space results completed!')

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

## ==============================================================================
# > ## Reproducibility information
#   > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-06-28 10:42:40 EDT"
# > proc.time()
# user   system  elapsed
# 49.177    2.336 1223.126
# > options(width = 120)
# > session_info()
# [39m CRAN (R 4.3.2)
#   benchmarkme              1.0.8       2022-06-12 [2] CRAN (R 4.3.2)
#   benchmarkmeData          1.0.4       2020-04-23 [2] CRAN (R 4.3.2)
#   Biobase                * 2.62.0      2023-10-24 [2] Bioconductor
#   BiocFileCache            2.10.1      2023-10-26 [2] Bioconductor
#   BiocGenerics           * 0.48.1      2023-11-01 [2] Bioconductor
#   BiocIO                   1.12.0      2023-10-24 [2] Bioconductor
#   BiocManager              1.30.22     2023-08-08 [2] CRAN (R 4.3.2)
#   BiocNeighbors            1.20.2      2024-01-07 [2] Bioconductor 3.18 (R 4.3.2)
#   BiocParallel             1.36.0      2023-10-24 [2] Bioconductor
#   BiocSingular             1.18.0      2023-10-24 [2] Bioconductor
#   BiocVersion              3.18.1      2023-11-15 [2] Bioconductor
#   Biostrings               2.70.2      2024-01-28 [2] Bioconductor 3.18 (R 4.3.2)
#   bit                      4.0.5       2022-11-15 [2] CRAN (R 4.3.2)
#   bit64                    4.0.5       2020-08-30 [2] CRAN (R 4.3.2)
#   bitops                   1.0-7       2021-04-24 [2] CRAN (R 4.3.2)
#   blob                     1.2.4       2023-03-17 [2] CRAN (R 4.3.2)
#   bslib                    0.6.1       2023-11-28 [2] CRAN (R 4.3.2)
#   cachem                   1.0.8       2023-05-01 [2] CRAN (R 4.3.2)
#   cli                      3.6.2       2023-12-11 [2] CRAN (R 4.3.2)
#   codetools                0.2-19      2023-02-01 [3] CRAN (R 4.3.2)
#   colorspace               2.1-0       2023-01-23 [2] CRAN (R 4.3.2)
#   config                   0.3.2       2023-08-30 [2] CRAN (R 4.3.2)
#   cowplot                  1.1.3       2024-01-22 [2] CRAN (R 4.3.2)
#   crayon                   1.5.2       2022-09-29 [2] CRAN (R 4.3.2)
#   curl                     5.2.0       2023-12-08 [2] CRAN (R 4.3.2)
#   data.table               1.15.0      2024-01-30 [2] CRAN (R 4.3.2)
#   DBI                      1.2.1       2024-01-12 [2] CRAN (R 4.3.2)
#   dbplyr                   2.4.0       2023-10-26 [2] CRAN (R 4.3.2)
#   DelayedArray             0.28.0      2023-10-24 [2] Bioconductor
#   DelayedMatrixStats       1.24.0      2023-10-24 [2] Bioconductor
#   digest                   0.6.34      2024-01-11 [2] CRAN (R 4.3.2)
#   doParallel               1.0.17      2022-02-07 [2] CRAN (R 4.3.2)
#   dotCall64                1.1-1       2023-11-28 [2] CRAN (R 4.3.2)
#   dplyr                    1.1.4       2023-11-17 [2] CRAN (R 4.3.2)
#   DT                       0.31        2023-12-09 [2] CRAN (R 4.3.2)
#   edgeR                    4.0.14      2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
#   ellipsis                 0.3.2       2021-04-29 [2] CRAN (R 4.3.2)
#   ExperimentHub            2.10.0      2023-10-24 [2] Bioconductor
#   fansi                    1.0.6       2023-12-08 [2] CRAN (R 4.3.2)
#   farver                   2.1.1       2022-07-06 [2] CRAN (R 4.3.2)
#   fastmap                  1.1.1       2023-02-24 [2] CRAN (R 4.3.2)
#   fields                   15.2        2023-08-17 [2] CRAN (R 4.3.2)
#   filelock                 1.0.3       2023-12-11 [2] CRAN (R 4.3.2)
#   foreach                  1.5.2       2022-02-02 [2] CRAN (R 4.3.2)
#   generics                 0.1.3       2022-07-05 [2] CRAN (R 4.3.2)
#   GenomeInfoDb           * 1.38.5      2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
#   GenomeInfoDbData         1.2.11      2024-02-09 [2] Bioconductor
#   GenomicAlignments        1.38.2      2024-01-16 [2] Bioconductor 3.18 (R 4.3.2)
#   GenomicRanges          * 1.54.1      2023-10-29 [2] Bioconductor
#   ggbeeswarm               0.7.2       2023-04-29 [2] CRAN (R 4.3.2)
#   ggplot2                * 3.4.4       2023-10-12 [2] CRAN (R 4.3.2)
#   ggrepel                  0.9.5       2024-01-10 [2] CRAN (R 4.3.2)
#   glue                     1.7.0       2024-01-09 [2] CRAN (R 4.3.2)
#   golem                    0.4.1       2023-06-05 [2] CRAN (R 4.3.2)
#   gridExtra                2.3         2017-09-09 [2] CRAN (R 4.3.2)
#   gtable                   0.3.4       2023-08-21 [2] CRAN (R 4.3.2)
#   here                   * 1.0.1       2020-12-13 [2] CRAN (R 4.3.2)
#   htmltools                0.5.7       2023-11-03 [2] CRAN (R 4.3.2)
#   htmlwidgets              1.6.4       2023-12-06 [2] CRAN (R 4.3.2)
#   httpuv                   1.6.14      2024-01-26 [2] CRAN (R 4.3.2)
#   httr                     1.4.7       2023-08-15 [2] CRAN (R 4.3.2)
#   interactiveDisplayBase   1.40.0      2023-10-24 [2] Bioconductor
#   IRanges                * 2.36.0      2023-10-24 [2] Bioconductor
#   irlba                    2.3.5.1     2022-10-03 [2] CRAN (R 4.3.2)
#   iterators                1.0.14      2022-02-05 [2] CRAN (R 4.3.2)
#   jquerylib                0.1.4       2021-04-26 [2] CRAN (R 4.3.2)
#   jsonlite                 1.8.8       2023-12-04 [2] CRAN (R 4.3.2)
#   KEGGREST                 1.42.0      2023-10-24 [2] Bioconductor
#   labeling                 0.4.3       2023-08-29 [2] CRAN (R 4.3.2)
#   later                    1.3.2       2023-12-06 [2] CRAN (R 4.3.2)
#   lattice                  0.22-5      2023-10-24 [3] CRAN (R 4.3.2)
#   lazyeval                 0.2.2       2019-03-15 [2] CRAN (R 4.3.2)
#   lifecycle                1.0.4       2023-11-07 [2] CRAN (R 4.3.2)
#   limma                    3.58.1      2023-10-31 [2] Bioconductor
#   locfit                   1.5-9.8     2023-06-11 [2] CRAN (R 4.3.2)
#   magick                   2.8.2       2023-12-20 [2] CRAN (R 4.3.2)
#   magrittr                 2.0.3       2022-03-30 [2] CRAN (R 4.3.2)
#   maps                     3.4.2       2023-12-15 [2] CRAN (R 4.3.2)
#   Matrix                   1.6-5       2024-01-11 [3] CRAN (R 4.3.2)
#   MatrixGenerics         * 1.14.0      2023-10-24 [2] Bioconductor
#   matrixStats            * 1.2.0       2023-12-11 [2] CRAN (R 4.3.2)
#   memoise                  2.0.1       2021-11-26 [2] CRAN (R 4.3.2)
#   mime                     0.12        2021-09-28 [2] CRAN (R 4.3.2)
#   munsell                  0.5.0       2018-06-12 [2] CRAN (R 4.3.2)
#   paletteer                1.6.0       2024-01-21 [2] CRAN (R 4.3.2)
#   pillar                   1.9.0       2023-03-22 [2] CRAN (R 4.3.2)
#   pkgconfig                2.0.3       2019-09-22 [2] CRAN (R 4.3.2)
#   plotly                   4.10.4      2024-01-13 [2] CRAN (R 4.3.2)
#   png                      0.1-8       2022-11-29 [2] CRAN (R 4.3.2)
#   Polychrome             * 1.5.1       2022-04-29 [1] R-Forge (R 4.3.2)
#   promises                 1.2.1       2023-08-10 [2] CRAN (R 4.3.2)
#   purrr                    1.0.2       2023-08-10 [2] CRAN (R 4.3.2)
#   R6                       2.5.1       2021-08-19 [2] CRAN (R 4.3.2)
#   rappdirs                 0.3.3       2021-01-31 [2] CRAN (R 4.3.2)
#   RColorBrewer             1.1-3       2022-04-03 [2] CRAN (R 4.3.2)
#   Rcpp                     1.0.12      2024-01-09 [2] CRAN (R 4.3.2)
#   RCurl                    1.98-1.14   2024-01-09 [2] CRAN (R 4.3.2)
#   rematch2                 2.1.2       2020-05-01 [2] CRAN (R 4.3.2)
#   restfulr                 0.0.15      2022-06-16 [2] CRAN (R 4.3.2)
#   rjson                    0.2.21      2022-01-09 [2] CRAN (R 4.3.2)
#   rlang                    1.1.3       2024-01-10 [2] CRAN (R 4.3.2)
#   rprojroot                2.0.4       2023-11-05 [2] CRAN (R 4.3.2)
#   Rsamtools                2.18.0      2023-10-24 [2] Bioconductor
#   RSQLite                  2.3.5       2024-01-21 [2] CRAN (R 4.3.2)
#   rsvd                     1.0.5       2021-04-16 [2] CRAN (R 4.3.2)
#   rtracklayer              1.62.0      2023-10-24 [2] Bioconductor
#   S4Arrays                 1.2.0       2023-10-24 [2] Bioconductor
#   S4Vectors              * 0.40.2      2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
#   sass                     0.4.8       2023-12-06 [2] CRAN (R 4.3.2)
#   ScaledMatrix             1.10.0      2023-10-24 [2] Bioconductor
#   scales                   1.3.0       2023-11-28 [2] CRAN (R 4.3.2)
#   scater                 * 1.30.1      2023-11-16 [2] Bioconductor
#   scatterplot3d            0.3-44      2023-05-05 [1] R-Forge (R 4.3.2)
#   scuttle                * 1.12.0      2023-10-24 [2] Bioconductor
#   sessioninfo            * 1.2.2       2021-12-06 [2] CRAN (R 4.3.2)
#   shiny                    1.8.0       2023-11-17 [2] CRAN (R 4.3.2)
#   shinyWidgets             0.8.1       2024-01-10 [2] CRAN (R 4.3.2)
#   SingleCellExperiment   * 1.24.0      2023-10-24 [2] Bioconductor
#   spam                     2.10-0      2023-10-23 [2] CRAN (R 4.3.2)
#   SparseArray              1.2.3       2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
#   sparseMatrixStats        1.14.0      2023-10-24 [2] Bioconductor
#   SpatialExperiment      * 1.12.0      2023-10-24 [2] Bioconductor
#   spatialLIBD            * 1.15.4      2024-05-01 [1] Github (LieberInstitute/spatialLIBD@77a5303)
#   statmod                  1.5.0       2023-01-06 [2] CRAN (R 4.3.2)
#   SummarizedExperiment   * 1.32.0      2023-10-24 [2] Bioconductor
#   tibble                   3.2.1       2023-03-20 [2] CRAN (R 4.3.2)
#   tidyr                    1.3.1       2024-01-24 [2] CRAN (R 4.3.2)
#   tidyselect               1.2.0       2022-10-10 [2] CRAN (R 4.3.2)
#   utf8                     1.2.4       2023-10-22 [2] CRAN (R 4.3.2)
#   vctrs                    0.6.5       2023-12-01 [2] CRAN (R 4.3.2)
#   vipor                    0.4.7       2023-12-18 [2] CRAN (R 4.3.2)
#   viridis                  0.6.5       2024-01-29 [2] CRAN (R 4.3.2)
#   viridisLite              0.4.2       2023-05-02 [2] CRAN (R 4.3.2)
