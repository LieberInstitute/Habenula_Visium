library("here")
library("spatialLIBD")
library("DeconvoBuddies")
library("tidyverse")
library("sessioninfo")


## output directory
dir_rdata <- here("processed-data", "04_harmony_BayesSpace")

## Load the data
spe <- readRDS(file.path(dir_rdata, "spe_harmony.rds"))

## Import BayesSpace clusters
spe <- cluster_import(spe,
    cluster_dir = file.path(dir_rdata, "clusters_BayesSpace"),
    prefix = ""
)

BayesSpace_vars <-
    colnames(colData(spe))[grep("^BayesSpace_harmony_k", colnames(colData(spe)))]


BayesSpace_stats_list <-
    lapply(BayesSpace_vars, function(BayesSpace_current) {
        message(Sys.time(), " processing ", BayesSpace_current)
        markers_1vALL <- findMarkers_1vAll(spe,
            assay_name = "logcounts",
            cellType_col = "BayesSpace_harmony_k08",
            mod = NULL
        )

        markers_1vALL %>%
            reframe(
                t_stat = std.logFC,
                gene = gene
            ) %>%
            mutate(cellType.target = paste0("BayesSpace_harmony_k08", "_", cellType.target)) %>%
            pivot_wider(
                names_from = cellType.target,
                values_from = t_stat
            ) %>%
            as.data.frame()
    })
names(BayesSpace_stats_list) <- BayesSpace_vars

## Explore some results
dim(BayesSpace_stats_list$BayesSpace_harmony_k08)
head(BayesSpace_stats_list$BayesSpace_harmony_k08)

## This matches the format of 'stats' in
## https://research.libd.org/spatialLIBD/reference/layer_stat_cor.html

saveRDS(BayesSpace_stats_list,
    file = file.path(dir_rdata, "BayesSpace_stats_list.rds")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 (2023-10-31)
#  os       macOS Sonoma 14.3.1
#  system   aarch64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2024-03-01
#  rstudio  2023.12.1+402 Ocean Storm (desktop)
#  pandoc   3.1.12.1 @ /opt/homebrew/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package                * version     date (UTC) lib source
#  abind                    1.4-5       2016-07-21 [1] CRAN (R 4.3.0)
#  AnnotationDbi            1.64.1      2023-11-02 [1] Bioconductor
#  AnnotationHub            3.10.0      2023-10-26 [1] Bioconductor
#  attempt                  0.3.1       2020-05-03 [1] CRAN (R 4.3.0)
#  beachmat                 2.18.1      2024-02-17 [1] Bioconductor 3.18 (R 4.3.2)
#  beeswarm                 0.4.0       2021-06-01 [1] CRAN (R 4.3.0)
#  benchmarkme              1.0.8       2022-06-12 [1] CRAN (R 4.3.0)
#  benchmarkmeData          1.0.4       2020-04-23 [1] CRAN (R 4.3.0)
#  Biobase                * 2.62.0      2023-10-26 [1] Bioconductor
#  BiocFileCache            2.10.1      2023-10-26 [1] Bioconductor
#  BiocGenerics           * 0.48.1      2023-11-02 [1] Bioconductor
#  BiocIO                   1.12.0      2023-10-26 [1] Bioconductor
#  BiocManager              1.30.22     2023-08-08 [1] CRAN (R 4.3.0)
#  BiocNeighbors            1.20.2      2024-01-13 [1] Bioconductor 3.18 (R 4.3.2)
#  BiocParallel             1.36.0      2023-10-26 [1] Bioconductor
#  BiocSingular             1.18.0      2023-11-06 [1] Bioconductor
#  biocthis                 1.12.0      2023-10-26 [1] Bioconductor
#  BiocVersion              3.18.1      2023-11-18 [1] Bioconductor 3.18 (R 4.3.2)
#  Biostrings               2.70.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.2)
#  bit                      4.0.5       2022-11-15 [1] CRAN (R 4.3.0)
#  bit64                    4.0.5       2020-08-30 [1] CRAN (R 4.3.0)
#  bitops                   1.0-7       2021-04-24 [1] CRAN (R 4.3.0)
#  blob                     1.2.4       2023-03-17 [1] CRAN (R 4.3.0)
#  bluster                  1.12.0      2023-12-19 [1] Bioconductor 3.18 (R 4.3.2)
#  brio                     1.1.4       2023-12-10 [1] CRAN (R 4.3.1)
#  bslib                    0.6.1       2023-11-28 [1] CRAN (R 4.3.1)
#  cachem                   1.0.8       2023-05-01 [1] CRAN (R 4.3.0)
#  cli                      3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
#  cluster                  2.1.6       2023-12-01 [1] CRAN (R 4.3.1)
#  codetools                0.2-19      2023-02-01 [1] CRAN (R 4.3.2)
#  colorout                 1.3-0.2     2024-02-27 [1] Github (jalvesaq/colorout@c6113a2)
#  colorspace               2.1-0       2023-01-23 [1] CRAN (R 4.3.0)
#  config                   0.3.2       2023-08-30 [1] CRAN (R 4.3.0)
#  cowplot                  1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
#  crayon                   1.5.2       2022-09-29 [1] CRAN (R 4.3.0)
#  curl                     5.2.0       2023-12-08 [1] CRAN (R 4.3.1)
#  data.table               1.15.0      2024-01-30 [1] CRAN (R 4.3.1)
#  DBI                      1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
#  dbplyr                   2.4.0       2023-10-26 [1] CRAN (R 4.3.1)
#  DeconvoBuddies         * 0.99.0      2024-03-01 [1] Github (LieberInstitute/DeconvoBuddies@8efdf5b)
#  DelayedArray             0.28.0      2023-11-06 [1] Bioconductor
#  DelayedMatrixStats       1.24.0      2023-11-06 [1] Bioconductor
#  devtools               * 2.4.5       2022-10-11 [1] CRAN (R 4.3.0)
#  digest                   0.6.34      2024-01-11 [1] CRAN (R 4.3.1)
#  doParallel               1.0.17      2022-02-07 [1] CRAN (R 4.3.0)
#  dotCall64                1.1-1       2023-11-28 [1] CRAN (R 4.3.1)
#  dplyr                  * 1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
#  dqrng                    0.3.2       2023-11-29 [1] CRAN (R 4.3.1)
#  DT                       0.32        2024-02-19 [1] CRAN (R 4.3.1)
#  edgeR                    4.0.16      2024-02-20 [1] Bioconductor 3.18 (R 4.3.2)
#  ellipsis                 0.3.2       2021-04-29 [1] CRAN (R 4.3.0)
#  ExperimentHub            2.10.0      2023-10-26 [1] Bioconductor
#  fansi                    1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
#  fastmap                  1.1.1       2023-02-24 [1] CRAN (R 4.3.0)
#  fields                   15.2        2023-08-17 [1] CRAN (R 4.3.0)
#  filelock                 1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
#  forcats                * 1.0.0       2023-01-29 [1] CRAN (R 4.3.0)
#  foreach                  1.5.2       2022-02-02 [1] CRAN (R 4.3.0)
#  fs                       1.6.3       2023-07-20 [1] CRAN (R 4.3.0)
#  generics                 0.1.3       2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb           * 1.38.6      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
#  GenomeInfoDbData         1.2.11      2024-02-27 [1] Bioconductor
#  GenomicAlignments        1.38.2      2024-01-20 [1] Bioconductor 3.18 (R 4.3.2)
#  GenomicRanges          * 1.54.1      2023-10-30 [1] Bioconductor
#  ggbeeswarm               0.7.2       2023-04-29 [1] CRAN (R 4.3.0)
#  ggplot2                * 3.5.0       2024-02-23 [1] CRAN (R 4.3.1)
#  ggrepel                  0.9.5       2024-01-10 [1] CRAN (R 4.3.1)
#  glue                     1.7.0       2024-01-09 [1] CRAN (R 4.3.1)
#  golem                    0.4.1       2023-06-05 [1] CRAN (R 4.3.0)
#  gridExtra                2.3         2017-09-09 [1] CRAN (R 4.3.0)
#  gtable                   0.3.4       2023-08-21 [1] CRAN (R 4.3.0)
#  here                   * 1.0.1       2020-12-13 [1] CRAN (R 4.3.0)
#  hms                      1.1.3       2023-03-21 [1] CRAN (R 4.3.0)
#  htmltools                0.5.7       2023-11-03 [1] CRAN (R 4.3.1)
#  htmlwidgets              1.6.4       2023-12-06 [1] CRAN (R 4.3.1)
#  httpuv                   1.6.14      2024-01-26 [1] CRAN (R 4.3.1)
#  httr                     1.4.7       2023-08-15 [1] CRAN (R 4.3.0)
#  igraph                   2.0.2       2024-02-17 [1] CRAN (R 4.3.1)
#  interactiveDisplayBase   1.40.0      2023-10-26 [1] Bioconductor
#  IRanges                * 2.36.0      2023-10-26 [1] Bioconductor
#  irlba                    2.3.5.1     2022-10-03 [1] CRAN (R 4.3.2)
#  iterators                1.0.14      2022-02-05 [1] CRAN (R 4.3.0)
#  jquerylib                0.1.4       2021-04-26 [1] CRAN (R 4.3.0)
#  jsonlite                 1.8.8       2023-12-04 [1] CRAN (R 4.3.1)
#  KEGGREST                 1.42.0      2023-10-26 [1] Bioconductor
#  later                    1.3.2       2023-12-06 [1] CRAN (R 4.3.1)
#  lattice                  0.22-5      2023-10-24 [1] CRAN (R 4.3.1)
#  lazyeval                 0.2.2       2019-03-15 [1] CRAN (R 4.3.0)
#  lifecycle                1.0.4       2023-11-07 [1] CRAN (R 4.3.1)
#  limma                    3.58.1      2023-11-02 [1] Bioconductor
#  lobstr                   1.1.2       2022-06-22 [1] CRAN (R 4.3.0)
#  locfit                   1.5-9.8     2023-06-11 [1] CRAN (R 4.3.0)
#  lubridate              * 1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
#  magick                   2.8.3       2024-02-18 [1] CRAN (R 4.3.1)
#  magrittr                 2.0.3       2022-03-30 [1] CRAN (R 4.3.0)
#  maps                     3.4.2       2023-12-15 [1] CRAN (R 4.3.1)
#  Matrix                   1.6-5       2024-01-11 [1] CRAN (R 4.3.2)
#  MatrixGenerics         * 1.14.0      2023-10-26 [1] Bioconductor
#  matrixStats            * 1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
#  memoise                  2.0.1       2021-11-26 [1] CRAN (R 4.3.0)
#  metapod                  1.10.1      2023-12-23 [1] Bioconductor 3.18 (R 4.3.2)
#  mime                     0.12        2021-09-28 [1] CRAN (R 4.3.0)
#  miniUI                   0.1.1.1     2018-05-18 [1] CRAN (R 4.3.0)
#  munsell                  0.5.0       2018-06-12 [1] CRAN (R 4.3.0)
#  paletteer                1.6.0       2024-01-21 [1] CRAN (R 4.3.1)
#  pillar                   1.9.0       2023-03-22 [1] CRAN (R 4.3.0)
#  pkgbuild                 1.4.3       2023-12-10 [1] CRAN (R 4.3.1)
#  pkgconfig                2.0.3       2019-09-22 [1] CRAN (R 4.3.0)
#  pkgload                  1.3.4       2024-01-16 [1] CRAN (R 4.3.1)
#  plotly                   4.10.4      2024-01-13 [1] CRAN (R 4.3.1)
#  png                      0.1-8       2022-11-29 [1] CRAN (R 4.3.0)
#  profvis                  0.3.8       2023-05-02 [1] CRAN (R 4.3.0)
#  promises                 1.2.1       2023-08-10 [1] CRAN (R 4.3.0)
#  prompt                   1.0.2.9000  2024-02-27 [1] Github (gaborcsardi/prompt@17bd0e1)
#  purrr                  * 1.0.2       2023-08-10 [1] CRAN (R 4.3.0)
#  R.cache                  0.16.0      2022-07-21 [1] CRAN (R 4.3.0)
#  R.methodsS3              1.8.2       2022-06-13 [1] CRAN (R 4.3.0)
#  R.oo                     1.26.0      2024-01-24 [1] CRAN (R 4.3.1)
#  R.utils                  2.12.3      2023-11-18 [1] CRAN (R 4.3.1)
#  R6                       2.5.1       2021-08-19 [1] CRAN (R 4.3.0)
#  rafalib                  1.0.0       2015-08-09 [1] CRAN (R 4.3.0)
#  rappdirs                 0.3.3       2021-01-31 [1] CRAN (R 4.3.0)
#  RColorBrewer             1.1-3       2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                     1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
#  RCurl                    1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
#  readr                  * 2.1.5       2024-01-10 [1] CRAN (R 4.3.1)
#  rematch2                 2.1.2       2020-05-01 [1] CRAN (R 4.3.0)
#  remotes                  2.4.2.1     2023-07-18 [1] CRAN (R 4.3.0)
#  restfulr                 0.0.15      2022-06-16 [1] CRAN (R 4.3.0)
#  rjson                    0.2.21      2022-01-09 [1] CRAN (R 4.3.0)
#  rlang                    1.1.3       2024-01-10 [1] CRAN (R 4.3.1)
#  rprojroot                2.0.4       2023-11-05 [1] CRAN (R 4.3.1)
#  Rsamtools                2.18.0      2023-10-26 [1] Bioconductor
#  RSQLite                  2.3.5       2024-01-21 [1] CRAN (R 4.3.1)
#  rsthemes                 0.4.0       2024-02-27 [1] Github (gadenbuie/rsthemes@34a55a4)
#  rstudioapi               0.15.0      2023-07-07 [1] CRAN (R 4.3.0)
#  rsvd                     1.0.5       2021-04-16 [1] CRAN (R 4.3.0)
#  rtracklayer              1.62.0      2023-10-26 [1] Bioconductor
#  S4Arrays                 1.2.0       2023-10-26 [1] Bioconductor
#  S4Vectors              * 0.40.2      2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
#  sass                     0.4.8.9000  2024-02-27 [1] Github (rstudio/sass@ae93a9a)
#  ScaledMatrix             1.10.0      2023-11-06 [1] Bioconductor
#  scales                   1.3.0       2023-11-28 [1] CRAN (R 4.3.1)
#  scater                   1.30.1      2023-11-16 [1] Bioconductor
#  scran                    1.30.2      2024-01-23 [1] Bioconductor 3.18 (R 4.3.2)
#  scuttle                  1.12.0      2023-11-06 [1] Bioconductor
#  sessioninfo            * 1.2.2       2021-12-06 [1] CRAN (R 4.3.0)
#  shiny                    1.8.0       2023-11-17 [1] CRAN (R 4.3.1)
#  shinyWidgets             0.8.1       2024-01-10 [1] CRAN (R 4.3.1)
#  SingleCellExperiment   * 1.24.0      2023-11-06 [1] Bioconductor
#  spam                     2.10-0      2023-10-23 [1] CRAN (R 4.3.1)
#  SparseArray              1.2.4       2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
#  sparseMatrixStats        1.14.0      2023-10-26 [1] Bioconductor
#  SpatialExperiment      * 1.12.0      2023-10-26 [1] Bioconductor
#  spatialLIBD            * 1.14.1      2023-11-30 [1] Bioconductor 3.18 (R 4.3.2)
#  statmod                  1.5.0       2023-01-06 [1] CRAN (R 4.3.0)
#  stringi                  1.8.3       2023-12-11 [1] CRAN (R 4.3.1)
#  stringr                * 1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
#  styler                   1.10.2      2023-08-29 [1] CRAN (R 4.3.0)
#  SummarizedExperiment   * 1.32.0      2023-11-06 [1] Bioconductor
#  suncalc                  0.5.1       2022-09-29 [1] CRAN (R 4.3.0)
#  testthat               * 3.2.1       2023-12-02 [1] CRAN (R 4.3.1)
#  tibble                 * 3.2.1       2023-03-20 [1] CRAN (R 4.3.0)
#  tidyr                  * 1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
#  tidyselect               1.2.0       2022-10-10 [1] CRAN (R 4.3.0)
#  tidyverse              * 2.0.0       2023-02-22 [1] CRAN (R 4.3.0)
#  timechange               0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
#  tzdb                     0.4.0       2023-05-12 [1] CRAN (R 4.3.0)
#  urlchecker               1.0.1       2021-11-30 [1] CRAN (R 4.3.0)
#  usethis                * 2.2.3       2024-02-19 [1] CRAN (R 4.3.1)
#  utf8                     1.2.4       2023-10-22 [1] CRAN (R 4.3.1)
#  vctrs                    0.6.5       2023-12-01 [1] CRAN (R 4.3.1)
#  vipor                    0.4.7       2023-12-18 [1] CRAN (R 4.3.1)
#  viridis                  0.6.5       2024-01-29 [1] CRAN (R 4.3.1)
#  viridisLite              0.4.2       2023-05-02 [1] CRAN (R 4.3.0)
#  withr                    3.0.0       2024-01-16 [1] CRAN (R 4.3.1)
#  XML                      3.99-0.16.1 2024-01-22 [1] CRAN (R 4.3.1)
#  xtable                   1.8-4       2019-04-21 [1] CRAN (R 4.3.0)
#  XVector                  0.42.0      2023-10-26 [1] Bioconductor
#  yaml                     2.3.8       2023-12-11 [1] CRAN (R 4.3.1)
#  zlibbioc                 1.48.0      2023-10-26 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
