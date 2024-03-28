library("here")
library("SingleCellExperiment")
library("spatialLIBD")
library("sessioninfo")

## Create output directories
dir_rdata <- here("processed-data", "05_snRNA-seq_model_stats")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load Habenula snRNA-seq data from the
## https://github.com/LieberInstitute/Habenula_Pilot
## project.

## If we were inside that project
## https://github.com/LieberInstitute/Habenula_Pilot/blob/85bfdf36505c829523358a98328a3e2b22bf2042/code/99_paper_figs/12_longer_Heatmap.R#L17-L18
##

habenula_pilot_path <- "~/Dropbox/Code/Habenula_Pilot"
if (!file.exists(habenula_pilot_path)) {
    habenula_pilot_path <-
        "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula"
}
input_sce_path <-
    file.path(
        habenula_pilot_path,
        "processed-data",
        "04_snRNA-seq",
        "sce_objects",
        "sce_final.Rdata"
    )
load(input_sce_path, verbose = TRUE)

table(sce_final$final_Annotations)
# Astrocyte       Endo Excit.Thal Inhib.Thal      LHb.1      LHb.2
#       538         38       1800       7612        201        266
#     LHb.3      LHb.4      LHb.5      LHb.6      LHb.7      MHb.1
#       134        477         83         39       1014        152
#     MHb.2      MHb.3  Microglia      Oligo        OPC
#       540         18        145       2178       1796

sce_final$final_Annotations_broad <-
    gsub("\\.[0-9]+", "", sce_final$final_Annotations)
table(sce_final$final_Annotations_broad, useNA = "ifany")
# Astrocyte       Endo Excit.Thal Inhib.Thal        LHb        MHb
#       538         38       1800       7612       2214        710
# Microglia      Oligo        OPC
#       145       2178       1796

registration_vars <-
    c("final_Annotations", "final_Annotations_broad")

enrichment_tstats <- lapply(registration_vars, function(current_var) {
    message(Sys.time(), " processing ", current_var)
    sce_pseudo <-
        registration_pseudobulk(
            sce_final,
            var_registration = current_var,
            var_sample_id = "RealSample",
            pseudobulk_rds_file = file.path(dir_rdata, paste0("pseudobulk_", current_var, ".rds"))
        )

    registration_mod <-
        registration_model(sce_pseudo, covars = "Run")

    block_cor <-
        registration_block_cor(sce_pseudo, registration_model = registration_mod)

    results_enrichment <-
        registration_stats_enrichment(
            sce_pseudo,
            block_cor = block_cor,
            covars = "Run",
            gene_ensembl = "ID",
            gene_name = "Symbol"
        )

    saveRDS(results_enrichment, file = file.path(dir_rdata, paste0("enrichment_", current_var, ".rds")))

    return(results_enrichment)
})
names(enrichment_tstats) <- registration_vars

head(enrichment_tstats$final_Annotations_broad)
#                 t_stat_Astrocyte t_stat_Endo t_stat_Excit.Thal
# ENSG00000238009        0.4941558   0.9298525        1.06714338
# ENSG00000241860        0.7344461  -2.5144058        1.34551148
# ENSG00000237491       -0.3227681  -0.1655309        0.66012110
# ENSG00000228794       -1.3250693  -2.4485440        1.01537446
# ENSG00000225880        1.1542415   0.2481395        0.02484108
# ENSG00000230368        0.1713202   0.2497699        1.75613918
#                 t_stat_Inhib.Thal t_stat_LHb t_stat_MHb t_stat_Microglia
# ENSG00000238009         1.2926161 -0.6457847  -1.972686       -0.5768099
# ENSG00000241860         1.7087987 -1.3389169  -1.756039        0.8147626
# ENSG00000237491         1.5402734  0.5331491  -2.065161       -2.0884211
# ENSG00000228794         0.9112540  1.0247560  -1.352135        0.3081748
# ENSG00000225880         0.3622976 -0.7895245  -1.578354        0.6886194
# ENSG00000230368         1.0299177 -0.9241120  -1.919964       -1.6361879
#                 t_stat_Oligo  t_stat_OPC p_value_Astrocyte p_value_Endo
# ENSG00000238009   -0.6875536  0.15240352         0.6238008   0.35782581
# ENSG00000241860   -1.4088542  1.38287453         0.4668052   0.01589105
# ENSG00000237491   -0.2446255  1.14172024         0.7484909   0.86933005
# ENSG00000228794   -0.1951413  0.57213103         0.1923952   0.01865473
# ENSG00000225880    0.1278139  0.22600230         0.2550025   0.80525040
# ENSG00000230368    0.9856077 -0.02113169         0.8648043   0.80399772
#                 p_value_Excit.Thal p_value_Inhib.Thal p_value_LHb
# ENSG00000238009         0.29206259         0.20328569   0.5219690
# ENSG00000241860         0.18575023         0.09494879   0.1878762
# ENSG00000237491         0.51281562         0.13107542   0.5967704
# ENSG00000228794         0.31579807         0.36741798   0.3114067
# ENSG00000225880         0.98030053         0.71896584   0.4342874
# ENSG00000230368         0.08642671         0.30900658   0.3607631
#                 p_value_MHb p_value_Microglia p_value_Oligo p_value_OPC
# ENSG00000238009  0.05521193        0.56717844     0.4955587   0.8796056
# ENSG00000241860  0.08645330        0.41985147     0.1663270   0.1740825
# ENSG00000237491  0.04518580        0.04292434     0.8079519   0.2601005
# ENSG00000228794  0.18364849        0.75948825     0.8462346   0.5703126
# ENSG00000225880  0.12207045        0.49488810     0.8989136   0.8223065
# ENSG00000230368  0.06174729        0.10935016     0.3300372   0.9832416
#                 fdr_Astrocyte   fdr_Endo fdr_Excit.Thal fdr_Inhib.Thal
# ENSG00000238009     0.8637648 0.48210598      0.8088415      0.5970315
# ENSG00000241860     0.7806029 0.04995283      0.7711246      0.5016208
# ENSG00000237491     0.9169061 0.90904042      0.8610136      0.5482505
# ENSG00000228794     0.5555164 0.05595557      0.8154249      0.6731578
# ENSG00000225880     0.6245792 0.86195194      0.9938470      0.8529653
# ENSG00000230368     0.9603161 0.86098997      0.6768724      0.6449839
#                   fdr_LHb   fdr_MHb fdr_Microglia fdr_Oligo   fdr_OPC
# ENSG00000238009 0.9518003 0.5097506    0.65387831 0.8410937 0.9553117
# ENSG00000241860 0.9509080 0.5770869    0.52043485 0.6518411 0.8513704
# ENSG00000237491 0.9555829 0.4820231    0.09619046 0.9452218 0.8513704
# ENSG00000228794 0.9509080 0.6954191    0.81683188 0.9539368 0.8533436
# ENSG00000225880 0.9509080 0.6332252    0.58817081 0.9726442 0.9291777
# ENSG00000230368 0.9509080 0.5241986    0.18930031 0.7680382 0.9942273
#                 logFC_Astrocyte logFC_Endo logFC_Excit.Thal
# ENSG00000238009       0.5169832  1.1725697       0.81335916
# ENSG00000241860       0.8443182 -3.2925692       1.12174363
# ENSG00000237491      -0.3198813 -0.1993164       0.47987575
# ENSG00000228794      -1.3037550 -2.7909634       0.74193769
# ENSG00000225880       1.2065456  0.3197624       0.01942206
# ENSG00000230368       0.1779030  0.3147986       1.29600976
#                 logFC_Inhib.Thal  logFC_LHb logFC_MHb logFC_Microglia
# ENSG00000238009        0.9786517 -0.4940665 -1.573655      -0.6010286
# ENSG00000241860        1.4062677 -1.1111834 -1.558357       0.9328809
# ENSG00000237491        1.0940957  0.3864335 -1.551804      -1.9650995
# ENSG00000228794        0.6670497  0.7450721 -1.056991       0.3082264
# ENSG00000225880        0.2826911 -0.6098200 -1.294387       0.7250262
# ENSG00000230368        0.7774936 -0.6963536 -1.519774      -1.6424652
#                 logFC_Oligo   logFC_OPC         ensembl       gene
# ENSG00000238009  -0.6272347  0.11714778 ENSG00000238009 AL627309.1
# ENSG00000241860  -1.3923262  1.14599996 ENSG00000241860 AL627309.5
# ENSG00000237491  -0.2121026  0.81767947 ENSG00000237491  LINC01409
# ENSG00000228794  -0.1713219  0.41956817 ENSG00000228794  LINC01128
# ENSG00000225880   0.1186566  0.17575396 ENSG00000225880  LINC00115
# ENSG00000230368   0.8849279 -0.01608617 ENSG00000230368     FAM41C

## Format matches the input we need for running
## https://research.libd.org/spatialLIBD/reference/layer_stat_cor.html
modeling_results <- fetch_data(type = "modeling_results")
head(modeling_results$enrichment)
#    t_stat_WM t_stat_Layer1 t_stat_Layer2 t_stat_Layer3 t_stat_Layer4
# 1 -0.6344143    -1.0321320    0.17815008   -0.72835965    1.56703859
# 2 -2.4758891     1.2232062   -0.87337451    1.93793650    1.33150141
# 3 -3.0079360    -0.8564572    2.13358520    0.48741121    0.35212807
# 4 -1.2916584    -0.9494234   -0.94854397    0.56378302   -0.11206713
# 5  2.3175897     0.6156900    0.11274780   -0.09907566   -0.03376771
# 6 -2.2686017    -0.6536163   -0.08615251    1.84786166    0.77710957
#   t_stat_Layer5 t_stat_Layer6  p_value_WM p_value_Layer1 p_value_Layer2
# 1    -0.2202707     0.7438713 0.527700348      0.3052551     0.85907467
# 2     0.4773214    -1.6152865 0.015497447      0.2249981     0.38518446
# 3     0.6071363     0.3832779 0.003557367      0.3944138     0.03607273
# 4     1.1536114     1.2634726 0.200358650      0.3453884     0.34583092
# 5    -0.3434730    -2.4629827 0.023143222      0.5399225     0.91052499
# 6     0.0355838     0.1945122 0.026108214      0.5153146     0.93156963
#   p_value_Layer3 p_value_Layer4 p_value_Layer5 p_value_Layer6     fdr_WM
# 1     0.46861049      0.1212213      0.8262455     0.45922628 0.63711486
# 2     0.05630996      0.1869669      0.6344905     0.11035355 0.03959651
# 3     0.62735659      0.7257077      0.5455535     0.70257361 0.01107944
# 4     0.57454588      0.9110629      0.2522414     0.21024494 0.31550730
# 5     0.92133658      0.9731501      0.7321823     0.01601949 0.05551142
# 6     0.06847596      0.4394838      0.9717067     0.84628893 0.06148729
#   fdr_Layer1 fdr_Layer2 fdr_Layer3 fdr_Layer4 fdr_Layer5 fdr_Layer6
# 1  0.5644497  0.9418694  0.8284720  0.4139767  0.9596911  0.8115244
# 2  0.4828399  0.6944277  0.3831776  0.4938944  0.9051814  0.5481106
# 3  0.6380674  0.2127681  0.9022535  0.8830293  0.8698217  0.9150742
# 4  0.6001298  0.6683198  0.8845131  0.9674045  0.7017599  0.6698502
# 5  0.7356770  0.9635319  0.9848034  0.9910492  0.9343088  0.2431891
# 6  0.7192630  0.9733708  0.4139515  0.7028656  0.9936069  0.9633326
#           ensembl        gene
# 1 ENSG00000243485 MIR1302-2HG
# 2 ENSG00000238009  AL627309.1
# 3 ENSG00000237491  AL669831.5
# 4 ENSG00000177757      FAM87B
# 5 ENSG00000225880   LINC00115
# 6 ENSG00000230368      FAM41C

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
#  date     2024-03-05
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
#  BiocVersion              3.18.1      2023-11-18 [1] Bioconductor 3.18 (R 4.3.2)
#  Biostrings               2.70.2      2024-01-30 [1] Bioconductor 3.18 (R 4.3.2)
#  bit                      4.0.5       2022-11-15 [1] CRAN (R 4.3.0)
#  bit64                    4.0.5       2020-08-30 [1] CRAN (R 4.3.0)
#  bitops                   1.0-7       2021-04-24 [1] CRAN (R 4.3.0)
#  blob                     1.2.4       2023-03-17 [1] CRAN (R 4.3.0)
#  brio                     1.1.4       2023-12-10 [1] CRAN (R 4.3.1)
#  bslib                    0.6.1.9001  2024-03-02 [1] Github (rstudio/bslib@243499a)
#  cachem                   1.0.8       2023-05-01 [1] CRAN (R 4.3.0)
#  cli                      3.6.2       2023-12-11 [1] CRAN (R 4.3.1)
#  codetools                0.2-19      2023-02-01 [1] CRAN (R 4.3.2)
#  colorout               * 1.3-0.2     2024-02-27 [1] Github (jalvesaq/colorout@c6113a2)
#  colorspace               2.1-0       2023-01-23 [1] CRAN (R 4.3.0)
#  config                   0.3.2       2023-08-30 [1] CRAN (R 4.3.0)
#  cowplot                  1.1.3       2024-01-22 [1] CRAN (R 4.3.1)
#  crayon                   1.5.2       2022-09-29 [1] CRAN (R 4.3.0)
#  curl                     5.2.0       2023-12-08 [1] CRAN (R 4.3.1)
#  data.table               1.15.0      2024-01-30 [1] CRAN (R 4.3.1)
#  DBI                      1.2.2       2024-02-16 [1] CRAN (R 4.3.1)
#  dbplyr                   2.4.0       2023-10-26 [1] CRAN (R 4.3.1)
#  DelayedArray             0.28.0      2023-11-06 [1] Bioconductor
#  DelayedMatrixStats       1.24.0      2023-11-06 [1] Bioconductor
#  devtools               * 2.4.5       2022-10-11 [1] CRAN (R 4.3.0)
#  digest                   0.6.34      2024-01-11 [1] CRAN (R 4.3.1)
#  doParallel               1.0.17      2022-02-07 [1] CRAN (R 4.3.0)
#  dotCall64                1.1-1       2023-11-28 [1] CRAN (R 4.3.1)
#  dplyr                    1.1.4       2023-11-17 [1] CRAN (R 4.3.1)
#  DT                       0.32        2024-02-19 [1] CRAN (R 4.3.1)
#  edgeR                    4.0.16      2024-02-20 [1] Bioconductor 3.18 (R 4.3.2)
#  ellipsis                 0.3.2       2021-04-29 [1] CRAN (R 4.3.0)
#  ExperimentHub            2.10.0      2023-10-26 [1] Bioconductor
#  fansi                    1.0.6       2023-12-08 [1] CRAN (R 4.3.1)
#  fastmap                  1.1.1       2023-02-24 [1] CRAN (R 4.3.0)
#  fields                   15.2        2023-08-17 [1] CRAN (R 4.3.0)
#  filelock                 1.0.3       2023-12-11 [1] CRAN (R 4.3.1)
#  foreach                  1.5.2       2022-02-02 [1] CRAN (R 4.3.0)
#  fs                       1.6.3       2023-07-20 [1] CRAN (R 4.3.0)
#  generics                 0.1.3       2022-07-05 [1] CRAN (R 4.3.0)
#  GenomeInfoDb           * 1.38.6      2024-02-10 [1] Bioconductor 3.18 (R 4.3.2)
#  GenomeInfoDbData         1.2.11      2024-02-27 [1] Bioconductor
#  GenomicAlignments        1.38.2      2024-01-20 [1] Bioconductor 3.18 (R 4.3.2)
#  GenomicRanges          * 1.54.1      2023-10-30 [1] Bioconductor
#  ggbeeswarm               0.7.2       2023-04-29 [1] CRAN (R 4.3.0)
#  ggplot2                  3.5.0       2024-02-23 [1] CRAN (R 4.3.1)
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
#  lubridate                1.9.3       2023-09-27 [1] CRAN (R 4.3.1)
#  magick                   2.8.3       2024-02-18 [1] CRAN (R 4.3.1)
#  magrittr                 2.0.3       2022-03-30 [1] CRAN (R 4.3.0)
#  maps                     3.4.2       2023-12-15 [1] CRAN (R 4.3.1)
#  Matrix                   1.6-5       2024-01-11 [1] CRAN (R 4.3.2)
#  MatrixGenerics         * 1.14.0      2023-10-26 [1] Bioconductor
#  matrixStats            * 1.2.0       2023-12-11 [1] CRAN (R 4.3.1)
#  memoise                  2.0.1       2021-11-26 [1] CRAN (R 4.3.0)
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
#  purrr                    1.0.2       2023-08-10 [1] CRAN (R 4.3.0)
#  R6                       2.5.1       2021-08-19 [1] CRAN (R 4.3.0)
#  rappdirs                 0.3.3       2021-01-31 [1] CRAN (R 4.3.0)
#  RColorBrewer             1.1-3       2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                     1.0.12      2024-01-09 [1] CRAN (R 4.3.1)
#  RCurl                    1.98-1.14   2024-01-09 [1] CRAN (R 4.3.1)
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
#  scuttle                  1.12.0      2023-11-06 [1] Bioconductor
#  sessioninfo              1.2.2       2021-12-06 [1] CRAN (R 4.3.0)
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
#  stringr                  1.5.1       2023-11-14 [1] CRAN (R 4.3.1)
#  SummarizedExperiment   * 1.32.0      2023-11-06 [1] Bioconductor
#  suncalc                  0.5.1       2022-09-29 [1] CRAN (R 4.3.0)
#  testthat               * 3.2.1       2023-12-02 [1] CRAN (R 4.3.1)
#  tibble                   3.2.1       2023-03-20 [1] CRAN (R 4.3.0)
#  tidyr                    1.3.1       2024-01-24 [1] CRAN (R 4.3.1)
#  tidyselect               1.2.0       2022-10-10 [1] CRAN (R 4.3.0)
#  timechange               0.3.0       2024-01-18 [1] CRAN (R 4.3.1)
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
