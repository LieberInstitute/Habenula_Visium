################################################################################
## Compute Spatial-Registration for both Fine and Broad snRNAseq vs Multiome snRNAseq (CSC)
##
## Notes:
## For 60 to 80k spots: $srun --pty --mem=60GB --x11 bash
##
## Authors. Implementation CSC
##
#################### BayesSpace vs Multiome-snRNAseq ###########################

library("here")
library("purrr")
library("spatialLIBD")
library("ComplexHeatmap")
library("grid") # need to print the plot, otherwise is clipped by internal function of layer_stat_cor_plot()
library("sessioninfo")

# spatialLIBD / * 1.21.5 / 2025-05-16 [1] Github (LieberInstitute/spatialLIBD@aff00db)


## Input dir
dir_input <- here("processed-data", "05_brain_area_differential_expression", "modeling_results_BS")

## Create output directories
dir_rdata <- here("processed-data", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
## Plot dir 
dir_plot <- here("plots", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

## specify the number of BayesSpace k to use 
#k=13
k <- seq(2,28)

## Load Visium Registration Results 
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
# Testing reproducibility with k=3 - OLD results before remove Br.6522:
head(registration_t_stats[[1]])
#                   Sp03D01    Sp03D02    Sp03D03
# ENSG00000228794 -1.6049998  0.3178627  1.2189367
# ENSG00000223764 -1.2436968  3.1261408 -1.3020605
# ENSG00000187634  0.2270575  2.6447516 -3.0780594
# ENSG00000188976 -0.7268748 -0.7570369  1.5731080
# ENSG00000187961 -2.6701112  1.4897980  0.8193468
# ENSG00000272512 -0.1668633 -0.9080208  1.0890365

# current results after removing Br.6522:
#                   Sp03D01     Sp03D02    Sp03D03
# ENSG00000187634 -3.633848 -0.07005073  3.8709752
# ENSG00000188976 -2.656709  0.86672109  1.6543627
# ENSG00000188290 -3.408015  0.53459515  2.7546149
# ENSG00000187608 -1.215121  1.63777544 -0.4070374
# ENSG00000188157 -1.953529  1.22393635  0.6704573
# ENSG00000078808  3.201552 -0.08369606 -3.2462649

registration_vars <-
    c("final_Annotations", "final_Annotations_broad")

message(' Processing Spatial Registration for BayesSpace k=', k)

compute_cor <- function(current_var) {
    # Load input snRNA-seq data
    # testing: current_var = "final_Annotations_broad"
    results_enrichment <-
        readRDS(here(
            "processed-data",
            "05_snRNA-seq_model_stats",
            paste0("enrichment_", current_var, ".rds")
        ))
    colnames(results_enrichment)
    modeling_res_enrichment <- list("enrichment" = results_enrichment)
    results_enrichment[1:5,]

    lapply(
      registration_t_stats,
      layer_stat_cor, #results_enrichment
      modeling_results = modeling_res_enrichment,
      top_n = 100
    )
}

cor_fine <- compute_cor("final_Annotations")
head(cor_fine[[1]])

cor_broad <- compute_cor("final_Annotations_broad")
head(cor_broad[[1]])

## Annotate clusters / classify by layer confidence classes (good/poor)
annotated_clusters_broad <-
  lapply(cor_broad, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)
# git show c8a68c62df15af418ff3a78dc29bfd7fbedc31a0
# annotated_clusters_broad <-
#     lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
annotated_clusters_broad[[1]]

## With default confidence and cutoff_merge_ratio 
annotated_clusters_fine <-
  lapply(cor_fine, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)
annotated_clusters_fine[[1]]

## Use annotation labels on the correlation matrices
# cor_fine <- mapply(function(cor, label_data) {
#     rownames(cor) <- paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
#     return(cor)
# }, cor_fine, annotated_clusters_fine, SIMPLIFY = FALSE)
stopifnot(is.list(cor_fine))
head(cor_fine[[1]])
# Check if the row names and the cluster names match
data.frame(
    "broad-corr" = sort(rownames(cor_fine[[1]])),
    "broad-ann" = sort(annotated_clusters_fine[[1]]$cluster)
)

# cor_broad <- mapply(function(cor, label_data) {
#     rownames(cor) <- paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
#     return(cor)
# }, cor_broad, annotated_clusters_broad, SIMPLIFY = FALSE)
stopifnot(is.list(cor_broad))
head(cor_broad[[1]])
## Verify if levels match
data.frame(
    "broad-corr" = sort(rownames(
        cor_broad[[1]]
    )),
    "broad-ann" = sort((annotated_clusters_broad[[1]]$cluster))
)


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
    file = file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq_top100.Rdata")
)

message("Make heatmaps broad res")

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_broadRes.pdf"), 
    width = 10, height = 10)

for (i in seq_len(length(cor_broad))) {
    hm <- (layer_stat_cor_plot(
          cor_broad[[i]], 
          annotation = annotated_clusters_broad[[i]],
          heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
          column_names_gp = gpar(fontsize = 16),
          row_names_gp = gpar(fontsize = 16),
          cluster_rows = FALSE)
     )
    draw(
        hm,
        column_title = "Spatial-Registration: Visium vs snRNA (Broad res)",
        column_title_gp = gpar(fontsize = 20, fontface = "bold")
    )
}

dev.off()

message("Make heatmaps fine res")

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_fineRes.pdf"), width = 10, height = 10)

for (i in seq_len(length(cor_fine))) {
  hm <- (layer_stat_cor_plot(
      cor_fine[[i]], 
      annotation = annotated_clusters_fine[[i]],
      heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
      column_names_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      cluster_rows = FALSE)
  )
    draw(
        hm,
        column_title = "Spatial-Registration: Visium vs snRNA (Fine res)",
        column_title_gp = gpar(fontsize = 20, fontface = "bold")
    )
}

dev.off()


################################################################################
##  Compute correlation for FINE cluster annotations
##  Subset Hb cell-types of interest
################################################################################

message("Make heatmaps fine res for Hb-subset")

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_fineRes_HbSubset.pdf"), width = 10, height = 10)

for (i in seq_len(length(cor_fine))) {
    # testing: i=2
    
    ## Prepare matrix 
    cor_fine_subset <- as.data.frame(cor_fine[[i]])
    # subset columns that contain "LHb", "MHb", or "Thal"? in the matrix 
    colnames(cor_fine_subset)
    rownames(cor_fine_subset)
    cor_fine_subset <- cor_fine_subset[, grep("LHb|MHb|Thal", colnames(cor_fine_subset))]
    colnames(cor_fine_subset)
    # check
    head(cor_fine_subset)
    
    # extract Hb annotations of interest from 'annotated_clusters_fine' if the SpD(s) are annotated
    # annotated_clusters_fine_subset <- annotated_clusters_fine[[i]][
    #     grepl("LHb|MHb|Thal", annotated_clusters_fine[[i]]$cluster), ]
    # if not:
    annotated_clusters_fine_subset <- annotated_clusters_fine[[i]]
    head(annotated_clusters_fine_subset)
    
    # prepare heatmap
    hm <- (layer_stat_cor_plot(
        as.matrix(cor_fine_subset), 
        annotation = annotated_clusters_fine_subset,
        heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
        column_names_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 14),
        cluster_rows = FALSE  # <-- turn off row clustering
    )
    )
    draw(
        hm,
        column_title = "Spatial-Registration: Hb Visium vs snRNA (Fine res)",
        column_title_gp = gpar(fontsize = 16, fontface = "bold")
    )
    
}

dev.off()


################################################################################
##  Compute correlation for BROAD cluster annotations
##  Subset Hb cell-types of interest
################################################################################

message("Make heatmaps broad res for Hb-subset")

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_broadRes_HbSubset.pdf"), 
    width = 10, height = 10)

for (i in seq_len(length(cor_fine))) {
    # testing: i=2
    
    ## Prepare matrix 
    cor_broad_subset <- as.data.frame(cor_broad[[i]])
    cor_broad_subset <- cor_broad_subset[, grep("LHb|MHb", colnames(cor_broad_subset))]
    colnames(cor_broad_subset)
    # check
    head(cor_broad_subset)
    # if not annotated:
    annotated_clusters_broad_subset <- annotated_clusters_broad[[i]]
    head(annotated_clusters_broad_subset)
    
    # prepare heatmap
    hm <- (layer_stat_cor_plot(
        as.matrix(cor_broad_subset), 
        annotation = annotated_clusters_broad_subset,
        heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
        column_names_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 14),
        cluster_rows = FALSE  # <-- turn off row clustering
    )
    )
    draw(
        hm,
        column_title = "Spatial-Registration: Hb Visium vs snRNA (Broad res)",
        column_title_gp = gpar(fontsize = 16, fontface = "bold")
    )
    
}

dev.off()


################################################################################
##  Compute correlation for BROAD Hb MERGE cluster annotations
##  Subset Hb cell-types of interest
################################################################################

message("Make heatmaps broad res for Hb-subset merged")

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_broadRes_HbSubset_merge.pdf"), 
    width = 10, height = 10)

for (i in seq_len(length(cor_fine))) {
    # testing: i=2
    
    ## Prepare matrix 
    cor_broad_subset <- as.data.frame(cor_broad[[i]])
    # subset columns that contain "LHb", "MHb", or "Thal"? in the matrix 
    cor_broad_subset <- cor_broad_subset[, grep("LHb|MHb|Thal", colnames(cor_broad_subset))]
    # if not annotated:
    annotated_clusters_broad_subset <- annotated_clusters_broad[[i]]
    head(annotated_clusters_broad_subset)
    
    ## Merge MHb and LHb clusters
    cor_broad_subset$Hb_merge <- rowMeans(cor_broad_subset[, c("LHb", "MHb")])
    cor_broad_subset <- cor_broad_subset[, !(colnames(cor_broad_subset) %in% c("LHb", "MHb"))]
    cor_broad_subset <- cor_broad_subset[, c("Excit.Thal", "Inhib.Thal", "Hb_merge")]
    # prepare heatmap
    hm_merge <- (layer_stat_cor_plot(
        as.matrix(cor_broad_subset), 
        annotation = annotated_clusters_broad_subset,
        heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
        column_names_gp = gpar(fontsize = 14),
        row_names_gp = gpar(fontsize = 14),
        cluster_rows = FALSE  # <-- turn off row clustering
    )
    )
    draw(
        hm_merge,
        column_title = "Spatial-Registration: Hb Visium vs snRNA (Broad res)",
        column_title_gp = gpar(fontsize = 16, fontface = "bold")
    )
    
}

dev.off()

message("Spatial Registration DONE!!!")



# library(slurmjobs)
# slurmjobs::job_single('01_compute_cor', create_shell = TRUE, memory = '20G', command = "01_compute_cor.R")
# 
# To submit the job use: sbatch 01_compute_cor.sh


## Reproducibility information
# print("Reproducibility information:")
# Sys.time()
# proc.time()
# options(width = 120)
# session_info()
# 
# > Sys.time()
# [1] "2025-05-26 14:17:23 EDT"
# > proc.time()
# user   system  elapsed 
# 20.314    2.126 7554.341 
# > options(width = 120)
# > session_info()
# warm               0.4.0     2021-06-01 [2] CRAN (R 4.4.0)
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
# Cairo                  1.6-2     2023-11-28 [2] CRAN (R 4.4.0)
# circlize               0.4.16    2024-02-20 [2] CRAN (R 4.4.0)
# cli                    3.6.5     2025-04-23 [2] CRAN (R 4.4.3)
# clue                   0.3-66    2024-11-13 [2] CRAN (R 4.4.2)
# cluster                2.1.8     2024-12-11 [3] CRAN (R 4.4.3)
# codetools              0.2-20    2024-03-31 [3] CRAN (R 4.4.3)
# colorspace             2.1-1     2024-07-26 [2] CRAN (R 4.4.1)
# ComplexHeatmap       * 2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
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
# purrr                * 1.0.4     2025-02-05 [2] CRAN (R 4.4.2)
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
# rstudioapi             0.17.1    2024-10-22 [2] CRAN (R 4.4.2)
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
