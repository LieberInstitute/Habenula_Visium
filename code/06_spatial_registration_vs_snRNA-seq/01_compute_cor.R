################################################################################
## Compute Spatial-Registration for both Fine and Broad snRNAseq vs Multiome snRNAseq (CSC)
##
## Notes:
## For 60 to 80k spots: $srun --pty --mem=60GB --x11 bash
##
## Authors. CSC
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
# Old version: 
dir_input <- here("processed-data", "05_layer_differential_expression", "modeling_results_BS_old")
#dir_input <- here("processed-data", "05_layer_differential_expression", "modeling_results_BS")


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
# Testing reproducibility - OLD results:
head(registration_t_stats[[1]])
#                   Sp03D01    Sp03D02    Sp03D03
# ENSG00000228794 -1.6049998  0.3178627  1.2189367
# ENSG00000223764 -1.2436968  3.1261408 -1.3020605
# ENSG00000187634  0.2270575  2.6447516 -3.0780594
# ENSG00000188976 -0.7268748 -0.7570369  1.5731080
# ENSG00000187961 -2.6701112  1.4897980  0.8193468
# ENSG00000272512 -0.1668633 -0.9080208  1.0890365

# Testing reproducibility - ReRun results:
#                   Sp03D01    Sp03D02    Sp03D03
# ENSG00000187634 -3.772719 -0.0582146  3.9975042
# ENSG00000188976 -2.853437  0.9240683  1.7720477
# ENSG00000188290 -3.035521  0.5067032  2.4564039
# ENSG00000187608 -1.277774  1.6266120 -0.3385263
# ENSG00000188157 -1.891043  1.1974406  0.6448430
# ENSG00000078808  3.130198 -0.0788782 -3.1740205

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
# annotated_clusters_fine <-
#     lapply(cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
# annotated_clusters_fine[[2]]

# ## relaxed merging threshold 0.1 
# annotated_clusters_broad <-
#     lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
# # head(annotated_clusters_broad)

annotated_clusters_broad <-
  lapply(cor_broad, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)
# git show c8a68c62df15af418ff3a78dc29bfd7fbedc31a0
# annotated_clusters_broad <-
#     lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
annotated_clusters_broad[[1]]

# Testing OLD results:
#   cluster layer_confidence layer_label
# 1 Sp03D03             good  Inhib.Thal
# 2 Sp03D01             good       Oligo
# 3 Sp03D02             good         MHb

# Testing reproducibility - ReRun results:
# cluster layer_confidence layer_label
# 1 Sp03D03             good         MHb
# 2 Sp03D01             good       Oligo
# 3 Sp03D02             good   Astrocyte

## With default confidence and cutoff_merge_ratio 
annotated_clusters_fine <-
  lapply(cor_fine, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)
annotated_clusters_fine[[1]]

# Testing OLD results:
#   cluster layer_confidence       layer_label
# 1 Sp03D03             good        Inhib.Thal
# 2 Sp03D01             good             Oligo
# 3 Sp03D02             good MHb.1/MHb.2/LHb.6

# Testing reproducibility - ReRun results:
#   cluster layer_confidence             layer_label
# 1 Sp03D03             good MHb.2/MHb.1/LHb.6/LHb.2
# 2 Sp03D01             good                   Oligo
# 3 Sp03D02             good               Astrocyte

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
    file = file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq_top100.Rdata")
)

##   Make heatmaps broad res

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_broadRes.pdf"), 
    width = 10, height = 10)

for (i in seq_len(length(cor_broad))) {
    hm <- (layer_stat_cor_plot(
          cor_broad[[i]], annotation = annotated_clusters_broad[[i]],
          heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
          )
     )
    draw(
        hm,
        column_title = "Spatial-Registration: Visium vs snRNA (Broad res)",
        column_title_gp = gpar(fontsize = 16, fontface = "bold")
    )
}

dev.off()

##   Make heatmaps fine res

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_fineRes.pdf"), width = 10, height = 10)

for (i in seq_len(length(cor_fine))) {
  hm <- (layer_stat_cor_plot(
      cor_fine[[i]], annotation = annotated_clusters_fine[[i]],
      heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
    )
  )
    draw(
        hm,
        column_title = "Spatial-Registration: Visium vs snRNA (Fine res)",
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
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
