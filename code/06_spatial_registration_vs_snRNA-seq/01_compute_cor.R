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
# dir_input <- here("processed-data", "05_layer_differential_expression", "modeling_results_BS")
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
# head(registration_t_stats[[1]])
# Sp02D01    Sp02D02
# ENSG00000237491 -1.6855500  1.6855500
# ENSG00000228794 -1.6513162  1.6513162
# ENSG00000223764 -0.8991561  0.8991561

registration_vars <-
    c("final_Annotations", "final_Annotations_broad")

message(' Processing Spatial Registration for BayesSpace k=', k)

compute_cor <- function(current_var) {
    # Load input snRNA-seq data
    # testing: current_var = "final_Annotations"
    results_enrichment <-
        readRDS(here(
            "processed-data",
            "05_snRNA-seq_model_stats",
            paste0("enrichment_", current_var, ".rds")
        ))
    colnames(results_enrichment)
    modeling_res_enrichment <- list("enrichment" = results_enrichment)
    ## check out table
    # results_enrichment[1:5, 1:5]
    #                 t_stat_Sp02D01 t_stat_Sp02D02 p_value_Sp02D01 p_value_Sp02D02 fdr_Sp02D01
    # ENSG00000237491     -1.6855500      1.6855500       0.1263326       0.1263326   0.5432782
    # ENSG00000228794     -1.6513162      1.6513162       0.1332383       0.1332383   0.5432782
    # ENSG00000223764     -0.8991561      0.8991561       0.3921055       0.3921055   0.6688699
    # ENSG00000187634      0.5440598     -0.5440598       0.5996823       0.5996823   0.8063719
    # ENSG00000188976     -1.2522290      1.2522290       0.2422019       0.2422019   0.5724077

    lapply(
      registration_t_stats,
      layer_stat_cor, #results_enrichment
      modeling_results = modeling_res_enrichment,
      top_n = 100
    )
}

cor_fine <- compute_cor("final_Annotations")
# head(cor_fine[[1]])
#         Astrocyte       Endo Excit.Thal Inhib.Thal       LHb.1       LHb.2       LHb.3       LHb.4       LHb.5
# Sp02D01 -0.2038447  0.1249809 -0.1806295  -0.228601 -0.04530781  0.03818429  0.06463135  0.04138618  0.03397194
# Sp02D02  0.2038447 -0.1249809  0.1806295   0.228601  0.04530781 -0.03818429 -0.06463135 -0.04138618 -0.03397194

cor_broad <- compute_cor("final_Annotations_broad")
# head(cor_broad[[1]])
#         Astrocyte      Endo Excit.Thal Inhib.Thal        LHb          MHb    Microglia      Oligo         OPC
# Sp02D01 -0.2094872  0.123983 -0.1993256 -0.3184479  0.1457124 -0.009116221  0.003706224  0.5103645 -0.05427858
# Sp02D02  0.2094872 -0.123983  0.1993256  0.3184479 -0.1457124  0.009116221 -0.003706224 -0.5103645  0.05427858

## Annotate clusters / classify by layer confidence classes (good/poor)
# annotated_clusters_fine <-
#     lapply(cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
# annotated_clusters_fine[[2]]
# cluster layer_confidence     layer_label
# 1 Sp06D04             good     MHb.2/MHb.1
# 2 Sp06D02             good            Endo
# 3 Sp06D06             good           MHb.1
# 4 Sp06D01             good           Oligo
# 5 Sp06D03             poor Astrocyte/Endo*
# 6 Sp06D05             good       Astrocyte

# ## relaxed merging threshold 0.1 
# annotated_clusters_broad <-
#     lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
# # head(annotated_clusters_broad)

annotated_clusters_broad <-
  lapply(cor_broad, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)

## With default confidence and cutoff_merge_ratio 
annotated_clusters_fine <-
  lapply(cor_fine, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)

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

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_broadRes.pdf"), width = 10, height = 10)

for (i in seq_len(length(cor_broad))) {
    hm <- (layer_stat_cor_plot(
          cor_broad[[i]], annotation = annotated_clusters_broad[[i]],
          heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
          )
     )
    # Draw the heatmap with title
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
    # Draw the heatmap with title
    draw(
        hm,
        column_title = "Spatial-Registration: Visium vs snRNA (Fine res)",
        column_title_gp = gpar(fontsize = 16, fontface = "bold")
    )
  
}

dev.off()

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
