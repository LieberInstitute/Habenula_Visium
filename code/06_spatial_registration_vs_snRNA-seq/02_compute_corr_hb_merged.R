
######### Compute Spatial-Registration for both Fine and Broad snRNAseq - Habenula Merged Clusters (CSC) ########

library("here")
library("purrr")
library("spatialLIBD")
library("sessioninfo")


## Input dir
dir_input <- here("processed-data", "05_layer_differential_expression", "modeling_results_BS")

## Create output directories
dir_rdata <- here("processed-data", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
## Plot dir 
dir_plot <- here("plots", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

## specify the number of BayesSpace k to use 
k <- seq(2,28)

## Load Registration Results 
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

message(' Processing Spatial Registration for ', length(k), ' BayesSpace k')

# load the enrichment results for Hb - merged version
compute_cor <- function(current_var) {
  # Load input snRNA-seq data
  # testing: current_var = "final_Annotations"
  results_enrichment <-
    readRDS(here(
      "processed-data",
      "05_snRNA-seq_model_stats",
      paste0("enrichment_Hb_merged_", current_var, ".rds")
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
# # head(cor_fine[[1]])
# #         Astrocyte       Endo Excit.Thal Inhib.Thal       LHb.1       LHb.2       LHb.3       LHb.4       LHb.5
# # Sp02D01 -0.2038447  0.1249809 -0.1806295  -0.228601 -0.04530781  0.03818429  0.06463135  0.04138618  0.03397194
# # Sp02D02  0.2038447 -0.1249809  0.1806295   0.228601  0.04530781 -0.03818429 -0.06463135 -0.04138618 -0.03397194

cor_broad <- compute_cor("final_Annotations_broad")
# head(cor_broad[[1]])
#         Astrocyte      Endo Excit.Thal Inhib.Thal        LHb          MHb    Microglia      Oligo         OPC
# Sp02D01 -0.2094872  0.123983 -0.1993256 -0.3184479  0.1457124 -0.009116221  0.003706224  0.5103645 -0.05427858
# Sp02D02  0.2094872 -0.123983  0.1993256  0.3184479 -0.1457124  0.009116221 -0.003706224 -0.5103645  0.05427858

## Annotate clusters / classify by layer confidence classes (good/poor)
# annotated_clusters_fine <-
#     lapply(cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
# annotated_clusters_fine[[5]]
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


# Use annotation labels on the correlation matrices / restricted to 10 characters

cor_fine <- mapply(function(cor, label_data) {
    rownames(cor) <-
      substr(
        paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)]), 1, 20
      )
    return(cor)
}, cor_fine, annotated_clusters_fine)

cor_broad <- mapply(function(cor, label_data) {
  rownames(cor) <- # paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
    substr(
      paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)]), 1, 20
    )
  return(cor)
}, cor_broad, annotated_clusters_broad)

head(cor_fine)
head(cor_broad)

stopifnot(is.list(cor_fine))
head(cor_fine[[1]])
stopifnot(is.list(cor_broad))
head(cor_broad[[1]])
rownames(cor_broad)

## Confidence marks "x" need to be re-loaded ?

annotated_clusters_broad <-
  lapply(cor_broad, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)
 
## With default confidence and cutoff_merge_ratio 
annotated_clusters_fine <-
  lapply(cor_fine, annotate_registered_clusters, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25)

## add customized labels to rownames to match on heatmap

data.frame(
  "broad" = sort(rownames(substr(
    cor_broad[[8]])), 1, 20),
  "fine" = sort(rownames(substr(
    cor_fine[[8]])),1 , 20)
)

save(cor_fine,
     cor_broad,
     file = file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq_top100_Hb_merged.Rdata")
)

##   Make heatmaps broad res
cor_broad[[2]]
annotated_clusters_broad[[2]]

pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_broadRes_Hb_merged.pdf"))
for (i in seq_len(length(cor_broad))) {
  # i = 9
  print(
    layer_stat_cor_plot(
      cor_broad[[i]], annotation = annotated_clusters_broad[[i]],
      heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
    )
  )
}
dev.off()

# ##   Make heatmaps fine res
# 
# pdf(here(dir_plot, "cor_top100_visium_snRNAseq_registration_fineRes.pdf"))
# for (i in seq_len(length(cor_fine))) {
#   print(
#     layer_stat_cor_plot(
#       cor_fine[[i]], annotation = annotated_clusters_fine[[i]],
#       heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1))
#     )
#   )
# }
# dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
