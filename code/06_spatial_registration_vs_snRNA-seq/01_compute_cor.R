# library(slurmjobs)
# slurmjobs::job_single('01_compute_cor', create_shell = TRUE, memory = '20G', command = "01_compute_cor.R")

# To submit the job use: sbatch 01_compute_cor.sh

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

## For testing
if (is.na(k)) {
  k <- 2
}

library("here")
library("spatialLIBD")
library("sessioninfo")

# # ==================== Here my enrichment datasets
# load(here("processed-data", "05_layer_differential_expression", "modeling_results_BS", "modeling_results_BayesSpace_k02.Rdata"))
# ## Quick inspection to the BS models 
# results_enrichment <- modeling_results$enrichment
# # ==================== /

## Create output directories
dir_rdata <-
    here("processed-data", "06_spatial_registration_vs_snRNA-seq")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
## Input dir
dir_input <- 
  here("processed-data", "05_layer_differential_expression", "modeling_results_BS")

## Load BayesSpace enrichment t-stats
# BayesSpace_stats_list <-
#     readRDS(here(
#         "processed-data",
#         "04_harmony_BayesSpace",
#         "BayesSpace_stats_list.rds"
#     ))
# # head(BayesSpace_stats_list)
# # BayesSpace_stats_list$BayesSpace_harmony_k08$k08_1[1:10]
# # head(rownames(BayesSpace_stats_list$BayesSpace_harmony_k08))
# 
# ## Reformat to match expect format for spatialLIBD::layer_stat_cor
# ## Aka, match the format from
# ## head(spatialLIBD::tstats_Human_DLPFC_snRNAseq_Nguyen_topLayer)
# 
# BayesSpace_stats_list <- lapply(BayesSpace_stats_list, function(x) {
#     rownames(x) <- x$gene
#     x$gene <- NULL
#     colnames(x) <- gsub("BayesSpace_harmony_", "", colnames(x))
#     return(x)
# })
# # # names(BayesSpace_stats_list)
# # [1] "BayesSpace_harmony_k02" "BayesSpace_harmony_k03" "BayesSpace_harmony_k04" "BayesSpace_harmony_k05"
# # [5] "BayesSpace_harmony_k06" "BayesSpace_harmony_k07" "BayesSpace_harmony_k08" "BayesSpace_harmony_k09"
# # [9] "BayesSpace_harmony_k10" "BayesSpace_harmony_k11" "BayesSpace_harmony_k12" "BayesSpace_harmony_k13"

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
bayesSpace_registration <-
  lapply(bayesSpace_registration_fn, function(x) {
    get(load(x))
  })


## Select t-stats from the registration enrichment data
registration_t_stats <-
  map(bayesSpace_registration, function(data) {
    x <- data$enrichment
    t_stats <- x[, grep("^t_stat_", colnames(x))]
    colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))
    return(t_stats)
  })

registration_vars <-
    c("final_Annotations", "final_Annotations_broad")

compute_cor <- function(current_var) {
    # Load input snRNA-seq data
    results_enrichment <-
        readRDS(here(
            "processed-data",
            "05_snRNA-seq_model_stats",
            paste0("enrichment_", current_var, ".rds")
        ))
    # compute the correlation 
    lapply(
        BayesSpace_stats_list,
        layer_stat_cor,
        modeling_results = "enrichment",
        top_n = 100
    )
}

cor_fine <- compute_cor("final_Annotations")
head(cor_fine)
cor_broad <- compute_cor("final_Annotations_broad")

## Explore the correlation matrix in BS K8
#head(cor_broad)
head(cor_broad$BayesSpace_harmony_k08[, seq_len(3)])
summary(cor_broad$BayesSpace_harmony_k08)

## Annotate clusters / classify by layer confidence classes (good/poor)
annotated_clusters_fine <-
    lapply(cor_fine, annotate_registered_clusters, cutoff_merge_ratio = 0.1)
head(annotated_clusters_fine)
annotated_clusters_broad <-
    lapply(cor_broad, annotate_registered_clusters, cutoff_merge_ratio = 0.1)

## Use annotation labels on the correlation matrices
cor_fine <- mapply(function(cor, label_data) {
    rownames(cor) <-
        paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
    return(cor)
}, cor_fine, annotated_clusters_fine)

cor_broad <- mapply(function(cor, label_data) {
    rownames(cor) <-
        paste0(rownames(cor), " ~ ", label_data$layer_label[match(rownames(cor), label_data$cluster)])
    return(cor)
}, cor_broad, annotated_clusters_broad)

save(cor_fine,
    cor_broad,
    file = file.path(dir_rdata, "cor_BayesSpace_vs_snRNA-seq.Rdata")
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
