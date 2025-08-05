# copied from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/08_spatial_registration/

################################################################################
## Compute Spatial-Registration between Hb Visium and Hb Multiome WNN clusters
##
## Notes:
## For 60 to 80k spots: $srun --pty --mem=60GB --x11 bash
##
## Authors. CSC
##
#################### BayesSpace vs Multiome-snRNAseq ###########################

library("spatialLIBD")
library("tidyverse")
library("jaffelab")
library("ComplexHeatmap")
library("grid") # need to print the plot, otherwise is clipped by internal function of layer_stat_cor_plot()
library("here")
library("sessioninfo")

## Input dir
dir_input <- here(
  "processed-data",
  "05_brain_area_differential_expression",
  "modeling_results_BS"
)

## Set up plotting
plot_dir <- here("plots", 
    "07_spatial_registration_vs_multiome_snRNA-seq"
)
data_dir <- here(
  "processed-data",
  "07_spatial_registration_vs_multiome_snRNA-seq"
)
plt_sufix <- "v5" 

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
}
if (!dir.exists(data_dir)) {
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
}

## Load enrichment data from multiome
rds_input <- here("processed-data", "05_snRNA-seq_model_stats", paste0("enrichment_snRNA-multiome_", plt_sufix, ".rds"))


sn_multiome_data <- readRDS(rds_input)
head(sn_multiome_data$enrichment[5:10])
#                   t_stat_C.05.DD_LHb t_stat_C.06 t_stat_C.07.DD_MHb t_stat_C.08
# ENSG00000238009          1.2234449    1.759594         -0.5585619  2.26415183
# ENSG00000241860          0.2151696    1.234082         -0.2209549  0.07914740
# ENSG00000237491          0.4194972    1.935467         -0.7272817  1.36389851
# ENSG00000228794          1.1454835    1.268358         -0.7307091  1.83122297
# ENSG00000225880          1.0900785    1.994628          0.2073577 -0.06977143
# ENSG00000230368          0.5377206    1.413340         -0.4024732  0.61983127

grep("t_stat", colnames(sn_multiome_data$enrichment), value = TRUE)


## extract only enrichment stats and sorted the t-stats by hb and no hb clusters

multiome_t_stats_sorted <- function(sn_data) {  
  
      sn_data <- sn_multiome_data
      x <- sn_data$enrichment
      # extract t-stats
      t_stats <- x[, grep("^t_stat_", colnames(x)), drop = FALSE]

      newname_clusters <- colnames(t_stats)
      ## split hb-clusters from not hb clusters
      # logical vector of matches
      is_habenula <- grepl("MHb|LHb", newname_clusters)
      habenula_clusters <- newname_clusters[is_habenula]
      non_habenula_clusters <- newname_clusters[!is_habenula]
      desired_order <- c(habenula_clusters, non_habenula_clusters)
      # Reorder t_stats 
      t_stats <- t_stats[, desired_order]
      
      # add ensembl and gene columns
      t_stats$ensembl <- rownames(t_stats)
      if ("gene" %in% colnames(x)) {
          t_stats$gene <- x$gene
      } else {
          t_stats$gene <- t_stats$ensembl
      }
      
      # reorder columns: ensembl, gene, then t-statistics
      t_stats <- t_stats[, c("ensembl", "gene", desired_order)]
      
      # replace sn_data$enrichment
      sn_data$enrichment <- t_stats
      colnames(t_stats)
      
    return(sn_data)
  
}



## Load Visium Registration Results 
k_list <- c(2:28) 
#k_list <- c(3,9,13,21,26) # --> previous version
#k_list <- c(3,11,15,20,24,28)

# Use naming convention
names(k_list) <- paste0("k", sprintf("%02d", k_list)) 
names(k_list)

## Load Visium Registration Results 

bayesSpace_registration_fn <-
  map(k_list, ~ here(
    dir_input,
    paste0(
      "modeling_results_BayesSpace_k",
      sprintf("%02d", .x),
      ".Rdata"
    )
  ))
bayesSpace_registration_fn

## Load the 3 model results (anova, enrichment, pairwise) for each domain in the BS
bayesSpace_registration <-
  lapply(bayesSpace_registration_fn, function(x) {
    get(load(x))
  })

stopifnot(is.list(bayesSpace_registration))
names(bayesSpace_registration[[1]])
# [1] "anova"      "enrichment" "pairwise"  

## Select t-stats from the registration enrichment data
registration_t_stats <-
  map(bayesSpace_registration, function(data) {
    x <- data$enrichment
    t_stats <- x[, grep("^t_stat_", colnames(x))]
    colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))
    return(t_stats)
  })
stopifnot(is.list(registration_t_stats))
#str(registration_t_stats)
#head(registration_t_stats[[1]])
                  # Sp02D01    Sp02D02
# ENSG00000237491 -1.6855500  1.6855500
# ENSG00000228794 -1.6513162  1.6513162
# ENSG00000223764 -0.8991561  0.8991561

map(registration_t_stats, jaffelab::corner)


## compute t-statistic and plot correlations for selected BayesSpace k(s)
plt_corr_snmultiome <- function(k_lst, 
                                tstats, tstats_multiome,
                                plt_name) {

    pdf(plt_name, width = 10, height = 10)

    for (kl in names(k_lst)) {
    # kl = names(k_lst[1])
      
        k = names(k_lst[kl])
        message("Processing Spatial-Registration for BayesSpace ", k)
        # extract t-stats for the specific BS k
        tstats_k <- tstats[[k]]
        #head(tstats_k)
        
        cor_layer <- layer_stat_cor(
          stats = tstats_k,
          modeling_results = tstats_multiome,
          model_type = "enrichment",
          top_n = 100
        )
        #print(head(cor_layer))
        message("Correlation computed ...")
        
        annotated_clusters <- annotate_registered_clusters(cor_layer, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)
        #head(annotated_clusters)

        # cor_layer <-
        #   rownames(cor_layer) <- paste0(rownames(cor_layer), " ~ ", annotated_clusters[match(rownames(cor_layer), annotated_clusters$cluster)])
        
        rdata_name <- paste0("bayesSpace_cor_top100_Visium_Multiome_", k, ".Rdata")
        save(cor_layer, file = here(data_dir, rdata_name))
        
        ## make spatial-registration heatmap for specific k
        p1 <- layer_stat_cor_plot(
            cor_layer, annotation = annotated_clusters,
            heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
            column_names_gp = gpar(fontsize = 16),
            row_names_gp = gpar(fontsize = 16),
            cluster_rows = FALSE
            )
        draw(
            p1,
            column_title = "Spatial-Registration: Visium vs Multiome (Fine res)",
            column_title_gp = gpar(fontsize = 20, fontface = "bold")
            )
    
    }
   
    dev.off()
    message("Plot done!")

}


## call function to compute SpatialRegistration of Visium vs snMultiome
plt_name <- here(plot_dir, paste0("cor_top100_registration_Visium_snMultiome_all_clusters_", plt_sufix, ".pdf"))
plt_corr_snmultiome(
    k_list,
    registration_t_stats, sn_multiome_data,
    plt_name
)

plt_name <- here(plot_dir, paste0("cor_top100_registration_Visium_snMultiome_all_clusters_sorted_", plt_sufix, ".pdf"))
sn_multiome_sorted <- multiome_t_stats_sorted(
    sn_multiome_data
)

plt_corr_snmultiome(
    k_list, 
    registration_t_stats, sn_multiome_sorted,
    plt_name
)

## =============================================================================

# ## Prepare Hb subset for cell-types of interest: "LHb", "MHb" and "Thal"
# sn_multiome_Hb_subset <- sn_multiome_data
# 
# #colnames(sn_multiome_Hb_subset$enrichment)
# #rownames(sn_multiome_Hb_subset$enrichment)
# sn_multiome_Hb_subset$enrichment <- sn_multiome_Hb_subset$enrichment[, grep("LHb|MHb|Thal", 
#                                                                             colnames(sn_multiome_Hb_subset$enrichment))]
# 
# plt_name <- here(plot_dir, paste0("cor_top100_registration_Visium_snMultiome_Hb_clusters_", plt_sufix,".pdf"))
# 
# plt_corr_snmultiome(
#     k_list,
#     registration_t_stats, 
#     sn_multiome_Hb_subset,
# plt_name)

message("Spatial registration Visium vs Multiome, Done!")


#library("slurmjobs")
## A regular job with 10 cores on the 'imaginary' partition
#job_single("02_compute_cor_visium_multiomeRnaseq", cores = 2, partition = "katun", create_shell = TRUE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


