# copied from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/08_spatial_registration/

################################################################################
## Compute Spatial-Registration heatmpas between Visium cell-types and Multiome WNN clusters
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
library("here")
library("sessioninfo")

## Input dir
dir_input <- here(
  "processed-data",
  #"05_layer_differential_expression",
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
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
}
if (!dir.exists(data_dir)) {
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
}
## Load data
# load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

### testing layer_cor on rna-multiome reference 

## Load enrichment data from multiome

# suffix_reference = "v1" # have all the clusters containing Hb, even those with mix cells. This is an EDA version
# rds_input <- here("processed-data", "05_snRNA-seq_model_stats", "enrichment_snRNA-multiome.rds") 

suffix_reference = "v2" # have Hb clusters containing Hb in 50%+. This is a polished version
rds_input <- here("processed-data", "05_snRNA-seq_model_stats", "enrichment_snRNA-multiome_v2.rds")


sn_multiome_data <- readRDS(rds_input)
head(sn_multiome_data$enrichment[5:10])
#                   t_stat_C.05.DD_LHb t_stat_C.06 t_stat_C.07.DD_MHb t_stat_C.08
# ENSG00000238009          1.2234449    1.759594         -0.5585619  2.26415183
# ENSG00000241860          0.2151696    1.234082         -0.2209549  0.07914740
# ENSG00000237491          0.4194972    1.935467         -0.7272817  1.36389851
# ENSG00000228794          1.1454835    1.268358         -0.7307091  1.83122297
# ENSG00000225880          1.0900785    1.994628          0.2073577 -0.06977143
# ENSG00000230368          0.5377206    1.413340         -0.4024732  0.61983127

head(colnames(unique(sn_multiome_data$enrichment)))
# # ======= clusters on v2 (polished version)
# [1] "t_stat_C.01"         "t_stat_C.02"         "t_stat_C.03"        
# [4] "t_stat_C.04"         "t_stat_C.05.DD_LHb"  "t_stat_C.06"        
# [7] "t_stat_C.07.DD_MHb"  "t_stat_C.08"         "t_stat_C.09"        
# [10] "t_stat_C.10.DD_MHb"  "t_stat_C.11.DD_MHb"  "t_stat_C.12"      
# ======= clusters on v1
# [1] "t_stat_C.01.DD_LHb"  "t_stat_C.02"         "t_stat_C.03"        
# [4] "t_stat_C.04.DD_LHb"  "t_stat_C.05.DD_LHb"  "t_stat_C.06"        
# [7] "t_stat_C.07.DD_MHb"  "t_stat_C.08.DD_LHb"  "t_stat_C.09"        
# [10] "t_stat_C.10.DD_MHb"  "t_stat_C.11.DD_MHb"  "t_stat_C.12.DD_LHb" 


## extract only enrichment stats and sorted the t-stats by hb and no hb clusters

snRNA_t_stats_sorted <- function(sn_data) {  
  
  sn_data <- sn_multiome_data
  x <- sn_data$enrichment
  t_stats <- x[, grep("^t_stat_", colnames(x))]
  colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))
  newname_clusters <- colnames(t_stats)
  
  ## split hb-clusters from not hb clusters
  no_hb_clust = list()
  hb_clusters = list()
  desired_order = list()
  for (idx in newname_clusters) { if (nchar(idx) <= 4) { no_hb_clust <- append(no_hb_clust, idx) } }
  no_hb_clust <- sort(c(unlist(no_hb_clust)))
  hb_clusters <- sort(newname_clusters[! newname_clusters %in% c(no_hb_clust)])
  desired_order <- c(hb_clusters, no_hb_clust)
  
  ## Set new class order on the enrichment data set
  t_stats <- t_stats[, desired_order]
  colnames(t_stats)
  sn_data$enrichment <- t_stats
  colnames(sn_data$enrichment)
  return(sn_data)
  
}


## Load Registration Results 

#k_list <- c(2:28) --> old version all clusters included
k_list <- c(3,9,13,21,26)
# testing k_list=13
names(k_list) <- paste0("k", sprintf("%02d", k_list)) ## Use naming convention
# [1] "k03" "k09" "k13" "k21" "k26"

## Load Registration Results 

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
# bayesSpace_registration_fn[12]
# [1] "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/05_brain_area_differential_expression/modeling_results_BS/modeling_results_BayesSpace_k11.Rdata"

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
head(registration_t_stats[[1]])
                  # Sp02D01    Sp02D02
# ENSG00000237491 -1.6855500  1.6855500
# ENSG00000228794 -1.6513162  1.6513162
# ENSG00000223764 -0.8991561  0.8991561

map(registration_t_stats, jaffelab::corner)


#### Calculate Correlation Matrix ####

## get layer data
# layer_modeling_results <- fetch_data(type = "modeling_results")


## CSC. Added from https://research.libd.org/spatialLIBD/articles/guide_to_spatial_registration.html#correlate-statsics-with-layer-reference
#str(registration_t_stats$k09)

#colnames(bayesSpace_registration$k09$enrichment)

## correlate registration bayes space registration vs snRNA multiome enrichment

# sn_multiome_data <-  snRNA_t_stats_sorted(sn_multiome_data)
# colnames(sn_multiome_data$enrichment)

## loop on the k(s) 

plt_corr_snmultiome <- function(
        suffix_name,
        plt_name) {

  #plt_name <- paste0("cor_top100_spatial_registration_snMultiome_", suffix_name,".pdf")
  pdf(here(plot_dir, plt_name))

  for (k in names(k_list)) {
    # k = "13"
    message("Processing Spatial-Registration for BayesSpace k", k)
    # bayesSpace_registration[['k13']]$enrichment
    k = paste0("k", k)
    bayesSpace_registration_k <- bayesSpace_registration[[k]]$enrichment
    # head(bayesSpace_registration_k)
    cor_layer <- layer_stat_cor(
      stats = bayesSpace_registration_k, #bayesSpace_registration$k13$enrichment,
      modeling_results = sn_multiome_data,
      model_type = "enrichment",
      top_n = 100
    )
    head(cor_layer)
    #             C.01       C.02       C.03       C.04 C.05.DD_LHb       C.06
    # Sp02D01 -0.2444767  0.5483361 -0.1990111 -0.1114921   -0.150348 -0.1617148
    # Sp02D02  0.2444767 -0.5483361  0.1990111  0.1114921    0.150348  0.1617148
    # colnames(cor_layer)

    annotated_clusters <- annotate_registered_clusters(cor_layer, confidence_threshold = 0.25, cutoff_merge_ratio = 0.1)
    head(annotated_clusters)
    # cluster     layer_confidence layer_label
    # 1 Sp02D01             good   C.02/C.22
    # 2 Sp02D02             good        C.37
    
    # cor_layer <- 
    #   rownames(cor_layer) <- paste0(rownames(cor_layer), " ~ ", annotated_clusters[match(rownames(cor_layer), annotated_clusters$cluster)])

    rdata_name <- paste0("bayesSpace_cor_top100_", k,"_", suffix_name, ".Rdata")
    save(cor_layer, file = here(data_dir, rdata_name))
    
    ## print layer correlation plot for specific k
    # plt1 <- layer_stat_cor_plot(cor_layer)
    
    print(
      layer_stat_cor_plot(
        cor_layer, annotation = annotated_clusters,
        heatmap_legend_param = list(title = "Cor", at = c(-1, 0, 1)),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10)
      )
    )
    
    # print(plt1)
  }
  
  dev.off()

}

#plt_corr_snmultiome(suffix_reference)
#"v2"
f_name <- paste0("testing_spatial_registration_snMultiome_", suffix_name,".pdf")
plt_corr_snmultiome(suffix_reference, f_name)


message("Spatial correlation vs snRNAseq multiome data done!")



#library("slurmjobs")
## A regular job with 10 cores on the 'imaginary' partition
#job_single("02_compute_cor_visium_multiomeRnaseq", cores = 2, partition = "katun", create_shell = TRUE)



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


