# copied from https://github.com/LieberInstitute/spatialDLPFC/blob/bd93c980d7653579f81ff1c91c309cea0c7474a6/code/analysis/08_spatial_registration/

library("spatialLIBD")
library("tidyverse")
#library("xlsx")
library("jaffelab")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

## Input dir
dir_input <- here("processed-data", "05_layer_differential_expression", "modeling_results_BS")

## Set up plotting
plot_dir <- here("plots", "07_spatial_registration_vs_multiome_snRNA-seq")
data_dir <- here("processed-data", "07_spatial_registration_vs_multiome_snRNA-seq")
if (!dir.exists(plot_dir)) { dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE) }
if (!dir.exists(data_dir)) { dir.create(data_dir, showWarnings = FALSE, recursive = TRUE) }
## Load data
# load(here("processed-data", "rdata","spe", "01_build_spe", "spe_filtered_final_with_clusters.Rdata"))

### testing layer_cor on rna-multiome reference 

## Load enrichment data from multiome

rds_input <- here("processed-data", "05_snRNA-seq_model_stats", "enrichment_snRNA-multiome.rds") 

sn_multiome_data <- readRDS(rds_input)
head(sn_multiome_data$enrichment[1:4])
#                     t_stat_C.01.DD_LHb t_stat_C.02 t_stat_C.03 t_stat_C.04.DD_LHb
# ENSG00000238009          1.4486161  -1.5812290   1.9542309          1.9335732
# ENSG00000241860          1.6045860  -1.2039448   1.2864945          2.1393666
# ENSG00000237491          2.6221958  -0.7505787  -0.1465673          1.5780745
# ENSG00000228794          1.1455886  -3.3670794   0.4186114          1.0890923
# ENSG00000225880          0.6636483  -0.2710459   0.3707127          0.5807029
# ENSG00000230368          0.2200720   0.8031336   1.1624824          2.2581398

## colnames(unique(sn_multiome_data$enrichment))
# [1] "t_stat_C.01.DD_LHb"  "t_stat_C.02"         "t_stat_C.03"        
# [4] "t_stat_C.04.DD_LHb"  "t_stat_C.05.DD_LHb"  "t_stat_C.06"        
# [7] "t_stat_C.07.DD_MHb"  "t_stat_C.08.DD_LHb"  "t_stat_C.09"        
# [10] "t_stat_C.10.DD_MHb"  "t_stat_C.11.DD_MHb"  "t_stat_C.12.DD_LHb" 

## Load Registration Results 

# k_list <- c(9, 16, 28)
k_list <- c(2:28)
names(k_list) <-
  paste0("k", sprintf("%02d", k_list)) ## Use paper naming convention

bayesSpace_registration_fn <-
  map(k_list, ~ here(
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

# registration_t_stats <-
#   map(bayesSpace_registration, function(data) {
#     x <- data$enrichment
#     t_stats <- x[, grep("^t_stat_", colnames(x))]
#     colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))
#     return(t_stats)
#   })
# 
# map(registration_t_stats, jaffelab::corner)


#### Calculate Correlation Matrix ####

## get layer data
layer_modeling_results <- fetch_data(type = "modeling_results")


## CSC. Added from https://research.libd.org/spatialLIBD/articles/guide_to_spatial_registration.html#correlate-statsics-with-layer-reference
#str(registration_t_stats$k09)

#colnames(bayesSpace_registration$k09$enrichment)

## correlate registration bayes space registration vs snRNA multiome enrichment

## loop on the k(s) 
pdf(here(plot_dir, "cor_top100_spatial_registration_snMultiome.pdf"))

for (k in names(k_list)) {
  # k = "k02"
  print(k)
  bayesSpace_registration_k <- bayesSpace_registration[[k]]$enrichment
  # head(bayesSpace_registration_k)
  cor_layer <- layer_stat_cor(
    stats = bayesSpace_registration_k, #bayesSpace_registration$k09$enrichment,
    modeling_results = sn_multiome_data,
    model_type = "enrichment",
    top_n = 100
  )
  # colnames(cor_layer)
  save(cor_layer,
       file = here(data_dir, paste0("bayesSpacce_layer_cor_top100_",k,".Rdata"))
  )
  ## print layer correlation plot for specific k
  plt1 <- layer_stat_cor_plot(cor_layer)
  print(plt1)
}

dev.off()

# head(bayesSpace_registration$k09$enrichment)
# head(bayesSpace_registration_k)
#
# cor_layer <- layer_stat_cor(
#   stats = bayesSpace_registration$k09$enrichment,
#   modeling_results = sn_multiome_data,
#   model_type = "enrichment",
#   top_n = 100
# )
# layer_stat_cor_plot(cor_layer)



#### Correlate with modeling results ####
# cor_top100 <- map(
#   registration_t_stats,
#   ~ layer_stat_cor(
#     .x,
#     layer_modeling_results,
#     model_type = "enrichment",
#     reverse = FALSE,
#     top_n = 100
#   )
# )


# save(cor_top100,
#      file = here(data_dir, "bayesSpacce_layer_cor_top100.Rdata")
# )

## Plot all for portability
pdf(here(plot_dir, "cor_top100_spatial_registration.pdf"))
#map(cor_top100, layer_stat_cor_plot, max = 1) # unused argument (max = 1)
map(cor_top100, layer_stat_cor_plot)
dev.off()

