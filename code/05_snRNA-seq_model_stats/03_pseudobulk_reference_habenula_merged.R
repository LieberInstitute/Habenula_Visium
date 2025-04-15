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

## Access cluster assignments to merge "LHb" and "MHb" in one cluster as "Habenula"

sce_final$final_Annotations_broad <-
  gsub("^[M|L]+", "", sce_final$final_Annotations_broad)

# merged_clusters <- as.character(colData(sce_final)$final_Annotations_broad)
# merged_clusters[merged_clusters %in% c("LHb", "MHb")] <- "Habenula"
# # Add new labels to colData
# colData(sce_final)$final_Annotations_broad <- factor(merged_clusters)
# # verify results
table(colData(sce_final)$final_Annotations_broad)
# Astrocyte       Endo Excit.Thal   Habenula Inhib.Thal  Microglia      Oligo 
# 538         38       1800       2924       7612        145       2178 
# OPC 
# 1796 

registration_vars <-
  c("final_Annotations", "final_Annotations_broad")

enrichment_tstats <- lapply(registration_vars, function(current_var) {
  
  message(Sys.time(), " processing ", current_var)
  sce_pseudo <-
    registration_pseudobulk(
      sce_final,
      var_registration = current_var,
      var_sample_id = "RealSample",
      pseudobulk_rds_file = file.path(dir_rdata, paste0("pseudobulk_Hb_merged_", current_var, ".rds"))
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
  
  file_name = file.path(dir_rdata, paste0("enrichment_Hb_merged_", current_var, ".rds"))
  print(paste0('Enrichment model RDS object', file_name))
  saveRDS(results_enrichment, file = file_name) # file.path(dir_rdata, paste0("enrichment_", current_var, ".rds"))
  
  return(results_enrichment)

})

names(enrichment_tstats) <- registration_vars
colnames(head(enrichment_tstats$final_Annotations_broad))
head(enrichment_tstats$final_Annotations_broad)
#               t_stat_Astrocyte t_stat_Endo t_stat_Excit.Thal    t_stat_Hb
# ENSG00000238009       0.29608900   0.8356582         0.9129169 -1.120325695
# ENSG00000241860       0.57969314  -2.8229767         1.2448011 -1.745431694
# ENSG00000237491      -0.71530375  -0.1077206         0.6790003 -0.003513703
# ENSG00000228794      -1.82896551  -3.0418040         1.1934486  0.709096568
# ENSG00000225880       1.03412913   0.1606375        -0.1271377 -1.164305495
# ENSG00000230368      -0.02740123   0.0909449         1.6500947 -1.409937399
# t_stat_icroglia t_stat_Inhib.Thal t_stat_Oligo  t_stat_OPC
# ENSG00000238009     -0.92501220         1.1531130   -1.0366136 -0.11231046
# ENSG00000241860      0.51275865         1.5895079   -1.8167641  1.16194748
# ENSG00000237491     -3.53388899         1.7795774   -0.6940170  1.04341778
# ENSG00000228794      0.03397147         1.0317376   -0.4978867  0.43480180
# ENSG00000225880      0.43486327         0.2143562   -0.1508804  0.01707407
# ENSG00000230368     -2.09274247         0.8711746    0.7472908 -0.27985199
# p_value_Astrocyte p_value_Endo p_value_Excit.Thal p_value_Hb
# ENSG00000238009        0.76887222  0.408892478          0.3673806 0.27002349
# ENSG00000241860        0.56575063  0.007722469          0.2212751 0.08947554
# ENSG00000237491        0.47905884  0.914819017          0.5014948 0.99721595
# ENSG00000228794        0.07574142  0.004384197          0.2405261 0.48285027
# ENSG00000225880        0.30800842  0.873282385          0.8995413 0.25198925
# ENSG00000230368        0.97829183  0.928043563          0.1076529 0.16717438

## Format matches the input we need for running
## https://research.libd.org/spatialLIBD/reference/layer_stat_cor.html
modeling_results <- fetch_data(type = "modeling_results")
colnames(head(modeling_results$enrichment))
head(modeling_results$enrichment)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

