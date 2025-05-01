# Load required libraries
library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("edgeR")
library("purrr")
library("IRanges")
library("S4Vectors")
# library("limma")
# library("BayesSpace")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("EnhancedVolcano")
library("dplyr")
library("here")
library("sessioninfo")

#### Set up dirs ####
# input_dir <- here("processed-data", "06_differential_expression", "01_pseudobulk_data")
input_dir <- here("processed-data", "05_brain_area_differential_expression")

#data_dir <- here("processed-data", "06_differential_expression", "03_pseudoBulkDGE")
data_dir <- here("processed-data", "05_brain_area_differential_expression")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

#plot_dir <- here("plots", "06_differential_expression", "03_pseudoBulkDGE")
plot_dir <- here("plots", "05_brain_area_differential_expression", "05_pseudobulk_DEG")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Define k values to iterate over
# k_values <- c(3, 9, 17) --- previous k(s) selected
k_values <- c(3, 13, 21, 26) # new k(s) selected based on the Habenula reference merged to only one Habenula class


##########functions for pseudobulk, saving results and volcano plots##########

## Function to perform pseudoBulkDGE and return results
run_pseudoBulkDGE <- function(data, design, coef, method) {
  #  runs quasi-likelihood F-tests using the edgeR pipeline
  de_results <- scran::pseudoBulkDGE(
    data,
    label = data$BayesSpace,
    design = design,
    coef = coef,                    # "brain_area2G2" "brain_area2G3" "brain_area2G4"
    condition = data$brain_area2,
    row.data = rowData(data),
    method = method
  )

}

# ## Function to filter and save results based on p-value < 0.05
# save_filtered_results <- function(results, model_name, output_dir) {
#   filtered_list <- list()
#
#   for (domain in names(results)) {
#     filtered_genes <- as.data.frame(results[[domain]]) |>
#       filter(!is.na(P.Value) & !is.na(logFC) & P.Value < 0.05) |>
#       select(gene_name, P.Value, logFC, gene_id, AveExpr)
#
#     filtered_list[[domain]] <- filtered_genes
#   }
#
#   combined_results <- bind_rows(filtered_list, .id = "BayesSpace_Domain")
#   write.csv(combined_results, file = file.path(output_dir, paste0("combined_", model_name, ".csv")), row.names = FALSE)
#
#   return(combined_results)
# }

## Function to create and save Enhanced Volcano plots
create_volcano_plots <- function(results, model_name, output_dir, brain_area) {
  
  # create directory for specific model  
  if (!dir.exists(model_name)) dir.create(here(plot_dir, model_name), recursive = TRUE, showWarnings = FALSE)
  # results = de_results_1
  # model_name = paste0("model1_k", k_nice)
  # output_dir = here(plot_dir, paste0("model1_k", k_nice))
  for (domain in names(results)) {
    
    # domain = "Sp13D11 ~ Habenula"
    # domain = "Sp13D01 ~ Oligo"
    domain_results <- results[[domain]]

    # replace tilde " ~ " from file-name with "-"
    f_name = paste0("volcano_", domain, "_", model_name, ".pdf")
    f_name = gsub(" ~ ", "-", f_name)
    f_name = gsub("/", "", f_name)
    
    pdf_file_path <- file.path(output_dir, f_name)
    print(pdf_file_path)
    pdf(file = pdf_file_path, width = 8, height = 8)
    
    # Keep only rows where both logFC and FDR (adj.P.Val) are not NA
    domain_results <- domain_results[!is.na(domain_results$logFC) & !is.na(domain_results$FDR), ]
    sum(is.na(domain_results$FDR))
    # Ensure adj.P.Val is truly numeric
    domain_results$FDR <- as.numeric(domain_results$FDR)
    
    lab = domain_results$gene_name
    top_genes <- head(domain_results$gene_name[order(domain_results$FDR)], 20)
    
    plot(EnhancedVolcano(domain_results,
                         lab = domain_results$gene_name,
                         selectLab = top_genes,     # Only label these
                         x = 'logFC',
                         y =  "FDR",
                         title = paste("BayesSpace cluster", domain),
                         subtitle = paste(brain_area, " - ", model_name),
                         pCutoff = 0.05,            # Adjust as needed
                         FCcutoff = 0.5,            # Adjust log2 fold change threshold
                         pointSize = 2.0,
                         labSize = 4.0,
                         drawConnectors = TRUE,   # Optional: lines from points to labels
                         widthConnectors = 0.5
                         # col=c('black', 'black', 'black', 'red3'),
                         # max.overlaps = inf        # Helps manage overcrowding
    ))
    
    dev.off()
    

  }
}

# Iterate over each k value and process the pseudoBulkDGE analysis


for (k in k_values) {
  
  # k = 13
  k_nice <- sprintf("%02d", k)  # Format k
  
  message("Processing BayesSpace k=", k_nice)
  
  # Load data for the current k
  data_file <- file.path(input_dir, paste0("sce_pseudo_PCA_brain_area_k", k, ".rds"))
  data <- readRDS(data_file)
  #colnames(colData(data))
  # [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"    
  # [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
  # [9] "nspots"          "pmi"             "rin"             "sample_id"      
  # [13] "sex"             "sum_umi"   
  #table(data$brain_area2)
  #ncol(data) 
  
  # Remove level "G0" from the data, as we only have one sample for this group
  levels(data$brain_area2)
  table(data$brain_area2, data$brain_id)
  coldata_tbl <- as_tibble(colData(data))
  filtered_coldata <- coldata_tbl |>
    filter(brain_area2 != "G0")
  # Get the matching sample names
  keep_samples <- filtered_coldata$sample_id # <- replace with your actual column name
  # Subset the SpatialExperiment object
  data <- data[, colnames(data) %in% keep_samples]
  # drop unused levels
  colData(data)$brain_area2 <- droplevels(
    colData(data)$brain_area2
  )
  levels(data$brain_area2)
  table(data$brain_area2, data$brain_id)
  # Br8518 Br9037 Br9090
  # G1     12     12     12
  # G2     12     12     11
  # G3     12     12     11
  # G4      0     12     12
  
  ######### 1: Dx Naive Model #############
  
  message("Model 1")
  message("Brain-Area: Naive Model")
  # model1 = model.matrix(~brain_area2 + sample_id + brain_id, colData(data))
  mtx_model <- model.matrix(~ brain_area2 + brain_id, colData(data))
  colnames(mtx_model)
  # [1] "(Intercept)"    "brain_area2G1"  "brain_area2G2"  "brain_area2G3" 
  # [5] "brain_area2G4"  "brain_idBr9037" "brain_idBr9090" 
  # head(mtx_model)
  #                 (Intercept) brain_area2G2 brain_area2G3 brain_area2G4
  # V13B23-280_A1           1             0             0             0
  # V13B23-280_B1           1             1             0             0
  # V13B23-280_C1           1             0             1             0
  # V13B23-280_D1           1             0             0             1
  # V13B23-285_B1           1             0             0             0
  # V13B23-285_C1           1             1             0             0
  
  ## extract the brain areas to compare
  ba_to_compare <- colnames(mtx_model)[grepl("^brain_area2G[1-4]$", colnames(mtx_model))]
  print(ba_to_compare)
  # [1] "brain_area2G2" "brain_area2G3" "brain_area2G4"
  
  ## compute pseudoBulkDGE for all brain_areas, plus build volcano plots by SpD
  
  for (ba in ba_to_compare) {
    # ba = "brain_area2G2"
    
    de_results_1 <- run_pseudoBulkDGE(data, ~ brain_area2 + brain_id, ba, "edgeR")
    ## quick inspection
    names(de_results_1)
    head(de_results_1[[1]])[c("logFC", "logCPM", "F", "PValue", "FDR")]
    map(names(de_results_1), ~ (de_results_1[[.x]][c("logFC", "logCPM", "F", "PValue", "FDR")]))
    #summary(pvals)
    map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["PValue"]])))
    map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["FDR"]])))
    
    # saveRDS(de_results_1, file = here(data_dir, paste0("de_results_1_k", k_nice, ".rds")))
    
    # create_volcano_plots(de_results_1, paste0("model1_k", k_nice), here(plot_dir, paste0("model1_k", k_nice)), "brain_area2G2")
    subdir_name <- paste0("model1_k", k_nice,"-", ba)
    create_volcano_plots(de_results_1, subdir_name, here(plot_dir, subdir_name), ba)

  }
  
}








