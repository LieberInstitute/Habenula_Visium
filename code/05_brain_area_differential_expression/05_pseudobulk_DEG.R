# Load required libraries
library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("edgeR")
library("purrr")
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
  de_results <- pseudoBulkDGE(
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
    domain_results <- results[[domain]]

    # replace tilde " ~ " from file-name with "-"
    f_name = paste0("volcano_", domain, "_", model_name, ".pdf")
    f_name = gsub(" ~ ", "-", f_name)
    f_name = gsub("/", "", f_name)
    
    pdf_file_path <- file.path(output_dir, f_name)
    print(pdf_file_path)
    pdf(file = pdf_file_path, width = 8, height = 8)
    
    # remove NA's
    sum(is.na(domain_results$PValue))
    # Keep only rows where both logFC and adj.P.Val are not NA
    domain_results <- domain_results[!is.na(domain_results$logFC) & !is.na(domain_results$PValue), ]
    sum(is.na(domain_results$PValue))
    
    lab = domain_results$gene_name
    top_genes <- head(domain_results$gene_name[order(domain_results$PValue)], 20)
    
    plot(EnhancedVolcano(domain_results,
                         lab = domain_results$gene_name,
                         selectLab = top_genes,     # Only label these
                         x = 'logFC',
                         y =  "PValue",  #'adj.P.Val',
                         title = paste("BayesSpace cluster", domain),
                         subtitle = paste(brain_area, " - ", model_name),
                         # pCutoff = 0.05,          # Adjust as needed
                         # FCcutoff = 1,            # Adjust log2 fold change threshold
                         # pointSize = 2.0,
                         # labSize = 4.0,
                         # drawConnectors = TRUE,   # Optional: lines from points to labels
                         # widthConnectors = 0.5,
                         max.overlaps = inf        # Helps manage overcrowding
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
  head(mtx_model)
  ba_to_compare <- colnames(mtx_model)[grepl("^brain_area2G[1-4]$", colnames(mtx_model))]
  # [1] "brain_area2G2" "brain_area2G3" "brain_area2G4"
  # run pseudoBulkDGE for brain_area2G2
  de_results_1 <- run_pseudoBulkDGE(data, ~ brain_area2 + brain_id, "brain_area2G2", "edgeR")
  
  # fast verification of results 
  map(names(de_results_1), ~ (de_results_1[[.x]][c("logFC", "logCPM", "F", "PValue", "FDR")]))
  #pvals <- as.vector(de_results_1[["Sp13D11 ~ Habenula"]][["PValue"]])
  #summary(pvals)
  map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["PValue"]])))
  map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["FDR"]])))
  
  # saveRDS(de_results_1, file = here(data_dir, paste0("de_results_1_k", k_nice, ".rds")))
  
  
  create_volcano_plots(de_results_1, paste0("model1_k", k_nice), here(plot_dir, paste0("model1_k", k_nice)), "brain_area2G2")
  
  # Remove NA values from each element in the SimpleList
  de_results <- de_results_1
  # Function to remove NA values safely
  library(IRanges)
  library(S4Vectors)
  
  # Function to remove NA rows from each DFrame
  remove_na_dframe <- function(df) {
    if (is(df, "DFrame")) {
      return(df[complete.cases(df), , drop = FALSE])  # Remove rows with NA
    }
    return(df)  # Return unchanged if not a DFrame
  }
  
  # Apply function to each DFrame inside the SimpleDFrameList
  de_results@listData <- lapply(de_results@listData, remove_na_dframe)
  head(de_results[[1]], n=3) 
  # DataFrame with 5418 rows and 5 columns
  # logFC    logCPM         F    PValue       FDR
  # <numeric> <numeric> <numeric> <numeric> <numeric>
  # ENSG00000188290  0.0296771   6.20166 0.0128958 0.9105208  0.967854
  # ENSG00000187608  0.1583498   6.36799 0.3202612 0.5766405  0.819354
  # ENSG00000078808 -0.5644011   6.83174 4.6106839 0.0419523  0.264309
  library(tidyverse)
  ## Quick inspection
  map(names(de_results), ~ summary(as.vector(de_results[[.x]][["FDR"]])))
  df <- as.data.frame(de_results[[1]]) |> filter(FDR < 0.05)
  head(df, n=3)
  # >   as.data.frame(de_results$Sp07D05) |> filter(FDR < 0.05)
  # logFC    logCPM        F       PValue         FDR
  # ENSG00000117614  1.0465828  6.728571 17.33725 3.399336e-04 0.036112950
  # ENSG00000171812 -1.1657048  6.031604 15.00548 7.102045e-04 0.046855044
  # ENSG00000198162  1.3927120  6.737752 21.72077 9.539201e-05 0.022426378

  # Save results
  saveRDS(de_results, file = here(data_dir, paste0("DE_brain_are_results_k", k_nice, ".rds")))
  
  
  # Convert SimpleDFrameList to a single data frame to make on Volcano plot
  de_df <- do.call(rbind, de_results@listData)
  # Rename columns if needed
  colnames(de_df) <- tolower(colnames(de_df))  # Convert to lowercase for consistency
  # Ensure it is a data.frame
  de_df <- as.data.frame(de_df)
  
  # Ensure that de_df has logFC and pval (or padj for adjusted p-values):
  if (!all(c("logfc", "pval") %in% colnames(de_df))) {
    stop("Error: Required columns 'logFC' and 'pval' are missing!")
  }
  
  # Add -log10(p-value) for better visualization
  library(dplyr)
  colnames(de_df)
  de_df <- de_df |> 
    mutate(
      neg_log10_pval = -log10(pvalue),
      significance = case_when(
        pvalue < 0.05 & abs(logfc) > 1 ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  # Check column names
  colnames(de_df)
  head(de_df, n=3)
  
  volcanoPlt <- ggplot(de_df, aes(x = logfc, y = neg_log10_pval, color = significance)) +
    geom_point(alpha = 0.6, size = 2) +  # Scatter plot
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +  # LogFC threshold
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # p-value threshold
    labs(
      title = paste0("Brain-Area model - BS:", k_nice),
      x = "Log2 Fold Change (logFC)",
      y = "-log10(P-value)",
      color = "Significance"
    ) +
    theme_minimal()
  print(volcanoPlt)  
  # # create_volcano_plots <- function(results, model_name, output_dir)
  # create_volcano_plots(de_results, 
  #                      paste0("DE_brain_are_results_k", k_nice), 
  #                      #here(plot_dir, paste0("DE_brain_are_results_k", k_nice)))
  #                      here(plot_dir))
  # Error in h(simpleError(msg, call)) : 
  #   error in evaluating the argument 'x' in selecting a method for function 'plot': adj.P.Val is not numeric!

}








