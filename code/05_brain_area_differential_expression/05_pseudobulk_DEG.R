# Load required libraries
library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("edgeR")
# library("limma")
library("BayesSpace")
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
plot_dir <- here("plots", "05_brain_area_differential_expression")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

k_values <- c(3, 9, 17)


##########functions for pseudobulk, saving results and volcano plots##########

## Function to perform pseudoBulkDGE and return results
run_pseudoBulkDGE <- function(data, design, coef, method) {
  de_results <- pseudoBulkDGE(
    data,
    label = data$BayesSpace,
    design = design,
    coef = coef,
    condition = data$brain_area_DEG,
    # condition = data$diagnosis,
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
create_volcano_plots <- function(results, model_name, output_dir) {
  # results = de_results
  for (domain in names(results)) {
    # domain_results = "Sp07D01"
    domain_results <- results[[domain]]
    pdf_file_path <- file.path(output_dir, paste0("volcano_", domain, "_", model_name, ".pdf"))
    
    pdf(file = pdf_file_path, width = 8, height = 8)
    plot(EnhancedVolcano(domain_results,
                         lab = domain_results$gene_name,
                         x = 'logFC',
                         y = 'adj.P.Val',
                         title = paste("BayesSpace cluster", domain),
                         subtitle = paste("Brain Area -", model_name)
    ))
    dev.off()
  }
}

# Iterate over each k value and process the pseudoBulkDGE analysis


for (k in k_values) {
  
  # k = 9
  k_nice <- sprintf("%02d", k)  # Format k
  
  message("Processing BayesSpace k=", k_nice)
  
  # Load data for the current k
  # data_file <- file.path(input_dir, paste0("round_test_summed_k", k_nice, ".rds"))
  data_file <- file.path(input_dir, paste0("sce_pseudo_PCA_brain_area_k", sprintf("%02d", k), ".rds"))
  data <- readRDS(data_file)
  colnames(colData(data))
  # [1] "age"            "BayesSpace"     "brain_area_DEG" "brain_id"      
  # [5] "diagnosis"      "nspots"         "sample_id"      "sex"
  
  table(data$brain_area_DEG)
  # ncol(data) 
  
  message("Brain Area Model")
  model1 = model.matrix(~brain_area_DEG + sample_id + brain_id, colData(data))
  #de_results_1 <- run_pseudoBulkDGE(data, ~ brain_area_DEG, "brain_area_DEG", "edgeR")
  #de_results_1 <- run_pseudoBulkDGE(data, ~ diagnosis + Visium_Reagent, "diagnosisAutism", "edgeR")

  # creates a design (or model) matrix, e.g., by expanding factors to a set of dummy variables
  # (depending on the contrasts) and expanding interactions similarly
  mtx_model <- model.matrix(~brain_area_DEG, colData(data))
  colnames(mtx_model)
  
  colnames(mtx_model)
  
  ## Wrapper function around edgeR's quasi-likelihood methods to conveniently perform DE analyses on pseudo-bulk 
  
  de_results <- pseudoBulkDGE(
    data,
    label = data$BayesSpace,         # specifying the cluster or cell type assignment for each column
    condition = data$brain_area_DEG, # specifying the experimental condition
    design = ~brain_area_DEG,        # represents the null hypothesis
    # design = model.matrix(~brain_area_DEG + sex + age, colData(data))
    coef="brain_area_DEGPosterior"  
  )
  
  # Genes that are filtered out will still show up in the DataFrame with all statistics set to NA
  #de_results$Sp09D04
  de_results$Sp07D05
  # DataFrame with 26610 rows and 5 columns
  # logFC    logCPM         F    PValue       FDR
  # <numeric> <numeric> <numeric> <numeric> <numeric>
  # ENSG00000243485        NA        NA        NA        NA        NA
  # ENSG00000238009        NA        NA        NA        NA        NA
  # ENSG00000241860        NA        NA        NA        NA        NA
  
  class(de_results)
  # [1] "SimpleList"
  # attr(,"package")
  # [1] "S4Vectors"
  #str(de_results)
  library("S4Vectors")

  # Remove NA values from each element in the SimpleList
  
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
  head(de_results$Sp07D05, n=3) 
  # DataFrame with 5418 rows and 5 columns
  # logFC    logCPM         F    PValue       FDR
  # <numeric> <numeric> <numeric> <numeric> <numeric>
  # ENSG00000188290  0.0296771   6.20166 0.0128958 0.9105208  0.967854
  # ENSG00000187608  0.1583498   6.36799 0.3202612 0.5766405  0.819354
  # ENSG00000078808 -0.5644011   6.83174 4.6106839 0.0419523  0.264309
  library(tidyverse)
  ## Quick inspection
  tmp_x <- as.data.frame(de_results$Sp07D05) |> filter(FDR < 0.05)
  head(tmp_x, n=3)
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








