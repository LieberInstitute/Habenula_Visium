# Load required libraries
#library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("scuttle")
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
library("pheatmap")
library("dplyr")
library("here")
library("sessioninfo")


# Note that scran::pseudoBulkDGE() does not allow you to specify which contrasts or groups to test directly in the function call when you have multiple coefficients 
#   (e.g., 4 brain areas). In short:
# As this is a wrapper that helps run DGE (differential gene expression) analysis on aggregated pseudobulk data,  
#.  when we pass a design matrix with multiple coefficients (like 4 brain areas), it:
# - By default, return all model coefficients (not contrasts).
# - It does not interpret or apply contrasts for you, and 
# - It runs a likelihood ratio test (LRT) or Wald test, depending on the method used (edgeR vs DESeq2).


#### Set up dirs ####
input_dir <- here("processed-data", "05_brain_area_differential_expression")

data_dir <- here("processed-data", "05_brain_area_differential_expression")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "05_brain_area_differential_expression", "05_pseudobulk_DEG")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Define k values to iterate over
# k_values <- c(3, 9, 17) --- previous k(s) selected
k_values <- c(3, 13, 21, 26) # new k(s) selected based on the Habenula reference merged to only one Habenula class


########## functions for pseudobulk, saving results and volcano plots##########

run_pseudobulk_linear_DE <- function(data) {

    # Aggregate 
    aggregated <- scuttle::aggregateAcrossCells(
        x = data,
        ids = colData(data)$pseudo_sample,
        statistics = "sum",
        use.assay.type = "counts"
    )
    
    # Parse pseudo-sample names into metadata - splot sample-name from bayesspace annotated domain
    sample_names <- colnames(aggregated)
    first_parts <- sub("^(([^_]+_[^_]+))_.*", "\\1", sample_names)
    second_parts <- sub("^[^_]+_[^_]+_(.*)", "\\1", sample_names)
    agg_coldata <- data.frame(cluster = second_parts, sample_id = first_parts)
    colnames(agg_coldata) <- c("cluster", "sample_id")
    
    # Map brain_area2_numeric to each pseudo-bulk sample
    map_df <- unique(colData(data)[, c("sample_id", "brain_area2_numeric")])
    agg_coldata$brain_area2_numeric <- map_df$brain_area2_numeric[
        match(agg_coldata$sample_id, map_df$sample_id)
    ]
    
    # Create design matrix
    design <- model.matrix(~ brain_area2_numeric, data = agg_coldata)
    
    # Step 5: Run EdgeR differential expression
    dge <- DGEList(counts = aggregated, assay_type)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, coef = 2)  # brain_area2_numeric
    
    # Return results
    return(list(
        lrt = lrt,
        top_genes = topTags(lrt),
        fit = fit,
        design = design,
        coldata = agg_coldata
    ))
    
}


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

## Iterate over each k value and process the pseudoBulkDGE analysis

for (k in k_values) {
  
  # k = 13
  k_nice <- sprintf("%02d", k)  # Format k
  
  message("Processing BayesSpace k=", k_nice)
  
  # Load data for the current k
  data_file <- file.path(input_dir, paste0("sce_pseudo_PCA_brain_area_k", k, ".rds"))
  spe_data <- readRDS(data_file)
  #colnames(colData(spe_data))
  #table(data$brain_area2, spe_data$brain_id)
  
  ## Compute DEG using brain-area2 as a linear predictor
  ## Encode brain_area2 numerically assuming a natural ordering (e.g., G0 < G1 < G2 < G3 < G4)

  # Create numeric variable
  colData(spe_data)$brain_area2_numeric <- as.numeric(factor(spe_data$brain_area2, 
                                                         levels = c("G0", "G1", "G2", "G3", "G4"))) - 1
  
  # Create unique pseudo-bulk sample labels
  colData(spe_data)$BayesSpace <- gsub("_$", "", colData(spe_data)$BayesSpace) # remove '-' added during ann
  colData(spe_data)$BayesSpace <- gsub("_", ".", colData(spe_data)$BayesSpace) # remove '_' added during ann
  colData(spe_data)$BayesSpace == spe_data$BayesSpace

  # creates unique factor levels for every combination of brain_id and BayesSpace and avoids accidental duplicates.
  colData(spe_data)$pseudo_sample <- paste0(colData(spe_data)$sample_id, "_", colData(spe_data)$BayesSpace)
  
  # inspect uniqueness
  table(duplicated(colData(spe_data)$pseudo_sample))
  #dups <- colData(spe_data)$pseudo_sample[duplicated(colData(spe_data)$pseudo_sample)]
    #head(table(colData(spe_data)$pseudo_sample))  
  length(unique(colData(spe_data)$brain_id)) * length(unique(colData(spe_data)$BayesSpace))
    # = 3 × 12 = 36
  # How many pseudo-bulk groups?
  length(unique(colData(spe_data)$pseudo_sample))
    
  
  # Check for duplicates (important)
  #stopifnot(!any(duplicated(colData(spe_data)$pseudo_sample)))
  
  # Optional: inspect Habenula samples
  grep("Habenula", unique(colData(spe_data)$pseudo_sample), value = TRUE)
  
  de_results <- run_pseudobulk_linear_DE(spe_data)

  names(de_results)
  #head(de_results[[1]])[c("logFC", "logCPM", "F", "PValue", "FDR")]

  res_table <- edgeR::topTags(de_results$lrt, n = Inf)$table
  
  # Create a volcano data frame
  res_table$Gene <- rownames(res_table)
  res_table$Significant <- res_table$FDR < 0.05
  
  ggplot(res_table, aes(x = logFC, y = -log10(FDR), color = Significant)) +
      geom_point(alpha = 0.8, size = 1) +
      scale_color_manual(values = c("grey", "red")) +
      theme_minimal() +
      labs(
          title = "Volcano Plot",
          x = "log2 Fold Change",
          y = "-log10(FDR)"
      ) +
      theme(legend.position = "top")
  
  # Create a volcano data frame
  res_table$Gene <- rownames(res_table)
  res_table$Significant <- res_table$PValue < 0.05
  
  ggplot(res_table, aes(x = logFC, y = -log10(PValue), color = Significant)) +
      geom_point(alpha = 0.8, size = 1) +
      scale_color_manual(values = c("grey", "red")) +
      theme_minimal() +
      labs(
          title = "Volcano Plot",
          x = "log2 Fold Change",
          y = "-log10(PValue)"
      ) +
      theme(legend.position = "top")
  
  # saveRDS(de_results_1, file = here(data_dir, paste0("de_results_1_k", k_nice, ".rds")))
  
  # Top 50 DE genes by FDR
  top_genes <- rownames(topTags(de_results$lrt, n = 50)$table)
  dge <- de_results$fit
  counts <- cpm(dge$counts, log = TRUE)  # log2 CPM
  
  # Subset to top genes
  heatmap_matrix <- counts[top_genes, ]
  # Z-score normalize by gene (row-wise)
  heatmap_matrix_z <- t(scale(t(heatmap_matrix)))
  
  # Add annotation for brain_area2_numeric (for columns)
  sample_metadata <- de_results$coldata
  rownames(sample_metadata) <- colnames(heatmap_matrix)
  
  annotation_col <- data.frame(
      brain_area2_numeric = sample_metadata$brain_area2_numeric
  )
  rownames(annotation_col) <- colnames(heatmap_matrix)
  

  ## prepare heatmaps with gene-expr data
  
  pheatmap(
      heatmap_matrix_z,
      annotation_col = annotation_col,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = FALSE,
      fontsize_row = 6,
      main = "Top DE Genes Heatmap"
  )
  
  
  # Top 50 genes (already from previous steps)
  top_genes <- rownames(edgeR::topTags(de_results$lrt, n = 50)$table)
  
  # LogCPM from DGEList
  dge <- de_results$fit
  logCPM <- edgeR::cpm(dge$counts, log = TRUE)
  
  # Subset to top genes
  heatmap_matrix <- logCPM[top_genes, ]
  
  # Z-score normalization (gene-wise)
  heatmap_matrix_z <- t(scale(t(heatmap_matrix)))
  
  # Column annotation
  sample_metadata <- de_results$coldata
  rownames(sample_metadata) <- colnames(heatmap_matrix_z)
  
  annotation_col <- data.frame(
      brain_area2 = sample_metadata$brain_area2_numeric
  )
  rownames(annotation_col) <- colnames(heatmap_matrix_z)
  
  # Plot heatmap with clustering enabled
  pheatmap(
      heatmap_matrix_z,
      annotation_col = annotation_col,
      cluster_rows = TRUE,     # cluster genes
      cluster_cols = TRUE,     # cluster samples
      show_rownames = TRUE,
      show_colnames = FALSE,
      cutree_cols = 3,
      fontsize_row = 6,
      main = "Top 50 DE Genes (Clustered)"
  )
  
  
}

  

#============= Alternative analysis ============================================
#============= Prepare data to compute DGE with linear model on BayesSpace k=20 

k_merge = 20
k_merge <- sprintf("%02d", k_merge)

message("Processing BayesSpace k=", k_merge)

## Load data for the current k
data_file <- file.path(input_dir, paste0("sce_pseudo_PCA_brain_area_k", k_merge, ".rds"))
spe_data <- readRDS(data_file)

## inspect data
colnames(colData(spe_data))
# [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"    
# [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
# [9] "nspots"          "pmi"             "rin"             "sample_id"      
# [13] "sex"             "sum_umi" 
table(spe_data$BayesSpace)
# Sp20D01 Sp20D02 Sp20D03 Sp20D04 Sp20D06 Sp20D07 Sp20D08 Sp20D09 Sp20D10 Sp20D11 
#   9      12      12      12      10      12      11      12      12      11 
# Sp20D12 Sp20D13 Sp20D14 Sp20D15 Sp20D16 Sp20D17 Sp20D18 Sp20D19 Sp20D20 
#   10      11      12      12      11      12      12      11      12 
table(spe_data$brain_id, spe_data$BayesSpace)
#           Sp20D18 Sp20D19 Sp20D20
# Br8518       4       3       4
# Br9037       4       4       4
# Br9090       4       4       4
table(spe_data$sample_id, spe_data$BayesSpace)

## merge Hb SpD(s) based on SpatialRegistration "Fine" resolution
levels(spe_data$BayesSpace)
# [1] "Sp20D01" "Sp20D02" "Sp20D03" "Sp20D04" "Sp20D06" "Sp20D07" "Sp20D08"
# [8] "Sp20D09" "Sp20D10" "Sp20D11" "Sp20D12" "Sp20D13" "Sp20D14" "Sp20D15"
# [15] "Sp20D16" "Sp20D17" "Sp20D18" "Sp20D19" "Sp20D20"

k20_SpD_to_merge <- "^Sp20D06|^Sp20D08|^Sp20D16|^Sp20D19"

spe_data$hb_BSk20_merged <-
    gsub(k20_SpD_to_merge, "hb_BSk20_merged", spe_data$BayesSpace)
spe_data$hb_BSk20_merged <- as.factor(spe_data$hb_BSk20_merged)
levels(spe_data$hb_BSk20_merged)
#table(spe_data$sample_id, spe_data$BayesSpace)
#table(spe_data$sample_id, spe_data$hb_BSk20_merged)

## create new variable describing 'pseudo_brain_area' to test Habenula A-P DGE
# Br8518 = V13B23-285 defined
# Br9037 = V13B23-280
# Br9090 = V14F07-340

colData(spe_data)$pseudo_brain_area <- case_when(
    colData(spe_data)$sample_id == "V13B23-285_A1" ~ 0,
    
    colData(spe_data)$sample_id == "V13B23-285_B1" |
        colData(spe_data)$sample_id == "V14F07-340_D1" |
        colData(spe_data)$sample_id == "V13B23-280_D1" ~ 1,
    
    colData(spe_data)$sample_id == "V13B23-285_C1" |
        colData(spe_data)$sample_id == "V14F07-340_C1" |
        colData(spe_data)$sample_id == "V13B23-280_C1" ~ 2,
    
    colData(spe_data)$sample_id == "V13B23-285_D1" |
        colData(spe_data)$sample_id == "V14F07-340_B1" |
        colData(spe_data)$sample_id == "V14F07-340_A1" |
        colData(spe_data)$sample_id == "V13B23-280_B1" ~ 3,
    
    colData(spe_data)$sample_id == "V13B23-280_A1" ~ 4
)
# Make a factor
spe_data$pseudo_brain_area <- factor(spe_data$pseudo_brain_area)
table(spe_data$BayesSpace, spe_data$pseudo_brain_area)
table(spe_data$hb_BSk20_merged, spe_data$pseudo_brain_area)
table(spe_data$brain_id, spe_data$pseudo_brain_area)
#         0  1  2  3  4
# Br8518 15 18 19 19  0
# Br9037  0 18 19 19 19
# Br9090  0 17 17 36  0

#===============================================================================



  #### Compute model against brain_area linear variable == this model doesn't have the correct desing, need to use manual DEG using contrast
  
  # # Remove level "G0" from the data, as we only have one sample for this group
  # levels(data$brain_area2)
  # table(data$brain_area2, data$brain_id)
  # coldata_tbl <- as_tibble(colData(data))
  # filtered_coldata <- coldata_tbl |>
  #   filter(brain_area2 != "G0")
  # # Get the matching sample names
  # keep_samples <- filtered_coldata$sample_id # <- replace with your actual column name
  # # Subset the SpatialExperiment object
  # data <- data[, colnames(data) %in% keep_samples]
  # # drop unused levels
  # colData(data)$brain_area2 <- droplevels(
  #   colData(data)$brain_area2
  # )
  # levels(data$brain_area2)
  # table(data$brain_area2, data$brain_id)
  # # Br8518 Br9037 Br9090
  # # G1     12     12     12
  # # G2     12     12     11
  # # G3     12     12     11
  # # G4      0     12     12
  # 
  # ######### 1: Brain Naive Model #############
  # 
  # message("Model 1")
  # message("Brain-Area: Naive Model")
  # mtx_model <- model.matrix(~ brain_area2 + brain_id, colData(data))
  # colnames(mtx_model)
  # # [1] "(Intercept)"    "brain_area2G1"  "brain_area2G2"  "brain_area2G3" 
  # # [5] "brain_area2G4"  "brain_idBr9037" "brain_idBr9090" 
  # # head(mtx_model)
  # #                 (Intercept) brain_area2G2 brain_area2G3 brain_area2G4
  # # V13B23-280_A1           1             0             0             0
  # # V13B23-280_B1           1             1             0             0
  # # V13B23-280_C1           1             0             1             0
  # 
  # ## extract the brain areas to compare
  # ba_to_compare <- colnames(mtx_model)[grepl("^brain_area2G[1-4]$", colnames(mtx_model))]
  # print(ba_to_compare)
  # # [1] "brain_area2G2" "brain_area2G3" "brain_area2G4"
  # 
  # ## compute pseudoBulkDGE for all brain_areas, plus build volcano plots by SpD
  # 
  # for (ba in ba_to_compare) {
  #   # ba = "brain_area2G2"
  #   
  #   de_results_1 <- run_pseudoBulkDGE(data, ~ brain_area2 + brain_id, ba, "edgeR")
  #   names(de_results_1)
  #   head(de_results_1[[1]])[c("logFC", "logCPM", "F", "PValue", "FDR")]
  # 
  #   # Quick verification
  #   map(names(de_results_1), ~ (de_results_1[[.x]][c("logFC", "logCPM", "F", "PValue", "FDR")]))
  #   map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["PValue"]])))
  #   map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["FDR"]])))
  #   
  #   # saveRDS(de_results_1, file = here(data_dir, paste0("de_results_1_k", k_nice, ".rds")))
  #   
  #   # create_volcano_plots(de_results_1, paste0("model1_k", k_nice), here(plot_dir, paste0("model1_k", k_nice)), "brain_area2G2")
  #   subdir_name <- paste0("model1_k", k_nice,"-", ba)
  #   create_volcano_plots(de_results_1, subdir_name, here(plot_dir, subdir_name), ba)
  # 
  # }








