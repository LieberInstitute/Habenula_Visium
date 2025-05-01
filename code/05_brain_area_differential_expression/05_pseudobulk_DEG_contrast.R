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


## Function to create and save Enhanced Volcano plots
create_volcano_plots <- function(results, 
                                 model_name,
                                 output_dir,
                                 brain_area) {
  
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


# Function to apply a specific contrast
runPseudobulkEdgeR <- function(
    pb,                  # pseudobulk data
    group_var,           # column name in colData(sce) to use for group comparison (e.g. "brain_area2")
    contrast_levels,     # vector of two levels to contrast (e.g., c("G2", "G1"))
    min_genes = 10
) {
  library(scuttle)
  library(edgeR)
  
  # Step 1: Aggregate counts across cells - We already have pseudobulk data
  col_data <- colData(pb)
  
  # Step 2: Create DGEList and filter
  y <- DGEList(counts = counts)
  
  # Drop low-expression genes
  keep <- filterByExpr(y, group = col_data[[group_var]])
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  if (nrow(y) < min_genes) {
    stop("Too few genes left after filtering.")
  }
  
  # Normalize
  y <- calcNormFactors(y)
  
  # Step 3: Design matrix with no intercept
  design <- model.matrix(~ 0 + col_data[[group_var]])
  group_names <- levels(factor(col_data[[group_var]]))
  colnames(design) <- group_names
  
  # Step 4: Fit model
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  
  # Step 5: Define contrast
  stopifnot(all(contrast_levels %in% group_names))
  contrast_str <- paste0(contrast_levels[1], " - ", contrast_levels[2])
  contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)
  
  # Step 6: QL test
  qlf <- glmQLFTest(fit, contrast = contrast_mat)
  
  # Step 7: Return results
  result <- topTags(qlf, n = Inf)$table
  result$gene <- rownames(result)
  result$contrast <- contrast_str
  return(result)
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
  
  
  ## The new scran::pseudoBulkDGE() no longer lets you run the full model without coef — which means it doesn’t expose the fit object needed for custom contrasts.
  ## So, I use the fitted model to manually define and test contrasts using the DGE package
  ## I have to manually aggregate pseudobulk counts and then run the model using standard edgeR workflow
  ## Note: scran::pseudoBulkDGE() works in two modes: (1) If provided a design matrix, so I DO NOT need `coef` 
  ##      and I can get `dge_result$fit`, `dge_result$design`, and run custom contrasts for the brain-areas
  
  
  # for (ba in ba_to_compare) {
  #   # ba = "brain_area2G2"
    
    ######### Model 1: Brain Naive Model #############
    
    message("Model 1")
    message("Brain-Area: Naive Model")

    res_G2vsG1 <- runPseudobulkEdgeR(
      pb = data,
      group_var = "brain_area2",
      contrast_levels = c("G2", "G1")
    )
    
    head(res_G2vsG1)
    #                     logFC   logCPM         F      PValue       FDR
    # ENSG00000074657  0.3563410 5.786637 10.001880 0.001884399 0.9999613
    # ENSG00000157593  0.5712547 5.386361  9.274807 0.002712389 0.9999613
    # ENSG00000100225  0.3704179 6.019549  8.567266 0.003958746 0.9999613
    # ENSG00000183513  0.4371085 5.865209  8.357765 0.004402942 0.9999613
    # ENSG00000175265  0.4104273 5.624887  7.529694 0.006773970 0.9999613
    # ENSG00000114805 -0.6765182 4.774126  6.885194 0.009451898 0.9999613
    #                               gene contrast
    # ENSG00000074657 ENSG00000074657  G2 - G1
    # ENSG00000157593 ENSG00000157593  G2 - G1
    # ENSG00000100225 ENSG00000100225  G2 - G1
    # ENSG00000183513 ENSG00000183513  G2 - G1
    # ENSG00000175265 ENSG00000175265  G2 - G1
    # ENSG00000114805 ENSG00000114805  G2 - G1
    
    # ## (2) Prepare edgeR object
    # library(edgeR)
    # counts <- assay(data, "counts")
    # col_data <- colData(data)
    # 
    # # Step 3a: Create edgeR object and design matrix
    # y <- DGEList(counts = counts)
    # # Add gene filtering (optional but recommended)
    # keep <- filterByExpr(y, group = col_data$brain_area2)
    # y <- y[keep, , keep.lib.sizes = FALSE]
    # # Normalize
    # y <- calcNormFactors(y)
    # # Build design matrix
    # design <- model.matrix(~ 0 + brain_area2, data = col_data)
    # colnames(design) <- gsub("brain_area2", "", colnames(design))  # Cleaner names like G1, G2, G3, etc.
    # colnames(design)    
    # head(design)
    # 
    # # Estimate dispersions
    # y <- estimateDisp(y, design)
    # # Fit negative binomial GLM
    # fit <- glmQLFit(y, design)
    # 
    # # Define contrast
    # contrast_matrix <- makeContrasts(G2vsG1 = G2 - G1, levels = design)
    # # Perform QL test
    # qlf <- glmQLFTest(fit, contrast = contrast_matrix)

    
    # ## quick inspection
    # names(de_results_1)
    # head(de_results_1[[1]])[c("logFC", "logCPM", "F", "PValue", "FDR")]
    # map(names(de_results_1), ~ (de_results_1[[.x]][c("logFC", "logCPM", "F", "PValue", "FDR")]))
    # #summary(pvals)
    # map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["PValue"]])))
    # map(names(de_results_1), ~ summary(as.vector(de_results_1[[.x]][["FDR"]])))


    # saveRDS(de_results_1, file = here(data_dir, paste0("de_results_1_k", k_nice, ".rds")))
    
    subdir_name <- paste0("model1_k", k_nice,"-", ba)
    create_volcano_plots(de_results_1, subdir_name, here(plot_dir, subdir_name), ba)

  # }
  
}








