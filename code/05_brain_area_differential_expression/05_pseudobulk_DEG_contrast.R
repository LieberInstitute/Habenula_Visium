# Load required libraries
library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("edgeR")
library("scuttle")
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


## In pseudobulk DE using tools like edgeR or scran::pseudoBulkDGE(), the PValue and FDR (False Discovery Rate) are calculated:
#  Across all genes, per contrast, not per cluster
#  Here, I plot the results to build a base line and continue polishing the script


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


create_volcano_plots <- function(sce,
                                 contrast_grps, 
                                 output_dir,
                                 model_name) {
  
  ## Plot all contrast computed for specific Bayes-Space k
  
  # for testing
  # sce = data
  # contrast_grps = dge_results_df
  # model_name = paste0("model1_k", k_nice)
  
  # customize subdir name
  output_subdir = here(output_dir, paste0(model_name, "_contrast"))
  if (!dir.exists(output_subdir)) dir.create(output_subdir, , recursive = TRUE, showWarnings = FALSE)
  
  #for (cg in seq_along(contrast_grps)[-length(contrast_grps)]) {
  for (cg in seq_along(contrast_grps)) {
    
    # First contrast grp
    # cg = 1 
    cg_name <- gsub(" ", "", unique(contrast_grps[[cg]]["contrast"]))
    
    print(cg)
    message("Processing contrast group: ", cg_name)
    # Processing contrast group: G1-G2
    
    contrast_grp_results <- contrast_grps[[cg]] 
    
    ## Add gene names to your DE results. Create a named vector and map the gene_names
    gene_name_map <- rowData(data)$gene_name
    names(gene_name_map) <- rownames(data)
    contrast_grp_results$gene_name <- gene_name_map[rownames(contrast_grp_results)]
    
    # cutomize file name
    f_name = paste0("volcano_", cg_name, "_", model_name, ".pdf")
    # volcano_G1-G2_model1_k13.pdf
    
    pdf_file_path <- file.path(output_subdir, f_name)
    print(pdf_file_path)
    pdf(file = pdf_file_path, width = 8, height = 8)
    
    # Keep only rows where both logFC and FDR (adj.P.Val) are not NA
    contrast_grp_results <- contrast_grp_results[!is.na(contrast_grp_results$logFC) & !is.na(contrast_grp_results$FDR), ]
    sum(is.na(contrast_grp_results$FDR))
    # Ensure adj.P.Val is truly numeric
    contrast_grp_results$FDR <- as.numeric(contrast_grp_results$FDR)
    
    lab = contrast_grp_results$gene_name
    top_genes <- head(contrast_grp_results$gene_name[order(contrast_grp_results$FDR)], 20)
    
    plot(EnhancedVolcano(contrast_grp_results,
                         lab = contrast_grp_results$gene_name,
                         selectLab = top_genes,     # Only label these
                         x = 'logFC',
                         y =  "FDR",
                         title = paste("Contrast grp: ", cg_name),
                         subtitle = model_name,
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


## Function to compare multiple contrast
runPseudobulkEdgeR <- function(
    sce,                 # pseudobulk data
    group_var,           # column name in colData(sce) to use for group comparison (e.g. "brain_area2")
    contrast_levels,     # vector of two levels to contrast (e.g., c("G2", "G1"))
    covariates = c("brain_id"),
    assay_type = "counts",
    min_genes = 50
) {
  
  # Step 1: Aggregate counts across cells - We already have pseudobulk data, so we ommit this step
  counts <- assay(sce, assay_type)
  col_data <- colData(sce)
  
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
  # '~ 0' is used to specify NO Intercept; explicit group means without a baseline group (e.g., for clean contrasts)
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

## In pseudobulk DE using tools like edgeR or scran::pseudoBulkDGE(), the PValue and FDR (False Discovery Rate) are calculated:
#  Across all genes, per contrast, not per cluster
#  Here, I plot the results to build a base line and continue polishing the script

## Iterate over each k value and process the pseudoBulkDGE analysis

for (k in k_values) {
  
  # k = 13
  k_nice <- sprintf("%02d", k)  # Format k
  
  message("Processing BayesSpace k=", k_nice)
  
  ## Load data for the current k
  data_file <- file.path(input_dir, paste0("sce_pseudo_PCA_brain_area_k", k, ".rds"))
  data <- readRDS(data_file)
  #colnames(colData(data))
  # [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"    
  # [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
  # [9] "nspots"          "pmi"             "rin"             "sample_id"      
  # [13] "sex"             "sum_umi"   
  #table(data$brain_area2)
  #ncol(data) 
  #rowData(data)$gene_name  # or $symbol, or $GeneSymbol depending on source
  
  ## Remove level "G0" from the data, as we only have one sample for this group
  
  levels(data$brain_area2)
  table(data$brain_area2, data$brain_id)
  # Br8518 Br9037 Br9090
  # G0     11      0      0
  # G1     12     12     12
  # G2     12     12     11
  # G3     12     12     11
  # G4      0     12     12
  
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

  
  data <- data[, data$BayesSpace == habenula_levels]
  # Drop unused levels in brain_area2
  data$brain_area2 <- droplevels(data$brain_area2)
  # Optional check brain-area grp in the subseted data
  table(data$brain_area2)
  # G1 G2 G3 G4 
  # 3  3  3  2 
  
  ## The new scran::pseudoBulkDGE() no longer lets you run the full model without coef — which means it doesn’t expose the fit object needed for custom contrasts.
  ## So, I use the fitted model to manually define and test contrasts using the DGE package
  ## I have to manually aggregate pseudobulk counts and then run the model using standard edgeR workflow
  ## Note: scran::pseudoBulkDGE() works in two modes: (1) If provided a design matrix, so I DO NOT need `coef` 
  ##      and I can get `dge_result$fit`, `dge_result$design`, and run custom contrasts for the brain-areas
  
  
  ## Get all pairwise combinations of group levels
  group_levels <- levels(factor(colData(data)$brain_area2))
  # [1] "G1" "G2" "G3" "G4"
  pairwise_contrasts <- combn(group_levels, 2, simplify = FALSE)
  head(pairwise_contrasts)
  
  ######### Model 1: Brain-Area Naive Model #############
  
  message("Model 1")
  message("Brain-Area: Naive Model")
  
  # Compute DE results across every pair of group comparisons (e.g. G2 vs G1, G3 vs G1, G3 vs G2, etc.)
  
  ## Store all results in a list
  dge_results_list <- lapply(pairwise_contrasts, function(contrast_pair) {
    
    cat("Running contrast:", paste(contrast_pair, collapse = " vs "), "\n")
    
    res <- runPseudobulkEdgeR(
      sce = data,
      group_var = "brain_area2",
      contrast_levels = contrast_pair
    )

    return(res)
    
  })
  
  ######### Model 2: Brain-Area + 1 co-variable #############
  
  #use coef = "brain_area2G2:sexM" or create a contrast.
  #design = ~ brain_area2 * brain_id
  
  # e.g. First contrast group
  #head(dge_results_list[[1]])
  #                       logFC   logCPM         F      PValue       FDR
  # ENSG00000074657 -0.3563410 5.786637 10.001880 0.001884399 0.9999613
  # ENSG00000157593 -0.5712547 5.386361  9.274807 0.002712389 0.9999613
  # ENSG00000100225 -0.3704179 6.019549  8.567266 0.003958746 0.9999613
  # ENSG00000183513 -0.4371085 5.865209  8.357765 0.004402942 0.9999613
  # ENSG00000175265 -0.4104273 5.624887  7.529694 0.006773970 0.9999613
  # ENSG00000114805  0.6765182 4.774126  6.885194 0.009451898 0.9999613
  # gene contrast
  # ENSG00000074657 ENSG00000074657  G1 - G2
  # ENSG00000157593 ENSG00000157593  G1 - G2
  # ENSG00000100225 ENSG00000100225  G1 - G2
  # ENSG00000183513 ENSG00000183513  G1 - G2
  # ENSG00000175265 ENSG00000175265  G1 - G2
  # ENSG00000114805 ENSG00000114805  G1 - G2
  
  # Combine all into one data.frame
  
  dge_results_df <- dge_results_list # testing
  # dge_results_df <- bind_rows(dge_results_list)
  # head(dge_results_df)
  # table(dge_results_df$contrast)
  # # G1 - G2 G1 - G3 G1 - G4 G2 - G3 G2 - G4 G3 - G4 
  # # 6352    6352    6352    6352    6352    6352 
  
  # ## Filter data to prepare for plots 
  # subset(dge_results_df, FDR < 0.05 & abs(logFC) > 1)

  # saveRDS(de_results_1, file = here(data_dir, paste0("de_results_1_k", k_nice, ".rds")))

  ## all contrast over all clusters
  #create_volcano_plots(data, dge_results_df, plot_dir, paste0("model1_k", k_nice))
  ## all contrast over all clusters only in the Hb domain(s) 
  create_volcano_plots(data, dge_results_df, plot_dir, paste0("model1_k", k_nice, "_HbSpD_ba_bID"))
  
}








