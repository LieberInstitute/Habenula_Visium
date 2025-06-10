library("SingleCellExperiment")
library("scuttle")
library("edgeR")
library("purrr")
library("ggplot2")
library("EnhancedVolcano")
#library("pheatmap")
library("dplyr")
library("here")
# library("sessioninfo")


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

plot_dir <- here("plots", "05_brain_area_differential_expression", "10_pseudobulk_DEG_linear_model")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Define k values to iterate over
# k_values <- c(3, 9, 17) --- previous k(s) selected
k_values <- c(3, 13, 21, 26) # new k(s) selected based on the Habenula reference merged to only one Habenula class


########## functions for pseudobulk, saving results and volcano plots##########

run_pseudobulk_linear_DE <- function(spe_data, assay_type = "counts") {
    # spe_data = spe_data
    # assay_type = "counts"
    message("Aggregating pseudo-bulk samples...")
    
    # Aggregate counts across pseudo-samples
    aggregated <- scuttle::aggregateAcrossCells(
        x = spe_data,
        ids = colData(spe_data)$pseudo_sample_id,
        statistics = "sum",
        use.assay.type = assay_type
    )
    
    # Extract cluster label from pseudo-sample name (e.g., SpD20_merged)
    agg_coldata <- data.frame(
        sample_id = colnames(aggregated),  # this is the full pseudo_sample_id
        cluster = sub(".*_", "", colnames(aggregated))  # extract SpD label
    )
    # Now extract true sample base for mapping metadata
    agg_coldata$base_id <- sub("^(([^_]+_[^_]+))_.*", "\\1", agg_coldata$sample_id)
    
    # Map numeric pseudo_brain_area and donor info
    map_meta <- as.data.frame(colData(spe_data)) |>
        dplyr::distinct(sample_id, .keep_all = TRUE)
    
    agg_coldata$pseudo_brain_area_numeric <- map_pseudo$pseudo_brain_area_numeric[
        match(agg_coldata$base_id, map_pseudo$sample_id)]
    
    agg_coldata$donor <- map_pseudo$brain_id[
        match(agg_coldata$base_id, map_pseudo$sample_id)]
    agg_coldata$donor <- factor(agg_coldata$donor)
    
    if (anyNA(agg_coldata$pseudo_brain_area_numeric)) {
        warning("Missing values found in pseudo_brain_area_numeric. Check mapping.")
    }
    if (anyNA(agg_coldata$donor)) {
        warning("Missing donor information for some samples.")
    }
    
    message("Donors:")
    print(table(agg_coldata$donor))
    
    # Create design matrix
    design <- model.matrix(~ pseudo_brain_area_numeric + donor, data = agg_coldata)
    message("Design matrix columns:")
    print(colnames(design))
    # Run EdgeR differential expression
    dge <- edgeR::DGEList(counts = assay(aggregated, assay_type))
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, coef = 2)  # Coef 2 = pseudo_brain_area_numeric
    
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
create_volcano_plots <- function(results, model_name, plt_dir, brain_area) {
  
  # create directory for specific model  
  if (!dir.exists(model_name)) dir.create(here(plot_dir, model_name), recursive = TRUE, showWarnings = FALSE)
  # results = return_dge
  # model_name = "K20_brainID_model_naive"
  # plt_dir = plot_dir
  # brain_area = "hb_BSk20_merged"
  
  for (domain in names(results)) {
    
    # domain = "hb_BSk20_merged"
    # domain = "Sp13D01 ~ Oligo"
    
    domain_results <- results[domain]

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
#table(spe_data$sample_id, spe_data$BayesSpace)

#===============================================================================
## merge Hb SpD(s) based on SpatialRegistration "Fine" resolution
levels(spe_data$BayesSpace)
# [1] "Sp20D01" "Sp20D02" "Sp20D03" "Sp20D04" "Sp20D06" "Sp20D07" "Sp20D08"
# [8] "Sp20D09" "Sp20D10" "Sp20D11" "Sp20D12" "Sp20D13" "Sp20D14" "Sp20D15"
# [15] "Sp20D16" "Sp20D17" "Sp20D18" "Sp20D19" "Sp20D20"

SpD20_merged <- "^Sp20D06|^Sp20D08|^Sp20D16|^Sp20D19"
spe_data$SpD20_merged <-gsub(SpD20_merged, "SpD20_Hb_merged", spe_data$BayesSpace)

# ensure all unmatched levels remain unchanged
spe_data$SpD20_merged <- ifelse(grepl(SpD20_merged, spe_data$BayesSpace),
                                "SpD20_Hb_merged", as.character(spe_data$BayesSpace))
spe_data$SpD20_merged <- factor(spe_data$SpD20_merged)
levels(spe_data$SpD20_merged)
levels(spe_data$BayesSpace)

head(table(spe_data$sample_id, spe_data$SpD20_merged))
# ...
#               SpD20_Hb_merged
# V13B23-280_A1               4
# V13B23-280_B1               4
# V13B23-280_C1               4
# V13B23-280_D1               4
# V13B23-285_A1               0
# V13B23-285_B1               3

#===============================================================================
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
# make sure all sample IDs are included
table(is.na(spe_data$pseudo_brain_area))  # should be FALSE

# Make a factor
spe_data$pseudo_brain_area <- factor(spe_data$pseudo_brain_area)
table(spe_data$BayesSpace, spe_data$pseudo_brain_area)
table(spe_data$brain_id, spe_data$pseudo_brain_area)
#         0  1  2  3  4
# Br8518 15 18 19 19  0
# Br9037  0 18 19 19 19
# Br9090  0 17 17 36  0

# creates unique factor levels for every combination of brain_id and BayesSpace and avoids accidental duplicates.
colData(spe_data)$pseudo_sample_id <- paste0(colData(spe_data)$sample_id, "_", colData(spe_data)$SpD20_merged)
head(spe_data$pseudo_sample_id)
# [1] "V13B23-280_A1_Sp20D01" "V13B23-280_B1_Sp20D01" "V13B23-280_C1_Sp20D01"
# [4] "V13B23-280_D1_Sp20D01" "V13B23-285_A1_Sp20D01" "V13B23-285_B1_Sp20D01"

# inspect uniqueness
table(duplicated(colData(spe_data)$pseudo_sample_id))
dups <- colData(spe_data)$pseudo_sample_id[duplicated(colData(spe_data)$pseudo_sample_id)]
head(table(colData(spe_data)$pseudo_sample_id))  
length(unique(colData(spe_data)$brain_id)) * length(unique(colData(spe_data)$SpD20_merged))
# = 3 × 12 = 36
# How many pseudo-bulk groups?
length(unique(colData(spe_data)$pseudo_sample_id))


#===============================================================================
## compute DGE across colData(spe_data)$pseudo_brain_area

aggregated <- scuttle::aggregateAcrossCells(
    x = spe_data,
    ids = colData(spe_data)$pseudo_sample_id,
    statistics = "sum",
    use.assay.type = "counts"
)

agg_coldata <- data.frame(
    sample_id = colnames(aggregated),  # this is the full pseudo_sample_id
    cluster = sub(".*_", "", colnames(aggregated))  # extract SpD label
)

# Now extract true sample base for mapping metadata
agg_coldata$base_id <- sub("^(([^_]+_[^_]+))_.*", "\\1", agg_coldata$sample_id)

# Add metadata (e.g. pseudo_brain_area and brain_id/donor)
map_pseudo <- as.data.frame(colData(spe_data)) |> distinct(sample_id, .keep_all = TRUE)

agg_coldata$pseudo_brain_area <- map_pseudo$pseudo_brain_area[
    match(agg_coldata$base_id, map_pseudo$sample_id)
]
agg_coldata$donor <- map_pseudo$brain_id[
    match(agg_coldata$base_id, map_pseudo$sample_id)
]
agg_coldata$donor <- factor(agg_coldata$donor)

table(agg_coldata$donor)
# Br8518 Br9037 Br9090 
# 63     63     58


#===============================================================================
## Create design matrix for contrasts specific brain areas 

# This sets up a linear model where pseudo_brain_area0 is the reference. It coefficient= 2 compares against it.
# - pseudo_brain_area is a factor with 5 levels (0 to 4)
# run EdgeR differential expression

design <- model.matrix(~ pseudo_brain_area + donor, data = agg_coldata)
colnames(design)
# [1] "(Intercept)"        "pseudo_brain_area1" "pseudo_brain_area2"
# [4] "pseudo_brain_area3" "pseudo_brain_area4"

dge <- DGEList(counts = assay(aggregated, "counts"))
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)  # area2 vs area0 

# DGE return results
return_dge <- list(
    lrt = lrt,
    top_genes = topTags(lrt),
    fit = fit,
    design = design,
    coldata = agg_coldata
)

#===============================================================================
# compute a differential gene expression (DGE) test using a single linear model (GLM) per cluster
# -> with intercept and coefficients for each level of pseudo_brain_area compared to the reference level 0
# - Subsets the DGE input data (return_dge$coldata) to include only pseudo-bulk samples belonging to that cluster 
# - highlight (filter) genes that are actually expressed in each cluster
# -> NOTE. It is not computing all pairwise contrasts, only the coefficient specified in: edgeR::glmLRT(fit_cluster, coef = 2)


plot_clusterwise_volcanos <- function(
        spe_data,
        return_dge,
        model_name,
        cluster_var = "SpD20_merged", # relative to SpD k=20 with Hb merged clusters
        expr_threshold = 1,           # only include genes with average CPM > 1 in the cluster (minimum expression level)
        FDR_thr = 0.05,
        output_dir = NULL,
        # assay_name = "counts",
        contrast_label,
        contrast_coef,                 # ge. coef=2 ~ area1 vs area0 or linear -> depends on the call
        expression_quantile = NULL
) {

    expr_cpm <- edgeR::cpm(return_dge$fit$counts)
    colnames(expr_cpm)
    # [1] "V13B23-280_A1_Sp20D01"         "V13B23-280_A1_Sp20D02"        
    # [3] "V13B23-280_A1_Sp20D03"         "V13B23-280_A1_Sp20D04"        
    # ... "V14F07-340_D1_SpD20_Hb_merged"
    
    # Add gene names
    gene_name_map <- rowData(spe_data)$gene_name
    names(gene_name_map) <- rownames(spe_data)
    
    # Map cluster annotation to return_dge$coldata; only samples mapped to the cluster (ge.Sp20D01)
    # This has one row per pseudo-bulk sample with its cluster
    map_cluster <- agg_coldata[, c("sample_id", "cluster")]
    head(map_cluster)
    #               sample_id cluster
    # 1 V13B23-280_A1_Sp20D01 Sp20D01
    # 2 V13B23-280_A1_Sp20D02 Sp20D02
    # 3 V13B23-280_A1_Sp20D03 Sp20D03  
    # Sanity check:
    identical(colnames(aggregated), map_cluster$sample_id)
    # TRUE
    
    # Direct match on full pseudo sample IDs
    return_dge$coldata$cluster <- map_cluster$cluster[
        match(return_dge$coldata$sample_id, map_cluster$sample_id)
    ]
    
    table(return_dge$coldata$cluster)
    # merged Sp20D01 Sp20D02 Sp20D03 Sp20D04 Sp20D07 Sp20D09 Sp20D10 Sp20D11 Sp20D12 
    # 11       9      12      12      12      12      12      12      11      10 
    # Sp20D13 Sp20D14 Sp20D15 Sp20D17 Sp20D18 Sp20D20 
    # 11      12      12      12      12      12 
    clusters <- unique(return_dge$coldata$cluster)
    clusters
    # [1] "Sp20D01" "Sp20D02" "Sp20D03" "Sp20D04" "Sp20D07" "Sp20D09" "Sp20D10"
    # [8] "Sp20D11" "Sp20D12" "Sp20D13" "Sp20D14" "Sp20D15" "Sp20D17" "Sp20D18"
    # [15] "Sp20D20" "merged"
    
    message("Generating volcano plots by cluster...")
    message("Processing ", length(clusters), " levels")
    
    # pdf_file <- file.path(output_dir, paste0("volcano_", model_name, "_FDR05.pdf"))
    pdf_file <- file.path(output_dir, paste0("volcano_", model_name, "_FDR", FDR_thr ,"_ExpQuantile", expression_quantile,".pdf"))
    pdf(pdf_file, width = 8, height = 8)

    for (clust in clusters) {
        
        # clust="merged"
        # clust="Sp20D01"
        message("Volcano plot for cluster: ", clust)
        
        # Identify samples in cluster
        table(return_dge$coldata$cluster)
        samples_in_cluster <- return_dge$coldata$sample_id[return_dge$coldata$cluster == clust]
        samples_in_cluster
        # [1] "V13B23-280_A1" "V13B23-280_B1" "V13B23-280_C1" "V13B23-280_D1"
        # [5] "V13B23-285_B1" "V13B23-285_C1" "V13B23-285_D1" "V14F07-340_A1"
        # [9] "V14F07-340_B1" "V14F07-340_C1" "V14F07-340_D1"

        # Filter samples in the cluster
        sample_mask <- grepl(paste(samples_in_cluster, collapse = "|"), colnames(aggregated))
        dge_cluster <- DGEList(counts = assay(aggregated, "counts")[, sample_mask])
        
        # Subset coldata + design, and drop levels that does not exist in that cluster to avoid error on estimateDisp()
        coldata_cluster <- droplevels(
            # return_dge$coldata[return_dge$coldata$sample_id %in% samples_in_cluster, ]
            return_dge$coldata[return_dge$coldata$cluster == clust, ]
        )
        # skip clusters that have fewer than 2 levels of pseudo_brain_area
        if (nlevels(coldata_cluster$pseudo_brain_area) < 2) {
            message("Skipping ", clust, ": not enough brain area levels.")
            next
        }
        
        design_cluster <- model.matrix(~ pseudo_brain_area + donor, data = coldata_cluster)
        dge_cluster <- calcNormFactors(dge_cluster)
        dge_cluster <- estimateDisp(dge_cluster, design_cluster)
        fit_cluster <- glmFit(dge_cluster, design_cluster)
        # likelihood ratio test (LRT)
        lrt_cluster <- glmLRT(fit_cluster, coef = contrast_coef)  
        
        res_cluster <- edgeR::topTags(lrt_cluster, n = Inf)$table
        res_cluster$gene_name <- gene_name_map[rownames(res_cluster)]
        
        ## Filter for low expressed genes in the current cluster
        if (is.null(expression_quantile)) {
            cluster_samples <- coldata_cluster$sample_id
            sample_ids_cpm <- colnames(expr_cpm)
            sample_mask_cpm <- grepl(paste(cluster_samples, collapse = "|"), sample_ids_cpm)
            expressed_genes <- rownames(expr_cpm)[rowMeans(expr_cpm[, sample_mask_cpm, drop = FALSE]) > expr_threshold]
        } else {
            ## filter for low expressed genes in the current cluster using CPM quantiles
            cluster_samples <- coldata_cluster$sample_id
            sample_mask_cpm <- grepl(paste(cluster_samples, collapse = "|"), colnames(expr_cpm))
            avg_expr <- rowMeans(expr_cpm[, sample_mask_cpm, drop = FALSE])
            nonzero_expr <- avg_expr[avg_expr > 0]
            # dynamically set threshold; g.e 25th / 50th percentile threshold
            cluster_expr_threshold <- quantile(nonzero_expr, probs = expression_quantile, na.rm = TRUE)
            message("Percentile threshold set (", (expression_quantile*100), "th): ", cluster_expr_threshold)
            expressed_genes <- names(avg_expr)[avg_expr > cluster_expr_threshold]
        }        
        
        res_filtered <- res_cluster[rownames(res_cluster) %in% expressed_genes, ]
        if (nrow(res_filtered) == 0) {
            message("No expressed genes found for ", clust)
            next
        }
        
        top_genes <- head(res_cluster$gene_name[order(res_cluster$FDR)], 20)
        
        p1 <- EnhancedVolcano::EnhancedVolcano(
            res_filtered,
            lab = res_filtered$gene_name,
            selectLab = top_genes,
            x = "logFC",
            y = "FDR",
            title = paste("Habenula AP", contrast_label, "- Cluster:", clust),
            subtitle = paste0(model_name, " (FDR=", FDR_thr ,"; Expr.Quantile=", expression_quantile ,")"),
            # pCutoff = 0.05, # statistical significance threshold (FDR ≤ 0.05)
            # FCcutoff = 0.5, # biological effect size threshold (log2FC > ±0.5)
            pCutoff = FDR_thr,     # even 0.2 for discovery / Default 0.05
            FCcutoff = 0.25,       # for smaller effect genes
            pointSize = 2.0,
            labSize = 4.0,
            max.overlaps = 50, 
            drawConnectors = TRUE
        )
        print(p1)
        
    }
    dev.off()
    
}


#===============================================================================

# LHb++: Column 1 (coef = 0): The intercept (baseline expression for pseudo_brain_area0)
# LHb+: Column 2 (coef = 1): log fold change of pseudo_brain_area1 relative to baseline
# LHb: Column 3 (coef = 2): logFC for pseudo_brain_area2 vs. baseline
# MHb+: Column 4 (coef = 3): logFC for pseudo_brain_area3 vs. baseline
# MHb++: Column 5 (coef = 4): logFC for pseudo_brain_area4 vs. baseline

# model1:  from LHb+ to LHb++
# - covar: donor

colnames(design)
model_name = paste0("BSk", k_merge, "_model1_AP1-0")
model_name
#[1] "BSk20_AP1-0"

#===============================================================================
## Compare pseudo_brain_area1 vs pseudo_brain_area4

plot_clusterwise_volcanos(
    spe_data = spe_data,
    return_dge = return_dge,
    model_name = model_name,
    cluster_var = "SpD20_merged",
    output_dir = plot_dir,
    FDR_thr = 0.05, 
    contrast_label = " From LHb++ to LHb+",
    contrast_coef = 2, # area1 vs area0
    expression_quantile = 0.97 
)


# model1:  from LHb+ to MHb+

model_name = paste0("BSk", k_merge, "_model1_AP3-0")
model_name
# [1] "BSk20_AP3-0"
plot_clusterwise_volcanos(
    spe_data = spe_data,
    return_dge = return_dge,
    model_name = model_name,
    cluster_var = "SpD20_merged",
    output_dir = plot_dir,
    contrast_label = "AP3-0",
    contrast_coef = 4 # area3 vs area0
)


#===============================================================================
## Compare pseudo_brain_area1 vs pseudo_brain_area3

model_name = paste0("BSk", k_merge, "_model1_AP1-3")
model_name
# [1] "BSk20_model1_AP1-3"

# Set pseudo_brain_area3 as the reference level to compare pseudo_brain_area1 vs pseudo_brain_area3
# model1:  from LHb+ to MHb+
coldata_cluster$pseudo_brain_area <- factor(coldata_cluster$pseudo_brain_area)
coldata_cluster$pseudo_brain_area <- relevel(coldata_cluster$pseudo_brain_area, ref = "3")
# re-create the design matrix
design <- model.matrix(~ pseudo_brain_area + donor, data = coldata_cluster)
colnames(design)
# [1] "(Intercept)"        "pseudo_brain_area0" "pseudo_brain_area1"
# [4] "pseudo_brain_area2" "pseudo_brain_area4" "donorBr9037"       
# [7] "donorBr9090"  

# FDR at 5%, with filtering by gene expr at 50% percentile, with relaxed FCcutoff for smaller effect genes
plot_clusterwise_volcanos(
    spe_data = spe_data,
    return_dge = return_dge,
    model_name = model_name,
    cluster_var = "SpD20_merged",
    output_dir = plot_dir,
    FDR_thr = 0.05, 
    contrast_label = " From LHb to MHb", # From LHb to MHb
    contrast_coef = 3, 
    expression_quantile = 0.50 # Relaxed; allows moderately expressed genes
)

# FDR at 5%, with filtering by gene expr at 75% percentile, with relaxed FCcutoff for smaller effect genes
plot_clusterwise_volcanos(
    spe_data = spe_data,
    return_dge = return_dge,
    model_name = model_name,
    cluster_var = "SpD20_merged",
    output_dir = plot_dir,
    FDR_thr = 0.05, 
    contrast_label = " From LHb to MHb", # From LHb to MHb
    contrast_coef = 3, 
    expression_quantile = 0.75 # More strict expressed genes
)

# FDR at 5%, with filtering by gene expr at 97% percentile, with relaxed FCcutoff for smaller effect genes
plot_clusterwise_volcanos(
    spe_data = spe_data,
    return_dge = return_dge,
    model_name = model_name,
    cluster_var = "SpD20_merged",
    output_dir = plot_dir,
    FDR_thr = 0.05, # more permissive
    contrast_label = " From LHb to MHb", # From LHb to MHb
    contrast_coef = 3, 
    expression_quantile = 0.97 # More strict, ideally should be >1
)

#===============================================================================
## FDR at 10%, with filtering by gene expr at 97% percentile, with  strict expressed genes
# plot_clusterwise_volcanos(
#     spe_data = spe_data,
#     return_dge = return_dge,
#     model_name = model_name,
#     cluster_var = "SpD20_merged",
#     output_dir = plot_dir,
#     FDR_thr = 0.1, # more permissive
#     contrast_label = " From LHb to MHb", # From LHb to MHb
#     contrast_coef = 3, 
#     expression_quantile = 0.95 # More strict expressed genes
# )


#===============================================================================
# NOTE: pseudo_brain_area is modeled as a linear numeric variable,
# enabling tests of trend along the AP axis (e.g., increasing index)
#===============================================================================

colData(spe_data)$pseudo_brain_area_numeric <- case_when(
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
# make sure all sample IDs are included
table(is.na(spe_data$pseudo_brain_area_numeric))  # should be FALSE
# Make a factor
#spe_data$pseudo_brain_area <- factor(spe_data$pseudo_brain_area_numeric)
table(spe_data$BayesSpace, spe_data$pseudo_brain_area_numeric)
table(spe_data$brain_id, spe_data$pseudo_brain_area_numeric)
# creates unique factor levels for every combination of brain_id and BayesSpace and avoids accidental duplicates.
#colData(spe_data)$pseudo_sample_id <- paste0(colData(spe_data)$sample_id, "_", colData(spe_data)$SpD20_merged)
head(spe_data$pseudo_sample_id)
# inspect uniqueness
table(duplicated(colData(spe_data)$pseudo_sample_id))
head(table(colData(spe_data)$pseudo_sample_id))  
length(unique(colData(spe_data)$brain_id)) * length(unique(colData(spe_data)$SpD20_merged))
# = 3 × 12 = 36
# How many pseudo-bulk groups?
length(unique(colData(spe_data)$pseudo_sample_id))

return_dge_linear <- run_pseudobulk_linear_DE(spe_data)

# FDR at 5%, with filtering by gene expr at 50th percentile, with relaxed FCcutoff for smaller effect genes
plot_clusterwise_volcanos(
    spe_data = spe_data,
    return_dge = return_dge_linear,
    model_name = "BSk20_AP_linear",
    cluster_var = "SpD20_merged",
    output_dir = plot_dir,
    contrast_label = "linearAP",
    contrast_coef = 2,  # pseudo_brain_area linear slope
    FDR_thr = 0.05,     # relaxed for EDA
    expression_quantile = 0.25
)








# top_genes <- rownames(topTags(return_dge$lrt, n = 50)$table)
# dge <- return_dge$fit
# counts <- cpm(dge$counts, log = TRUE)  # log2 CPM
# 
# # Subset to top genes
# heatmap_matrix <- counts[top_genes, ]
# # Z-score normalize by gene (row-wise)
# heatmap_matrix_z <- t(scale(t(heatmap_matrix)))
# 
# # Add annotation for brain_area2_numeric (for columns)
# sample_metadata <- return_dge$coldata
# rownames(sample_metadata) <- colnames(heatmap_matrix)
# 
# annotation_col <- data.frame(
#     brain_area2_numeric = sample_metadata$pseudo_brain_area
# )
# rownames(annotation_col) <- colnames(heatmap_matrix)
# 
# 
# ## prepare heatmaps with gene-expr data
# 
# pheatmap(
#     heatmap_matrix_z,
#     annotation_col = annotation_col,
#     cluster_rows = TRUE,
#     cluster_cols = TRUE,
#     show_rownames = TRUE,
#     show_colnames = FALSE,
#     fontsize_row = 6,
#     main = "Top DE Genes Heatmap"
# )
# 
# 
# # Top 50 genes (already from previous steps)
# top_genes <- rownames(edgeR::topTags(return_dge$lrt, n = 50)$table)
# 
# # LogCPM from DGEList
# dge <- return_dge$fit
# logCPM <- edgeR::cpm(dge$counts, log = TRUE)
# 
# # Subset to top genes
# heatmap_matrix <- logCPM[top_genes, ]
# 
# # Z-score normalization (gene-wise)
# heatmap_matrix_z <- t(scale(t(heatmap_matrix)))
# 
# # Column annotation
# sample_metadata <- de_results$coldata
# rownames(sample_metadata) <- colnames(heatmap_matrix_z)
# 
# annotation_col <- data.frame(
#     pseudo_brain_area = sample_metadata$pseudo_brain_area
# )
# rownames(annotation_col) <- colnames(heatmap_matrix_z)
# 
# # Plot heatmap with clustering enabled
# pheatmap(
#     heatmap_matrix_z,
#     annotation_col = annotation_col,
#     cluster_rows = TRUE,     # cluster genes
#     cluster_cols = TRUE,     # cluster samples
#     show_rownames = TRUE,
#     show_colnames = FALSE,
#     cutree_cols = 3,
#     fontsize_row = 6,
#     main = "Top 50 DE Genes (Clustered)"
# )




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








