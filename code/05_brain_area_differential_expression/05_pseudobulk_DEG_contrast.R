# Load required libraries
library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("edgeR")
library("scuttle")
library("purrr")
library("IRanges")
library("S4Vectors")
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
input_dir <- here("processed-data", "05_brain_area_differential_expression")
data_dir <- here("processed-data", "05_brain_area_differential_expression")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
plot_dir <- here("plots", "05_brain_area_differential_expression", "05_pseudobulk_DEG_contrast")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# # Define k values to iterate over
# k_values <- c(3, 9, 17) --- previous k(s) selected
# k_values <- c(3, 13, 21, 26) # new k(s) selected based on the Habenula reference merged to only one Habenula class

args = commandArgs(trailingOnly = TRUE)
k <- as.integer(args[2])

## For testing
if (is.na(k)) {
    k <- 2
}

## Function plot volcanos by FDR to compare multiple contrast
create_volcano_plots <- function(sce_p,
                                 contrast_grps,
                                 output_dir,
                                 model_name) {

  ## Plot all contrast computed for specific Bayes-Space k - By FDR

  # for testing
  # sce = spe_pseudo
  # contrast_grps = dge_results_df
  # model_name = paste0("model1_k", k_nice)

    f_name = paste0("volcano_FDR_", model_name, ".pdf")
    pdf_file_path <- file.path(output_dir, f_name)
    print(pdf_file_path)
    pdf(file = pdf_file_path, width = 8, height = 8)

  for (cg in seq_along(contrast_grps)) {

    # First contrast grp
    # cg = 1
    cg_name <- gsub(" ", "", unique(contrast_grps[[cg]]["contrast"]))

    print(cg)
    message("Processing contrast group: ", cg_name)
    # Processing contrast group: G1-G2

    contrast_grp_results <- contrast_grps[[cg]]

    ## Add gene names to your DE results. Create a named vector and map the gene_names
    gene_name_map <- rowData(sce_p)$gene_name
    names(gene_name_map) <- rownames(sce_p)
    contrast_grp_results$gene_name <- gene_name_map[rownames(contrast_grp_results)]

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

  }
  dev.off()

}

## Function plot volcanos by PValue just to check flat-volcanos
create_volcano_plots_pValue <- function(sce_p,
                                        contrast_grps,
                                        output_dir,
                                        model_name) {

    ## Plot all contrast computed for specific Bayes-Space k
    ## plot Volcano with PValues just to check

    # Testing:
    # sce = spe_pseudo
    # contrast_grps = dge_results_df
    # model_name = paste0("model1_k", k_nice)

    f_name = paste0("volcano_PValue_", model_name, ".pdf")
    pdf_file_path <- file.path(output_dir, f_name)
    print(pdf_file_path)
    pdf(file = pdf_file_path, width = 8, height = 8)

    for (cg in seq_along(contrast_grps)) {

        # First contrast grp
        # cg = 1
        cg_name <- gsub(" ", "", unique(contrast_grps[[cg]]["contrast"]))

        print(cg)
        message("Processing contrast group: ", cg_name)
        # Processing contrast group: G1-G2

        contrast_grp_results <- contrast_grps[[cg]]

        ## Add gene names to your DE results. Create a named vector and map the gene_names
        gene_name_map <- rowData(sce_p)$gene_name
        names(gene_name_map) <- rownames(sce_p)
        contrast_grp_results$gene_name <- gene_name_map[rownames(contrast_grp_results)]

        # Keep only rows where both logFC and FDR (adj.P.Val) are not NA
        contrast_grp_results <- contrast_grp_results[!is.na(contrast_grp_results$logFC) & !is.na(contrast_grp_results$PValue), ]
        sum(is.na(contrast_grp_results$PValue))
        # Ensure adj.P.Val is truly numeric
        contrast_grp_results$PValue <- as.numeric(contrast_grp_results$PValue)

        lab = contrast_grp_results$gene_name
        top_genes <- head(contrast_grp_results$gene_name[order(contrast_grp_results$PValue)], 20)

        plot(EnhancedVolcano(contrast_grp_results,
                             lab = contrast_grp_results$gene_name,
                             selectLab = top_genes,     # Only label these
                             x = 'logFC',
                             y =  "PValue",
                             title = paste("Contrast grp: ", cg_name),
                             subtitle = model_name,
                             pCutoff = 0.05,            # Adjust as needed
                             FCcutoff = 0.5,            # Adjust log2 fold change threshold
                             pointSize = 2.0,
                             labSize = 4.0,
                             drawConnectors = TRUE,   # Optional: lines from points to labels
                             widthConnectors = 0.5
        ))

    }

    dev.off()

}

## Function to compare multiple contrast
runPseudobulkEdgeR <- function(
    sce_p,               # pseudobulk data
    group_var,           # column name in colData(sce) to use for group comparison (e.g. "brain_area2")
    contrast_levels,     # vector of two levels to contrast (e.g., c("G2", "G1"))
    covariates = covars, # c("sample_id", "brain_id"),
    model_name = model_name,
    min_genes = 200
) {

    # testing:
    # sce_p = spe_pseudo
    # group_var = "brain_area2"
    # contrast_levels = pairwise_contrasts[[1]]

  # I already have pseudobulk data, so I just pull the assay. I use counts, because edgeR and limma-trend/voom expect integer counts
  counts <- assay(sce_p, "counts")
  col_data <- colData(sce_p)

  # Create DGEList and filter
  # samples = colData(current) argument is optional. Provided as sample-level metadata (e.g., brain_area, brain_id, cluster, etc.)
  sample_metadata <- as.data.frame(colData(sce_p))
  rownames(sample_metadata) <- colnames(sce_p)
  # verify that data match
  stopifnot(all(rownames(sample_metadata) == colnames(counts)))

  # y <- DGEList(counts = counts)
  y <- DGEList(counts = counts, samples = sample_metadata)

  # Drop low-expression genes that are not expressed above a log-CPM threshold in a minimum number of samples
  keep <- filterByExpr(y, group = col_data[[group_var]])
  message("Genes retained after filtering:")
  print(summary(keep))
  # Subset DGEList to expressed genes + recalc lib size
  y <- y[keep, , keep.lib.sizes = FALSE]

  if (nrow(y) < min_genes) {
    stop("Too few genes left after filtering.")
  }

  # Normalize
  y <- calcNormFactors(y)

  # # Plot mean-difference (MD) plot for each normalized pseudo-bulk profile
  # # It should exhibit a trumpet shape centered at zero indicating that the normalization successfully removed systematic bias between profiles
  # f_name <- here(plt_subdir, paste0(" mean-difference-", model_name, ".pdf"))
  # pdf(f_name, width = 7, height = 5)
  # par(mfrow=c(2,3))
  # for (i in seq_len(ncol(y))) {
  #     plotMD(y, column=i)
  # }
  # dev.off()
  # message("Plot MD done")
  
  # Design matrix with no intercept
  # '~ 0' is used to specify NO Intercept; explicit group means without a baseline group (e.g., for clean contrasts)
  #design <- model.matrix(~ 0 + col_data[[group_var]])
  design <- model.matrix(~ 0 + col_data[[group_var]], y$sample_id)
  group_names <- levels(factor(col_data[[group_var]]))
  colnames(design) <- group_names

  # Fit model. Estimate 'negative binomial (NB) dispersions' and 'quasi-likelihood dispersions' 
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)

  # Define contrast
  stopifnot(all(contrast_levels %in% group_names))
  contrast_str <- paste0(contrast_levels[1], " - ", contrast_levels[2]) # g.e. "G1 - G2"
  contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)

  # QL test. Test for differences in expression
  qlf <- glmQLFTest(fit, contrast = contrast_mat)

  # Return results
  result <- topTags(qlf, n = Inf)$table
  result$gene <- rownames(result)
  result$contrast <- contrast_str

  return(result)

}

## In pseudobulk DE using tools like edgeR or scran::pseudoBulkDGE(), the PValue and FDR (False Discovery Rate) are calculated:
#  Across all genes, per contrast, not per cluster
#  Here, I plot the results to build a base line and continue polishing the script


# k = 13
k_nice <- sprintf("%02d", k) # Format k

message("Processing DEG for BayesSpace k=", k_nice)

## Load data for the current k
data_file <- file.path(
  input_dir,
  paste0("sce_pseudo_PCA_brain_area_k", k_nice, ".rds")
)
spe_pseudo <- readRDS(data_file)
#colnames(colData(spe_pseudo))
# [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"
# [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
# [9] "nspots"          "pmi"             "rin"             "sample_id"
# [13] "sex"             "sum_umi"
#table(spe_pseudo$brain_area2)
#ncol(spe_pseudo)
#rowData(spe_pseudo)$gene_name  # $symbol or $GeneSymbol depending on source

## Remove level "G0" from the spe_pseudo, as we only have one sample for this group

levels(spe_pseudo$brain_area2)
table(spe_pseudo$brain_area2, spe_pseudo$brain_id)
# Br8518 Br9037 Br9090
# G0     11      0      0
# G1     12     12     12
# G2     12     12     11
# G3     12     12     11
# G4      0     12     12

spe_pseudo <- spe_pseudo[, spe_pseudo$brain_area2 != "G0"]
spe_pseudo$brain_area2 <- droplevels(spe_pseudo$brain_area2)
levels(spe_pseudo$brain_area2)
# [1] "G1" "G2" "G3" "G4"
table(spe_pseudo$brain_area2, spe_pseudo$brain_id)
# Br8518 Br9037 Br9090
# G1     12     12     12
# G2     12     12     11
# G3     12     12     11
# G4      0     12     12

## Subset all Habenula domains from Bayes-Space - g.e. For BS k13 = SpD11

levels_vec <- levels(colData(spe_pseudo)$BayesSpace)
habenula_levels <- levels_vec[grepl("Habenula", levels_vec)]
habenula_levels
# Sp03D02.Habenula

## Iterate on every Habenula cluster in the specific BS (k)

for (hab_level in habenula_levels) {
    
    # hab_levels = habenula_levels[1]
    message("Processing DEG for Hb domain: ", hab_level)

      spe_pseudo_subset <- spe_pseudo[, spe_pseudo$BayesSpace == hab_level]
      spe_pseudo_subset$BayesSpace <- droplevels(spe_pseudo_subset$BayesSpace)
      spe_pseudo_subset$BayesSpace

      ## quick inspection: check how many genes expressed by cluster we have after subset Hb domains from the pseudobulk data

      #expr_mat <- assay(spe_pseudo_subset, "logcounts")
      expr_mat <- assay(spe_pseudo_subset, "counts")
      clusters <- spe_pseudo_subset$BayesSpace
      unique_clusters <- unique(clusters)
      # Count expressed genes per cluster
      genes_per_cluster <- sapply(unique_clusters, function(clust) {
        cluster_expr <- expr_mat[, clusters == clust]
        # Count genes with expression > 0 in at least one pseudobulked sample
        sum(rowSums(cluster_expr > 0) > 0)
      })
      gene_counts_df <- data.frame(
        Cluster = unique_clusters,
        Num_Expressed_Genes = genes_per_cluster,
        row.names = NULL
      )
      print("Number of genes after subset:")
      print(gene_counts_df)
      #                 Cluster Num_Expressed_Genes
      # 1      Sp03D01.Oligo               18798
      # 2   Sp03D02.Habenula               26223
      # 3 Sp03D03.Inhib.Thal               21454

      ## The new scran::pseudoBulkDGE() no longer lets you run the full model without coef — which means it doesn’t expose the fit object needed for custom contrasts.
      ## So, I use the fitted model to manually define and test contrasts using the DGE package
      ## I have to manually aggregate pseudobulk counts and then run the model using standard edgeR workflow
      ## Note: scran::pseudoBulkDGE() works in two modes: (1) If provided a design matrix, so I DO NOT need `coef`
      ##      and I can get `dge_result$fit`, `dge_result$design`, and run custom contrasts for the brain-areas

      ## Get all pairwise combinations of group levels
      group_levels <- levels(factor(colData(spe_pseudo_subset)$brain_area2))
      # [1] "G1" "G2" "G3" "G4"
      pairwise_contrasts <- combn(group_levels, 2, simplify = FALSE)
      head(pairwise_contrasts)
      # [[1]]
      # [1] "G1" "G2"

      ######### Model 1: Brain-Area Naive Model ####################################

      message("Model 1")
      message("Brain-Area: Naive Model")
      model_name = paste0(hab_level, "_model1__", "sampleID-brainID")
      plt_subdir <- here(plot_dir, "sampleID-brainID")
      if (!dir.exists(plt_subdir)) dir.create(plt_subdir, recursive = TRUE)

      # Compute DE results across every pair of group comparisons (e.g. G2 vs G1, G3 vs G1, G3 vs G2, etc.)

      ## Store all results in a list
      dge_results_list <- lapply(pairwise_contrasts, function(contrast_pair) {
        cat("Running contrast:", paste(contrast_pair, collapse = " vs "), "\n")

        res <- runPseudobulkEdgeR(
          sce = spe_pseudo_subset,
          group_var = "brain_area2",
          contrast_levels = contrast_pair,
          covariates = c("sample_id", "brain_id"),
          model_name = model_name
        )

        return(res)
      })

      # Plot volcano with all contrast using p-value (just to check):
      create_volcano_plots_pValue(spe_pseudo_subset, dge_results_list, plt_subdir, model_name)
      # Plot volcano with all contrast using FDR:
      create_volcano_plots(spe_pseudo_subset, dge_results_list, plt_subdir, model_name)
      ## Plot raw p-values to check if we have a flat or U-shaped (not enriched near 0), suggesting no strong DE
      f_name <- here(plt_subdir, paste0("histograms_of_pvalues-", model_name, ".pdf"))
      pdf(f_name, width = 7, height = 5)
      map(dge_results_list, ~ hist(.x$PValue, breaks = 50, main = "Histogram of raw p-values"))
      dev.off()

      ######### Model 2: Brain-Area: sample_id nspots ###########################

      message("Model 1")
      message("Brain-Area: Naive Model")
      model_name = paste0(hab_level, "_model1__", "_sampleID-nspots")
      plt_subdir <- here(plot_dir, "sampleID-nspots")
      if (!dir.exists(plt_subdir)) dir.create(plt_subdir, recursive = TRUE)

      # Compute DE results across every pair of group comparisons (e.g. G2 vs G1, G3 vs G1, G3 vs G2, etc.)

      ## Store all results in a list
      dge_results_list <- lapply(pairwise_contrasts, function(contrast_pair) {
          cat("Running contrast:", paste(contrast_pair, collapse = " vs "), "\n")

          res <- runPseudobulkEdgeR(
              sce = spe_pseudo_subset,
              group_var = "brain_area2",
              contrast_levels = contrast_pair,
              covariates = c("sample_id", "nspots"),
          )

          return(res)
      })

      # Plot volcano with all contrast using p-value (just to check):
      create_volcano_plots_pValue(spe_pseudo_subset, dge_results_list, plt_subdir, model_name)
      # Plot volcano with all contrast using FDR:
      create_volcano_plots(spe_pseudo_subset, dge_results_list, plt_subdir, model_name)
      ## Plot raw p-values to check if we have a flat or U-shaped (not enriched near 0), suggesting no strong DE
      f_name <- here(plt_subdir, paste0("histograms_of_pvalues-", model_name, ".pdf"))
      pdf(f_name, width = 7, height = 5)
      map(dge_results_list, ~ hist(.x$PValue, breaks = 50, main = "Histogram of raw p-values"))
      dev.off()

      ## checking DE results ####################################################

      # message("Genes tested: ", nrow(dge_results_list[[1]]))
      # if (nrow(dge_results_list[[1]]) < 200) { message("Small number of genes tested. Likely inflated FDRs. Tested: ", nrow(dge_results_list[[1]])) }
      # summary(dge_results_list[[1]])
      # table(dge_results_list[[1]]$PValue < 0.05)
      # table(dge_results_list[[1]]$FDR < 0.05)
      #

      # Combine all into one data.frame

      # dge_results_df <- bind_rows(dge_results_list)
      # head(dge_results_df)
      # table(dge_results_df$contrast)
      # # G1 - G2 G1 - G3 G1 - G4 G2 - G3 G2 - G4 G3 - G4
      # # 6352    6352    6352    6352    6352    6352

      # ## Filter data to prepare for plots
      # subset(dge_results_df, FDR < 0.05 & abs(logFC) > 1)

      saveRDS(dge_results_list, file = here(data_dir, paste0("de_results_k", k_nice, "_", model_name,  ".rds")))

}



message("Processing DEG for BayesSpace k=", k_nice, " done!!")

# library("slurmjobs")
# slurmjobs::job_single('05_pseudobulk_DEG_contrast',
#                       create_shell = TRUE, memory = '30G',
#                       command = "05_pseudobulk_DEG_contrast.R",
#                       partition = "katun")



