########################################################################
## Compute build a boxplot with principal components (PCs) on the y-axis and brain area on the x-axis 
## Input: BayesSpace pseudobulk data
## Output: Box-Plots
## Authors. CSC
## Data: Apr24, 2025
## For 60 to 80k spots: $srun --pty --mem=60GB --x11 bash
########################################################################

library("SingleCellExperiment")
library("variancePartition")
library("pheatmap")
library("grid")
library("ggplot2")
library("tidyr")
library("here")
library("sessioninfo")

#### Set up input dirs
data_dir <- here("processed-data", "05_brain_area_differential_expression")

#### Set up output dirs 
plot_dir <- here(
  "plots",
  "05_brain_area_differential_expression",
  "04_CCA_variance_partition"
)

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Define k values to iterate over
k_values <- c(3, 13, 21, 26)


for (k in k_values) {
  # k = 3
  k_nice <- sprintf("%02d", k) # Formatting k as two digits
  
  print(paste0("Plotting CCA for BayesSpace k = ", k_nice))
  ## load spe_pseudo data
  spe_pseudo <-
    readRDS(
      file.path(
        data_dir,
        paste0("sce_pseudo_PCA_brain_area_k", k_nice, ".rds")
      )
    )
  
  # Access meta-data
  
  #colnames(colData(spe_pseudo))
  # [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"    
  # [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
  # [9] "nspots"          "pmi"             "rin"             "sample_id"      
  # [13] "sex"             "sum_umi"         "PC1"             
  
  ## Prior to the variance partition analysis, I evaluate the correlation between sample variables
  # - highly correlated variables can produce unstable estimates of the variance fractions and 
  # - impede the identification of the variables that really contribute to the expression variation
  
  ## First, I plot heatmap of correlations
  
  ## Use expressed genes only (i.e. that passed the filtering step)
  rse_gene_filt <- spe_pseudo[
    rowData(spe_pseudo)$gene_id,
  ]
  
  ## variables to examine
  formula <- ~ BayesSpace + brain_id  + brain_area2 + sex  + expr_chrM_ratio + nspots + pmi + rin
  ## Measure correlations
  CCA <- canCorPairs(formula, colData(rse_gene_filt))
  
  plot_name <- paste0("CCA_pseudo_BS_k", k_nice, ".pdf")
  pdf(file = here(plot_dir, plot_name), width = 5, height = 5)
  
  ## Heatmap
  pheatmap(
    CCA, 
    color = hcl.colors(50, "YlOrRd", rev = TRUE), ## color scale
    fontsize = 8, ## text size
    border_color = "black", ## border color for heatmap cells
    cellwidth = unit(0.4, "cm"), ## height of cells
    cellheight = unit(0.4, "cm") ## width of cells
  )
  # Add title manually using grid
  grid::grid.text(paste0("CCA_pseudo_BS_k", k_nice), x = 0.5, y = 0.9, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  dev.off()
  
  # Access PC1 reduced dimensions
  
  #pcs <- reducedDim(spe_pseudo, "PCA")[, 1:20]
  #colnames(pcs)
  # [1] "PC01" "PC02" "PC03" "PC04" "PC05" "PC06" "PC07" "PC08" "PC09" "PC10"
  # [11] "PC11" "PC12" "PC13" "PC14" "PC15" "PC16" "PC17" "PC18" "PC19" "PC20"
  # Add PC to the colData() for plotting
  # colData(spe_pseudo)$PC1 <- pcs[, 1]
  # colData(spe_pseudo)$brain_area_ln <- factor(colData(spe_pseudo)$brain_area2)  
  # ggplot(as.data.frame(colData(spe_pseudo)), aes(x = brain_area_ln, y = PC1, fill = brain_area_ln)) +
  #   geom_boxplot() +
  #   theme_bw() +
  #   labs(x = "Brain Area", y = "PC1", title = "PC1 by Brain Area")
  
  
  ## Then, build Box-Plots of PCs vs Brain-Areas to check variance among Brain-Areas  
  
  ## Access PCA reduced dimensions
  
  df_pseudo <- as.data.frame(colData(spe_pseudo))
  pcs_df <- as.data.frame(reducedDim(spe_pseudo, "PCA")[, 1:20])
  
  df_long <- cbind(df_pseudo, pcs_df) |>
    pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "value")
  ## forces the PC levels to follow "PC1" to "PC20" in numeric order (avoid unsorted by string value on facet_wrap)
  df_long$PC <- factor(df_long$PC, levels = paste0("PC", 1:20))
  head(df_long)
  
  pdf(file = file.path(plot_dir, paste0("BoxPlot_pseudo_PCs-BrainArea_k", k_nice, ".pdf")), width = 8, height = 8)
  
  plt1 <- ggplot(df_long, aes(x = brain_area2, y = value, fill = brain_area2)) +
    geom_boxplot() +
    # Tendency line: mean value per brain_area2
    stat_summary(
      fun = mean,
      geom = "line",
      aes(group = 1),
      color = "black",
      linewidth = 1
    ) +
    facet_wrap(~ PC, scales = "free_y") +
    ggtitle(paste("BayesSpace Domain:", k_nice)) +
    theme_bw() + 
    labs(x = NULL, y = NULL)
  print(plt1)
  
  dev.off()
  
  rm("spe_pseudo")
  
}

message("CCA done!")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
