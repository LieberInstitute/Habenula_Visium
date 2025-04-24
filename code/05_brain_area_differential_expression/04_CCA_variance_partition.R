
library("SingleCellExperiment")
library("BayesSpace")
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
  # k = 13
  k_nice <- sprintf("%02d", k) # Formatting k as two digits
  
  ## load spe_pseudo data
  
  spe_pseudo <-
    readRDS(
      file.path(
        data_dir,
        paste0("sce_pseudo_PCA_brain_area_k", sprintf("%02d", k), ".rds")
      )
    )
  
  # Access meta-data
  #colnames(colData(spe_pseudo))
  # [1] "age"             "BayesSpace"      "brain_area"      "brain_area2"    
  # [5] "brain_id"        "diagnosis"       "expr_chrM"       "expr_chrM_ratio"
  # [9] "nspots"          "pmi"             "rin"             "sample_id"      
  # [13] "sex"             "sum_umi"         "PC1"             "brain_area_ln"  
  
  # Access PCA reduced dimensions
  pcs <- reducedDim(spe_pseudo, "PCA")
  #colnames(pcs)
  # [1] "PC01" "PC02" "PC03" "PC04" "PC05" "PC06" "PC07" "PC08" "PC09" "PC10"
  # [11] "PC11" "PC12" "PC13" "PC14" "PC15" "PC16" "PC17" "PC18" "PC19" "PC20"
  
  # Add PC to the colData() for plotting
  colData(spe_pseudo)$PC1 <- pcs[, 1]
  colData(spe_pseudo)$brain_area_ln <- factor(colData(spe_pseudo)$brain_area2)  
  
  
  # ggplot(as.data.frame(colData(spe_pseudo)), aes(x = brain_area_ln, y = PC1, fill = brain_area_ln)) +
  #   geom_boxplot() +
  #   theme_bw() +
  #   labs(x = "Brain Area", y = "PC1", title = "PC1 by Brain Area")
  
  df <- as.data.frame(colData(spe_pseudo))
  pcs_df <- as.data.frame(reducedDim(spe_pseudo, "PCA"))
  
  df_long <- cbind(df, pcs_df) %>%
    pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "value")
  
  pdf(file = file.path(dir_plots, paste0("sce_pseudo_PCs_k", sprintf("%02d", k), ".pdf")), width = 8, height = 8)
  fontSize = 6
  
  ggplot(df_long, aes(x = brain_area_ln, y = value, fill = brain_area_ln)) +
    geom_boxplot() +
    facet_wrap(~ PC, scales = "free_y") +
    theme_bw()
  
  dev.off()
  
}


