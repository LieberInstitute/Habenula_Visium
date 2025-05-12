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
library("dplyr")
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


## Function to plot PCA by BayesSpace and specific SPD across the 4 brain-area groups

plot_PC_linear_brainArea <- function(
        spe_p,
        k_nice,
        total_sps,
        spD
        ) {
    # # for testing specific SpD
    # spe_p = spe_pseudo_subset
    # k_nice
    # total_sps = total_spots
    # spD = "xx"

    ## ectract PCs
    df_pseudo <- as.data.frame(colData(spe_p))
    pcs_df <- as.data.frame(reducedDim(spe_p, "PCA")[, 1:20])
    df_long <- cbind(df_pseudo, pcs_df) |>
        pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "value")
    ## forces the PC levels to follow "PC1" to "PC20" in numeric order (avoid unsorted by string value on facet_wrap)
    df_long$PC <- factor(df_long$PC, levels = paste0("PC", 1:20))
    head(df_long)

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
        ggtitle(paste("BayesSpace K:", k_nice, " - SpD ", spD), subtitle = paste("Number of spots:", total_sps)) +
        theme_bw() +
        labs(x = NULL, y = NULL)

    return(plt1)
}


## plot CCA Hetmap for all BS domains
## plot BoxPlots for all BS domains and Hb specific SpD

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


  ## Then, build Box-Plots of PCs vs Brain-Areas to check variance across Brain-Areas
  ## Access PCA reduced dimensions

  ## quick inspection
  table(spe_pseudo$BayesSpace, spe_pseudo$brain_area2)
  ## summarize total spots by brain-area2
  spot_summary <- as.data.frame(colData(spe_pseudo)) |>
      group_by(brain_area2) |>
      summarise(total_spots = sum(nspots))
  spot_summary
  total_spots <- sum(spot_summary$total_spots)
  ## plot the PC(s) by all SpD in the specific BS k
  pdf(file = file.path(plot_dir, paste0("BoxPlot_PC-BrainArea_BSk", k_nice, "_All.pdf")), width = 8, height = 8)
  plt1 <- plot_PC_linear_brainArea(spe_pseudo,
                           k_nice,
                           total_spots,
                           "ALL")
  print(plt1)
  dev.off()

  ## make it more specific extracting Hb SpD(s)

  levels_vec <- levels(colData(spe_pseudo)$BayesSpace)
  habenula_levels <- levels_vec[grepl("Habenula", levels_vec)]
  habenula_levels
  # "Sp13D11.Habenula"

  message("Habenula SpD found: ", length(habenula_levels))

  ## Iterate on every Habenula cluster in the specific BS (k)

  for (hab_level in habenula_levels) {

      message("Plotting PC(s) across Hb SpD: ", hab_level)

      spe_pseudo_subset <- spe_pseudo[, spe_pseudo$BayesSpace == hab_level]
      spe_pseudo_subset$BayesSpace <- droplevels(spe_pseudo_subset$BayesSpace)
      ## quick inspection
      table(spe_pseudo_subset$BayesSpace, spe_pseudo_subset$brain_area2)
      ## summarize total spots by brain-area2
      spot_summary <- as.data.frame(colData(spe_pseudo_subset)) |>
          group_by(brain_area2) |>
          summarise(total_spots = sum(nspots))
      spot_summary
      total_spots <- sum(spot_summary$total_spots)
      message("Total_spots in the SpD: ", total_spots )

      ## plot specific SpD
      f_name <- paste0("BoxPlot_PC-BrainArea_BSk", k_nice, "_Hb_", gsub(".Habenula","", hab_level), ".pdf")
      pdf(file = file.path(plot_dir, f_name),
          width = 8, height = 8)
      plt1 <- plot_PC_linear_brainArea(spe_pseudo_subset,
                               k_nice,
                               total_spots,
                               hab_level)
      print(plt1)
      dev.off()
  }

}


message("CCA done!")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
