library("spatialLIBD")
library("SingleCellExperiment")
library("scran")
library("BayesSpace")
library("ggplot2")
library("ggpubr")
library("ggrepel")
library("dplyr")
library("here")
library("sessioninfo")

#### Set up dirs ####
data_dir <- here("processed-data", "05_brain_area_differential_expression")

#### Set up dirs ####
plot_dir <- here(
  "plots",
  "05_brain_area_differential_expression",
  "03_covariate_analysis"
)
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Define k values to iterate over
# k_values <- c(3, 9, 17)
k_values <- c(3, 13, 21, 27)

for (k in k_values) {
  # k = 9
  k_nice <- sprintf("%02d", k) # Formatting k as two digits

  # Load the data from the .rds file for current value of k
  data_file <- file.path(
    data_dir,
    paste0("sce_pseudo_PCA_brain_area_k", sprintf("%02d", k), ".rds")
  )

  if (file.exists(data_file)) {
    # Load the RDS file and extract relevant data
    data <- readRDS(data_file)

    # Remove level "G0" from the data, as we only have one sample here
    levels(data$brain_area2)
    table(data$brain_area2, data$brain_id)
    # Br8518 Br9037 Br9090
    # G0      8      0      0
    # G1      8      8      8
    # G2      8      8      8
    # G3      8      8      7
    # G4      0      8      8
    coldata_tbl <- as_tibble(colData(data))
    filtered_coldata <- coldata_tbl |>
      filter(brain_area2 != "G0")
    # Get the matching sample names (assuming there's a column with sample/barcode IDs)
    keep_samples <- filtered_coldata$sample_id # <- replace with your actual column name
    # Subset the SpatialExperiment object
    spe_filtered <- data[, colnames(data) %in% keep_samples]
    # drop unused levels
    colData(spe_filtered)$brain_area2 <- droplevels(
      colData(spe_filtered)$brain_area2
    )
    levels(spe_filtered$brain_area2)
    table(spe_filtered$brain_area2, spe_filtered$brain_id)
    # Br8518 Br9037 Br9090
    # G1      8      8      8
    # G2      8      8      8
    # G3      8      8      7
    # G4      0      8      8

    # Extract the relevant columns from the filtered spe
    plot_data <- as.data.frame(colData(spe_filtered)) |>
      select(
        sample_id,
        brain_area,
        brain_area2,
        BayesSpace,
        nspots,
        sum_umi,
        expr_chrM_ratio,
        pmi,
        age
      )
    
  } else {
    
    message("*****************************************************************")
    warning(paste("Data file not found for k =", k_nice, ". Skipping this k."))
    message("*****************************************************************")
    next # Skip this k if the file is not found
    
  }

  # Generate nspots boxplots for each BayesSpace domain
  plot_name <- paste0("nspots_by_bayesspace_k", k_nice, ".pdf")
  pdf(here(plot_dir, plot_name), width = 8, height = 6)

  for (domain in unique(plot_data$BayesSpace)) {
    # test: domain <- "Sp09D08 ~ Hb"
    # test2: domain <- "Sp03D02 ~ Hb"

    y_max_nspots <- max(plot_data$nspots) + 20
    domain_data <- plot_data |> filter(BayesSpace == domain)
    plot <- ggboxplot(
      domain_data,
      x = "brain_area2",
      y = "nspots",
      color = "brain_area2",
      palette = c("blue", "red", "black", "purple"),
      add = "jitter",
      shape = 19,
      xlab = "Brain Area",
      ylab = "Number of Spots"
    ) +
      geom_text(
        aes(label = sample_id),
        position = position_jitter(width = 0.1, height = 0.2),
        hjust = 0.5,
        vjust = 1,
        size = 3
      ) +
      ggtitle(paste("BayesSpace Domain:", domain)) +
      theme_bw() +
      theme(legend.position = "bottom", axis.text.x = element_text(size = 10)) +
      stat_compare_means(
        aes(group = brain_area2, label = paste0("p = ", after_stat(p.format))),
        label.y = y_max_nspots
      ) #, method = "t.test"
    print(plot)
  }

  dev.off()

  ## nspots boxplots with all spatial domains in 1 plot

  p_values_nspots <- compare_means(
    nspots ~ brain_area2,
    data = plot_data,
    group.by = "BayesSpace"
  ) # method = "t.test"
  p_values_nspots$fdr <- p.adjust(p_values_nspots$p, method = "fdr")
  y_max_nspots <- max(plot_data$nspots) + 50
  plot <- ggboxplot(
    plot_data,
    x = "BayesSpace",
    y = "nspots",
    color = "brain_area2",
    palette = c("blue", "red", "black", "purple"),
    add = "jitter",
    shape = 19,
    xlab = "BayesSpace Domain",
    ylab = "Number of Spots"
  ) +
    geom_text(
      aes(label = sample_id, color = brain_area2),
      position = position_jitter(width = 0.1, height = 0.2),
      size = 1.75,
      hjust = 0.5,
      vjust = 1
    ) +
    theme_bw() +
    scale_x_discrete(labels = scales::label_wrap(10)) + # Wrap labels every 10 characters
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate and justify
    theme(legend.position = "bottom") +
    stat_compare_means(
      aes(group = brain_area2, label = paste0("p = ", after_stat(p.format))),
      label.y = y_max_nspots
    ) + #, method = "t.test"
    geom_text(
      data = p_values_nspots,
      aes(
        x = BayesSpace,
        y = y_max_nspots - 20,
        label = paste0("FDR = ", signif(fdr, 3))
      ),
      size = 3,
      color = "black"
    )
  print(plot)

  plot_name <- paste0("nspots_boxplot_k", k_nice, ".png")
  ggsave(
    filename = here(plot_dir, plot_name),
    plot = plot,
    width = 12,
    height = 8
  )

  # Generate sum_umi boxplots across BayesSpace domains

  p_values_sum_umi <- compare_means(
    nspots ~ brain_area2,
    data = plot_data,
    group.by = "BayesSpace",
    method = "wilcox.test"
  )
  p_values_sum_umi$fdr <- p.adjust(p_values_sum_umi$p, method = "fdr")
  y_max_sum_umi <- max(plot_data$sum_umi) + 600000
  plot <- ggboxplot(
    plot_data,
    x = "BayesSpace",
    y = "sum_umi",
    color = "brain_area2",
    palette = c("blue", "red", "black", "purple"),
    add = "jitter",
    shape = 19,
    xlab = "BayesSpace Domain",
    ylab = "sum_umi"
  ) +
    geom_text(
      aes(label = sample_id, color = brain_area2),
      position = position_jitter(width = 0.1, height = 0.2),
      size = 1.75,
      hjust = 0.5,
      vjust = 1
    ) +
    theme_bw() +
    scale_x_discrete(labels = scales::label_wrap(10)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    stat_compare_means(
      aes(group = brain_area2, label = paste0("p = ", after_stat(p.format))),
      label.y = y_max_sum_umi
    ) + #, method = "t.test"
    geom_text(
      data = p_values_sum_umi,
      aes(
        x = BayesSpace,
        y = y_max_sum_umi - 200000,
        label = paste0("FDR = ", signif(fdr, 2))
      ),
      size = 3,
      color = "black"
    )
  print(plot)

  plot_name <- paste0("sum_umi_boxplot_k", k_nice, ".png")
  ggsave(
    filename = here(plot_dir, plot_name),
    plot = plot,
    width = 12,
    height = 8
  )

  # Generate "expr_chrM_ratio" box plots with all spatial domains in 1 plot

  p_values_mito <- compare_means(
    expr_chrM_ratio ~ brain_area2,
    data = plot_data,
    group.by = "BayesSpace"
  ) # method = "t.test"
  p_values_mito$fdr <- p.adjust(p_values_mito$p, method = "fdr")
  y_max_mito <- max(plot_data$expr_chrM_ratio) + 0.02
  plot <- ggboxplot(
    plot_data,
    x = "BayesSpace",
    y = "expr_chrM_ratio",
    color = "brain_area2",
    palette = c("blue", "red", "black", "purple"),
    add = "jitter",
    shape = 19,
    xlab = "BayesSpace Domain",
    ylab = "expr_chrM_ratio"
  ) +
    geom_text(
      aes(label = sample_id, color = brain_area2),
      position = position_jitter(width = 0, height = 0),
      size = 1.75,
      hjust = 0,
      vjust = 1
    ) +
    theme_bw() +
    scale_x_discrete(labels = scales::label_wrap(10)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    stat_compare_means(
      aes(group = brain_area2, label = paste0("p = ", after_stat(p.format))),
      label.y = y_max_mito
    ) + # , method = "t.test"
    geom_text(
      data = p_values_mito,
      aes(
        x = BayesSpace,
        y = y_max_mito - 0.005,
        label = paste0("FDR = ", signif(fdr, 2))
      ),
      size = 3,
      color = "black"
    )
  print(plot)

  plot_name <- paste0("expr_chrM_ratio_boxplot_k", k_nice, ".png")
  ggsave(
    filename = here(plot_dir, plot_name),
    plot = plot,
    width = 12,
    height = 8
  )
}


###### Additional analysis (since age and PMI remain constant) #############

# k_nice <- "9"

data_file <- file.path(
  data_dir,
  paste0("sce_pseudo_PCA_brain_area_k", k_nice, ".rds")
)
plot_data <- as.data.frame(colData(spe_filtered)) |>
  select(
    sample_id,
    brain_id,
    brain_area,
    brain_area2,
    BayesSpace,
    nspots,
    sum_umi,
    expr_chrM_ratio,
    pmi,
    age
  )

plot_data_summary <- plot_data |>
  select(sample_id, brain_id, age, brain_area2, expr_chrM_ratio, pmi) |>
  distinct(sample_id, .keep_all = TRUE)

# Function for creating and saving box plots

create_boxplot <- function(data, x, y, y_label, filename) {
  pdf(here(plot_dir, filename), width = 8, height = 6)
  y_max <- max(data[[y]]) + 5
  print(
    ggboxplot(
      data,
      x = x,
      y = y,
      color = "brain_area2",
      palette = c("blue", "red", "black", "purple"),
      add = "jitter",
      shape = 19,
      xlab = "brain_area2",
      ylab = y_label
    ) +
      geom_text_repel(
        aes(label = sample_id),
        size = 3,
        color = 'black',
        max.overlaps = Inf,
        force_pull = 2,
        show.legend = FALSE
      ) +
      theme_bw() +
      ggtitle(paste("BayesSpace Domain:", k_nice)) +
      theme(legend.position = "bottom") +
      stat_compare_means(
        aes(group = brain_area2, label = paste0("p = ", after_stat(p.format))),
        label.y = y_max
      ) #, method = "t.test"
  )
  dev.off()
}

# Generate boxplots
colnames(plot_data_summary)
create_boxplot(
  plot_data_summary,
  "brain_area2",
  "pmi",
  "PMI",
  paste0("full_pmi_distribution_k", k_nice, ".pdf")
)
create_boxplot(
  plot_data_summary,
  "brain_area2",
  "age",
  "Age",
  paste0("full_age_distribution_k", k_nice, ".pdf")
)


# Age vs. expr_chrM_ratio scatter plot

# Check and remove NA's
summary(plot_data_summary[, c("age", "expr_chrM_ratio", "brain_area2")])
colSums(is.na(plot_data_summary[, c("age", "expr_chrM_ratio", "brain_area2")]))
plot_data_summary$brain_area2 <- droplevels(plot_data_summary$brain_area2)
levels(plot_data_summary$brain_area2)

plot_data_summary |>
  group_by(age) |>
  tally()

# Note. brain_area2 has fewer than 2 non-missing observations to plot by age, which is causing geom_smooth(method = "lm") to fail.

pdf(here(plot_dir, "age_vs_expr_chrM_by_brain-area.pdf"), width = 8, height = 6)
ggplot(
  plot_data_summary,
  aes(x = age, y = expr_chrM_ratio, color = brain_area2)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = brain_area2), alpha = 0.2) +
  geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 15) +
  scale_color_manual(
    values = c("G1" = "blue", "G2" = "red", "G3" = "black", "G4" = "purple")
  ) +
  scale_fill_manual(
    values = c("G1" = "blue", "G2" = "red", "G3" = "black", "G4" = "purple")
  ) +
  labs(x = "Age", y = "expr_chrM_ratio", color = "brain_area2") +
  theme_bw() +
  theme(legend.position = "bottom")
dev.off()


# library("slurmjobs")
# slurmjobs::job_single('03_covariate_analysis', create_shell = TRUE, memory = '30G', command = "Rscript 03_covariate_analysis.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
