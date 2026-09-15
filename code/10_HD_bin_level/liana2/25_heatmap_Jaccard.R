# Similar heatmap where each element is the Jaccard index of ligand-receptor pairs for that source-target combo in cellular vs extracellular

library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(readr)
library(here)
library(tidyr)

df_cellular <- read_csv(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana2/table/significant_interactions_across_donors_5000.0.csv",
  show_col_types = FALSE
)

df_extracellular <- read_csv(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana2/extracellular/table/significant_interactions_across_donors_5000.0.csv",
  show_col_types = FALSE
)


output_dir <- here(
  "plots", "10_HD_bin_level", "no_secondary", "liana2", "cellular_vs_extracellular"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# calculate Jaccard similarity of LR pairs for each source-target combo in each donor, then average across donors for a heatmap of mean Jaccard similarity between cellular and extracellular analyses for each source-target combo

cell_lr_donor <- df_cellular %>%
  mutate(lr_pair = paste(ligand, receptor, sep = "-")) %>%
  distinct(donor_id, source, target, lr_pair)

extra_lr_donor <- df_extracellular %>%
  mutate(lr_pair = paste(ligand, receptor, sep = "-")) %>%
  distinct(donor_id, source, target, lr_pair)

all_lr_donor <- bind_rows(
  cell_lr_donor %>% mutate(compartment = "cellular"),
  extra_lr_donor %>% mutate(compartment = "extracellular")
)

jaccard_donor <- all_lr_donor %>%
  group_by(donor_id, source, target) %>%
  summarise(
    cellular_set = list(unique(lr_pair[compartment == "cellular"])),
    extracellular_set = list(unique(lr_pair[compartment == "extracellular"])),
    n_cellular = length(cellular_set[[1]]),
    n_extracellular = length(extracellular_set[[1]]),
    n_intersection = length(intersect(cellular_set[[1]], extracellular_set[[1]])),
    n_union = length(union(cellular_set[[1]], extracellular_set[[1]])),
    jaccard = ifelse(n_union == 0, NA_real_, n_intersection / n_union),
    .groups = "drop"
  )

jaccard_mean <- jaccard_donor %>%
  group_by(source, target) %>%
  summarise(
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    sd_jaccard = sd(jaccard, na.rm = TRUE),
    n_donor = sum(!is.na(jaccard)),
    .groups = "drop"
  ) %>%
  mutate(
    mean_jaccard = ifelse(is.nan(mean_jaccard), NA, mean_jaccard),
    label = ifelse(is.na(mean_jaccard), "", sprintf("%.2f", mean_jaccard))
  )


p <- ggplot(jaccard_mean, aes(x = target, y = source, fill = mean_jaccard)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    na.value = "grey90",
    limits = c(0, 1),
    name = "Mean\nJaccard"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Target cell type",
    y = "Source cell type",
    title = "Donor-averaged LR-pair similarity between cellular and extracellular analyses"
  )

ggsave(
  filename = file.path(output_dir, "cellular_vs_extracellular_jaccard_heatmap.png"),
  plot = p,
  width = 10,
  height = 8
)
