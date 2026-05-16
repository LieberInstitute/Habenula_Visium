# At the middle bandwidth: scatter plot among the 3 high-scoring cell-type combos from the average-score heatmap you presented: 
# x-axis is mean LR score, y-axis is mean ratio, color by cell-type combination. 
# Consider plotting top-5 or so LR pairs by mean score (just to not have too many points to see, 
# and so you could add annotation of which pairs are being plotted). For the "mean ratio", 
# from among the significant pairs, ratio of highest mean to to second-highest mean among source-target combinations

# read in the data
library(ggplot2)
library(dplyr)  
df = read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana2/table/significant_interactions_across_donors_5000.0.csv")

# -----------------------------
# 1. Prepare source-target and LR-pair labels
# -----------------------------
df <- df %>%
  mutate(
    source_target = paste(source, target, sep = " -> "),
    lr_pair = paste(ligand, receptor, sep = "-")
  )

# -----------------------------
# Average LR score for each source-target-LR pair
# This is useful if the same LR pair appears multiple times,
# for example across donors.
# -----------------------------
lr_by_donor <- df %>%
  group_by(donor_id, source_target, source, target, lr_pair, ligand, receptor) %>%
  summarise(
    donor_mean_score = mean(mean, na.rm = TRUE),
    .groups = "drop"
  )

lr_score <- lr_by_donor %>%
  group_by(source_target, source, target, lr_pair, ligand, receptor) %>%
  summarise(
    mean_in_combo = mean(donor_mean_score, na.rm = TRUE),
    n_donors = n_distinct(donor_id),
    .groups = "drop"
  ) %>% filter(n_donors == 8) # keep only those with at least 2 donors

# -----------------------------
# 2. Calculate mean ratio within each source-target combo
# -----------------------------
library(dplyr)
library(purrr)

mean_ratio_lr <- lr_score %>%
  group_by(lr_pair) %>%
  filter(n() >= 2) %>%
  mutate(
    n_combos_for_lr = n(),
    
    strongest_other_combo = map_chr(
      source_target,
      ~ source_target[which.max(ifelse(source_target != .x, mean_in_combo, -Inf))]
    ),
    
    mean_other_max = map_dbl(
      source_target,
      ~ max(mean_in_combo[source_target != .x], na.rm = TRUE)
    ),
    
    mean_ratio = mean_in_combo / mean_other_max
  ) %>%
  ungroup() %>%
  rename(mean_target = mean_in_combo) %>%
  filter(
    mean_other_max > 0,
    !is.na(mean_ratio),
    is.finite(mean_ratio)
  )

# -----------------------------
# 3. Find top 3 source-target combos and Select top 5 LR pairs within each top source-target combo
# -----------------------------

# top 3 source-target combos by average LR score
top_combos <- lr_score %>%
  group_by(source_target) %>%
  summarise(avg_combo_score = mean(mean_in_combo, na.rm = TRUE), .groups = "drop") %>%
  slice_max(avg_combo_score, n = 5)

# change it to what we expect: MHb -> MHb, LHb -> LHb, MHb -> LHb
top_combos <- lr_score %>%
  group_by(source_target) %>%
  summarise(avg_combo_score = mean(mean_in_combo, na.rm = TRUE), .groups = "drop") %>%
  filter(source_target %in% c("MHb.1 -> LHb.2.7", "MHb.2 -> LHb.2.7", "MHb.1 -> Excit_LHb"))

# within each top combo, select top 5 LR pairs by mean LR score
top_lr_for_plot <- mean_ratio_lr %>%
  semi_join(top_combos, by = "source_target") %>%
  group_by(source_target) %>%
  slice_max(mean_target, n = 5, with_ties = FALSE) %>%
  ungroup()

# -----------------------------
# 7. Scatter plot
# -----------------------------

library(ggplot2)
library(ggrepel)

p <- ggplot(top_lr_for_plot,
       aes(x = mean_target,
           y = mean_ratio,
           color = source_target,
           label = lr_pair)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    x = "Mean LR score",
    y = "Mean ratio (highest / second-highest)",
    color = "Cell-type combination"
  ) +
  theme_bw()

plot_dir = here::here("plots", "10_HD_bin_level","no_secondary", "liana2")
ggsave(p, filename = paste0(plot_dir, "/scatter_top_LR_pairs.png"), width = 8, height = 6)



top_lr_for_plot <- mean_ratio_lr %>%
  semi_join(top_combos, by = "source_target") %>%
  group_by(source_target) %>%
  slice_max(mean_ratio, n = 5, with_ties = FALSE) %>%
  ungroup()


library(ggplot2)
library(ggrepel)

p <- ggplot(top_lr_for_plot,
       aes(x = mean_target,
           y = mean_ratio,
           color = source_target,
           label = lr_pair)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = Inf) +
  labs(
    x = "Mean LR score",
    y = "Mean ratio (highest / second-highest)",
    color = "Cell-type combination"
  ) +
  theme_bw()

plot_dir = here::here("plots", "10_HD_bin_level","no_secondary", "liana2")
ggsave(p, filename = paste0(plot_dir, "/scatter_top_LR_pairs_mean_ratio.png"), width = 8, height = 6)

