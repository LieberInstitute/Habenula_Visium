# At the middle bandwidth: scatter plot among the 3 high-scoring cell-type combos from the average-score heatmap you presented: 
# x-axis is mean LR score, y-axis is mean ratio, color by cell-type combination. 
# Consider plotting top-5 or so LR pairs by mean score (just to not have too many points to see, 
# and so you could add annotation of which pairs are being plotted). For the "mean ratio", 
# from among the significant pairs, ratio of highest mean to to second-highest mean among source-target combinations
# read in the data
library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(readr)
library(here)

df <- read_csv(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana2/extracellular/table/significant_interactions_across_donors_5000.0.csv",
  show_col_types = FALSE
)

# -----------------------------
# Parameters
# -----------------------------
top_prop <- 0.20
global_keep_prop <- 0.50
top_n_per_combo <- 5
min_n_donors <- 8

combos_of_interest <- c(
  "MHb.1 -> LHb.2.7",
  "MHb.2 -> LHb.2.7",
  "MHb.1 -> Excit_LHb",
  "MHb.2 -> Excit_LHb",
  "Astrocyte -> Excit_LHb",
  "Astrocyte -> LHb.2.7",
  "Astrocyte -> MHb.1",
  "Astrocyte -> MHb.2",
  "Astrocyte -> Subependymal",
  "Astrocyte -> LHb.4",
  "Astrocyte -> Ependymal",
  "Ependymal -> Subependymal",
  "Ependymal -> MHb.1",
  "Ependymal -> MHb.2",
  "MHb.1 -> Subependymal",
  "MHb.2 -> Subependymal",
  "MHb.1 -> Ependymal",
  "MHb.2 -> Ependymal",
  "Oligo -> Excit_LHb",
  "Oligo -> LHb.2.7",
  "Oligo -> LHb.4",
  "Oligo -> MHb.1",
  "Oligo -> MHb.2",
  "Oligo -> Subependymal",
  "Oligo -> Ependymal",
  "Oligo -> Excit.Thal/Inhib_LHb_4.2",
  "Oligo -> Astrocyte",
  "Subependymal -> MHb.1",
  "Subependymal -> MHb.2",
  "Subependymal -> Ependymal"
)

# -----------------------------
# 1. Prepare source-target and LR-pair labels
# -----------------------------
df <- df %>%
  mutate(
    source_target = paste(source, target, sep = " -> "),
    lr_pair = paste(ligand, receptor, sep = "-")
  )

# -----------------------------
# 2. Average LR absolute score within donor
# -----------------------------
lr_by_donor <- df %>%
  group_by(donor_id, source_target, source, target, lr_pair, ligand, receptor) %>%
  summarise(
    donor_abs_score = mean(mean, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 3. Average across donors
# -----------------------------
lr_score <- lr_by_donor %>%
  group_by(source_target, source, target, lr_pair, ligand, receptor) %>%
  summarise(
    abs_score = mean(donor_abs_score, na.rm = TRUE),
    n_donors = n_distinct(donor_id),
    .groups = "drop"
  ) %>%
  filter(n_donors >= min_n_donors)

# -----------------------------
# 4. Calculate mean ratio for each LR pair
#    ratio = abs_score in this source-target combo /
#            strongest abs_score in other source-target combos
# -----------------------------
mean_ratio_lr <- lr_score %>%
  group_by(lr_pair) %>%
  filter(n() >= 2) %>%
  mutate(
    n_combos_for_lr = n(),

    strongest_other_combo = map_chr(
      source_target,
      ~ source_target[
        which.max(ifelse(source_target != .x, abs_score, -Inf))
      ]
    ),

    mean_other_max = map_dbl(
      source_target,
      ~ max(abs_score[source_target != .x], na.rm = TRUE)
    ),

    mean_ratio = abs_score / mean_other_max
  ) %>%
  ungroup() %>%
  filter(
    mean_other_max > 0,
    !is.na(mean_ratio),
    is.finite(mean_ratio)
  )

# -----------------------------
# 5. Keep selected source-target combos
# -----------------------------
mean_ratio_lr_selected <- mean_ratio_lr %>%
  filter(source_target %in% combos_of_interest)

# -----------------------------
# 6. Apply filtering strategy across all selected source-target combos:
#    A. Within each selected combo, keep top 20% by abs_score
#    B. Among those retained LR interactions, remove bottom 50% globally by abs_score
#    C. For each LR pair, keep only the source-target combo with the highest abs_score
#    D. Show only the top N interactions overall
# -----------------------------

top_n_total <- 30

# First only keep combos of interest
mean_ratio_lr_selected <- mean_ratio_lr %>%
  filter(source_target %in% combos_of_interest)

# A. Top 20% within each selected source-target combo
top20_within_combo <- mean_ratio_lr_selected %>%
  group_by(source, target, source_target) %>%
  mutate(
    abs_score_cutoff = quantile(
      abs_score,
      probs = 1 - top_prop,
      na.rm = TRUE
    ),
    n_lr_before_filter = n()
  ) %>%
  filter(abs_score >= abs_score_cutoff) %>%
  ungroup()

# B. Remove bottom 50% globally after the top 20% filter
global_abs_score_cutoff_after_top20 <- quantile(
  top20_within_combo$abs_score,
  probs = 1 - global_keep_prop,
  na.rm = TRUE
)

top20_global50 <- top20_within_combo %>%
  mutate(
    global_abs_score_cutoff_after_top20 = global_abs_score_cutoff_after_top20
  ) %>%
  filter(abs_score >= global_abs_score_cutoff_after_top20)

# C. For the same LR pair, keep only the source-target combo with the highest abs_score
# D. Then show only top N interactions overall
top_lr_for_plot <- top20_global50 %>%
  group_by(lr_pair) %>%
  slice_max(
    order_by = abs_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(desc(abs_score)) %>%
  slice_head(n = top_n_total) %>%
  mutate(
    ct_pair = source_target,
    lr_pair_label = paste(ligand, receptor, sep = " -> ")
  )

# -----------------------------
# 7. Scatter plot
# -----------------------------
p <- ggplot(
  top_lr_for_plot,
  aes(
    x = abs_score,
    y = mean_ratio,
    color = ct_pair,
    label = lr_pair_label
  )
) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    x = "Mean absolute LR score",
    y = "Mean ratio: this combo / strongest other combo",
    color = "Cell-type pair"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

plot_dir <- here::here("plots", "10_HD_bin_level", "no_secondary", "liana2", "extracellular")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  p,
  filename = file.path(plot_dir, "scatter_top20_global50_uniqueLR_top20.png"),
  width = 9,
  height = 6,
  dpi = 300
)

cat("Number of LR interactions plotted:", nrow(top_lr_for_plot), "\n")
cat("Global abs_score cutoff after within-combo top 20%:", global_abs_score_cutoff_after_top20, "\n")