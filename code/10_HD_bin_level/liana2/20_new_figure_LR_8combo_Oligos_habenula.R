# ==============================================================================
# LIANA LR scatter plot
#
# Goal:
#   Identify ligand-receptor interactions that are:
#
#   1. reproducibly observed across donors
#   2. enriched in selected cell-type pairs
#   3. among the strongest interactions within each selected cell-type pair
#
# Plot:
#   x-axis = mean LR score across donors (absolute strength)
#   y-axis = LR enrichment relative to mean across all mapped cell-type pairs
#   color  = cell-type pair
#
# Selection:
#   1. enrichment > enrichment_cutoff
#   2. within each selected cell-type pair,
#      keep top N interactions by mean LR score
# ==============================================================================


# ==============================================================================
# 0. Libraries
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(here)


# ==============================================================================
# 1. Input files
# ==============================================================================

data_file <- paste0(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
  "Habenula_Visium/processed-data/10_HD_bin_level/",
  "no_secondary/liana2/table/",
  "significant_interactions_across_donors_5000.0.csv"
)


CELL_TYPE_MAP_FILE <- paste0(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
  "Habenula_Visium/raw-data/sample_info/",
  "hd_cell_type_map.csv"
)


# ==============================================================================
# 2. Parameters
# ==============================================================================

# Minimum number of donors in which an LR x cell-pair interaction
# must be present.
min_n_donors <- 4


# Relative enrichment threshold.
#
# enrichment =
#
#   mean LR score in current cell-type pair
#   ---------------------------------------
#   mean LR score across all mapped cell-type pairs
#
# enrichment = 1:
#   average level for this LR
#
# enrichment = 2:
#   twice the average level for this LR
#
# enrichment > 2:
#   preferential enrichment in this cell-type pair
#
enrichment_cutoff <- 3


# Maximum number of LR interactions plotted
# for each selected cell-type pair.
#
# IMPORTANT:
# These are selected by abs_score AFTER requiring
# enrichment > enrichment_cutoff.
#
top_n_per_combo <- 4


# Selected NEW / mapped cell-type combinations.
combos_of_interest <- c(
  "Oligo -> MHb_A",
  "Oligo -> MHb_B",
  "Oligo -> Excit_LHb",
  "Oligo -> LHb_A",
  "Oligo -> LHb_C"
)

analysis_label <- "Oligo_to_Habenula"

# ==============================================================================
# 3. Read data
# ==============================================================================

df <- read_csv(
  data_file,
  show_col_types = FALSE
)


cell_type_map <- read_csv(
  CELL_TYPE_MAP_FILE,
  show_col_types = FALSE
)


cat(
  "\n============================================================\n",
  "INPUT DATA\n",
  "============================================================\n",
  sep = ""
)

cat(
  "LIANA file:",
  data_file,
  "\n"
)

cat(
  "Rows in LIANA input:",
  nrow(df),
  "\n"
)

cat(
  "Unique donors:",
  n_distinct(df$donor_id),
  "\n"
)


# ==============================================================================
# 4. Validate input
# ==============================================================================

required_data_columns <- c(
  "donor_id",
  "source",
  "target",
  "ligand",
  "receptor",
  "mean"
)


missing_data_columns <- setdiff(
  required_data_columns,
  colnames(df)
)


if (length(missing_data_columns) > 0) {

  stop(
    paste0(
      "Input LIANA file is missing required columns: ",
      paste(
        missing_data_columns,
        collapse = ", "
      )
    )
  )
}


required_mapping_columns <- c(
  "old_cell_type",
  "new_cell_type",
  "color"
)


missing_mapping_columns <- setdiff(
  required_mapping_columns,
  colnames(cell_type_map)
)


if (length(missing_mapping_columns) > 0) {

  stop(
    paste0(
      "Cell-type mapping file is missing required columns: ",
      paste(
        missing_mapping_columns,
        collapse = ", "
      )
    )
  )
}


# Clean mapping strings.
cell_type_map <- cell_type_map %>%
  mutate(
    old_cell_type = trimws(old_cell_type),
    new_cell_type = trimws(new_cell_type),
    color = trimws(color)
  )


# Check duplicated old names.
duplicated_old_names <- cell_type_map %>%
  filter(
    duplicated(old_cell_type) |
      duplicated(
        old_cell_type,
        fromLast = TRUE
      )
  )


if (nrow(duplicated_old_names) > 0) {

  cat(
    "\nDuplicated old_cell_type entries detected:\n"
  )

  print(
    duplicated_old_names,
    n = Inf
  )

  stop(
    "Each old_cell_type must map to only one new_cell_type."
  )
}


# ==============================================================================
# 5. Check donor threshold
# ==============================================================================

n_total_donors <- n_distinct(
  df$donor_id
)


if (min_n_donors > n_total_donors) {

  stop(
    paste0(
      "min_n_donors = ",
      min_n_donors,
      ", but the input contains only ",
      n_total_donors,
      " unique donors. ",
      "Please lower min_n_donors before running the analysis."
    )
  )
}


# ==============================================================================
# 6. Map ORIGINAL cell types to NEW cell types
# ==============================================================================

source_map <- cell_type_map %>%
  select(
    source_original = old_cell_type,
    source_new = new_cell_type,
    source_color = color
  )


target_map <- cell_type_map %>%
  select(
    target_original = old_cell_type,
    target_new = new_cell_type,
    target_color = color
  )


df <- df %>%

  rename(
    source_original = source,
    target_original = target
  ) %>%

  left_join(
    source_map,
    by = "source_original"
  ) %>%

  left_join(
    target_map,
    by = "target_original"
  ) %>%

  mutate(

    # If mapping is unavailable, preserve original name.
    source = coalesce(
      source_new,
      source_original
    ),

    target = coalesce(
      target_new,
      target_original
    ),

    # NEW mapped cell-type pair.
    source_target = paste(
      source,
      target,
      sep = " -> "
    ),

    # LR identifier.
    lr_pair = paste(
      ligand,
      receptor,
      sep = "-"
    )
  )


# ==============================================================================
# 7. Mapping audit
# ==============================================================================

cat(
  "\n============================================================\n",
  "CELL-TYPE MAPPING\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nOriginal -> mapped SOURCE cell types:\n"
)

print(
  df %>%
    distinct(
      source_original,
      source
    ) %>%
    arrange(source),
  n = Inf
)


cat(
  "\nOriginal -> mapped TARGET cell types:\n"
)

print(
  df %>%
    distinct(
      target_original,
      target
    ) %>%
    arrange(target),
  n = Inf
)


cat(
  "\nAvailable NEW cell-type combinations:\n"
)

print(
  df %>%
    distinct(source_target) %>%
    arrange(source_target),
  n = Inf
)


# Check requested combinations.
missing_combos <- setdiff(
  combos_of_interest,
  unique(df$source_target)
)


if (length(missing_combos) > 0) {

  cat(
    "\nWARNING: These requested combinations were not found ",
    "after cell-type mapping:\n",
    sep = ""
  )

  print(
    missing_combos
  )
}


# ==============================================================================
# 8. Average LR score WITHIN donor
#
# This occurs AFTER cell-type mapping.
#
# If several original source-target combinations collapse onto
# the same NEW source-target combination, their scores are first
# averaged within each donor.
# ==============================================================================

lr_by_donor <- df %>%

  group_by(
    donor_id,
    source_target,
    source,
    target,
    lr_pair,
    ligand,
    receptor
  ) %>%

  summarise(

    donor_abs_score = mean(
      mean,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  # Remove pathological NaN / infinite results.
  filter(
    !is.na(donor_abs_score),
    is.finite(donor_abs_score)
  )


cat(
  "\n============================================================\n",
  "WITHIN-DONOR AGGREGATION\n",
  "============================================================\n",
  sep = ""
)

cat(
  "LR x mapped-cell-pair x donor observations:",
  nrow(lr_by_donor),
  "\n"
)


# ==============================================================================
# 9. Average LR score ACROSS donors
# ==============================================================================

lr_score <- lr_by_donor %>%

  group_by(
    source_target,
    source,
    target,
    lr_pair,
    ligand,
    receptor
  ) %>%

  summarise(

    # Mean LR interaction score across donors.
    abs_score = mean(
      donor_abs_score,
      na.rm = TRUE
    ),

    # Number of donors contributing to this interaction.
    n_donors = n_distinct(
      donor_id
    ),

    .groups = "drop"
  ) %>%

  filter(
    n_donors >= min_n_donors,
    !is.na(abs_score),
    is.finite(abs_score)
  )


cat(
  "\n============================================================\n",
  "ACROSS-DONOR AGGREGATION\n",
  "============================================================\n",
  sep = ""
)


cat(
  "Minimum donor requirement:",
  min_n_donors,
  "\n"
)


cat(
  "Mapped cell-type combinations after donor filter:",
  n_distinct(lr_score$source_target),
  "\n"
)


cat(
  "LR x cell-pair observations after donor filter:",
  nrow(lr_score),
  "\n"
)


if (nrow(lr_score) == 0) {

  stop(
    paste0(
      "No LR interactions remain after requiring n_donors >= ",
      min_n_donors,
      "."
    )
  )
}


# ==============================================================================
# 10. Calculate LR enrichment across ALL mapped cell-type combinations
#
# For EACH LR:
#
#                         score in current cell-type pair
# enrichment = ----------------------------------------------------
#              mean score of this LR across all mapped cell pairs
#
#
# Interpretation:
#
# enrichment = 1
#   current cell pair has the average score for this LR
#
# enrichment = 2
#   current cell pair has twice the average score for this LR
#
# enrichment > 2
#   LR is preferentially enriched in this cell-type pair
#
#
# IMPORTANT:
#
# This calculation is performed BEFORE restricting to
# combos_of_interest.
#
# Therefore, the denominator uses ALL mapped cell-type
# combinations where this LR survives the donor filter.
# ==============================================================================

lr_enrichment <- lr_score %>%

  group_by(
    lr_pair
  ) %>%

  # Need at least two cell-type combinations to define
  # meaningful relative enrichment.
  filter(
    n() >= 2
  ) %>%

  mutate(

    # Number of mapped cell-type combinations available
    # for this LR.
    n_combos_for_lr = n(),

    # Average score for this LR across ALL mapped
    # source-target combinations.
    mean_score_all_combos = mean(
      abs_score,
      na.rm = TRUE
    ),

    # Relative enrichment.
    enrichment =
      abs_score /
      mean_score_all_combos
  ) %>%

  ungroup() %>%

  filter(
    mean_score_all_combos > 0,
    !is.na(enrichment),
    is.finite(enrichment)
  )


# ==============================================================================
# 11. Inspect enrichment distribution
# ==============================================================================

cat(
  "\n============================================================\n",
  "LR ENRICHMENT\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nDistribution of LR enrichment across all mapped cell pairs:\n"
)


print(
  summary(
    lr_enrichment$enrichment
  )
)


cat(
  "\nNumber of available mapped cell-type combinations per LR:\n"
)


print(
  lr_enrichment %>%

    distinct(
      lr_pair,
      n_combos_for_lr
    ) %>%

    count(
      n_combos_for_lr,
      name = "n_lr_pairs"
    ) %>%

    arrange(
      n_combos_for_lr
    ),

  n = Inf
)


# ==============================================================================
# 12. Restrict to selected mapped source-target combinations
#
# IMPORTANT:
#
# Restriction happens AFTER calculating enrichment.
# The enrichment denominator therefore still represents
# ALL available mapped cell-type combinations.
# ==============================================================================

lr_enrichment_selected <- lr_enrichment %>%

  filter(
    source_target %in% combos_of_interest
  )


cat(
  "\n============================================================\n",
  "SELECTED CELL-TYPE PAIRS\n",
  "============================================================\n",
  sep = ""
)


cat(
  "LR x cell-pair observations in selected combinations:",
  nrow(lr_enrichment_selected),
  "\n"
)


# ==============================================================================
# 13. Apply specificity threshold
#
# Keep LR interactions with enrichment > enrichment_cutoff.
#
# This is the SPECIFICITY requirement.
# ==============================================================================

enriched_lr <- lr_enrichment_selected %>%

  filter(
    enrichment > enrichment_cutoff
  )


cat(
  "\n============================================================\n",
  "ENRICHED LR INTERACTIONS\n",
  "============================================================\n",
  sep = ""
)


cat(
  "Number of interactions with enrichment > ",
  enrichment_cutoff,
  ": ",
  nrow(enriched_lr),
  "\n",
  sep = ""
)


cat(
  "\nNumber of enriched interactions per selected cell-type pair ",
  "BEFORE top-N selection:\n",
  sep = ""
)


print(
  enriched_lr %>%

    count(
      source_target,
      name = "n_enriched"
    ) %>%

    mutate(
      source_target = factor(
        source_target,
        levels = combos_of_interest
      )
    ) %>%

    arrange(
      source_target
    ),

  n = Inf
)


if (nrow(enriched_lr) == 0) {

  stop(
    paste0(
      "No selected LR interactions have enrichment > ",
      enrichment_cutoff,
      "."
    )
  )
}


# ==============================================================================
# 14. Select TOP LR interactions for plotting
#
# NEW LOGIC:
#
# Step 1:
#   Require enrichment > enrichment_cutoff
#
# Step 2:
#   Within each selected cell-type pair,
#   select the TOP N interactions by abs_score
#
#
# Therefore:
#
#   enrichment = specificity threshold
#   abs_score   = ranking criterion
#
#
# This prioritizes interactions that are:
#
#   sufficiently cell-pair-enriched
#              +
#   strong in absolute LR score
# ==============================================================================

top_lr_for_plot <- enriched_lr %>%

  group_by(
    source,
    target,
    source_target
  ) %>%

  slice_max(
    order_by = abs_score,
    n = top_n_per_combo,
    with_ties = FALSE
  ) %>%

  ungroup() %>%

  mutate(

    # Preserve desired display order.
    ct_pair = factor(
      source_target,
      levels = combos_of_interest
    ),

    # Plot label.
    lr_pair_label = paste(
      ligand,
      receptor,
      sep = " -> "
    )
  )


# ==============================================================================
# 15. Inspect selected interactions
# ==============================================================================

cat(
  "\n============================================================\n",
  "TOP INTERACTIONS FOR PLOTTING\n",
  "============================================================\n",
  sep = ""
)


cat(
  "Number of LR interactions selected for plotting:",
  nrow(top_lr_for_plot),
  "\n"
)


cat(
  "\nNumber selected per cell-type pair:\n"
)


print(
  top_lr_for_plot %>%

    count(
      source_target,
      name = "n_selected"
    ) %>%

    mutate(
      source_target = factor(
        source_target,
        levels = combos_of_interest
      )
    ) %>%

    arrange(
      source_target
    ),

  n = Inf
)


cat(
  "\nInteractions selected for plotting:\n"
)


print(
  top_lr_for_plot %>%

    select(
      source_target,
      ligand,
      receptor,
      abs_score,
      enrichment,
      mean_score_all_combos,
      n_combos_for_lr,
      n_donors
    ) %>%

    mutate(
      source_target = factor(
        source_target,
        levels = combos_of_interest
      )
    ) %>%

    arrange(
      source_target,
      desc(abs_score)
    ),

  n = Inf
)


# ==============================================================================
# 16. Check whether the same LR appears in multiple selected cell-type pairs
#
# This is allowed.
# ==============================================================================

cat(
  "\nLR pairs appearing in more than one selected cell-type pair:\n"
)


repeated_lr <- top_lr_for_plot %>%

  count(
    lr_pair,
    name = "n_selected_combos"
  ) %>%

  filter(
    n_selected_combos > 1
  ) %>%

  arrange(
    desc(n_selected_combos),
    lr_pair
  )


print(
  repeated_lr,
  n = Inf
)


# ==============================================================================
# 17. Scatter plot
#
# X-axis:
#   Mean LR score across donors
#   -> absolute interaction strength
#
# Y-axis:
#   LR enrichment relative to this LR's mean score
#   across all mapped cell-type pairs
#   -> relative cell-pair specificity
#
# Therefore:
#
#   farther RIGHT = stronger
#   farther UP    = more specifically enriched
#
# The upper-right region contains interactions that
# are both strong and relatively cell-pair-specific.
# ==============================================================================

p <- ggplot(
  top_lr_for_plot,
  aes(
    x = abs_score,
    y = enrichment,
    color = ct_pair,
    label = lr_pair_label
  )
) +

  # Enrichment threshold used during selection.
  geom_hline(
    yintercept = enrichment_cutoff,
    linetype = "dashed"
  ) +

  geom_point(
    size = 3
  ) +

  geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    seed = 123
  ) +

  labs(

    x = "Mean LR score",

    y = paste0(
      "LR enrichment relative to mean across cell-type pairs"
    ),

    color = "Cell-type pair"
  ) +

  theme_bw(
    base_size = 14
  ) +

  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )


print(p)


# ==============================================================================
# 18. Save selected LR data
# ==============================================================================

plot_dir <- here::here(
  "plots",
  "10_HD_bin_level",
  "no_secondary",
  "liana2"
)


dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


selected_data_file <- file.path(
  plot_dir,
  paste0(
    "scatter_",
    analysis_label,
    "_LRenrichmentGT",
    enrichment_cutoff,
    "_top",
    top_n_per_combo,
    "ByMeanScorePerCombo_data.csv"
  )
)


write_csv(
  top_lr_for_plot %>%
    arrange(
      ct_pair,
      desc(abs_score)
    ) %>%
    select(
      source,
      target,
      source_target,
      ligand,
      receptor,
      lr_pair,
      abs_score,
      enrichment,
      mean_score_all_combos,
      n_combos_for_lr,
      n_donors
    ),
  selected_data_file
)


# ==============================================================================
# 19. Save plot as PDF
# ==============================================================================

output_pdf <- file.path(
  plot_dir,
  paste0(
    "scatter_",
    analysis_label,
    "_LRenrichmentGT",
    enrichment_cutoff,
    "_top",
    top_n_per_combo,
    "ByMeanScorePerCombo.pdf"
  )
)


ggsave(
  filename = output_pdf,
  plot = p,
  width = 9,
  height = 6,
  device = "pdf"
)


# ==============================================================================
# 20. Final summary
# ==============================================================================

cat(
  "\n============================================================\n",
  "ANALYSIS COMPLETE\n",
  "============================================================\n",
  sep = ""
)


cat(
  "Number of LR interactions plotted:",
  nrow(top_lr_for_plot),
  "\n"
)


cat(
  "Selection rule:\n",
  "  1. enrichment > ",
  enrichment_cutoff,
  "\n",
  "  2. top ",
  top_n_per_combo,
  " per cell-type pair by mean LR score\n",
  sep = ""
)


cat(
  "\nSaved selected interaction table:\n",
  selected_data_file,
  "\n"
)


cat(
  "\nSaved PDF:\n",
  output_pdf,
  "\n"
)