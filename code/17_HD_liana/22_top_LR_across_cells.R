# Generate a CSV file to share with the team (likely on this channel) of top significant LR interactions. 
# 1. For the cell-type pair, take LR pairs in the top 20% (of LR interactions within that cell-type pair) by absolute score
# 2. Drop LR pairs in the bottom 50% of all interactions across cell-type pairs
# 3. Sort by mean ratio (it's ok to have < 1)
# 4. For the CSV, let's just show the top 5 per cell-type pair

# just make these supplementary tables
library(readr)
library(dplyr)
library(purrr)

input_file <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana2/table/significant_interactions_across_donors_5000.0.csv"

output_file <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana2/table/top_significant_LR_interactions_5000.0.csv"

top_prop <- 0.20
top_n_per_combo <- 5
global_keep_prop <- 0.50

# If you want LR pairs present in all donors, use 8.
min_n_donors <- 4

df <- read_csv(input_file, show_col_types = FALSE)

# -----------------------------
# 1. Prepare labels
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
    donor_mean_score = mean(mean, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 3. Average across donors
# -----------------------------
lr_score <- lr_by_donor %>%
  group_by(source_target, source, target, lr_pair, ligand, receptor) %>%
  summarise(
    abs_score = mean(donor_mean_score, na.rm = TRUE),
    n_donors = n_distinct(donor_id),
    donor_ids = paste(sort(unique(donor_id)), collapse = ";"),
    min_pval = min(min_pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_donors >= min_n_donors)

# -----------------------------
# 4. Calculate mean ratio for each LR pair
#    ratio = score in this source-target combo /
#            strongest score in other source-target combos
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
# 5. Within each source-target combo:
#    A. top 20% within each source-target combo by abs_score
#    B. remove bottom 50% globally by abs_score
#    C. take top 5 per combo by mean_ratio
# -----------------------------

top_LR <- mean_ratio_lr %>%
  group_by(source, target, source_target) %>%
  mutate(
    abs_score_cutoff = quantile(
      abs_score,
      probs = 1 - top_prop,
      na.rm = TRUE
    ),
    n_LR_in_source_target = n()
  ) %>%
  filter(abs_score >= abs_score_cutoff) %>%
  ungroup() %>%
  mutate(
    global_abs_score_cutoff = quantile(
      abs_score,
      probs = 1 - global_keep_prop,
      na.rm = TRUE
    )
  ) %>%
  filter(abs_score >= global_abs_score_cutoff) %>%
  group_by(source, target, source_target) %>%
  slice_max(
    order_by = mean_ratio,
    n = top_n_per_combo,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(source, target, desc(mean_ratio)) %>%
  select(
    source,
    target,
    ligand,
    receptor,
    lr_pair,
    abs_score,
    mean_ratio
  )

# Add the MDD and SUD results into the table
# read in the MDD and SUD results
library(here)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana2',
    'open_targets'
)
mdd_file <- read.csv(file.path(out_path, "MDD_LorR_interactions.csv"))
sud_file <- read.csv(file.path(out_path, "Substance_dependence_LorR_interactions.csv"))

library(dplyr)

# 1. Create a combined "ligand-receptor" pair column in both the MDD/SUD files and the top_LR table for easy matching. Then, create a set of these pairs from the MDD and SUD files for quick lookup.
mdd_pairs <- mdd_file %>%
  mutate(lr_pair_check = paste(genesymbol_intercell_source,
                               genesymbol_intercell_target,
                               sep = "-")) %>%
  pull(lr_pair_check) %>%
  unique()

sud_pairs <- sud_file %>%
  mutate(lr_pair_check = paste(genesymbol_intercell_source,
                               genesymbol_intercell_target,
                               sep = "-")) %>%
  pull(lr_pair_check) %>%
  unique()

# 2. Annotate the top_LR table with new columns indicating whether each LR pair is a risk factor for MDD and/or SUD based on the presence of the pair in the respective sets created above.
top_LR_annotated <- top_LR %>%
  mutate(
    lr_pair_check = paste(ligand, receptor, sep = "-"),
    risk_MDD = ifelse(lr_pair_check %in% mdd_pairs, 1, 0),
    risk_SUD = ifelse(lr_pair_check %in% sud_pairs, 1, 0)
  ) %>%
  select(-lr_pair_check)

head(top_LR_annotated)

write_csv(top_LR_annotated, output_file)

cat("Saved file to:\n", output_file, "\n")
cat("Number of selected LR interactions:", nrow(top_LR_annotated), "\n")
