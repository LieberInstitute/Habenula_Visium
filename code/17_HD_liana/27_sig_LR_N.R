# ---------------------------------------------------------
# Libraries
# ---------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------

input_file <- paste0(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
  "Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/",
  "liana2/table/significant_interactions_across_donors_5000.0.csv"
)

output_dir <- paste0(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
  "Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/",
  "liana2/table"
)

CELL_TYPE_MAP_FILE <- paste0(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
  "Habenula_Visium/raw-data/sample_info/",
  "hd_cell_type_map.csv"
)


# ---------------------------------------------------------
# Read data
# ---------------------------------------------------------

df <- read_csv(
  input_file,
  show_col_types = FALSE
)

cell_type_map <- read_csv(
  CELL_TYPE_MAP_FILE,
  show_col_types = FALSE
)

cat(
  "Cell type map columns:",
  paste(names(cell_type_map), collapse = ", "),
  "\n"
)


# ---------------------------------------------------------
# Clean sample IDs
#
# donor_id:
#   Br8433_2_5000.0
#
# becomes:
#   Br8433_2
# ---------------------------------------------------------

df <- df %>%
  rename(
    sample_id = donor_id
  ) %>%
  mutate(
    sample_id = sub(
      "_5000\\.0$",
      "",
      sample_id
    )
  )


# ---------------------------------------------------------
# Map cell type names
#
# Assumes:
#   first column  = old cell type name
#   second column = new cell type name
#
# If your map file uses a different structure,
# change these two lines.
# ---------------------------------------------------------

old_cell_type_col <- names(cell_type_map)[1]
new_cell_type_col <- names(cell_type_map)[2]

cell_type_lookup <- cell_type_map %>%
  transmute(
    old_cell_type = .data[[old_cell_type_col]],
    new_cell_type = .data[[new_cell_type_col]]
  ) %>%
  distinct()


# Map source cell types
df <- df %>%
  left_join(
    cell_type_lookup %>%
      rename(
        source = old_cell_type,
        source_new = new_cell_type
      ),
    by = "source"
  ) %>%
  mutate(
    source = coalesce(source_new, source)
  ) %>%
  select(
    -source_new
  )


# Map target cell types
df <- df %>%
  left_join(
    cell_type_lookup %>%
      rename(
        target = old_cell_type,
        target_new = new_cell_type
      ),
    by = "target"
  ) %>%
  mutate(
    target = coalesce(target_new, target)
  ) %>%
  select(
    -target_new
  )


# ---------------------------------------------------------
# Check mapped cell types
# ---------------------------------------------------------

cat("\nSource cell types after mapping:\n")
print(
  sort(unique(df$source))
)

cat("\nTarget cell types after mapping:\n")
print(
  sort(unique(df$target))
)


# ---------------------------------------------------------
# 1. Count sample support for each LR interaction
#
# Unique interaction:
# source + target + ligand + receptor
# ---------------------------------------------------------

interaction_support <- df %>%
  distinct(
    sample_id,
    source,
    target,
    ligand,
    receptor
  ) %>%
  group_by(
    source,
    target,
    ligand,
    receptor
  ) %>%
  summarise(
    n_samples_significant = n_distinct(sample_id),

    sample_ids = paste(
      sort(unique(sample_id)),
      collapse = ";"
    ),

    .groups = "drop"
  ) %>%
  mutate(
    source_target = paste(
      source,
      target,
      sep = " -> "
    ),

    lr_pair = paste(
      ligand,
      receptor,
      sep = "-"
    )
  )


n_total_samples <- n_distinct(
  df$sample_id
)

cat(
  "Total number of samples:",
  n_total_samples,
  "\n"
)

cat(
  "Total number of unique interactions:",
  nrow(interaction_support),
  "\n"
)


# ---------------------------------------------------------
# 2. Count interactions present in AT LEAST k samples
# ---------------------------------------------------------

cutoff_summary <- tibble(
  min_samples_required = seq_len(n_total_samples)
) %>%
  rowwise() %>%
  mutate(
    n_significant_interactions = sum(
      interaction_support$n_samples_significant >=
        min_samples_required
    )
  ) %>%
  ungroup()

print(cutoff_summary)


# ---------------------------------------------------------
# 3. Plot
# ---------------------------------------------------------

p <- ggplot(
  cutoff_summary,
  aes(
    x = min_samples_required,
    y = n_significant_interactions
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_x_continuous(
    breaks = seq_len(n_total_samples)
  ) +
  labs(
    x = "Minimum number of samples in which the pair is significant",
    y = "Number of unique significant pairs",
    title = "Significant LR pairs across sample cutoffs",
    subtitle = paste0(
      "Each interaction is defined by source, target, ligand, and receptor; ",
      n_total_samples,
      " samples total"
    )
  ) +
  theme_classic(
    base_size = 13
  )

print(p)


# ---------------------------------------------------------
# 4. Save cutoff results
# ---------------------------------------------------------

write_csv(
  cutoff_summary,
  file.path(
    output_dir,
    "LR_pair_counts_by_sample_cutoff_5000.0.csv"
  )
)

ggsave(
  file.path(
    output_dir,
    "LR_pair_counts_by_sample_cutoff_5000.0.pdf"
  ),
  plot = p,
  width = 7,
  height = 5
)


# ---------------------------------------------------------
# 5. Distribution of exact sample support
# ---------------------------------------------------------

support_distribution <- interaction_support %>%
  count(
    n_samples_significant,
    name = "n_unique_interactions"
  ) %>%
  complete(
    n_samples_significant = seq_len(n_total_samples),
    fill = list(
      n_unique_interactions = 0
    )
  )

print(
  support_distribution
)


p_support <- ggplot(
  support_distribution,
  aes(
    x = n_samples_significant,
    y = n_unique_interactions
  )
) +
  geom_col() +
  scale_x_continuous(
    breaks = seq_len(n_total_samples)
  ) +
  labs(
    x = "Exact number of samples in which the interaction is significant",
    y = "Number of unique source–target–LR interactions",
    title = "Cross-sample reproducibility of significant LR interactions"
  ) +
  theme_classic(
    base_size = 13
  )

print(
  p_support
)

ggsave(
  file.path(
    output_dir,
    "LR_pair_support_distribution_5000.0.pdf"
  ),
  plot = p_support,
  width = 7,
  height = 5
)


# ---------------------------------------------------------
# 6. Identify interactions significant in >=4 samples
# ---------------------------------------------------------

sig_interactions_ge4 <- interaction_support %>%
  filter(
    n_samples_significant >= 4
  ) %>%
  arrange(
    desc(n_samples_significant),
    source,
    target,
    ligand,
    receptor
  )

cat(
  "Number of unique interactions significant in >=4 samples:",
  nrow(sig_interactions_ge4),
  "\n"
)


# ---------------------------------------------------------
# 7. Calculate global communication score
#
# For each retained interaction:
# mean of the "mean" communication score across samples
# in which that interaction was significant.
# ---------------------------------------------------------

sig_interactions_ge4_summary <- df %>%
  semi_join(
    sig_interactions_ge4,
    by = c(
      "source",
      "target",
      "ligand",
      "receptor"
    )
  ) %>%
  group_by(
    source,
    target,
    ligand,
    receptor
  ) %>%
  summarise(
    global_communication_score = mean(
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(
    sig_interactions_ge4,
    by = c(
      "source",
      "target",
      "ligand",
      "receptor"
    )
  ) %>%
  select(
    source,
    target,
    ligand,
    receptor,
    lr_pair,
    source_target,
    n_samples_significant,
    sample_ids,
    global_communication_score
  ) %>%
  arrange(
    desc(n_samples_significant),
    desc(global_communication_score)
  )

print(
  sig_interactions_ge4_summary
)


# ---------------------------------------------------------
# 8. Save final interaction table
# ---------------------------------------------------------

write_csv(
  sig_interactions_ge4_summary,
  file.path(
    output_dir,
    "significant_LR_interactions_ge4_samples_5000.0.csv"
  )
)