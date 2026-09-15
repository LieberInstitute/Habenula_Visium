library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(readr)
library(here)
library(tidyr)
library(ggforce)
library(stringr)

# =========================
# 0. User parameters
# =========================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) {
    return(default)
  }
  if (hit == length(args)) {
    stop("Missing value after ", flag)
  }
  args[hit + 1]
}

# Usage:
# Rscript this_script.R --min_donors 4
# Rscript this_script.R --min_donors 5
n_required_donors <- as.integer(get_arg("--min_donors", 4))

if (is.na(n_required_donors) || n_required_donors < 1) {
  stop("n_required_donors must be a positive integer.")
}

alpha <- 0.05
min_pairs_for_lm <- 10

donor_suffix <- paste0("_min", n_required_donors, "donors")

cat("Minimum donors required for repeated donor-level pattern:", n_required_donors, "\n")


# =========================
# 1. Read data
# =========================

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


# =========================
# 2. Use common donors only
# =========================

cellular_donors <- sort(unique(df_cellular$donor_id))
extracellular_donors <- sort(unique(df_extracellular$donor_id))

common_donors <- intersect(cellular_donors, extracellular_donors)

cat("Cellular donors:", length(cellular_donors), "\n")
cat("Extracellular donors:", length(extracellular_donors), "\n")
cat("Common donors:", length(common_donors), "\n")
cat("Required donors for repeated pattern:", n_required_donors, "\n")
print(common_donors)

df_cellular_use <- df_cellular %>%
  filter(donor_id %in% common_donors)

df_extracellular_use <- df_extracellular %>%
  filter(donor_id %in% common_donors)


# =========================
# 3. Summarize cellular at donor level
# =========================

cellular_donor <- df_cellular_use %>%
  mutate(
    lr_pair = paste(ligand, receptor, sep = "-")
  ) %>%
  group_by(donor_id, source, target, ligand, receptor, lr_pair) %>%
  summarise(
    mean_cellular = mean(mean, na.rm = TRUE),
    pval_cellular_min = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cellular_present = TRUE
  )


# =========================
# 4. Summarize extracellular at donor level
# =========================

extracellular_donor <- df_extracellular_use %>%
  mutate(
    lr_pair = paste(ligand, receptor, sep = "-")
  ) %>%
  group_by(donor_id, source, target, ligand, receptor, lr_pair) %>%
  summarise(
    mean_extracellular = mean(mean, na.rm = TRUE),
    pval_extracellular_min = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    extracellular_present = TRUE
  )


# =========================
# 5. Full join at donor level
#    Missing side is set to 0
# =========================

lr_by_donor_full <- cellular_donor %>%
  full_join(
    extracellular_donor,
    by = c("donor_id", "source", "target", "ligand", "receptor", "lr_pair")
  ) %>%
  mutate(
    cellular_present = replace_na(cellular_present, FALSE),
    extracellular_present = replace_na(extracellular_present, FALSE),

    mean_cellular = replace_na(mean_cellular, 0),
    mean_extracellular = replace_na(mean_extracellular, 0),

    donor_pattern = case_when(
      cellular_present & extracellular_present ~ "shared",
      cellular_present & !extracellular_present ~ "cellular_only",
      !cellular_present & extracellular_present ~ "extracellular_only",
      TRUE ~ "neither"
    )
  )

cat("Donor-level source-target-LR rows after full join:", nrow(lr_by_donor_full), "\n")
cat("Donor-level presence patterns:\n")
print(table(lr_by_donor_full$donor_pattern))


# =========================
# 6. Count reproducible donor-level patterns
# =========================

pattern_summary <- lr_by_donor_full %>%
  group_by(source, target, ligand, receptor, lr_pair) %>%
  summarise(
    n_shared_donor = sum(donor_pattern == "shared"),
    n_cellular_only_donor = sum(donor_pattern == "cellular_only"),
    n_extracellular_only_donor = sum(donor_pattern == "extracellular_only"),

    n_cellular_donor = sum(cellular_present),
    n_extracellular_donor = sum(extracellular_present),
    n_either_donor = n_distinct(donor_id),

    pattern_pass =
      n_shared_donor >= n_required_donors |
      n_cellular_only_donor >= n_required_donors |
      n_extracellular_only_donor >= n_required_donors,

    dominant_pattern = case_when(
      n_shared_donor >= n_cellular_only_donor &
        n_shared_donor >= n_extracellular_only_donor ~ "shared",

      n_cellular_only_donor >= n_shared_donor &
        n_cellular_only_donor >= n_extracellular_only_donor ~ "cellular_only",

      n_extracellular_only_donor >= n_shared_donor &
        n_extracellular_only_donor >= n_cellular_only_donor ~ "extracellular_only",

      TRUE ~ NA_character_
    ),

    dominant_pattern_n = pmax(
      n_shared_donor,
      n_cellular_only_donor,
      n_extracellular_only_donor
    ),

    .groups = "drop"
  )

cat("LR-source-target combinations by pattern-pass status:\n")
print(table(pattern_summary$pattern_pass))

cat("Dominant patterns before filtering:\n")
print(table(pattern_summary$dominant_pattern))


# =========================
# 7. Final LR-level table for regression/scatter
#    Means are averaged after donor-level zero filling
# =========================

df_lr_compare <- lr_by_donor_full %>%
  group_by(source, target, ligand, receptor, lr_pair) %>%
  summarise(
    mean_cellular = mean(mean_cellular, na.rm = TRUE),
    mean_extracellular = mean(mean_extracellular, na.rm = TRUE),

    pval_cellular_min = if (all(is.na(pval_cellular_min))) {
      NA_real_
    } else {
      min(pval_cellular_min, na.rm = TRUE)
    },

    pval_extracellular_min = if (all(is.na(pval_extracellular_min))) {
      NA_real_
    } else {
      min(pval_extracellular_min, na.rm = TRUE)
    },

    .groups = "drop"
  ) %>%
  left_join(
    pattern_summary,
    by = c("source", "target", "ligand", "receptor", "lr_pair")
  ) %>%
  filter(
    pattern_pass,
    is.finite(mean_cellular),
    is.finite(mean_extracellular)
  ) %>%
  mutate(
    presence_group = dominant_pattern
  )

cat("Source-target-LR pairs after donor-level pattern filter:", nrow(df_lr_compare), "\n")
cat("Dominant pattern groups after filter:\n")
print(table(df_lr_compare$presence_group))


# =========================
# 8. Log transform
# =========================

positive_values <- c(
  df_lr_compare$mean_cellular,
  df_lr_compare$mean_extracellular
)

eps <- min(positive_values[positive_values > 0], na.rm = TRUE) / 2

cat("Pseudo-count eps for log transform:", eps, "\n")

df_lr_compare_log <- df_lr_compare %>%
  mutate(
    log_cellular = log10(mean_cellular + eps),
    log_extracellular = log10(mean_extracellular + eps)
  )


# =========================
# 9. Regression outlier diagnostics
# =========================

run_response_outlier_diagnostics <- function(dat, alpha = 0.05, min_pairs = 10) {

  dat <- dat %>%
    filter(
      is.finite(log_cellular),
      is.finite(log_extracellular)
    )

  n <- nrow(dat)

  if (
    n < min_pairs ||
    n_distinct(dat$log_cellular) < 2 ||
    n_distinct(dat$log_extracellular) < 2
  ) {
    return(
      dat %>%
        mutate(
          fitted_log_extracellular = NA_real_,
          residual = NA_real_,
          studentized_residual = NA_real_,
          studentized_p = NA_real_,
          lm_intercept = NA_real_,
          lm_slope = NA_real_,
          lm_r2 = NA_real_,
          n_pairs_in_combo = n,
          outlier_direction = NA_character_
        )
    )
  }

  fit <- lm(log_extracellular ~ log_cellular, data = dat)
  sm <- summary(fit)

  stud_res <- rstudent(fit)

  p <- length(coef(fit))
  df_resid_ext <- n - p - 1

  stud_p <- 2 * pt(
    abs(stud_res),
    df = df_resid_ext,
    lower.tail = FALSE
  )

  dat %>%
    mutate(
      fitted_log_extracellular = fitted(fit),
      residual = resid(fit),
      studentized_residual = as.numeric(stud_res),
      studentized_p = as.numeric(stud_p),
      lm_intercept = unname(coef(fit)[1]),
      lm_slope = unname(coef(fit)[2]),
      lm_r2 = sm$r.squared,
      n_pairs_in_combo = n,

      outlier_direction = case_when(
        studentized_residual > 0 ~ "extracellular higher than expected",
        studentized_residual < 0 ~ "extracellular lower than expected",
        TRUE ~ NA_character_
      )
    )
}


df_lr_response_diag <- df_lr_compare_log %>%
  group_by(source, target) %>%
  group_modify(
    ~ run_response_outlier_diagnostics(
      dat = .x,
      alpha = alpha,
      min_pairs = min_pairs_for_lm
    )
  ) %>%
  ungroup()


# =========================
# 10. Within source-target multiple testing correction
# =========================

df_lr_response_diag <- df_lr_response_diag %>%
  group_by(source, target) %>%
  mutate(
    studentized_p_fdr_within_combo = if_else(
      !is.na(studentized_p),
      p.adjust(studentized_p, method = "BH"),
      NA_real_
    ),

    studentized_p_bonf_within_combo = if_else(
      !is.na(studentized_p),
      p.adjust(studentized_p, method = "bonferroni"),
      NA_real_
    ),

    is_response_outlier_fdr_within_combo =
      studentized_p_fdr_within_combo < alpha,

    is_response_outlier_bonf_within_combo =
      studentized_p_bonf_within_combo < alpha
  ) %>%
  ungroup()


# =========================
# 11. Count outliers
# =========================

outlier_count_check <- df_lr_response_diag %>%
  summarise(
    n_total_lr_combo_tests = sum(!is.na(studentized_p)),
    n_within_combo_fdr = sum(is_response_outlier_fdr_within_combo, na.rm = TRUE),
    n_within_combo_bonf = sum(is_response_outlier_bonf_within_combo, na.rm = TRUE)
  )

print(outlier_count_check)

outlier_count_by_presence <- df_lr_response_diag %>%
  filter(!is.na(studentized_p)) %>%
  group_by(presence_group) %>%
  summarise(
    n_tests = n(),
    n_outlier_fdr = sum(is_response_outlier_fdr_within_combo, na.rm = TRUE),
    n_outlier_bonf = sum(is_response_outlier_bonf_within_combo, na.rm = TRUE),
    .groups = "drop"
  )

print(outlier_count_by_presence)


# =========================
# 12. Final within-combo FDR outlier table
# =========================

response_outliers_main <- df_lr_response_diag %>%
  filter(is_response_outlier_fdr_within_combo) %>%
  mutate(
    confidence_tier = case_when(
      is_response_outlier_bonf_within_combo ~ "within-combo Bonferroni",
      TRUE ~ "within-combo FDR"
    )
  ) %>%
  arrange(
    source,
    target,
    desc(is_response_outlier_bonf_within_combo),
    studentized_p_fdr_within_combo,
    desc(abs(studentized_residual))
  ) %>%
  select(
    confidence_tier,
    presence_group,
    n_shared_donor,
    n_cellular_only_donor,
    n_extracellular_only_donor,
    dominant_pattern_n,
    source, target,
    ligand, receptor, lr_pair,
    mean_cellular, mean_extracellular,
    log_cellular, log_extracellular,
    fitted_log_extracellular,
    residual,
    studentized_residual,
    studentized_p,
    studentized_p_fdr_within_combo,
    studentized_p_bonf_within_combo,
    outlier_direction,
    lm_slope,
    lm_intercept,
    lm_r2,
    n_pairs_in_combo,
    pval_cellular_min,
    pval_extracellular_min,
    n_cellular_donor,
    n_extracellular_donor,
    n_either_donor
  )

print(head(response_outliers_main, 30))

write_csv(
  response_outliers_main,
  file.path(
    output_dir,
    paste0(
      "cellular_extracellular_response_outliers_within_combo_FDR_donor_pattern_full_join",
      donor_suffix,
      ".csv"
    )
  )
)

write_csv(
  df_lr_response_diag,
  file.path(
    output_dir,
    paste0(
      "cellular_extracellular_response_outlier_diagnostics_within_combo_FDR_donor_pattern_full_join",
      donor_suffix,
      ".csv"
    )
  )
)

write_csv(
  pattern_summary,
  file.path(
    output_dir,
    paste0(
      "cellular_extracellular_donor_level_presence_pattern_summary_full_join",
      donor_suffix,
      ".csv"
    )
  )
)

write_csv(
  df_lr_compare,
  file.path(
    output_dir,
    paste0(
      "cellular_extracellular_LR_compare_after_donor_pattern_filter_full_join",
      donor_suffix,
      ".csv"
    )
  )
)


# =========================
# 13. Facet plots: all source-target combos with within-combo FDR outliers
# =========================

plot_all_combo_dat <- df_lr_response_diag %>%
  filter(!is.na(studentized_p)) %>%
  group_by(source, target) %>%
  mutate(
    n_outliers_combo = sum(is_response_outlier_fdr_within_combo, na.rm = TRUE),
    n_points_combo = n()
  ) %>%
  ungroup() %>%
  filter(n_outliers_combo > 0) %>%
  mutate(
    combo = paste(source, "->", target)
  )

if (nrow(plot_all_combo_dat) > 0) {

  combo_order <- plot_all_combo_dat %>%
    distinct(source, target, combo, n_outliers_combo, n_points_combo) %>%
    arrange(desc(n_outliers_combo), desc(n_points_combo), combo) %>%
    pull(combo)

  plot_all_combo_dat <- plot_all_combo_dat %>%
    mutate(
      combo = factor(combo, levels = combo_order)
    )

  combo_label_df <- plot_all_combo_dat %>%
    distinct(source, target, combo, n_outliers_combo, n_points_combo) %>%
    mutate(
      combo_label = paste0(
        combo,
        "\n(n=", n_points_combo, ", outliers=", n_outliers_combo, ")"
      )
    )

  plot_all_combo_dat <- plot_all_combo_dat %>%
    left_join(
      combo_label_df %>% select(combo, combo_label),
      by = "combo"
    )

  combo_label_order <- combo_label_df %>%
    arrange(desc(n_outliers_combo), desc(n_points_combo), combo) %>%
    pull(combo_label)

  plot_all_combo_dat <- plot_all_combo_dat %>%
    mutate(
      combo_label = factor(combo_label, levels = combo_label_order)
    )

  ncol_plot <- 2
  nrow_plot <- 5
  panels_per_page <- ncol_plot * nrow_plot

  n_pages <- ceiling(length(combo_label_order) / panels_per_page)

  cat("Number of source-target combos with within-combo FDR outliers:", length(combo_label_order), "\n")
  cat("Total figure pages:", n_pages, "\n")


  make_facet_page <- function(page_num) {

    ggplot(
      plot_all_combo_dat,
      aes(x = log_cellular, y = log_extracellular)
    ) +
      geom_point(
        data = plot_all_combo_dat %>%
          filter(!is_response_outlier_fdr_within_combo),
        aes(shape = presence_group),
        color = "grey75",
        alpha = 0.6,
        size = 1.2
      ) +
      geom_point(
        data = plot_all_combo_dat %>%
          filter(is_response_outlier_fdr_within_combo),
        aes(shape = presence_group),
        color = "red",
        alpha = 0.9,
        size = 1.8
      ) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        color = "black",
        linewidth = 0.5
      ) +
      ggrepel::geom_text_repel(
        data = plot_all_combo_dat %>%
          filter(is_response_outlier_fdr_within_combo),
        aes(label = lr_pair),
        color = "red",
        size = 2.5,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.15,
        segment.size = 0.2,
        show.legend = FALSE
      ) +
      ggforce::facet_wrap_paginate(
        ~ combo_label,
        ncol = ncol_plot,
        nrow = nrow_plot,
        page = page_num,
        scales = "free"
      ) +
      labs(
        title = paste0(
          "Cellular vs extracellular LR scores by source-target combo (Page ",
          page_num, "/", n_pages, ")"
        ),
        subtitle = paste0(
          "Donor-level full join: missing side set to 0; kept patterns repeated in at least ",
          n_required_donors,
          " donors. Red points are within-combo FDR response outliers."
        ),
        x = "log10 cellular mean score",
        y = "log10 extracellular mean score",
        shape = "Donor-level pattern"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 8),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )
  }


  pdf(
    file = file.path(
      output_dir,
      paste0(
        "all_combo_response_outliers_within_combo_FDR_facets_donor_pattern_full_join",
        donor_suffix,
        ".pdf"
      )
    ),
    width = 10,
    height = 16
  )

  for (i in seq_len(n_pages)) {
    print(make_facet_page(i))
  }

  dev.off()

} else {

  cat("No within-combo FDR outliers found. Facet PDF was not generated.\n")
}