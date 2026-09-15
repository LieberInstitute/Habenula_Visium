# Compare extracellular and cellular heatmap raw or normalized scores in a scatter plot to identify pairs of cells types that communicate differently

library(ggplot2)
library(ggrepel)
library(dplyr)
library(purrr)
library(readr)
library(here)
library(tidyr)

df_cellular <- read_csv(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/plots/10_HD_bin_level/no_secondary/liana2/source_target_sum_mean_heatmap_data_5000.0.csv",
  show_col_types = FALSE
)

df_extracellular <- read_csv(
  "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/plots/10_HD_bin_level/no_secondary/liana2/extracellular/source_target_sum_mean_heatmap_data_5000.0.csv",
  show_col_types = FALSE
)

output_dir <- here(
  "plots", "10_HD_bin_level", "no_secondary", "liana2", "cellular_vs_extracellular"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

library(tidyverse)
library(ggrepel)

# cellular long
cell_long <- df_cellular %>%
  pivot_longer(
    cols = -target,
    names_to = "source",
    values_to = "cellular_score"
  )

# extracellular long
extra_long <- df_extracellular %>%
  pivot_longer(
    cols = -target,
    names_to = "source",
    values_to = "extracellular_score"
  )

# merge
plot_dat <- cell_long %>%
  inner_join(extra_long, by = c("source", "target")) %>%
  mutate(
    pair = paste(source, "→", target),
    diff = extracellular_score - cellular_score,
    abs_diff = abs(diff),
    log2_ratio = log2((extracellular_score + 1e-6) / (cellular_score + 1e-6))
  )

head(plot_dat)

p_raw <- ggplot(plot_dat, aes(x = cellular_score, y = extracellular_score)) +
  geom_point(aes(color = abs_diff), size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient(low = "grey70", high = "red") +
  theme_bw(base_size = 12) +
  labs(
    x = "Cellular score",
    y = "Extracellular score",
    color = "|Difference|",
    title = "Extracellular vs Cellular communication scores",
    subtitle = "Each point is one source-target cell-type pair"
  )

ggsave(
  filename = file.path(output_dir, "cellular_vs_extracellular_scatter_raw.png"),
  plot = p_raw,
  width = 8,
  height = 6
)

# normalize within each dataset separately
cell_range <- range(plot_dat$cellular_score, na.rm = TRUE)
extra_range <- range(plot_dat$extracellular_score, na.rm = TRUE)

plot_dat_norm <- plot_dat %>%
  mutate(
    cellular_norm = (cellular_score - cell_range[1]) / (cell_range[2] - cell_range[1]),
    extracellular_norm = (extracellular_score - extra_range[1]) / (extra_range[2] - extra_range[1]),
    norm_diff = extracellular_norm - cellular_norm,
    abs_norm_diff = abs(norm_diff)
  )

top_lab_norm <- plot_dat_norm %>%
  slice_max(order_by = abs_norm_diff, n = 15)

p_norm <- ggplot(plot_dat_norm, aes(x = cellular_norm, y = extracellular_norm)) +
  geom_point(aes(color = abs_norm_diff), size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = top_lab_norm,
    aes(label = pair),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_gradient(low = "grey70", high = "blue") +
  coord_equal() +
  theme_bw(base_size = 12) +
  labs(
    x = "Normalized cellular score",
    y = "Normalized extracellular score",
    color = "|Normalized\nDifference|",
    title = "Extracellular vs Cellular normalized communication scores",
    subtitle = "Points far from the diagonal indicate different communication patterns"
  )

ggsave(
  filename = file.path(output_dir, "cellular_vs_extracellular_scatter_normalized.png"),
  plot = p_norm,
  width = 8,
  height = 6
)

plot_dat_z <- plot_dat %>%
  mutate(
    cellular_z = as.numeric(scale(cellular_score)),
    extracellular_z = as.numeric(scale(extracellular_score)),
    z_diff = extracellular_z - cellular_z,
    abs_z_diff = abs(z_diff)
  )

top_lab_z <- plot_dat_z %>%
  slice_max(order_by = abs_z_diff, n = 15)

p_z <- ggplot(plot_dat_z, aes(x = cellular_z, y = extracellular_z)) +
  geom_point(aes(color = abs_z_diff), size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(
    data = top_lab_z,
    aes(label = pair),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_gradient(low = "grey70", high = "purple") +
  theme_bw(base_size = 12) +
  labs(
    x = "Cellular z-score",
    y = "Extracellular z-score",
    color = "|Z-score\nDifference|",
    title = "Extracellular vs Cellular communication (z-score normalized)"
  )

ggsave(
  filename = file.path(output_dir, "cellular_vs_extracellular_scatter_zscore.png"),
  plot = p_z,
  width = 8,
  height = 6
)

# -----------------------------

library(tidyverse)
library(ggrepel)

plot_dat <- cell_long %>%
  inner_join(extra_long, by = c("source", "target")) %>%
  mutate(
    pair = paste(source, "→", target),
    diff = extracellular_score - cellular_score,
    abs_diff = abs(diff),
    log2_ratio = log2((extracellular_score + 1e-6) / (cellular_score + 1e-6))
  )

# Fit regression: extracellular ~ cellular
fit_raw <- lm(
  extracellular_score ~ cellular_score,
  data = plot_dat,
  na.action = na.exclude
)

plot_dat_raw_resid <- plot_dat %>%
  mutate(
    pred_extracellular = predict(fit_raw),
    resid = residuals(fit_raw),
    std_resid = rstandard(fit_raw),
    abs_std_resid = abs(std_resid),
    residual_direction = ifelse(
      resid > 0,
      "Extracellular higher than expected",
      "Extracellular lower than expected"
    )
  )

top_lab_raw <- plot_dat_raw_resid %>%
  slice_max(order_by = abs_std_resid, n = 15)

p_raw <- ggplot(
  plot_dat_raw_resid,
  aes(x = cellular_score, y = extracellular_score)
) +
  geom_point(
    aes(color = abs_std_resid),
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = top_lab_raw,
    aes(label = pair),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_gradient(
    low = "grey70",
    high = "red",
    name = "|Standardized\nresidual|"
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Cellular score",
    y = "Extracellular score",
    title = "Extracellular vs Cellular communication scores",
    subtitle = "Labeled points are top outliers based on regression residuals"
  )

ggsave(
  filename = file.path(output_dir, "cellular_vs_extracellular_scatter_raw_regression_outliers.png"),
  plot = p_raw,
  width = 8,
  height = 6
)

cell_range <- range(plot_dat$cellular_score, na.rm = TRUE)
extra_range <- range(plot_dat$extracellular_score, na.rm = TRUE)

plot_dat_norm <- plot_dat %>%
  mutate(
    cellular_norm = (cellular_score - cell_range[1]) / 
      (cell_range[2] - cell_range[1]),
    extracellular_norm = (extracellular_score - extra_range[1]) / 
      (extra_range[2] - extra_range[1])
  )

fit_norm <- lm(
  extracellular_norm ~ cellular_norm,
  data = plot_dat_norm,
  na.action = na.exclude
)

plot_dat_norm_resid <- plot_dat_norm %>%
  mutate(
    pred_extracellular_norm = predict(fit_norm),
    resid_norm = residuals(fit_norm),
    std_resid_norm = rstandard(fit_norm),
    abs_std_resid_norm = abs(std_resid_norm),
    residual_direction = ifelse(
      resid_norm > 0,
      "Extracellular higher than expected",
      "Extracellular lower than expected"
    )
  )

top_lab_norm <- plot_dat_norm_resid %>%
  slice_max(order_by = abs_std_resid_norm, n = 15)

p_norm <- ggplot(
  plot_dat_norm_resid,
  aes(x = cellular_norm, y = extracellular_norm)
) +
  geom_point(
    aes(color = abs_std_resid_norm),
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = top_lab_norm,
    aes(label = pair),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_gradient(
    low = "grey70",
    high = "blue",
    name = "|Standardized\nresidual|"
  ) +
  coord_equal() +
  theme_bw(base_size = 12) +
  labs(
    x = "Normalized cellular score",
    y = "Normalized extracellular score",
    title = "Extracellular vs Cellular normalized communication scores",
    subtitle = "Outliers are defined by residuals from extracellular ~ cellular regression"
  )

ggsave(
  filename = file.path(output_dir, "cellular_vs_extracellular_scatter_normalized_regression_outliers.png"),
  plot = p_norm,
  width = 8,
  height = 6
)

plot_dat_z <- plot_dat %>%
  mutate(
    cellular_z = as.numeric(scale(cellular_score)),
    extracellular_z = as.numeric(scale(extracellular_score))
  )

fit_z <- lm(
  extracellular_z ~ cellular_z,
  data = plot_dat_z,
  na.action = na.exclude
)

plot_dat_z_resid <- plot_dat_z %>%
  mutate(
    pred_extracellular_z = predict(fit_z),
    resid_z = residuals(fit_z),
    std_resid_z = rstandard(fit_z),
    abs_std_resid_z = abs(std_resid_z),
    residual_direction = ifelse(
      resid_z > 0,
      "Extracellular higher than expected",
      "Extracellular lower than expected"
    )
  )

top_lab_z <- plot_dat_z_resid %>%
  slice_max(order_by = abs_std_resid_z, n = 15)

p_z <- ggplot(
  plot_dat_z_resid,
  aes(x = cellular_z, y = extracellular_z)
) +
  geom_point(
    aes(color = abs_std_resid_z),
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  ggrepel::geom_text_repel(
    data = top_lab_z,
    aes(label = pair),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_gradient(
    low = "grey70",
    high = "purple",
    name = "|Standardized\nresidual|"
  ) +
  theme_bw(base_size = 12) +
  labs(
    x = "Cellular z-score",
    y = "Extracellular z-score",
    title = "Extracellular vs Cellular communication scores",
    subtitle = "Z-score normalized; outliers are based on regression residuals"
  )

ggsave(
  filename = file.path(output_dir, "cellular_vs_extracellular_scatter_zscore_regression_outliers.png"),
  plot = p_z,
  width = 8,
  height = 6
)
