#   CleaningY vs. library-size-normalized FICTURE results comparison, for a
#   supplemental figure

library(here)
library(tidyverse)
library(sessioninfo)

cleany_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'deciding_k', 'multiome_cleaning_y.csv'
)
norm_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'deciding_k', 'multiome_normalized.csv'
)
plot_dir = here('plots', '14_supp_tables', 'ficture_comparison')
metrics = c(
    frac_unique_non_hb = 'Unique Non-Hb',
    frac_unique_hb = 'Unique Hb',
    frac_shared_hb = 'Shared Hb',
    final_score = 'Overall Score'
)

dir.create(plot_dir, showWarnings = FALSE)

cleany_df = read_csv(cleany_path, show_col_types = FALSE) |>
    mutate(dataset = 'cleaningY')
norm_df = read_csv(norm_path, show_col_types = FALSE) |>
    mutate(dataset = 'normalized')

stopifnot(setequal(cleany_df$k, norm_df$k))

metric_df = bind_rows(cleany_df, norm_df) |>
    select(k, dataset, all_of(names(metrics))) |>
    pivot_longer(
        cols = all_of(names(metrics)),
        names_to = 'metric',
        values_to = 'value'
    ) |>
    mutate(metric = metrics[metric])

top3 <- metric_df |>
    group_by(metric, dataset) |>
    slice_max(value, n = 3, with_ties = FALSE) |>
    ungroup() |>
    # globally unique key (metric included) so per-facet ordering works
    mutate(x_key = paste(metric, dataset, k, sep = "___"))

top3 <- top3 |>
    arrange(metric, desc(value), dataset) |>
    mutate(x_key = factor(x_key, levels = unique(x_key)))

p <- ggplot(top3, aes(x = x_key, y = value, fill = dataset)) +
    geom_col() +
    scale_x_discrete(labels = function(x) sub("^.*___.*___", "", x)) +
    facet_wrap(~ metric, scales = "free_x", ncol = 4) +
    labs(x = "FICTURE k Value", y = "Score", fill = "Dataset") +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_dir, "top3_by_metric.pdf"), p, width = 9, height = 5)

session_info()
