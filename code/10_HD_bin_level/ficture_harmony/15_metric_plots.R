library(here)
library(tidyverse)
library(sessioninfo)

in_paths = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'deciding_k', '%s_%s.csv'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'new_samples2', 'ficture_harmony', 'deciding_k'
)
ref_names = c('snRNAseq_fine', 'snRNAseq_broad')
input_methods = c('normalized', 'cleaning_y')

dir.create(plot_dir, showWarnings = FALSE)

score_df_list = list()
for (ref_name in ref_names) {
    for (input_method in input_methods) {
        score_df_list[[sprintf('%s_%s', ref_name, input_method)]] = sprintf(
                in_paths, ref_name, input_method
            ) |>
            read_csv(show_col_types = FALSE) |>
            mutate(
                ref_name = ref_name,
                input_method = input_method
            )
    }
}
score_df = bind_rows(score_df_list)

#   Final score vs k by FICTURE input method and reference
p = ggplot(
        score_df,
        aes(x = k, y = final_score, color = input_method, group = input_method)
    ) +
    geom_line() +
    geom_point() +
    facet_wrap(~ ref_name, nrow = 2) +
    labs(y = 'Final Score', color = 'Input Method') +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'score_vs_k_line.pdf'), width = 10)
print(p)
dev.off()

p = score_df |>
    select(-c(frac_ambig, num_clusters, final_score)) |>
    pivot_longer(
        cols = c(frac_unique_hb, frac_unique_non_hb, frac_shared_hb),
        names_to = "metric",
        values_to = "value"
    ) |>
    group_by(ref_name, input_method, metric) |>
    arrange(desc(value)) |>
    slice_head(n = 3) |>
    mutate(k = reorder(k, value, decreasing = TRUE)) |>
    ungroup() |>
    ggplot(aes(x = k, y = value, fill = input_method)) +
        geom_bar(stat = "identity") +
        facet_grid(ref_name ~ metric, scales = "free_x") +
        theme_bw(base_size = 15)
pdf(file.path(plot_dir, 'three_metrics_bar.pdf'), width = 10)
print(p)
dev.off()
