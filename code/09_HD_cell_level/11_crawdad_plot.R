library(here)
library(tidyverse)
library(crawdad)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
plot_path = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'crawdad', 'dot_plot_%s.pdf'
)
result_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad',
    '%s_results.csv'
)
plot_dir = here('plots', '09_HD_cell_level', 'probe_fix', 'crawdad')

################################################################################
#   Functions
################################################################################

custom_dotplot = function(result_df, filename) {
    p = ggplot(
            result_df, aes(x = reference, y = neighbor, color = Z, size = scale)
        ) +
        geom_point() +
        scale_color_gradient2(low = 'blue', mid = 'white', high = 'red') +
        scale_radius(
            trans = 'reverse',
            breaks = seq(
                min(result_df$scale), max(result_df$scale), length.out = 3
            ),
            range = c(2, 15)
        ) +
        coord_fixed() +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    pdf(file.path(plot_dir, filename), width = 9)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

sample_ids = readLines(sample_id_path)[1:3]

result_list = list()
for (sample_id in sample_ids) {
    result_list[[sample_id]] = sprintf(result_paths, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        mutate(sample_id = sample_id)

    result_list[[sample_id]]$Z_sig = correctZBonferroni(
        result_list[[sample_id]]
    )
}

result_df = do.call(rbind, result_list) |>
    filter(neighbor != reference) |>
    #   First average Z-scores across permutations
    group_by(sample_id, neighbor, scale, reference, Z_sig) |>
    summarize(Z = mean(Z)) |>
    ungroup() |>
    #   Then filter to the smallest spatial scale with significant Z-scores
    filter(abs(Z) >= Z_sig) |>
    group_by(sample_id, neighbor, reference) |>
    filter(scale == min(scale)) |>
    ungroup()

#   Create an order-agnostic pair identifier for each neighbor-reference
#   combination
result_df$pair = sapply(
    seq_len(nrow(result_df)),
    function(i) {
        paste(
            sort(c(result_df$neighbor[i], result_df$reference[i])),
            collapse = '_'
        )
    }
)

result_df = result_df |>
    #   Retain combinations where both directions of association are
    #   significant, and sign of Z scores agree across directions and samples
    group_by(pair) |>
    filter(n() == 2 * length(unique(sample_ids)), all(Z > 0) | all(Z < 0)) |>
    #   Take the maximum scale across samples and the min Z score at that scale
    group_by(neighbor, reference) |>
    summarize(scale = max(scale), Z = min(Z[which.max(scale)])) |>
    ungroup()

custom_dotplot(result_df, 'dot_plot_combined.pdf')

#   Somewhat arbitrarily increase stringency of effect size to narrow results
result_df |>
    filter(abs(Z) >= 8) |>
    custom_dotplot(filename = 'dot_plot_combined_strict.pdf')

session_info()
