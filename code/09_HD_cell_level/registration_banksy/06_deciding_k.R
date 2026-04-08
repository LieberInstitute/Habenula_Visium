#   This script provides a data-driven method for ranking spatial registration
#   results across all tested resolutions, ultimately deciding upon an optimal
#   clustering resolution from Banksy. It ultimately uses the multiome
#   reference data, favors covering many cell types (especially uniquely) and
#   penalizes ambiguously mapped clusters

library(here)
library(tidyverse)
library(sessioninfo)
library(spatialLIBD)
library(ComplexHeatmap)
library(viridis)

ref_name = "multiome_mid"

in_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    sprintf('cor_vs_%s.rds', ref_name)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'ranking', sprintf('%s.csv', ref_name)
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy'
)
resolution = c(seq_len(20) / 10, 4, 8)
max_ambig_clusters = 5

dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)
this_cor = readRDS(in_path)

#   Gather registration info across all resolutions in a single tibble
anno_df_list = list()
for (i in seq_len(length(this_cor))) {
    anno_df_list[[i]] = this_cor[[i]] |>
        annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
        as_tibble() |>
        mutate(res = resolution[i])
}
anno_df = bind_rows(anno_df_list)

#   Count fraction of habenula and non-habenula cell types covered uniquely by
#   at least one cluster (by resolution)
unique_df = anno_df |>
    filter(layer_confidence == 'good', !grepl('/', layer_label)) |>
    group_by(res) |>
    summarize(
        frac_unique_non_hb = length(
            unique(layer_label[!grepl("[ML]Hb", layer_label)])
        ) / 7,
        frac_unique_hb = length(
            unique(layer_label[grepl("[ML]Hb", layer_label)])
        ) / 9
    )

#   Count fraction of habenula and non-habenula cell types covered in any way by
#   at least one cluster (by resolution)
shared_df = anno_df |>
    filter(layer_confidence == 'good') |>
    separate_longer_delim(layer_label, delim = '/') |>
    group_by(res, cluster) |>
    filter(all(grepl('^MHb', layer_label)) | all(grepl('LHb', layer_label))) |>
    group_by(res) |>
    summarize(frac_shared_hb = length(unique(layer_label)) / 9)

#   Calculate fraction of ambiguous mappings, join with other metrics, and score
metric_df = anno_df |>
    group_by(res) |>
    summarize(
        num_clusters = n(),
        #   Fraction of clusters considered to map ambiguously
        frac_ambig = mean(
            (layer_confidence == 'poor') |
            !(
                #   Either it's split across 1+ MHb clusters
                grepl('^(MHb\\.[1-3]/*)+$', layer_label) |
                #   Or 1+ LHb clusters (had to interactively test this one)
                grepl('^((Inhib_)?LHb(\\.[1-7]|_4\\.[12])/*)+$', layer_label) |
                #   Or 1 cluster of any type
                !grepl('/', layer_label)
            )
        )
    ) |>
    left_join(unique_df, by = "res") |>
    left_join(shared_df, by = "res") |>
    mutate(
        #   Weight all 3 metrics equally, except don't even consider resolutions
        #   with too many ambiguous clusters
        final_score = ifelse(
            round(num_clusters * frac_ambig) > max_ambig_clusters,
            0,
            (frac_unique_non_hb + frac_unique_hb + frac_shared_hb) / 3
        )
    ) |>
    arrange(desc(final_score), frac_ambig)

write_csv(metric_df, out_path)

heatmap_mat = metric_df |>
    select(
        res, frac_unique_non_hb, frac_unique_hb, frac_shared_hb, final_score
    ) |>
    dplyr::rename(
        'Unique Non-Hb' = frac_unique_non_hb,
        'Unique Hb' = frac_unique_hb,
        'Shared Hb' = frac_shared_hb,
        'Final Score' = final_score
    ) |>
    column_to_rownames(var = 'res') |>
    as.matrix()

p = Heatmap(
    heatmap_mat,
    name = 'Score',
    row_title = 'Banksy Resolution',
    column_title = 'Metric',
    col = viridis(100),
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", heatmap_mat[i, j]), x, y,
        gp = gpar(fontsize = 10))
    }
)
pdf(file.path(plot_dir, sprintf('%s_top_results.pdf', ref_name)), width = 5)
draw(p)
dev.off()

session_info()
