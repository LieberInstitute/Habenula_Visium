#   This script provides a data-driven method for ranking spatial registration
#   results across all tested k values, ultimately deciding upon an optimal
#   clustering k value from FICTURE. It ultimately uses the Yalcinbas
#   snRNA-seq reference data, favors covering many cell types (especially
#   uniquely) and penalizes ambiguously mapped clusters

library(here)
library(tidyverse)
library(sessioninfo)
library(spatialLIBD)
library(ComplexHeatmap)
library(viridis)

input_method = c("normalized", "cleaning_y")[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]
ref_name = "multiome"

in_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'registration', 'cor_rds', input_method,
    sprintf('cor_vs_%s.rds', ref_name)
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'deciding_k', sprintf('%s_%s.csv', ref_name, input_method)
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'new_samples2', 'ficture_harmony', 'deciding_k'
)
k_values = c(seq(3, 40), 70, 100)
max_ambig_clusters = 5

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

this_cor = readRDS(in_path)

#   Gather registration info across all k values in a single tibble
anno_df_list = list()
for (i in seq_len(length(this_cor))) {
    anno_df_list[[i]] = this_cor[[i]] |>
        annotate_registered_clusters(cutoff_merge_ratio = 0.1) |>
        as_tibble() |>
        mutate(k = k_values[i])
}
anno_df = bind_rows(anno_df_list)

#   Count fraction of habenula and non-habenula cell types covered uniquely by
#   at least one cluster (by k)
unique_df = anno_df |>
    filter(layer_confidence == 'good', !grepl('/', layer_label)) |>
    group_by(k) |>
    summarize(
        frac_unique_non_hb = length(
            unique(layer_label[!grepl("[ML]Hb", layer_label)])
        ) / 7,
        frac_unique_hb = length(
            unique(layer_label[grepl("[ML]Hb", layer_label)])
        ) / 9
    )

#   Count fraction of habenula and non-habenula cell types covered in any way by
#   at least one cluster (by k)
shared_df = anno_df |>
    filter(layer_confidence == 'good') |>
    separate_longer_delim(layer_label, delim = '/') |>
    group_by(k, cluster) |>
    filter(all(grepl('^MHb', layer_label)) | all(grepl('LHb', layer_label))) |>
    group_by(k) |>
    summarize(frac_shared_hb = length(unique(layer_label)) / 9)

#   Calculate fraction of ambiguous mappings, join with other metrics, and score
metric_df = anno_df |>
    group_by(k) |>
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
    left_join(unique_df, by = "k") |>
    left_join(shared_df, by = "k") |>
    mutate(
        #   Weight all 3 metrics equally, except don't even consider k values
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
        k, frac_unique_non_hb, frac_unique_hb, frac_shared_hb, final_score
    ) |>
    dplyr::rename(
        'Unique Non-Hb' = frac_unique_non_hb,
        'Unique Hb' = frac_unique_hb,
        'Shared Hb' = frac_shared_hb,
        'Final Score' = final_score
    ) |>
    column_to_rownames(var = 'k') |>
    as.matrix()

p = Heatmap(
    heatmap_mat,
    name = 'Score',
    row_title = 'FICTURE k',
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
pdf(file.path(plot_dir, sprintf('%s__%s_top_results.pdf', ref_name, input_method)), width = 5)
draw(p)
dev.off()

session_info()
