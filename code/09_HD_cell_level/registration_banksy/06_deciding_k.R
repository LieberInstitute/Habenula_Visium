#   Since there are a huge number of spatial registration results to comb
#   through, this script intends to automate selecting interesting results.
#   Use Chunyu's code to rank Banksy results by 4 metrics, largely based on
#   "clean" matches against the snRNA-seq data

library(here)
library(tidyverse)
library(sessioninfo)
library(spatialLIBD)

ref_names = c('snRNAseq_fine', 'snRNAseq_broad')

#   Get the reference data for this task
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ref_name = ref_names[task_id]

in_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'lambda0_2', sprintf('cor_vs_%s_subset.rds', ref_name)
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'registration_banksy',
    'lambda0_2', 'ranking', sprintf('%s_subset.csv', ref_name)
)
resolution = c(seq_len(20) / 10, 4, 8)

dir.create(dirname(out_path), showWarnings = FALSE)

minmax_scale = function(x) (x - min(x)) / (max(x) - min(x))

this_cor = readRDS(in_path)

#   Annotate clusters
annotated_clusters = lapply(
    this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
)

# N of cell type covered
filtered_clusters = lapply(
    annotated_clusters,
    function(df) {
        # filter layer_confidence is good and multiple clusters matching one
        # cell type
        df_good <- df[df$layer_confidence == "good", ]
        df_good <- df_good[!grepl("/", df_good$layer_label), ]  
        return(df_good)
    }
)
filtered_cell_type_counts = sapply(
    filtered_clusters, function(df) length(unique(df$layer_label))
)

original_cluster_counts = sapply(
    annotated_clusters, function(df) length(unique(df$cluster))
)

# N of habenula cell type covered
mhb_lhb_counts = sapply(
    filtered_clusters,
    function(df) {
        length(unique(df$layer_label[grepl("MHb|LHb", df$layer_label)]))
    }
)

# N of one-to-one
filtered_clusters_one = lapply(
    annotated_clusters,
    function(df) {
        # filter layer_confidence is good
        df = df |> separate_rows(layer_label, sep = "/", convert = FALSE)
        df_good = df[df$layer_confidence == "good", ]

        # find duplicated layer_label and cluster
        dup_labels = df_good$layer_label[
            duplicated(df_good$layer_label) |
            duplicated(df_good$layer_label, fromLast = TRUE)
        ]
        dup_clusters = df_good$cluster[
            duplicated(df_good$cluster) |
            duplicated(df_good$cluster, fromLast = TRUE)
        ]

        # remove duplicated layer_label or cluster
        df_unique = df_good[
            !(df_good$layer_label %in% dup_labels) & 
            !(df_good$cluster %in% dup_clusters),
        ]

        return(df_unique)
    }
)

filtered_cell_type_counts_one = sapply(
    filtered_clusters_one, function(df) length(unique(df$layer_label))
)

# N of mhb_lhb_counts (clusters mapped to all mhb/lhb are also included)
filtered_clusters_pure = lapply(
    annotated_clusters,
    function(df) {
        # filter layer_confidence is good
        df = df |> separate_rows(layer_label, sep = "/", convert = FALSE)
        df_good = df[df$layer_confidence == "good", ]

        # for each cluster, only keep all the labels are MHb or LHb
        # specifically
        df_pure = df_good |>
            group_by(cluster) |>
            filter(
                all(grepl("MHb", layer_label)) |
                all(grepl("LHb", layer_label))
            ) |>
            ungroup()

        return(df_pure)
    }
)

filtered_cell_type_pure = sapply(
    filtered_clusters_pure,
    function(df) length(unique(df$layer_label))
)

# summary
cell_type_summary = tibble(
    resolution = resolution,
    original_cluster_count = original_cluster_counts,
    unique_cell_type_count = filtered_cell_type_counts,
    unique_mhb_lhb_counts = mhb_lhb_counts,
    one_to_one_cell_type_count = filtered_cell_type_counts_one,
    multiple_mhb_lhb_counts = filtered_cell_type_pure
)

score_unique <- minmax_scale(cell_type_summary$unique_cell_type_count)
score_mhb_lhb <- minmax_scale(cell_type_summary$unique_mhb_lhb_counts)
score_one_to_one <- minmax_scale(cell_type_summary$one_to_one_cell_type_count)
score_multiple_mhb_lhb <- minmax_scale(
    cell_type_summary$multiple_mhb_lhb_counts
)

cell_type_summary$score_total <- score_unique + score_mhb_lhb +
    score_one_to_one + score_multiple_mhb_lhb
cell_type_summary$score_total_scaled <- minmax_scale(
    cell_type_summary$score_total
)
cell_type_summary <- cell_type_summary[order(-cell_type_summary$score_total), ]

write_csv(cell_type_summary, out_path)

session_info()
