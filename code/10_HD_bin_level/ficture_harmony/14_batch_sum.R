library(here)
library(tidyverse)
library(SpatialExperiment)
library(sessioninfo)
library(spatialLIBD)
library(data.table)

ref_names = c(
    'snRNAseq_fine', 'snRNAseq_broad', 'multiome',
    sprintf('Visium_BayesSpace_k%02d', 2:28)
)

#   Get the reference data for this task
task_id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
ref_name = ref_names[task_id]

out_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'registration',"cor_rds","cleaning_y",sprintf('cor_vs_%s.rds', ref_name)
)
this_cor = readRDS(out_path)

#   Annotate clusters
annotated_clusters = lapply(
    this_cor, annotate_registered_clusters, cutoff_merge_ratio = 0.1
)

# N of cell type covered
filtered_clusters <- lapply(annotated_clusters, function(df) {
  # filter layer_confidence is good and multiple clusters matching one cell type
  df_good <- df[df$layer_confidence == "good", ]
  df_good <- df_good[!grepl("/", df_good$layer_label), ]  
  return(df_good)
})

original_cluster_counts <- sapply(annotated_clusters, function(df) {
  length(unique(df$cluster))
})

filtered_cell_type_counts <- sapply(filtered_clusters, function(df) {
  length(unique(df$layer_label))
})

if (task_id %in% 1:3){

# N of habenula cell type covered

mhb_lhb_counts <- sapply(filtered_clusters, function(df) {
  length(unique(df$layer_label[grepl("MHb|LHb", df$layer_label)]))
})

# N of one-to-one

filtered_clusters_one<- lapply(annotated_clusters, function(df) {
  # filter layer_confidence is good
  df <- df %>% separate_rows(layer_label, sep = "/", convert = FALSE)
  df_good <- df[df$layer_confidence == "good", ]

  #find duplicated layer_label and cluster
  dup_labels <- df_good$layer_label[duplicated(df_good$layer_label) | duplicated(df_good$layer_label, fromLast = TRUE)]
  dup_clusters <- df_good$cluster[duplicated(df_good$cluster)|duplicated(df_good$cluster, fromLast = TRUE)]

  #remove duplicated layer_label or cluster
  df_unique <- df_good[!(df_good$layer_label %in% dup_labels) & !(df_good$cluster %in% dup_clusters), ]
  
  return(df_unique)
})

filtered_cell_type_counts_one <- sapply(filtered_clusters_one, function(df) {
  length(unique(df$layer_label))
})

# N of mhb_lhb_counts (clusters mapped to all mhb/lhb are also included)

filtered_clusters_pure<- lapply(annotated_clusters, function(df) {
  # filter layer_confidence is good
  df <- df %>% separate_rows(layer_label, sep = "/", convert = FALSE)
  df_good <- df[df$layer_confidence == "good", ]

  # for each cluster, only keep all the labels are MHb or LHb specifically
  df_pure <- df_good %>%
    group_by(cluster) %>%
    filter(
      all(grepl("MHb", layer_label)) |
      all(grepl("LHb", layer_label))
    ) %>%
    ungroup()
  
  return(df_pure)
})

filtered_cell_type_pure <- sapply(filtered_clusters_pure, function(df) {
  length(unique(df$layer_label))
})

# summary
cell_type_summary <- data.frame(
  original_cluster_count = original_cluster_counts,
  unique_cell_type_count = filtered_cell_type_counts,
  unique_mhb_lhb_counts = mhb_lhb_counts,
  one_to_one_cell_type_count = filtered_cell_type_counts_one,
  multiple_mhb_lhb_counts = filtered_cell_type_pure
)

minmax_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

score_unique <- minmax_scale(cell_type_summary$unique_cell_type_count)
score_mhb_lhb <- minmax_scale(cell_type_summary$unique_mhb_lhb_counts)
score_one_to_one <- minmax_scale(cell_type_summary$one_to_one_cell_type_count)
score_multiple_mhb_lhb <- minmax_scale(cell_type_summary$multiple_mhb_lhb_counts)

cell_type_summary$score_total <- score_unique + score_mhb_lhb + score_one_to_one + score_multiple_mhb_lhb
cell_type_summary$score_total_scaled <- minmax_scale(cell_type_summary$score_total)
cell_type_summary <- cell_type_summary[order(-cell_type_summary$score_total), ]

write_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'registration',"sum_score","cleany", sprintf('heatmap_score_%s.csv', ref_name))

write.csv(as.data.frame(cell_type_summary),write_path, row.names = FALSE)
}

if (task_id > 3){
  
  filtered_clusters_one<- lapply(annotated_clusters, function(df) {
  # filter layer_confidence is good
  df_good <- df[df$layer_confidence == "good", ]
  df_good <- df_good[!grepl("/", df_good$layer_label), ]

  #find duplicated layer_label and cluster
  dup_labels <- df_good$layer_label[duplicated(df_good$layer_label) | duplicated(df_good$layer_label, fromLast = TRUE)]

  #remove duplicated layer_label or cluster
  df_unique <- df_good[!(df_good$layer_label %in% dup_labels), ]
  
  return(df_unique)
})

filtered_cell_type_counts_one <- sapply(filtered_clusters_one, function(df) {
  length(unique(df$layer_label))
})

# summary
cell_type_summary <- data.frame(
  original_cluster_count = original_cluster_counts,
  unique_cell_type_count = filtered_cell_type_counts,
  one_to_one_cell_type_count = filtered_cell_type_counts_one)

minmax_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

score_unique <- cell_type_summary$unique_cell_type_count/max(cell_type_summary$unique_cell_type_count)
score_one_to_one <- cell_type_summary$one_to_one_cell_type_count/max(cell_type_summary$unique_cell_type_count)

cell_type_summary$score_total <- score_unique + score_one_to_one
cell_type_summary$score_total_scaled <- minmax_scale(cell_type_summary$score_total)
cell_type_summary <- cell_type_summary[order(-cell_type_summary$score_total), ]

write_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'registration',"sum_score","cleany", sprintf('heatmap_score_%s.csv', ref_name))

write.csv(as.data.frame(cell_type_summary),write_path, row.names = FALSE)
}