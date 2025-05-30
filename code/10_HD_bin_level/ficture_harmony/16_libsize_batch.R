# Combine library size normalized and batch effect corrected results
library(here)
library(tidyverse)
library(SpatialExperiment)
library(HDF5Array)
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


cleany_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration',"sum_score","cleany", sprintf('heatmap_score_%s.csv', ref_name))
cleany = read.csv(cleany_path)

normalized_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration',"sum_score","normalized", sprintf('heatmap_score_%s.csv', ref_name))

normalized = read.csv(normalized_path)

merged_df <- merge(cleany, normalized, 
                   by = "original_cluster_count", 
                   suffixes = c("_cleany", "_normalized"))

merged_sorted <- merged_df[order(-merged_df$score_total_cleany), ]
merged_sorted$cleany_rank<-order(-merged_sorted$score_total_scaled_cleany)
merged_sorted$normalized_rank<-order(-merged_sorted$score_total_scaled_normalized)

cor=cor(merged_sorted$cleany_rank,merged_sorted$normalized_rank)

com_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'registration',"sum_score","cleany_normalized", sprintf('heatmap_score_%s_cor%s.csv', ref_name, round(cor,2)))


write.csv(as.data.frame(merged_sorted),com_path, row.names = FALSE)