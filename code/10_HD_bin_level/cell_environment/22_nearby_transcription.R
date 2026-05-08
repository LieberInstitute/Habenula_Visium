#   How do FICTURE clusters in the extracellular environment around cells agree
#   or disagree transcriptionally with those cells?

library(tidyverse)
library(here)
library(duckplyr)
library(sessioninfo)

array_task = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
k = c(3:10, 20)[array_task]

cor_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'registration', 'cor_vs_cell_types.rds'
)
extra_bin_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_bins.csv.gz'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'ficture_plotting', 'extracellular.parquet'
)
plot_path = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'nearby_transcription', sprintf('nearby_transcription_k%d.pdf', k)
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'nearby_transcription', sprintf('k%d.csv', k)
)
cell_type_levels = c(
    'MHb.1', 'MHb.2', 'Excit_LHb', 'LHb.2.7', 'LHb.4', 'LHb.4/Inhib_LHb_4.2',
    'Inhib_LHb_4.2', 'Excit.Thal/Inhib_LHb_4.2', 'Excit.Thal', 'Astrocyte',
    'Endo', 'Endo/microglia', 'Oligo', 'OPC', 'Ependymal', 'Subependymal'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(plot_path), showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

anno_df = read_csv_duckdb(anno_path, prudence = 'stingy') |>
    dplyr::rename(banksy = cluster, cell_type = fine_cell_type)

cluster_df = read_csv_duckdb(cluster_path, prudence = 'stingy') |>
    dplyr::rename(cell_key = key) |>
    left_join(anno_df, by = 'banksy') |>
    filter(cell_type != 'Drop') |>
    select(cell_key, cell_type)

ficture_df = read_parquet_duckdb(ficture_path, prudence = 'stingy') |>
    dplyr::rename(ficture_cluster = paste0('k', k)) |>
    filter(!is.na(ficture_cluster)) |>
    select(bin_key, ficture_cluster)

full_df = read_csv_duckdb(extra_bin_path, prudence = 'lavish') |>
    mutate(
        sample_id = str_extract(cell_key, '_(H1-.*)$', group = 1),
        bin_key = paste(bin_id, sample_id, sep = '_')
    ) |>
    select(cell_key, bin_key) |>
    inner_join(ficture_df, by = 'bin_key') |>
    inner_join(cluster_df, by = 'cell_key') |>
    collect()

prop_df = full_df |>
    group_by(ficture_cluster, cell_type) |>
    summarize(prop_bins = n()) |>
    group_by(cell_type) |>
    mutate(prop_bins = prop_bins / sum(prop_bins)) |>
    ungroup()

prop_df = readRDS(cor_path)[[array_task]] |>
    as.data.frame() |>
    rownames_to_column('ficture_cluster') |>
    pivot_longer(
        cols = -ficture_cluster, names_to = 'cell_type', values_to = 'cor_val'
    ) |>
    mutate(
        cell_type = case_when(
            cell_type == 'Endo.microglia' ~ 'Endo/microglia',
            cell_type == 'Excit.Thal.Inhib_LHb_4.2' ~ 'Excit.Thal/Inhib_LHb_4.2',
            cell_type == 'LHb.4.Inhib_LHb_4.2' ~ 'LHb.4/Inhib_LHb_4.2',
            TRUE ~ cell_type
        )
    ) |>
    left_join(
        prop_df |> mutate(ficture_cluster = as.character(ficture_cluster)),
        by = c('ficture_cluster', 'cell_type')
    ) |>
    mutate(cell_type = factor(cell_type, levels = cell_type_levels)) |>
    replace_na(list(prop_bins = 0))

write_csv(prop_df, out_path)

p = ggplot(
        prop_df,
        aes(x = ficture_cluster, y = cell_type, color = cor_val, size = prop_bins)
    ) +
    geom_point() +
    coord_fixed() +
    scale_color_gradient2() +
    scale_radius(range = c(1, 8)) +
    theme_bw(base_size = 16) +
    labs(
        x = 'FICTURE Cluster', y = 'Cell Type',
        color = 'Transcriptional\nCorrelation',
        size = 'Proportion of\nExtracellular\nBins'
    )
pdf(plot_path, width = k / 2 + 4)
print(p)
dev.off()

session_info()
