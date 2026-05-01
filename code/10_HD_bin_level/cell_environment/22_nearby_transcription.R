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
    'ficture_outputs', 'cleaningy', sprintf('k_%d', k), 'analysis',
    sprintf('nF%d.d_12', k), 'cleaningy_joined_input.tsv.gz'
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

ficture_df = read_csv_duckdb(ficture_path, prudence = 'stingy') |>
    distinct(barcode, sample_id, factor_K1) |>
    filter(factor_K1 != 'NA') |>
    dplyr::rename(bin_id = barcode)

#   Mapping from clusters to cell types
anno_df = read_csv(anno_path, show_col_types = FALSE)

cluster_df = read_csv_duckdb(extra_bin_path, prudence = 'lavish') |>
    mutate(sample_id = str_extract(cell_key, '_(H1-.*)$', group = 1)) |>
    inner_join(ficture_df, by = c('bin_id', 'sample_id')) |>
    #   Grab Banksy clusters for each cell
    left_join(
        read_csv_duckdb(cluster_path, prudence = 'stingy') |>
            dplyr::rename(cell_key = key),
        by = 'cell_key'
    ) |>
    filter(banksy != 16, !is.na(banksy)) |>
    #   Annotate cell types
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy), anno_df$cluster)
        ]
    ) |>
    collect()

prop_df = cluster_df |>
    group_by(factor_K1, cell_type) |>
    summarize(prop_bins = n()) |>
    group_by(cell_type) |>
    mutate(prop_bins = prop_bins / sum(prop_bins)) |>
    ungroup()

prop_df = readRDS(cor_path)[[array_task]] |>
    as.data.frame() |>
    rownames_to_column('factor_K1') |>
    pivot_longer(
        cols = -factor_K1, names_to = 'cell_type', values_to = 'cor_val'
    ) |>
    mutate(
        cell_type = case_when(
            cell_type == 'Endo.microglia' ~ 'Endo/microglia',
            cell_type == 'Excit.Thal.Inhib_LHb_4.2' ~ 'Excit.Thal/Inhib_LHb_4.2',
            cell_type == 'LHb.4.Inhib_LHb_4.2' ~ 'LHb.4/Inhib_LHb_4.2',
            TRUE ~ cell_type
        )
    ) |>
    left_join(prop_df, by = c('factor_K1', 'cell_type')) |>
    mutate(cell_type = factor(cell_type, levels = cell_type_levels)) |>
    replace_na(list(prop_bins = 0))

write_csv(prop_df, out_path)

p = ggplot(
        prop_df,
        aes(x = factor_K1, y = cell_type, color = cor_val, size = prop_bins)
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
