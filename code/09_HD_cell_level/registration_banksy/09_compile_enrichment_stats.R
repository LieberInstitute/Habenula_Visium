#   Compile enrichment t-stats for Yalcinbas et al broad + snRNA-seq, multiome
#   and Banksy Visium HD data into a single CSV for Brion to use

library(tidyverse)
library(here)
library(sessioninfo)

stat_paths = c(
    here(
        'processed-data', '09_HD_cell_level', 'new_samples2',
        'registration_banksy', 'modeling_results', 'lambda0_2', '1_7.rds'
    ),
    here(
        "processed-data", "05_snRNA-seq_model_stats",
        sprintf(
            "enrichment_%s.rds",
            c("final_Annotations", "final_Annotations_broad")
        )
    ),
    here(
        'processed-data', '05_snRNA-seq_model_stats',
        'enrichment_snRNA-multiome_v5.rds'
    )
)
names(stat_paths) = c("HD", "Yalcinbas_fine", "Yalcinbas_broad", "multiome")
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'compiled_enrichment_stats.csv.gz'
)

stat_df_list = list()
for (dataset in names(stat_paths)) {
    temp = readRDS(stat_paths[[dataset]])

    if (class(temp) == 'list') {
        stat_df = temp$enrichment
    } else {
        stat_df = temp
    }

    stat_df_list[[dataset]] = stat_df |>
        as_tibble() |>
        pivot_longer(
            cols = matches('^(t_stat|p_value|fdr|logFC)_'),
            names_to = c('.value', 'cluster'),
            names_pattern = '^(t_stat|p_value|fdr|logFC)_(.+)$'
        ) |>
        dplyr::rename(gene_id = ensembl, gene_name = gene) |>
        mutate(dataset = factor(dataset, levels = names(stat_paths)))
}

bind_rows(stat_df_list) |>
    write_csv(out_path)

session_info()
