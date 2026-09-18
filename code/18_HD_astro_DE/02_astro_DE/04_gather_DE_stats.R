library(tidyverse)
library(here)
library(sessioninfo)
library(duckplyr)

de_paths = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'DE_%d.parquet'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'aggregated_DE_stats.csv.gz'
)
num_permutations = 1000
min_permutations = 500

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

de_null_df_list = list()
for (i in seq_len(num_permutations)) {
    de_null_df_list[[i]] = read_parquet_duckdb(
        sprintf(de_paths, i), prudence = 'lavish'
    )
}
de_null_df = bind_rows(de_null_df_list) |>
    dplyr::rename(t_null = t) |>
    select(gene_id, t_null) |>
    collect()

de_true_df = read_parquet_duckdb(
        sprintf(de_paths, 0), prudence = 'stingy'
    ) |>
    dplyr::rename(fdr = adj.P.Val, p = P.Value) |>
    collect()

de_true_df |>
    select(gene_id, t) |>
    dplyr::rename(t_true = t) |>
    left_join(de_null_df, by = 'gene_id') |>
    group_by(gene_id) |>
    summarize(
        num_permutations = n(),
        #   How often does a t-stat sampled from the null distribution exceed
        #   the magnitude of the t-stat from the real-label test? Here the
        #   definition of "exceeed" depends on the sign of the observed t-stat
        p_empirical = ifelse(
            t_true[1] > 0,
            max(1 / num_permutations, mean(t_null >= t_true)),
            max(1 / num_permutations, mean(t_null <= t_true))
        )
    ) |>
    ungroup() |>
    #   Require the gene to pass expression filtering in a minimum number of
    #   permutations to achieve a reasonable p-value resolution
    filter(num_permutations >= min_permutations) |>
    #   FDR correction across all genes
    mutate(fdr_empirical = p.adjust(p_empirical, method = 'fdr')) |>
    select(gene_id, p_empirical, fdr_empirical) |>
    #   Join back with info like logFC and reorganize columns
    left_join(de_true_df, by = 'gene_id') |>
    select(gene_id, gene_name, t, logFC, p, fdr, p_empirical, fdr_empirical) |>
    write_csv(out_path)

session_info()
