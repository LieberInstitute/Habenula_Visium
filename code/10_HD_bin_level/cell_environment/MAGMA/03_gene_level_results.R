library(tidyverse)
library(here)
library(rtracklayer)
library(sessioninfo)

k = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

gwas_name_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/gwas_info.csv'
gwas_map = read_csv(gwas_name_path, show_col_types = FALSE) |>
    select(nickname, manuscript_name)
these_gwas_names = ifelse(
    gwas_map$nickname == 'MDD2019', 'MDD', gwas_map$nickname
)

gene_set_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', 'gene_sets', sprintf('k%d.tsv', k)
)
gene_stat_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', these_gwas_names, sprintf('%s.genes.out', these_gwas_names)
)
set_stat_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', '%s', '%s', sprintf('k%d.gsa.out', k)
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', sprintf('top_genes_k%d.csv', k)
)
out_low_genes_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'MAGMA', sprintf('low_gene_sets_k%d.csv', k)
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
names(gene_stat_paths) = these_gwas_names
sig_cutoff = 0.05
min_genes_per_set = 10

#   MAGMA set-level outputs have a variable amount of header lines. Auto-detect
#   the header length and read in dynamically
read_table_auto_skip = function(path, check_lines = 100) {
    n_skip = sum(grepl('^#', readLines(path, n = check_lines)))
    clean_df = read_table(path, skip = n_skip, show_col_types = FALSE)
    return(clean_df)
}

gene_df_list = list()
for (gwas in names(gene_stat_paths)) {
    gene_stat_df = read_table(
            gene_stat_paths[[gwas]], show_col_types = FALSE
        ) |>
        mutate(gwas = gwas) |>
        select(GENE, P, gwas)

    #   Read in set-level stats to identify significant sets
    set_df = sprintf(set_stat_paths, gwas, gwas) |>
        read_table_auto_skip() |>
        dplyr::rename(cell_type = VARIABLE) |>
        mutate(set_is_sig = P < sig_cutoff) |>
        select(cell_type, set_is_sig)

    #   Read in gene sets themselves and merge with gene- and set-level
    #   stats
    gene_df_list[[gwas]] = read_table(gene_set_path, show_col_types = FALSE) |>
        left_join(gene_stat_df, by = c('gene_id' = 'GENE')) |>
        dplyr::rename(cell_type = set_id) |>
        left_join(set_df, by = 'cell_type') |>
        dplyr::rename(p = P) |>
        select(gene_id, cell_type, gwas, p, set_is_sig)
}

gene_df = bind_rows(gene_df_list)

#   Warn about percentage of genes missing MAGMA stats
for (gwas in names(gene_stat_paths)) {
    message(
        sprintf(
            'Dropping %d%% of genes for %s GWAS missing MAGMA stats',
            round(100 * mean(is.na(gene_df$p[gene_df$gwas == gwas]))),
            gwas
        )
    )
}

gene_df = gene_df |>
    filter(!is.na(p)) |>
    mutate(gwas = ifelse(gwas == 'MDD', 'MDD2019', gwas)) |>
    left_join(gwas_map, by = c('gwas' = 'nickname'))

#   Do we have enough genes for meaningful testing? Export sets with too few
#   genes, so we can note them in the heatmaps
gene_df |>
    group_by(cell_type, gwas) |>
    filter(n() < min_genes_per_set) |>
    ungroup() |>
    distinct(cell_type, manuscript_name) |>
    write_csv(out_low_genes_path)

#   Read in the GTF to get gene symbols
gtf = import(reference_gtf)
gtf = gtf[gtf$type == 'gene'] |>
    as.data.frame() |>
    as_tibble() |>
    select(gene_id = gene_id, gene_name = gene_name)

#   Export final gene sets, only including genes where the set
#   as a whole was significant
gene_df |>
    filter(p < sig_cutoff, set_is_sig) |>
    #   Require sets to have a minimum number of genes (to accurately determine
    #   set-level significance)
    group_by(cell_type, gwas) |>
    filter(n() >= min_genes_per_set) |>
    ungroup() |>
    arrange(gwas, cell_type, p) |>
    left_join(gtf, by = 'gene_id') |>
    select(manuscript_name, cell_type, gene_id, gene_name, p) |>
    dplyr::rename(gwas = manuscript_name) |>
    write_csv(out_path)

session_info()
