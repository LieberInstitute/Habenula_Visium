library(tidyverse)
library(here)
library(spatialLIBD)
library(duckdb)
library(ggrepel)
library(sessioninfo)

ficture_k = 10
spe_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'spe_norm_filtered.rds'
)
extra_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'cell_environment',
    'extracellular_bins.csv.gz'
)
ficture_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'bin_level_clusters_batch.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'banksy', 'lambda0_2',
    'leiden_res1_7.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'registration_banksy',
    'cluster_annotation.csv'
)
model_path = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'ficture_harmony',
    'registration', 'modeling_results', 'cleaning_y',
    sprintf('%d.rds', ficture_k)
)
ficture_colnames = c('sample_id', 'barcode', sprintf('FICTURE_k%d', ficture_k))
plot_dir = here('plots', '10_HD_bin_level', 'new_samples2', 'cell_environment')
cell_type_colors = c(
    MHb.2 = '#942911',
    LHb.2.7 = '#0987AD',
    Other = '#848484'
)

dir.create(plot_dir, showWarnings = FALSE)

spe = readRDS(spe_path)

#   Merge in annotated Banksy results
anno_df = read_csv(ct_anno_path, show_col_types = FALSE)
spe$cell_type = tibble(key = spe$key) |>
    left_join(read_csv(banksy_path, show_col_types = FALSE), by = 'key') |>
    mutate(
        cell_type = anno_df$fine_cell_type[
            match(as.character(banksy_lambda0_2), anno_df$cluster)
        ]
    ) |>
    pull(cell_type)

con = dbConnect(duckdb())

col_data = tibble(
    cell_id = spe$key, cell_type = spe$cell_type, bins_per_cell = spe$bin_count
)
duckdb_register(con, "col_data", col_data)

#   Form a 2um bin-level tibble of extracellular bins with info about the
#   associated cell type, FICTURE cluster, and cellular size
sql_query = sprintf(
    "
    SELECT 
        extra.sample_id,
        extra.bin_id AS barcode,
        extra.cell_id || '_' || extra.sample_id AS cell_id,
        col_data.cell_type,
        col_data.bins_per_cell,
        ficture.FICTURE_k%d as ficture_cluster
    FROM read_csv_auto('%s') AS extra
    INNER JOIN col_data 
        ON extra.cell_id || '_' || extra.sample_id = col_data.cell_id
    LEFT JOIN read_csv_auto('%s') AS ficture
        ON extra.sample_id = ficture.sample_id AND extra.bin_id = ficture.barcode
    ",
    ficture_k, extra_path, ficture_path
)
bin_df = dbGetQuery(con, sql_query) |>
    as_tibble()
duckdb_register(con, "bin_df", bin_df)

#   By cell type, summarize average info about number of cellular and
#   extracellular bins
cell_df = dbGetQuery(
        con,
        "
        SELECT
            cell_id,
            cell_type,
            bins_per_cell AS num_cellular_bins,
            COUNT(*) AS num_extra_bins
        FROM bin_df
        GROUP BY cell_id, cell_type, bins_per_cell
        "
    ) |>
    as_tibble() |>
    group_by(cell_type) |>
    summarize(
        mean_cellular_bins = mean(num_cellular_bins),
        mean_extra_bins = mean(num_extra_bins)
    )

#   Examine how different cell types compare in terms of number of cellular and
#   extracellular bins. Done as a scatterplot with a regression line as a means
#   of indirectly showing which cell types are in relatively more-dense areas
#   (i.e. deviating downward from the regression line)
p = cell_df |>
    filter(cell_type != 'Ambig') |>
    ggplot(aes(x = mean_cellular_bins, y = mean_extra_bins, label = cell_type)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "blue") +
        geom_text_repel(size = 5) +
        labs(x = "Mean Cellular Bins", y = "Mean Extracellular Bins") +
        theme_bw(base_size = 20)
pdf(file.path(plot_dir, "cell_sizes_scatter.pdf"))
print(p)
dev.off()

#   For select cell types, examine proportions of FICTURE clusters present among
#   their extracellular bins
p = bin_df |>
    filter(ficture_cluster != "NA") |>
    mutate(
        cell_type = ifelse(
            cell_type %in% c('MHb.2', 'LHb.2.7'), cell_type, 'Other'
        )
    ) |>
    group_by(cell_type, ficture_cluster) |>
    summarize(num_bins = n()) |>
    group_by(cell_type) |>
    mutate(prop_bins = num_bins / sum(num_bins)) |>
    ungroup() |>
    mutate(
        ficture_cluster = factor(
            ficture_cluster,
            levels = as.character(sort(unique(as.integer(ficture_cluster))))
        )
    ) |>
    ggplot(aes(x = ficture_cluster, y = prop_bins, fill = cell_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = cell_type_colors) +
        labs(
            x = "FICTURE Cluster",
            y = "Proportion of Extracellular Bins",
            fill = "Cell Type"
        ) +
        theme_bw(base_size = 20)
pdf(file.path(plot_dir, "extracellular_profile_k10.pdf"), width = 8, height = 6)
print(p)
dev.off()

#   Mhb.2 is distinguished by FICTURE cluster 0, and LHb.2.7 by FICTURE cluster
#   3. Check top enriched genes for these two clusters (by significance)
message("Top enriched genes for FICTURE clusters 0 and 3:")
readRDS(model_path)$enrichment |>
    as_tibble() |>
    dplyr::rename(gene_id = ensembl, gene_name = gene) |>
    select(gene_id, gene_name, fdr_X0, fdr_X3) |>
    pivot_longer(
        cols = matches("^fdr_X[03]$"),
        names_to = "ficture_cluster",
        values_to = "fdr"
    ) |>
    mutate(ficture_cluster = str_remove(ficture_cluster, "^fdr_X")) |>
    group_by(ficture_cluster) |>
    arrange(fdr) |>
    slice_head(n = 5) |>
    print(n = 10)

session_info()
