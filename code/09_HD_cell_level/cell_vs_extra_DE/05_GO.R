library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

cell_map_path = here("raw-data", "sample_info", "hd_cell_type_map.csv")
de_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'main_results'
)
plot_dir = here(
    "plots", "09_HD_cell_level", "no_secondary", "cell_vs_extra_DE", "GO"
)
go_num_terms = 2

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

plot_go = function(plot_df, cell_map_df, plot_path) {
    #   Order GO terms by the first (highest-ranked) cell type they appear in
    term_order = plot_df |>
        mutate(cell_type = factor(cell_type, levels = cell_map_df$new_cell_type)) |>
        group_by(Description) |>
        slice_min(cell_type, n = 1, with_ties = FALSE) |>
        ungroup() |>
        arrange(cell_type) |>
        pull(Description)

    p = plot_df |>
        mutate(
            cell_type   = factor(cell_type, levels = cell_map_df$new_cell_type),
            Description = factor(Description, levels = term_order)
        ) |>
        ggplot(
            aes(x = cell_type, y = Description, color = log_fdr, size = gene_ratio)
        ) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(
            x = "Cell Type", y = "GO Term", color = "-log10(FDR)",
            size = "Gene Ratio"
        )

    pdf(plot_path, width = 3 + length(cell_types) / 3, height = 5)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

cell_map_df = read_csv(cell_map_path, show_col_types = FALSE) |>
    dplyr::rename(cell_type = old_cell_type) |>
    dplyr::select(cell_type, new_cell_type)

de_df = list.files(de_dir, full.names = TRUE) |>
    map_dfr(read_parquet_duckdb, prudence = 'lavish') |>
    filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
    dplyr::select(gene_id, cell_type, logFC) |>
    collect() |>
    #   Just drop cell types where we don't have enough DEGs
    group_by(cell_type, sign(logFC)) |>
    filter(n() > 10) |>
    ungroup() |>
    dplyr::rename(de_sign = `sign(logFC)`) |>
    left_join(cell_map_df, by = "cell_type") |>
    dplyr::select(gene_id, new_cell_type, de_sign) |>
    dplyr::rename(cell_type = new_cell_type)

#   Background universe is any gene that passed expression filters for the given
#   cell type just before running DE
universe_df = list.files(de_dir, full.names = TRUE) |>
    map_dfr(read_parquet_duckdb, prudence = 'lavish') |>
    dplyr::select(gene_id, cell_type) |>
    collect() |>
    left_join(cell_map_df, by = "cell_type") |>
    dplyr::select(gene_id, new_cell_type) |>
    dplyr::rename(cell_type = new_cell_type)

#   Run enrichGO per cell type for up- and down-regulated sets
ego_df_list = list()
for (this_cell_type in cell_map_df$new_cell_type) {
    universe_genes = universe_df |>
        filter(cell_type == this_cell_type) |>
        pull(gene_id)

    for (this_de_sign in c(-1, 1)) {
        de_sign_name = ifelse(this_de_sign == 1, "up", "down")

        gene_set = de_df |>
            filter(cell_type == this_cell_type, de_sign == this_de_sign) |>
            pull(gene_id)

        if(length(gene_set) == 0) {
            message(
                sprintf(
                    "Not enough %sregulated genes for cell type %s; dropping this combo",
                    de_sign_name, this_cell_type
                )
            )
            next
        }
  
        ego = enrichGO(
            gene          = gene_set,
            OrgDb         = org.Hs.eg.db,
            keyType       = "ENSEMBL",
            ont           = "BP",
            universe      = universe_genes,
            pAdjustMethod = "BH",
            pvalueCutoff  = 1,
            qvalueCutoff  = 1
        )
            
        ego_df_list[[this_cell_type]] = ego@result |>
            as_tibble() |>
            filter(p.adjust < 0.05) |>
            mutate(cell_type = this_cell_type, de_direction = de_sign_name)
    }
}

plot_df = bind_rows(ego_df_list) |>
    mutate(
        gene_ratio = Count / as.integer(str_extract(GeneRatio, "(?<=/)[0-9]+")),
        log_fdr = -log10(p.adjust)
    )

#   Custom dot plot by cell type for each DE direction
for (this_direction in c("up", "down")) {
    plot_go(
        plot_df |>
            filter(de_direction == this_direction),
        cell_map_df,
        file.path(plot_dir, sprintf("GO_%s.pdf", this_direction))
    )
}

session_info()
