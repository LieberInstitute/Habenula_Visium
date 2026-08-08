library(tidyverse)
library(here)
library(sessioninfo)
library(duckplyr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

de_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'main_results'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'check_results'
)
cell_marker_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA', 'gene_sets',
    'fine.tsv'
)
extra_marker_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'gene_sets', 'fine.tsv'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

################################################################################
#   Functions
################################################################################

custom_volcano <- function(
        data, FDR_cut = 0.05, FC_cut = 1, p_col = "P.Value",
        fdr_col = "adj.P.Val", lfc_col = "logFC", text = TRUE,
        highlight_genes = NULL
    ){
  
    # define colors
    signif_colors <- c("purple", "blue", "red")
    FDR_label = paste0("FDR<", FDR_cut)
    FC_label = paste0("abs(logFC)>", FC_cut)
    names(signif_colors) <- c("both", FDR_label, FC_label)
    
    volcano <- data |>
        mutate(
            DE_class = case_when(
                    !!sym(fdr_col) < FDR_cut & abs(!!sym(lfc_col)) > FC_cut ~ "both",
                    !!sym(fdr_col) < FDR_cut ~ FDR_label,
                    abs(!!sym(lfc_col)) > FC_cut ~ FC_label,
                    TRUE ~ "None"
                ) |>
                factor(levels = c(FDR_label, FC_label, "both", "None"))
        ) |>
        ggplot(
                aes(
                    x = !!sym(lfc_col), y = -log10(!!sym(p_col)),
                    color = DE_class
                )
            ) +
            geom_point(alpha = 0.5, size = 0.3) +
            scale_color_manual(values = signif_colors) +
            facet_wrap(~cell_type, scales = "free") +
            labs(x = "log(FC)", y = "-log10(P value)") +
            theme_bw(base_size = 10) +
            theme(legend.position = "right")
    
    if(text) {
        volcano <- volcano + 
            geom_text_repel(
                aes(label = ifelse(!!sym(fdr_col) < FDR_cut, gene_name, "")),
                size = 1.5
            ) 
    }
    
    if(!is.null(highlight_genes)){
        volcano <- volcano + 
            geom_text_repel(
                aes(
                    label = ifelse(
                        gene_name %in% highlight_genes, gene_name, ""
                    )
                ),
                size = 1.5
            )
    }
    
    return(volcano)
}

################################################################################
#   Main
################################################################################

dir.create(plot_dir, showWarnings = FALSE)

de_df = list.files(de_dir, full.names = TRUE) |>
    map_dfr(read_parquet_duckdb, prudence = 'lavish') |>
    collect()

p = custom_volcano(de_df, text = FALSE)
pdf(file.path(plot_dir, 'volcano_plots.pdf'))
print(p)
dev.off()

# Jaccard similarity between significant gene sets across cell types
sig_genes = de_df |>
    filter(adj.P.Val < 0.05) |>
    as.data.frame() |>
    split(~ cell_type) |>
    lapply(\(x) x$gene_id)

cell_types = names(sig_genes)
n = length(cell_types)

jaccard_mat = matrix(NA, n, n, dimnames = list(cell_types, cell_types))
for (i in seq_len(n)) {
    for (j in seq_len(n)) {
        a = sig_genes[[i]]; b = sig_genes[[j]]
        jaccard_mat[i, j] = length(intersect(a, b)) / length(union(a, b))
    }
}

col_fun = colorRamp2(c(0, 1), c("white", "#08306b"))

pdf(file.path(plot_dir, 'jaccard_heatmap.pdf'), width = 8, height = 7)
Heatmap(
    jaccard_mat,
    name = "Jaccard",
    col = col_fun,
    clustering_distance_rows = function(m) as.dist(1 - m),
    clustering_distance_columns = function(m) as.dist(1 - m),
    clustering_method_rows = "average",
    clustering_method_columns = "average",
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 90
)
dev.off()

#   Interestingly, DEGs up in extracellular are never cell-type markers in
#   either cellular (this part is not surprising) or extracellular data. We
#   might conclude that extracellular bins are not providing additional signal
#   for differentiating cell types on top of what the cellular bins provide
message("Among markers, which are DEGs and in which direction?")
rbind(
        read_tsv(cell_marker_path, show_col_types = FALSE) |>
            mutate(compartment = 'cell'),
        read_tsv(extra_marker_path, show_col_types = FALSE) |>
            mutate(compartment = 'extra')
    ) |>
    dplyr::rename(cell_type = set_id) |>
    left_join(
        de_df |>
            filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
            mutate(cell_type = str_replace_all(cell_type, '/', '.')) |>
            select(gene_id, cell_type, logFC),
        by = c("cell_type", "gene_id")
    ) |>
    group_by(cell_type, compartment) |>
    summarize(
        frac_up = mean(!is.na(logFC) & (logFC > 0)),
        frac_down = mean(!is.na(logFC) & (logFC < 0)),
        not_sig = mean(is.na(logFC))
    ) |>
    ungroup() |>
    print(n = 30)

marker_df = rbind(
        read_tsv(cell_marker_path, show_col_types = FALSE) |>
            mutate(compartment = 'cell'),
        read_tsv(extra_marker_path, show_col_types = FALSE) |>
            mutate(compartment = 'extra')
    ) |>
    dplyr::rename(cell_type = set_id) |>
    select(cell_type, gene_id, MeanRatio)

cell_df = de_df |>
    filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
    mutate(cell_type = str_replace_all(cell_type, '/', '.')) |>
    select(gene_id, cell_type, logFC) |>
    left_join(
        read_tsv(cell_marker_path, show_col_types = FALSE) |>
            mutate(compartment = 'cell') |>
            dplyr::rename(cell_type = set_id) |>
            select(cell_type, gene_id, MeanRatio),
        by = c("cell_type", "gene_id")
    ) |>
    group_by(sign(logFC)) |>
    summarize(prop_markers = mean(!is.na(MeanRatio))) |>
    mutate(marker_type = 'cell')
extra_df = de_df |>
    filter(adj.P.Val < 0.05, abs(logFC) > 1) |>
    mutate(cell_type = str_replace_all(cell_type, '/', '.')) |>
    select(gene_id, cell_type, logFC) |>
    left_join(
        read_tsv(extra_marker_path, show_col_types = FALSE) |>
            mutate(compartment = 'extra') |>
            dplyr::rename(cell_type = set_id) |>
            select(cell_type, gene_id, MeanRatio),
        by = c("cell_type", "gene_id")
    ) |>
    group_by(sign(logFC)) |>
    summarize(prop_markers = mean(!is.na(MeanRatio))) |>
    mutate(marker_type = 'extra')
p = rbind(cell_df, extra_df) |>
    mutate(
        DEG_type = ifelse(`sign(logFC)` > 0, 'extra DEG', 'cellular DEG'),
        prop_not_markers = 1 - prop_markers
    ) |>
    pivot_longer(
        cols = c(prop_markers, prop_not_markers), names_to = "type",
        values_to = "prop"
    ) |>
    mutate(type = ifelse(type == "prop_markers", "marker", "not marker")) |>
    ggplot(aes(x = marker_type, y = prop, fill = type)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ DEG_type) +
        labs(x = "Marker type", y = "Proportion of DEGs that are Markers") +
        theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'DEG_marker_barplot.pdf'))
print(p)
dev.off()

session_info()
