#!/usr/bin/env Rscript

# Visualize gene set enrichment results as two heatmaps:
# 1) Genes per cell type × disease
# 2) TFs per cell type × disease
# Fill = -log10(FDR), annotations = significance stars

library(here)
library(tidyverse)
library(pheatmap)

# -----------------------------
# paths
# -----------------------------
out_path <- here(
    "processed-data", "10_HD_bin_level", "no_secondary", "liana2", "open_targets"
)

plot_dir <- here(
    "plots", "10_HD_bin_level", "no_secondary", "liana2", "open_targets"
)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# load results
# -----------------------------
res <- read_csv(
    file.path(out_path, "geneset_enrichment_MDD_Substance_dependence.csv"),
    show_col_types = FALSE
)

# -----------------------------
# cell type order (subset of interest)
# -----------------------------
cell_type_order <- c(
    "MHb_A", "MHb_B", "MHb_C",
    "LHb_A", "LHb_B", "LHb_C",
    "GABA_LHb_C.1", "GABA_LHb_C.2"
)

# -----------------------------
# separate genes and TFs, extract cell type
# -----------------------------
res <- res |>
    mutate(
        set_type = case_when(
            str_ends(gene_set_name, "_genes") ~ "genes",
            str_ends(gene_set_name, "_TFs") ~ "TFs",
            TRUE ~ NA_character_
        ),
        cell_type = str_replace(gene_set_name, "_(genes|TFs)$", "")
    ) |>
    filter(!is.na(set_type), cell_type %in% cell_type_order) |>
    mutate(
        cell_type = factor(cell_type, levels = cell_type_order),
        sig_label = case_when(
            fdr <= 0.01 ~ "**",
            fdr <= 0.05 ~ "*",
            TRUE ~ ""
        )
    )

# -----------------------------
# helper: make heatmap for a given set_type
# -----------------------------
make_heatmap <- function(data, set_type_label, out_file) {
    # Pivot to matrix: rows = disease, cols = cell_type
    mat <- data |>
        select(disease, cell_type, neglog10_fdr) |>
        pivot_wider(names_from = cell_type, values_from = neglog10_fdr) |>
        column_to_rownames("disease") |>
        as.matrix()

    # Annotation matrix (stars)
    annot_mat <- data |>
        select(disease, cell_type, sig_label) |>
        pivot_wider(names_from = cell_type, values_from = sig_label) |>
        column_to_rownames("disease") |>
        as.matrix()

    # Ensure column order
    cols <- cell_type_order[cell_type_order %in% colnames(mat)]
    mat <- mat[, cols, drop = FALSE]
    annot_mat <- annot_mat[, cols, drop = FALSE]

    # Replace NA with 0
    mat[is.na(mat)] <- 0
    annot_mat[is.na(annot_mat)] <- ""

    # Plot
    pdf(out_file, width = 10, height = 4)
    pheatmap(
        mat,
        display_numbers = annot_mat,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red"))(100),
        border_color = "grey80",
        main = paste0(
            "Disease-risk Enrichment: ", set_type_label,
            "\n-log10(FDR); * FDR\u22640.05, ** FDR\u22640.01"
        ),
        fontsize_number = 14,
        number_color = "black",
        angle_col = 45
    )
    dev.off()
    cat("Saved:", out_file, "\n")
}

# -----------------------------
# heatmap 1: genes
# -----------------------------
genes_data <- res |> filter(set_type == "genes")

make_heatmap(
    genes_data,
    "Target Genes",
    file.path(plot_dir, "geneset_enrichment_genes_heatmap.pdf")
)

# -----------------------------
# heatmap 2: TFs
# -----------------------------
tfs_data <- res |> filter(set_type == "TFs")

make_heatmap(
    tfs_data,
    "Transcription Factors",
    file.path(plot_dir, "geneset_enrichment_TFs_heatmap.pdf")
)

# -----------------------------
# heatmaps with ALL cell types
# -----------------------------
cell_type_order_all <- c(
    "MHb_A", "MHb_B", "MHb_C",
    "LHb_A", "LHb_B", "LHb_C",
    "GABA_LHb_C.1", "GABA_LHb_C.2",
    "Excit.Thal", "Inhib.Thal",
    "Astrocyte", "Microglia", "Ependymal", "Oligo", "OPC"
)

# Re-read and re-process with all cell types
res_all <- read_csv(
    file.path(out_path, "geneset_enrichment_MDD_Substance_dependence.csv"),
    show_col_types = FALSE
) |>
    mutate(
        set_type = case_when(
            str_ends(gene_set_name, "_genes") ~ "genes",
            str_ends(gene_set_name, "_TFs") ~ "TFs",
            TRUE ~ NA_character_
        ),
        cell_type = str_replace(gene_set_name, "_(genes|TFs)$", "")
    ) |>
    filter(!is.na(set_type), cell_type %in% cell_type_order_all) |>
    mutate(
        cell_type = factor(cell_type, levels = cell_type_order_all),
        sig_label = case_when(
            fdr <= 0.01 ~ "**",
            fdr <= 0.05 ~ "*",
            TRUE ~ ""
        )
    )

make_heatmap_all <- function(data, set_type_label, out_file) {
    mat <- data |>
        select(disease, cell_type, neglog10_fdr) |>
        pivot_wider(names_from = cell_type, values_from = neglog10_fdr) |>
        column_to_rownames("disease") |>
        as.matrix()

    annot_mat <- data |>
        select(disease, cell_type, sig_label) |>
        pivot_wider(names_from = cell_type, values_from = sig_label) |>
        column_to_rownames("disease") |>
        as.matrix()

    cols <- cell_type_order_all[cell_type_order_all %in% colnames(mat)]
    mat <- mat[, cols, drop = FALSE]
    annot_mat <- annot_mat[, cols, drop = FALSE]

    mat[is.na(mat)] <- 0
    annot_mat[is.na(annot_mat)] <- ""

    pdf(out_file, width = 14, height = 4)
    pheatmap(
        mat,
        display_numbers = annot_mat,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        color = colorRampPalette(c("white", "orange", "red"))(100),
        border_color = "grey80",
        main = paste0(
            "Disease-risk Enrichment: ", set_type_label, " (all cell types)",
            "\n-log10(FDR); * FDR\u22640.05, ** FDR\u22640.01"
        ),
        fontsize_number = 14,
        number_color = "black",
        angle_col = 45
    )
    dev.off()
    cat("Saved:", out_file, "\n")
}

genes_data_all <- res_all |> filter(set_type == "genes")
make_heatmap_all(
    genes_data_all,
    "Target Genes",
    file.path(plot_dir, "geneset_enrichment_genes_heatmap_all.pdf")
)

tfs_data_all <- res_all |> filter(set_type == "TFs")
make_heatmap_all(
    tfs_data_all,
    "Transcription Factors",
    file.path(plot_dir, "geneset_enrichment_TFs_heatmap_all.pdf")
)

cat("Done.\n")
sessionInfo()
