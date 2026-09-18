library(tidyverse)
library(here)
library(sessioninfo)
library(ggrepel)
library(clusterProfiler)
library(rrvgo)
library(org.Hs.eg.db)

de_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'aggregated_DE_stats.csv.gz'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'astro_DE', 'results'
)
out_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE'
)
go_num_terms = 5
highlight_genes = c(
    'SLC6A11', 'SLC1A2', 'ENTPD2', 'LTF', 'CHRNA3', 'CALB2', 'SYT6', 'CACNA1E'
)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

custom_volcano <- function(
        data, FDR_cut = 0.05, FC_cut = 1, p_col = "p",
        fdr_col = "fdr", spatial_fdr_col = 'fdr_empirical', 
        lfc_col = "logFC", text = FALSE, highlight_genes = NULL
    ){
    volcano <- data |>
        mutate(DE_class = !!sym(fdr_col) < FDR_cut) |>
        ggplot(
                aes(
                    x = !!sym(lfc_col), y = -log10(!!sym(p_col)),
                    color = DE_class
                )
            ) +
            geom_point(alpha = 0.5, size = 0.5) +
            scale_color_manual(values = c('TRUE' = 'red', 'FALSE' = 'gray')) +
            labs(
                x = "log(FC)", y = "-log10(P value)",
                color = sprintf("FDR < %.2f?", FDR_cut)
            ) +
            theme_bw(base_size = 15) +
            theme(legend.position = "right")
    
    if(!is.null(highlight_genes)){
        volcano <- volcano + 
            geom_text_repel(
                aes(
                    label = ifelse(
                        gene_name %in% highlight_genes, gene_name, ""
                    )
                ),
                size = 3, max.overlaps = 50
            )
    }
    
    return(volcano)
}

plot_go = function(plot_df, plot_path) {
    #   Order GO terms by the first (highest-ranked) direction they appear in
    term_order = plot_df |>
        mutate(de_direction = factor(de_direction, levels = c('up', 'down'))) |>
        group_by(Description) |>
        slice_min(de_direction, n = 1, with_ties = FALSE) |>
        ungroup() |>
        arrange(de_direction) |>
        pull(Description)

    p = plot_df |>
        mutate(
            de_direction = factor(de_direction, levels = c('up', 'down')),
            Description = factor(Description, levels = term_order)
        ) |>
        ggplot(
            aes(
                x = de_direction, y = Description, color = log_fdr,
                size = gene_ratio
            )
        ) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(
            x = "DE Direction", y = "GO Term", color = "-log10(FDR)",
            size = "Gene Ratio"
        )

    pdf(plot_path, width = 5, height = 5)
    print(p)
    dev.off()
}

#   Tree maps using semantic similarity of GO terms
plot_treemap_go = function(go_df, plot_path) {
    sim_matrix = calculateSimMatrix(
        go_df$ID,
        orgdb = "org.Hs.eg.db",
        ont = "BP",
        method = "Rel"
    )

    scores = -log10(go_df$p.adjust)
    names(scores) = go_df$ID
    reduced_terms = reduceSimMatrix(
        sim_matrix,
        scores,
        threshold = 0.7,
        orgdb="org.Hs.eg.db"
    )

    pdf(plot_path)
    print(treemapPlot(reduced_terms))
    dev.off()
}

################################################################################
#   Main
################################################################################

de_df = read_csv(de_path, show_col_types = FALSE) |>
    #   I was asked to show up in MHb as positive
    mutate(logFC = -1 * logFC)

pdf(file.path(plot_dir, 'volcano.pdf'))
print(custom_volcano(de_df, highlight_genes = highlight_genes))
dev.off()

ego_df_list = list()
for (this_de_sign in c(-1, 1)) {
    de_sign_name = ifelse(this_de_sign == 1, "up", "down")

    gene_set = de_df |>
        filter(
            fdr_empirical < 0.05, fdr < 0.05, abs(logFC) > 1,
            sign(logFC) == this_de_sign
        ) |>
        pull(gene_id)

    ego = enrichGO(
        gene          = gene_set,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENSEMBL",
        ont           = "BP",
        universe      = unique(de_df$gene_id),
        pAdjustMethod = "BH",
        pvalueCutoff  = 1,
        qvalueCutoff  = 1
    )
        
    ego_df_list[[de_sign_name]] = as.data.frame(ego) |>
        as_tibble() |>
        filter(p.adjust < 0.05) |>
        mutate(de_direction = de_sign_name)
}

plot_df = bind_rows(ego_df_list) |>
    mutate(
        gene_ratio = Count / as.integer(str_extract(GeneRatio, "(?<=/)[0-9]+")),
        log_fdr = -log10(p.adjust)
    ) |>
    group_by(de_direction) |>
    slice_min(p.adjust, n = go_num_terms, with_ties = FALSE) |>
    ungroup()

#   Custom dot plot by cell type for each DE direction
plot_go(plot_df, file.path(plot_dir, "GO.pdf"))

for (this_de_direction in c("up", "down")) {
    this_plot_df = bind_rows(ego_df_list) |>
        filter(de_direction == this_de_direction)

    plot_treemap_go(
        this_plot_df,
        file.path(plot_dir, sprintf('GO_treemap_%s.pdf', this_de_direction))
    )
}

#   Record the GO terms and associated genes for each DE direction (for a supp
#   table)
bind_rows(ego_df_list) |>
    dplyr::rename(term_id = ID, term_description = Description, gene_id = geneID) |>
    mutate(DE_higher_in = ifelse(de_direction == "up", "MHb", "LHb")) |>
    dplyr::select(term_id, term_description, gene_id, DE_higher_in, p.adjust) |>
    separate_rows(gene_id, sep = "/") |>
    write_csv(file.path(out_dir, "GO_terms.csv"))

session_info()
