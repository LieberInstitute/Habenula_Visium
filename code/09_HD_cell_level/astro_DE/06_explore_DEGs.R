library(tidyverse)
library(here)
library(sessioninfo)
library(qs2)
library(jaffelab)
library(SpatialExperiment)
library(spatialLIBD)
library(edgeR)
library(limma)

de_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'aggregated_DE_stats.csv.gz'
)
dge_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'DGE.qs2'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'spe', 'cleaned.qs2'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'astro_DE', 'explore_DEGs'
)
de_formula = ~ astro_label + ncells + expr_chrM_ratio
cont_covariates = c('ncells', 'expr_chrM_ratio')

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

deg_boxplots = function(dge, genes, cleaningY = FALSE) {
    plot_df_list = list()
    for (gene in genes) {
        if (cleaningY) {
            gene_expr = as.numeric(
                cleaningY(
                    dge$E$E[gene, , drop = FALSE], mod = dge$design, P = 2
                )
            )
            y_label = "cleaningY logcounts"
        } else {
            y <- dge$EList$E[gene, ]
            w <- dge$EList$weights[match(gene, rownames(dge$EList$E)), ]
            X <- dge$design
            fit <- lm.wfit(X, y, w)

            #   Manually implement cleaningY with P = 2
            nuisance_fit <- X[, -(1:2), drop = FALSE] %*%
                fit$coefficients[-(1:2)]
            gene_expr <- y - nuisance_fit
            y_label = "Limma-cleaned logcounts"
        }
      
        plot_df_list[[gene]] = tibble(
            astro_label = str_extract(rownames(dge$design), '(lateral|medial)'),
            sample_id = str_extract(rownames(dge$design), '^Br[0-9]{4}_[12]'),
            gene_id = dge$genes[gene, ]$gene_name,
            expr = gene_expr
        )
    }
  
    p = bind_rows(plot_df_list) |>
        ggplot(aes(x = astro_label, y = expr)) +
            geom_boxplot() +
            geom_jitter(width = 0.1) +
            facet_wrap(~gene_id, scales = "free_y") +
            labs(x = "Astrocyte label", y = y_label) +
            theme_bw(base_size = 15)
    
    return(p)
}

################################################################################
#   Main
################################################################################

dge = qs_read(dge_path)

de_df = read_csv(de_path, show_col_types = FALSE)

message(
    sprintf(
        "Spearman corr between voom-lmFit FDR and empirical FDR: %.2f",
        cor(de_df$fdr, de_df$fdr_empirical, method = 'spearman')
    )
)

genes = de_df |>
    filter(logFC > 10) |>
    pull(gene_id)
all_genes = genes
pdf(file.path(plot_dir, "high_logFC_boxplots.pdf"))
print(deg_boxplots(dge, genes))
dev.off()

genes = de_df |>
    arrange(p_empirical, desc(logFC)) |>
    slice_head(n = 2) |>
    pull(gene_id)
all_genes = c(all_genes, genes)
pdf(file.path(plot_dir, "top_2_p_boxplots.pdf"))
print(deg_boxplots(dge, genes))
dev.off()

spe = qs_read(spe_path)
for (gene in all_genes) {
    gene_name = dge$genes[gene, ]$gene_name
  
    p = vis_gene(
        spe, sampleid = 'Br9090_1', geneid = gene, assay = 'counts',
        is_stitched = TRUE, point_size = 20, spatial = FALSE
    )
    png(
        file.path(plot_dir, sprintf('spatial_%s.png', gene_name)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

session_info()
