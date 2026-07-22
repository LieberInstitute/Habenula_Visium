library(tidyverse)
library(here)
library(sessioninfo)
library(qs2)
library(jaffelab)
library(SpatialExperiment)
library(spatialLIBD)

de_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'main_results', 'aggregated_DE_stats.csv.gz'
)
sce_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'astro_sce.qs2'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'contamination',
    'spe', 'raw.qs2'
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

deg_boxplots = function(sce, mod, genes, clean = TRUE) {
    plot_df_list = list()
    for (gene in genes) {
        if (clean) {
            gene_expr = as.numeric(
                cleaningY(
                    logcounts(sce)[gene, , drop = FALSE], mod = mod, P = 2
                )
            )
            y_label = "Cleaned logcounts"
        } else {
            gene_expr = unname(logcounts(sce)[gene, ])
            y_label = "Logcounts"
        }
      
        plot_df_list[[gene]] = tibble(
            astro_label = sce$astro_label,
            sample_id = sce$sample_id,
            gene_id = rowData(sce_pb[gene, ])$gene_name,
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

#   Recompute sce_pb as in the main DE
sce = qs_read(sce_path)
sce = sce[, sce$astro_label %in% c('medial', 'lateral')]
sce_pb = registration_pseudobulk(
    sce, var_registration = "astro_label", var_sample_id = "sample_id"
)
sce_pb$astro_label = factor(sce_pb$astro_label, levels = c('medial', 'lateral'))
  sce_pb$sum_umi = unname(colSums(counts(sce_pb)))
sce_pb$expr_chrM = colSums(
    counts(sce_pb)[which(seqnames(sce_pb) == "chrM"), , drop = FALSE]
)
sce_pb$expr_chrM_ratio = sce_pb$expr_chrM / sce_pb$sum_umi

#   Center and scale continuous covariates
for (this_covariate in cont_covariates) {
    sce_pb[[this_covariate]] = as.numeric(scale(sce_pb[[this_covariate]]))
}
  
mod = model.matrix(de_formula, colData(sce_pb))

de_df = read_csv(de_path, show_col_types = FALSE)

genes = de_df |>
    filter(logFC > 10) |>
    pull(gene_id)
all_genes = genes
pdf(file.path(plot_dir, "high_logFC_boxplots.pdf"))
print(deg_boxplots(sce_pb, mod, genes))
dev.off()

genes = de_df |>
    arrange(p_empirical, desc(logFC)) |>
    slice_head(n = 2) |>
    pull(gene_id)
all_genes = c(all_genes, genes)
pdf(file.path(plot_dir, "top_2_p_boxplots.pdf"))
print(deg_boxplots(sce_pb, mod, genes))
dev.off()

for (gene in all_genes) {
    gene_name = rowData(sce_pb[gene, ])$gene_name
  
    p = vis_clus(
            spe, sampleid = sample_id, geneid = gene,
            is_stitched = TRUE, point_size = 20, spatial = FALSE
        ) +
        guides(fill = guide_legend(override.aes = list(size = 8)))
    png(
        file.path(plot_dir, sprintf('spatial_%s.png', gene_name)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

session_info()
