library(tidyverse)
library(here)
library(sessioninfo)
library(duckplyr)
library(ComplexHeatmap)
library(circlize)

de_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'main_results'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'cell_vs_extra_DE',
    'check_results'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

################################################################################
#   Functions
################################################################################

custom_volcano <- function(data, 
                           clus, 
                           FDR_cut = 0.05, 
                           p_col = "vlmf_P.Value", 
                           fdr_col = "vlmf_adj.P.Val", 
                           lfc_col = "vlmf_logFC", 
                           model_name = "carrier",
                           save = TRUE,
                           save_size = 4,
                           text = TRUE,
                           highlight_genes = NULL){
  
  # define colors
  signif_colors <- c("purple", "blue", "red")
  names(signif_colors) <- c("both", paste("FDR<", FDR_cut) , "abs(logFC)>1" )
  
  volcano <- data |>
    filter(cluster == clus) |>
    mutate(DE_class = case_when(!!sym(fdr_col) < FDR_cut & abs(vlmf_logFC) > 1  ~ "both",
                                !!sym(fdr_col) < FDR_cut ~ paste("FDR<", FDR_cut),
                                abs(!!sym(lfc_col)) > 1 ~ "abs(logFC)>1",
                                TRUE ~ "None")) |>
    ggplot(aes(x = vlmf_logFC, y = -log10(!!sym(p_col)), color = DE_class)) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_color_manual(values = signif_colors) +
    labs(x = "log(FC)", y = "-log10(P value)") +
    theme_bw() +
    labs(title = clus) +
    theme(legend.position = "right")
  
  if(!text){
    ## dont add text
  } else if(nrow(data) > 1000){
    volcano <- volcano + 
      # geom_text_repel(aes(label = ifelse(DE_class != "None", gene_name, "")), size = 1.5) 
      geom_text_repel(aes(label = ifelse(!!sym(fdr_col) < FDR_cut, gene_name, "")), size = 1.5) 
  } else {
    volcano <- volcano + geom_text_repel(aes(label = gene_name), size = 1.5)
  }
  
  if(!is.null(highlight_genes)){
    volcano <- volcano + 
      geom_text_repel(aes(label = ifelse(gene_name %in% highlight_genes, gene_name, "")), size = 1.5)
  }
  
  ## Save or return
  if(save) ggsave(volcano, filename = here(plot_dir, sprintf("Volcano_%s_%s-%s.png", datatype, model_name, clus)), 
                  height = save_size + 0.5, 
                  width = save_size)
  else return(volcano)    
}

################################################################################
#   Main
################################################################################

dir.create(plot_dir, showWarnings = FALSE)

de_df = list.files(de_dir, full.names = TRUE) |>
    map_dfr(read_parquet_duckdb, prudence = 'lavish') |>
    collect()

p = de_df |>
    mutate(neg_log10_p = -log10(adj.P.Val), is_sig = adj.P.Val < 0.05) |>
    ggplot(aes(x = logFC, y = neg_log10_p, color = is_sig)) +
        geom_point(size = 0.3, alpha = 0.5) +
        scale_color_manual(
            values = c("FALSE" = "grey70", "TRUE" = "firebrick"), guide = "none"
        ) +
        facet_wrap(~ cell_type, ncol = 3) +
        labs(x = "log FC", y = expression(-log[10](adj.~p))) +
        theme_bw(base_size = 15)
pdf(file.path(plot_dir, 'volcano_plots.pdf'), width = 6, height = 10)
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

session_info()
