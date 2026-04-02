library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(sessioninfo)
library(ComplexHeatmap)
library(grid)
library(paletteer)

plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'registration_banksy'
)

# fixed 27-cluster model
model_path_27 = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/registration_banksy/modeling_results/1_8.rds"

# List all paths and names for reference data
ref_paths = here(
    "processed-data", "05_snRNA-seq_model_stats",
    "enrichment_final_Annotations.rds"
)

ref_name = "snRNAseq_fine"

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(
    here('processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy'),
    showWarnings = FALSE, recursive = TRUE
)

# read the fixed 27-cluster model once
t_stat = readRDS(model_path_27)$enrichment

#   Make heatmaps
cluster_anno <- read.csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/registration_banksy/cluster_annotation.csv")

fine_colors = c(
    OPC = '#d3c871',
    Oligo = '#4d5802',
    Microglia = '#222222',
    Astrocyte = '#8d363c',
    Endo = '#ee6c14',
    MHb.1 = '#FF00FF',
    MHb.2 = '#FAA0A0',
    LHb.2.7 = '#00a900',
    LHb.1.3.4 = '#004F2D',
    LHb.4 = '#84DCC6',
    Excit.Thal = '#9e4ad1'
)

cluster_palette = setNames(
    as.character(paletteer::paletteer_d("Polychrome::palette36", 27)),
    as.character(1:27)
)

message("Processing: ", ref_name)

out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    sprintf('cor_vs_%s.rds', ref_name)
)

# load reference data
results_enrichment = list(enrichment = readRDS(ref_paths))

# correlation
this_cor = layer_stat_cor(
    t_stat,
    modeling_results = results_enrichment,
    model_type = "enrichment",
    top_n = 100
)

rownames(this_cor) = sub('^X', '', rownames(this_cor))

annotated_clusters = annotate_registered_clusters(
    this_cor,
    cutoff_merge_ratio = 0.1
)

row_ids = rownames(this_cor)

fine_lookup = setNames(
    cluster_anno$fine_cell_type,
    as.character(cluster_anno$cluster)
)
fine_vec = fine_lookup[row_ids]

if (any(is.na(fine_vec))) {
    stop("Some row names in this_cor do not match cluster_anno.")
}

cluster_vec = factor(row_ids, levels = as.character(1:27))

ha_right = rowAnnotation(
    cluster_id = cluster_vec,
    fine_cell_type = fine_vec,
    col = list(
        cluster_id = cluster_palette,
        fine_cell_type = fine_colors
    ),
    show_annotation_name = FALSE,
    gp = gpar(col = NA),
    annotation_width = unit(c(3, 3), "mm")


p = layer_stat_cor_plot(
    this_cor,
    annotation = annotated_clusters,
    query_colors = NULL,  
    heatmap_legend_param = list(
        title = "Cor",
        at = c(-1, 0, 1)
    ),
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 8)
) + ha_right

pdf(
    file.path(plot_dir, sprintf("cor_vs_%s_27clusters.pdf", ref_name)),
    width = 8,
    height = 10
)
draw(p, merge_legends = TRUE)
dev.off()

# -----------------------------
# Right strips + numbers + cell type labels
# -----------------------------
cluster_palette = setNames(
    as.character(paletteer::paletteer_d("Polychrome::palette36", 27)),
    as.character(1:27)
)

row_ids = rownames(this_cor)

fine_lookup = setNames(
    cluster_anno$fine_cell_type,
    as.character(cluster_anno$cluster)
)
fine_vec = fine_lookup[row_ids]

if (any(is.na(fine_vec))) {
    stop("Some row names in this_cor do not match cluster_anno.")
}

cluster_vec = factor(row_ids, levels = as.character(1:27))

# 调这里
num_fontsize = 12
label_fontsize = 12
cluster_strip_width = unit(4, "mm")
fine_strip_width = unit(4, "mm")

ha_right = rowAnnotation(
    cluster_num = anno_text(
        row_ids,
        just = "left",
        location = 0,
        gp = gpar(fontsize = num_fontsize, col = "black")
    ),
    cluster_id = cluster_vec,
    fine_cell_type_strip = fine_vec,
    fine_cell_type_label = anno_text(
        fine_vec,
        just = "left",
        location = 0,
        gp = gpar(fontsize = label_fontsize, col = "black")
    ),
    col = list(
        cluster_id = cluster_palette,
        fine_cell_type_strip = fine_colors
    ),
    show_annotation_name = FALSE,
    show_legend = c(cluster_id = FALSE, fine_cell_type_strip = FALSE),
    gp = gpar(col = NA),
    annotation_width = unit.c(
        max_text_width(row_ids, gp = gpar(fontsize = num_fontsize)) + unit(1.5, "mm"),
        cluster_strip_width,
        fine_strip_width,
        max_text_width(fine_vec, gp = gpar(fontsize = label_fontsize)) + unit(2, "mm")
    )
)

p = layer_stat_cor_plot(
    this_cor,
    annotation = annotated_clusters,
    query_colors = NULL,
    heatmap_legend_param = list(
        title = "Cor",
        at = c(-1, 0, 1)
    )
) + ha_right



pdf(
    file.path(plot_dir, sprintf("cor_vs_%s_27clusters.pdf", ref_name)),
    width = 8,
    height = 8
)
draw(p, merge_legends = TRUE)
dev.off()

