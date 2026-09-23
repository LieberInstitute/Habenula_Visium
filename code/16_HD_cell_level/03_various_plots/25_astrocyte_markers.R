#   Are literature astrocyte markers expressed in the MHb?

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)
library(scran)
library(qs2)

spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
spe_nuclear_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'contamination','spe','cleaned.qs2'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'astrocyte_markers'
)
cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
sample_ids = c('Br9090_1', 'Br8433_1')
markers = c(
    "GFAP", "AQP4", "S100B", "SLC1A2", "SLC1A3", "APOE", "VIM"
)
cluster_colors = c(
    '6' = "#502419", '11' = '#19647E', '18' = '#F29559',
    'Other' = '#959595'
)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vis_clus_hd = function(
        spe, clustervar, sample_id, plot_path, colors = cluster_colors
    ) {
    #   Run twice to overcome a bug with different behavior on the first plot
    for (i in seq_len(2)) {
        p = vis_clus(
                spe, sampleid = sample_id, clustervar = clustervar,
                is_stitched = TRUE, point_size = 30, spatial = FALSE,
                colors = colors
            ) +
            guides(fill = guide_legend(override.aes = list(size = 8)))
    }
    png(plot_path, width = 1500, height = 1500)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)
spe_nuclear = qs_read(spe_nuclear_path)

#Spe_nuclear needs logcounts
spe_nuclear = computeLibraryFactors(spe_nuclear)
spe_nuclear = logNormCounts(spe_nuclear)


#Ensembl IDs for markers
stopifnot(all(markers %in% rowData(spe)$gene_name))
markers_ensembl = rowData(spe)$gene_id[match(markers, rowData(spe)$gene_name)]
#Clusters and cluster annotations for full cell segmentated data
anno_df = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    dplyr::rename(cluster = banksy) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster')

spe$banksy_cluster = anno_df$cluster
spe$astro_cluster = ifelse(
    anno_df$cluster %in% c(6, 11, 18), as.character(anno_df$cluster), 'Other'
)
spe$cell_type = anno_df$fine_cell_type
stopifnot(!any(is.na(spe$cell_type)))

#Clusters and cluster annotations for nuclear segmentated data
anno_nuclear_df = tibble(key = spe_nuclear$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    dplyr::rename(cluster = banksy) |>
    left_join(read_csv(ct_anno_path, show_col_types = FALSE), by = 'cluster')

spe_nuclear$banksy_cluster = anno_nuclear_df$cluster
spe_nuclear$astro_cluster = ifelse(
    anno_nuclear_df$cluster %in% c(6, 11, 18), as.character(anno_nuclear_df$cluster), 'Other'
)
spe_nuclear$cell_type = anno_nuclear_df$fine_cell_type
stopifnot(!any(is.na(spe_nuclear$cell_type)))

 


for (sample_id in sample_ids) {
    vis_clus_hd(
        spe, 'astro_cluster', sample_id,
        file.path(plot_dir, sprintf('Astro_banksy_%s.png', sample_id))
    )

    p = vis_gene(
        spe, sampleid = sample_id, geneid = markers_ensembl,
        is_stitched = TRUE, point_size = 30, spatial = FALSE,
        cap_percentile = 0.99
    )
    png(
        file.path(plot_dir, sprintf('Astro_markers_%s.png', sample_id)),
        width = 1500, height = 1500
    )
    print(p)
    dev.off()
}

#Bubble plot of astro and habenula markers for the astro and MHb fine resolution clusters

top_markers = markers = c(
    "GFAP", "AQP4", "S100B", "SLC1A2", "SLC1A3", "APOE", "VIM",
    'GPR151','POU4F1','TAC1'
)
top_markers_ensembl = rowData(spe)$gene_id[match(top_markers, rowData(spe)$gene_name)]



bubble_plot_spatial = function(spe, donor_id, clusters_of_interest, top_markers, top_markers_ensembl, title, group_order = NULL){

  donor_spe = spe[ , spe$sample_id == donor_id]
  cell_index = donor_spe$banksy_cluster %in% clusters_of_interest
  expr_data <- assay(donor_spe, 'logcounts')[top_markers_ensembl, cell_index ]
  metadata <- colData(donor_spe)[cell_index, ]
  
  # Convert to data frame for plotting
  plot_data <- as.data.frame(t(as.matrix(expr_data))) %>%
    tibble::rownames_to_column("cell_id") %>%
    cbind(group_var = metadata[['banksy_cluster']]) %>%
    tidyr::pivot_longer(cols = -c(cell_id, group_var), names_to = "gene", values_to = "expression")
  # Calculate mean expression and percent expressing per cluster
  summary_data <- plot_data %>% 
    group_by(gene, group_var) %>%
    summarise(
      mean_expression = mean(expression),
      pct_expressing = sum(expression > 0) / n() * 100,
      .groups = "drop"
    ) %>%
    # Calculate z-score of mean_expression per gene across clusters
    group_by(gene) %>%
    mutate(mean_expression_zscore = scale(mean_expression)[,1]) %>%
    ungroup()

    index = match(summary_data$gene, top_markers_ensembl)
    summary_data$gene_name = top_markers[index]
  # Set factor levels to control axis order
  summary_data$gene_name <- factor(summary_data$gene_name, levels = top_markers)
  if(is.null(group_order)){
    summary_data$group_var <- factor(summary_data$group_var, 
                                      levels = sort(unique(summary_data$group_var)))
  } else {
    summary_data$group_var <- factor(summary_data$group_var, levels = group_order)
  }

  p2 = ggplot(summary_data, aes(x = gene_name, y = group_var, size = pct_expressing, color = mean_expression_zscore)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = 'white', high = "red", midpoint = 0,
    #limits = c(-2, 2),
    #oob = scales::squish,
    name = "Mean Exp. z-score") +
    theme_minimal() + ggtitle(sprintf("%s: %s", title, donor_id)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Gene", y = 'Fine res spatial cluster', size = "% Expressing", color = "Mean Exp. z-score")
  
  return(p2)
}

#We don't have nuclear counts for the 9902 sample
all_samples = names(table(spe$sample_id))
all_samples = all_samples[!grepl('9902', all_samples)]

clusters_of_interest = c(6, 18,11, 8, 15)
group_order = clusters_of_interest

for(i in 1:length(all_samples)){
    donor_id = all_samples[i]

    p_cell = bubble_plot_spatial(spe = spe, donor_id = donor_id, 
        clusters_of_interest = clusters_of_interest, 
        top_markers = top_markers, top_markers_ensembl = top_markers_ensembl, 
        title = "Full cell", group_order = group_order)

    p_nuclear = bubble_plot_spatial(spe = spe_nuclear, donor_id = donor_id, 
        clusters_of_interest = clusters_of_interest, 
        top_markers = top_markers, top_markers_ensembl = top_markers_ensembl, 
        title = "Nuclear", group_order = group_order)
  
  ggsave(
    filename = file.path(plot_dir, sprintf("bubble_plot_%s.pdf", donor_id)),
    plot = cowplot::plot_grid(p_cell, p_nuclear, nrow = 1), device = 'pdf',
    width = 12, height = 6
  )
}

session_info()
