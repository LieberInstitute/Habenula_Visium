#   Plot requests for Jonathan

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)
library(cowplot)
library(viridis)

cluster_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'jonathan_plots'
)


dir.create(plot_dir, showWarnings = FALSE)

#Human markers
deconvo_marker_path = file.path(dirname(here()),'Hb_multiome','processed-data','05_03_annotation_adjustments','14_deconvoBuddies_markers')
marker_stats_MeanRatio = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_MeanRatio.rds'))
marker_stats_1vAll = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_1vAll.rds'))
marker_stats = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_combo.rds'))


################################################################################
#   Functions
################################################################################

# vis_clus_improved and plot_combo are for plotting annotated spatial clusters, with off-target cell-types in alpha
vis_clus_improved = function(spe, sampleid, clustervar,target_celltypes, colors, alpha_value = 1, flip = FALSE) {
    small_spe = spe[,spe$sample_id == sampleid]

    p = spatialCoords(small_spe) |>
        as_tibble() |>
        mutate(
            x = if (flip) max(pxl_col_in_fullres) - pxl_col_in_fullres else pxl_col_in_fullres,
            y = if (flip) pxl_row_in_fullres else max(pxl_row_in_fullres) - pxl_row_in_fullres,
            cluster = small_spe[[clustervar]],
            alpha_vec = ifelse(small_spe[[clustervar]] %in% target_celltypes, 1, alpha_value)
        ) |>
        ggplot(aes(x = x, y = y, color = cluster, alpha = alpha_vec)) +
            geom_point(size = 0.75) +
            scale_color_manual(values = colors) +
            coord_fixed() +
            labs(color = 'Cell Type', title = sampleid) +
            guides(color = guide_legend(override.aes = list(size = 5)),
                  alpha = 'none') +
            theme_void(base_size = 15)
    
    return(p)
}


plot_combo = function(spe, combo_name, region_name, target_celltypes, alpha_value = 1, flip = FALSE, is_pdf = FALSE) {
    stopifnot(!is.null(spe[[combo_name]]))

    spe[[combo_name]] = factor(
        spe[[combo_name]], levels = names(cluster_combos[[combo_name]])
    )

    all_samples = unique(spe$sample_id)

    for (this_sample_id in all_samples) {
        #   Run twice to overcome a bug with different behavior on the first plot
        for (i in seq_len(2)) {
            p = vis_clus_improved(
                spe, sampleid = this_sample_id, clustervar = combo_name,
                colors = cluster_combos[[combo_name]], target_celltypes, alpha_value, flip
            )
        }

        if (is_pdf) {
            pdf(
                file.path(plot_dir, sprintf('%s_%s_%s.pdf', combo_name, region_name, this_sample_id)),
                width = 8, height = 8
            )
        } else {
            png(
                file.path(plot_dir, sprintf('%s_%s_%s.png', combo_name, region_name, this_sample_id)),
                width = px_per_plot, height = px_per_plot, units = "px"
            )
        }

        print(p)
        dev.off()
    }
}


#Plotting specific genes, spatially
plot_marker_spatial = function(spe_bin, sample_id, genes, px_per_plot, plot_title, 
    flip = FALSE, alpha_value = 1, is_pdf = FALSE, max_exp = NULL){
    
    spe_bin = spe_bin[, spe_bin$sample_id == sample_id]
    grey_viridis <- c("grey90", rev(magma(255)))
    num_genes = length(genes)

    vg_exp = as.matrix(assays(spe_bin)$logcounts[genes,])

    # Gather expression and spatial coordinates into a tidy tibble
    exp_df = vg_exp |>
        t() |>
        as_tibble() |>
        cbind(spatialCoords(spe_bin)) |>
        as_tibble() |>
        mutate(
            x = if (flip) max(pxl_col_in_fullres) - pxl_col_in_fullres else pxl_col_in_fullres,
            y = if (flip) pxl_row_in_fullres else max(pxl_row_in_fullres) - pxl_row_in_fullres,
        )

    plot_list = list()
    for (i in seq_len(length(genes))) {
        #Build gene specific color scale, all maxxed out to the max_exp
        gene_vals = exp_df[[genes[i]]]
        data_max = max(gene_vals, na.rm = TRUE)
        
        if (!is.null(max_exp) && data_max > max_exp) {
            breakpoints = c(
                seq(0, max_exp / data_max, length.out = 256),
                1
            )
            gene_colors = c("grey90", rev(magma(254)), rev(magma(1)))
        } else {
            breakpoints = NULL
            gene_colors = c("grey90", rev(magma(255)))
        }
        
        plot_list[[i]] = ggplot(
                exp_df,
                aes(
                    x = x, y = y, color = !!sym(genes[i]),
                    fill = !!sym(genes[i]),
                    alpha = !!sym(genes[i])
                )
                
            ) +
            geom_point(aes(size = ifelse(!!sym(genes[i]) > 0, 1.5, 0.1))) +
            scale_size_identity() +
            scale_fill_gradientn(colors = gene_colors, values = breakpoints) +
            scale_color_gradientn(colors = gene_colors, values = breakpoints) +
            scale_alpha_continuous(range = c(alpha_value, 1)) +
            coord_fixed() +
            labs(
                fill = 'expr',
                title = rowData(spe_bin)$gene_name[match(genes[i], rownames(spe_bin))]
            ) +
            guides(color = "none", alpha = 'none') +
            theme_bw(base_size = 15) +
            #   Remove pretty much everything related to x- and y-axis labels
            theme(
                axis.title.x = element_blank(), axis.title.y = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                plot.title = element_text(size = 25),
                plot.margin = margin(0, 0, 0, 0, 'pt'),
                legend.key.size = unit(1, "cm"),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18)
            )
    }

    if (is_pdf) {
        pdf(
            file.path(plot_dir, sprintf('%s_marker_genes_%s.pdf', plot_title, sample_id)),
            width = 12, height = 8
        )
    } else {
        png(
            file.path(plot_dir, sprintf('%s_marker_genes_%s.png', plot_title, sample_id)),
            width = px_per_plot * num_genes *0.5, height = px_per_plot * 1
        )
    }
  
    print(plot_grid(plotlist = plot_list, nrow = 1, rel_widths = rep(1, length(plot_list))))
    dev.off()
  
}

#Plotting the z-scored value of aggregated log counts for a set of markers
plot_marker_signature_spatial = function(spe_bin, sample_id, genes, px_per_plot, plot_title, 
    flip = FALSE, alpha_value = 1, is_pdf = FALSE, max_exp = NULL){
    
    spe_bin = spe_bin[, spe_bin$sample_id == sample_id]
    grey_viridis <- c("grey90", rev(magma(255)))

    vg_exp = as.matrix(assays(spe_bin)$logcounts[genes,])
    
    # Sum across genes for each cell, then z-score
    summed_exp = colSums(vg_exp)
    max_norm_exp = summed_exp / max(summed_exp)
    #z_scored_exp = scale(summed_exp)[, 1]

    # Gather expression and spatial coordinates into a tidy tibble
    exp_df = tibble(
        #z_score = z_scored_exp
        max_norm_exp = max_norm_exp
    ) |>
        cbind(spatialCoords(spe_bin)) |>
        as_tibble() |>
        mutate(
            x = if (flip) max(pxl_col_in_fullres) - pxl_col_in_fullres else pxl_col_in_fullres,
            y = if (flip) pxl_row_in_fullres else max(pxl_row_in_fullres) - pxl_row_in_fullres,
            alpha_val = ifelse(max_norm_exp == 0 , 0.05, 1)
        )

    # Compute symmetric range for diverging scale
    #data_max = max(exp_df$z_score, na.rm = TRUE)
    #data_min = min(exp_df$z_score, na.rm = TRUE)
  
    p = ggplot(
            exp_df,
            aes(
                x = x, y = y, color = max_norm_exp,
                fill = max_norm_exp, alpha = alpha_val
            )
        ) +
        geom_point(size = .5) +
        scale_fill_gradient(low = "#D3D3D3", high = "#8B0000") +
        scale_color_gradient(low = "#D3D3D3", high = "#8B0000") +
        scale_alpha_identity() +
        coord_fixed() +
        labs(
            fill = 'max_norm_exp',
            title = plot_title
        ) +
        guides(color = "none", alpha = 'none') +
        theme_bw(base_size = 15) +
        theme(
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
            plot.title = element_text(size = 25),
            plot.margin = margin(0, 0, 0, 0, 'pt'),
            legend.key.size = unit(1, "cm"),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18)
        )

    if (is_pdf) {
        pdf(
            file.path(plot_dir, sprintf('%s_signature_%s.pdf', plot_title, sample_id)),
            width = 8, height = 8
        )
    } else {
        png(
            file.path(plot_dir, sprintf('%s_signature_%s.png', plot_title, sample_id)),
            width = px_per_plot, height = px_per_plot
        )
    }
  
    print(p)
    dev.off()
  
}


#And what about the z-scored values of the marker enrichments
plot_marker_enrichment = function(spe_bin, sample_id, marker_set, px_per_plot, flip = FALSE){
    
    spe_test = spe_bin[, spe_bin$sample_id == sample_id]
  
    #CPM normalization and swap gene names
    assay(spe_test, "cpm") = MetaMarkers::convert_to_cpm(assay(spe_test, "counts"))
    rownames(spe_test) = rowData(spe_test)$gene_name

    #Filter markers for those present in data
    top_current_markers = marker_set  %>%
    select(gene, cellType.target) %>% filter(gene %in% rownames(spe_test))

    colnames(top_current_markers) = c('gene', 'cell_type')
    top_current_markers$group = 'All'


    ct_scores = MetaMarkers::score_cells(log1p(cpm(spe_test)), top_current_markers)
    ct_enrichment = MetaMarkers::compute_marker_enrichment(ct_scores)
    scaled_enrichment = scale(t(ct_enrichment))[ ,]


    colnames(scaled_enrichment) = sapply(strsplit(colnames(scaled_enrichment), split = '|', fixed = TRUE), `[`, 2)
    all_celltypes = colnames(scaled_enrichment)

    # Gather expression and spatial coordinates into a tidy tibble
    exp_df = as.matrix(scaled_enrichment) |>
        as_tibble() |>
        cbind(spatialCoords(spe_test)) |>
        as_tibble() |>
        mutate(
            x = if (flip) max(pxl_col_in_fullres) - pxl_col_in_fullres else pxl_col_in_fullres,
            y = if (flip) pxl_row_in_fullres else max(pxl_row_in_fullres) - pxl_row_in_fullres
        )

    for(this_celltype in all_celltypes){

        p = ggplot(
                exp_df,
                aes(
                    x = x, y = y, color = !!sym(this_celltype)
                )
            ) +
            geom_point(size = .5) +
            scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                                limits = c(NA, 5), oob = scales::squish) +
            coord_fixed() +
            labs(
                color = 'scaled marker enrichment',
                title = sprintf('%s %s', sample_id, this_celltype)
            ) +
            theme_bw(base_size = 15) +
            theme(
                axis.title.x = element_blank(), axis.title.y = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(),
                axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
                plot.title = element_text(size = 25),
                plot.margin = margin(0, 0, 0, 0, 'pt'),
                legend.key.size = unit(1, "cm"),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 18)
            )

        png(
            file.path(plot_dir, sprintf('%s_%s_marker_enrichment.png', sample_id, this_celltype)),
            width = px_per_plot, height = px_per_plot
        )
        print(p)
        dev.off()
    }

}


################################################################################
#   Main
################################################################################

spe = readRDS(spe_path)

anno_df = read_csv(anno_path, show_col_types = FALSE)
spe$banksy = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    pull(banksy)
stopifnot(!any(is.na(spe$banksy)))

#-------------------------------------------------------------------------------
#   Define combinations of clusters to view based on Jonathan's request
#-------------------------------------------------------------------------------
other_color = '#d4d4d4'


my_colors_mid = c(
  Astrocyte = "#972f2f",
  OPC = "#829454",
  Oligo =  "#384a08",
  Microglia = "#141b02",
  Endo = "#f65a45",
  Ependymal = "#dbb369",
  
  Excit.Thal = "#2e6296",
  Inhib.Thal = "#8DADCA",
  LHb.4 = "#082844",
  Inhib_LHb_4.1 = "#9c66c0",
  Inhib_LHb_4.2 = "#5e0c56",
  LHb.1.3 = "#527BAA",
  LHb.1 = "#0C383E",
  LHb.2.7 = "#ee9630",
  LHb.1.3.4 = "#306171",
  MHb.1 = "#5e0c01",
  MHb.1.2 = "#f67104",
  MHb.2 = "#943f02",
  MHb.3 = "#f4d5ab"
) 

my_colors_class <- c(
    LHb = "#0269ae",
    MHb = "#ad1d8c",
    `Non-neurons` = "#532222", 
    Thalamus = "#5906df"
    
)




cluster_combos = list(
    MHb.1 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Thalamus = my_colors_mid[['Excit.Thal']],
        Other = other_color
    ),
    MHb.2 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Thalamus = my_colors_mid[['Excit.Thal']],
        Other = other_color
    ),
    LHb.2.7 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Thalamus = my_colors_mid[['Excit.Thal']],
        Other = other_color
    ),
    LHb.4 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Thalamus = my_colors_mid[['Excit.Thal']],
        Other = other_color
    ),
    Excit_LHb = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Thalamus = my_colors_mid[['Excit.Thal']],
        Other = other_color
    ),
    Inhib_LHb_4.2 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Thalamus = my_colors_mid[['Excit.Thal']],
        Other = other_color
    ),
    summary = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Other = other_color
    ),
    Excit.Thal = c(
        LHb = my_colors_class[["LHb"]],
        MHb = my_colors_class[["MHb"]], 
        Excit.Thal = my_colors_class[["Thalamus"]]     
    ),
    Ependymal = c(
        Ependymal = my_colors_mid[["Ependymal"]],
        Subependymal = "#8d9707",
        Other = other_color
    ),
    Astrocyte = c(
        LHb = my_colors_class[["LHb"]],
        MHb = my_colors_class[["MHb"]],
        OPC = my_colors_mid[["OPC"]],
        Oligo = my_colors_mid[["Oligo"]],
        Astrocyte = my_colors_mid[["Astrocyte"]],
        Endo = my_colors_mid[["Endo"]],
        Ependymal = my_colors_mid[["Ependymal"]],
        Subependymal = "#8d9707",
        Other = other_color
    ),
    Oligo = c(
        LHb = my_colors_class[["LHb"]],
        MHb = my_colors_class[["MHb"]],
        OPC = my_colors_mid[["OPC"]],
        Oligo = my_colors_mid[["Oligo"]],
        Astrocyte = my_colors_mid[["Astrocyte"]],
        Endo = my_colors_mid[["Endo"]],
        Ependymal = my_colors_mid[["Ependymal"]],
        Subependymal = "#8d9707",
        Other = other_color
    ),
    OPC = c(
        LHb = my_colors_class[["LHb"]],
        MHb = my_colors_class[["MHb"]],
        OPC = my_colors_mid[["OPC"]],
        Oligo = my_colors_mid[["Oligo"]],
        Astrocyte = my_colors_mid[["Astrocyte"]],
        Endo = my_colors_mid[["Endo"]],
        Ependymal = my_colors_mid[["Ependymal"]],
        Subependymal = "#8d9707",
        Other = other_color
    ),
    MHb.1_only = c(MHb.1 = my_colors_mid[["MHb.1"]], Other = other_color),
    MHb.2_only = c(MHb.2 = my_colors_mid[["MHb.2"]], Other = other_color),
    LHb.2.7_only = c(LHb.2.7 = my_colors_mid[["LHb.2.7"]], Other = other_color),
    LHb.4_only = c(LHb.4 = my_colors_mid[["LHb.4"]], Other = other_color),
    Excit_LHb_only = c(Excit_LHb = "#0587f9", Other = other_color),
    Astro_only = c(Astrocyte = my_colors_mid[["Astrocyte"]], Other = other_color),
    OPC_only = c(OPC = my_colors_mid[["OPC"]], Other = other_color),
    Oligo_only = c(Oligo = my_colors_mid[["Oligo"]], Other = other_color)

)


spe$MHb.1 = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy == 20 ~ "Thalamus",
    TRUE ~ "Other"
)
spe$MHb.2 = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)
spe$LHb.2.7 = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)
spe$LHb.4 = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)
spe$Excit_LHb = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)
spe$Inhib_LHb_4.2 = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
    TRUE ~ "Other"
)
spe$summary = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    TRUE ~ "Other"
)
spe$Excit.Thal = case_when(
    spe$banksy %in% c(7, 25, 23, 19) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 20 ~ "Excit.Thal"
)
spe$Ependymal = case_when(
    spe$banksy %in% c(12, 16, 24) ~ "Ependymal",
    spe$banksy == 13 ~ "Subependymal",
    TRUE ~ "Other"
)
spe$Astrocyte = case_when(
    spe$banksy %in% c(7, 25, 23, 19) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 1 ~ "OPC",
    spe$banksy %in% c(3,4,10,27) ~ "Oligo",
    spe$banksy %in% c(6, 11, 18) ~ "Astrocyte",
    spe$banksy %in% c(17, 21, 22, 26) ~ "Endo",
    spe$banksy %in% c(12, 16, 24) ~ "Ependymal",
    spe$banksy == 13 ~ "Subependymal",
    TRUE ~ "Other"   
)
spe$Oligo = case_when(
    spe$banksy %in% c(7, 25, 23, 19) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 1 ~ "OPC",
    spe$banksy %in% c(3,4,10,27) ~ "Oligo",
    spe$banksy %in% c(6, 11, 18) ~ "Astrocyte",
    spe$banksy %in% c(17, 21, 22, 26) ~ "Endo",
    spe$banksy %in% c(12, 16, 24) ~ "Ependymal",
    spe$banksy == 13 ~ "Subependymal",
    TRUE ~ "Other"  
)
spe$OPC = case_when(
    spe$banksy %in% c(7, 25, 23, 19) ~ "LHb",
    spe$banksy %in% c(8, 15) ~ "MHb",
    spe$banksy == 1 ~ "OPC",
    spe$banksy %in% c(3,4,10,27) ~ "Oligo",
    spe$banksy %in% c(6, 11, 18) ~ "Astrocyte",
    spe$banksy %in% c(17, 21, 22, 26) ~ "Endo",
    spe$banksy %in% c(12, 16, 24) ~ "Ependymal",
    spe$banksy == 13 ~ "Subependymal",
    TRUE ~ "Other"  
)

spe$MHb.1_only = case_when(
    spe$banksy == 8 ~ "MHb.1",
    TRUE ~ "Other"
)
spe$MHb.2_only = case_when(
    spe$banksy == 15 ~ "MHb.2",
    TRUE ~ "Other"
)
spe$LHb.2.7_only = case_when(
    spe$banksy == 19 ~ "LHb.2.7",
    TRUE ~ "Other"
)
spe$LHb.4_only = case_when(
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    TRUE ~ "Other"
)
spe$Excit_LHb_only = case_when(
    spe$banksy == 23 ~ "Excit_LHb",
    TRUE ~ "Other"
)
spe$Astro_only = case_when(
    spe$banksy %in% c(6, 11, 18) ~ "Astrocyte",
    TRUE ~ "Other"
)
spe$Oligo_only = case_when(
    spe$banksy %in% c(3,4,10,27) ~ "Oligo",
    TRUE ~ "Other"
)
spe$OPC_only = case_when(
    spe$banksy == 1 ~ "OPC",
    TRUE ~ "Other"
)


#-------------------------------------------------------------------------------
#   Plot combinations
#-------------------------------------------------------------------------------

#for (combo_name in names(cluster_combos)) {
#    plot_combo(spe, combo_name)
#}

px_per_plot = 600
for (combo_name in names(cluster_combos)) {
    #pngs
    if(combo_name == 'Ependymal'){
      plot_combo(spe, combo_name,region_name = 'Habenula', target_celltypes = c('Ependymal','Subependymal'), 
      alpha_value = .005, flip = FALSE)
    }else{
    plot_combo(spe, combo_name,region_name = 'Habenula', target_celltypes = combo_name, alpha_value = .005, flip = FALSE)
    }
      
      #pdfs
    #plot_combo(spe, combo_name,region_name = 'Habenula', target_celltypes = combo_name, alpha_value = .005, flip = FALSE, is_pdf = TRUE)
}

#And summary plots of the neurons

plot_combo(spe, combo_name = 'summary', region_name = 'Habenula', 
target_celltypes = c('MHb.1', 'MHb.2', 'LHb.2.7', 'LHb.4', 'Excit_LHb'), 
alpha_value = .005, flip = FALSE, is_pdf = FALSE)


#############
#Marker genes
#############

all_samples = unique(spe$sample_id)

genes = c('GFAP', 'AQP4', 'S100B','SLC1A2','SLC1A3','APOE','VIM','NFIA', 'NFIB')
genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
plot_marker_spatial(spe, sample_id = all_samples[1], genes, px_per_plot, plot_title = 'astrocyte_markers_cells', alpha_value = .1, flip = FALSE)

genes = c('FOXJ1', 'PIFO')
genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
plot_marker_spatial(spe, sample_id = all_samples[1], genes, px_per_plot, plot_title = 'ependymal_markers_cells', alpha_value = .1, flip = FALSE)

genes = c('KCNJ10', 'HCN1', 'SLC12A5', 'CACNA1G', 'CACNA1I')
genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
plot_marker_spatial(spe, sample_id = all_samples[5], genes, px_per_plot, plot_title = 'ephys_genes_of_interest', alpha_value = .1, flip = FALSE)


#MHb1 - TAC3
#MHb2 - GPR149
#LHb2.7, maybe GALR1, 

for(i in 1:length(all_samples)){

    genes = c('SLC32A1', 'GAD1', 'GAD2')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'inhibitory_markers_cells', alpha_value = .1, flip = FALSE)

  
    genes = c('OPRM1', 'TAC3', 'GPR149', 'COL25A1', 'HTR4', 'SLIT1')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'oprm1_markers_cells', alpha_value = .1, flip = FALSE)

    genes = c('TAC1', 'TAC3')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'MHb1_markers_cells', alpha_value = .1, flip = FALSE)


    genes = c('CHAT', 'GPR149')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'MHb2_markers_cells', alpha_value = .1, flip = FALSE)

    genes = c('MME', 'GFRA1', 'SSTR2', 'PDGFD', 'MCC')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'MHb12_markers_cells', alpha_value = .1, flip = FALSE)

    genes = c('ADH1B', 'EBF3', 'SMIM35', 'BHLHE22')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'MHb3_markers_cells', alpha_value = .1, flip = FALSE)

    genes = c('COL25A1', 'GALR1')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'LHb27_markers_cells', alpha_value = .1, flip = FALSE)

    genes = c('HTR4', 'GRIK4')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'LHb134_markers_cells', alpha_value = .1, flip = FALSE)

    genes = c('SLIT1', 'HS3ST4')
    genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
    plot_marker_spatial(spe, sample_id = all_samples[i], genes, px_per_plot, plot_title = 'LHb4_markers_cells', alpha_value = .1, flip = FALSE)

}


#Combinations of top marker genes, using the top 10 1vsAll markers as a starting point

top_1vsAll_marker_df = marker_stats_1vAll %>% group_by(cellType.target) %>% filter(std.logFC.rank <= 50)

#top_1vsAll_marker_df = marker_stats_MeanRatio %>% group_by(cellType.target) %>% filter(MeanRatio.rank <= 20)


dup_genes = top_1vsAll_marker_df$gene[which(duplicated(top_1vsAll_marker_df$gene))]
#For each duplicate, assign it to the cell-type with the better (minimum) rank
keep_dups = top_1vsAll_marker_df %>% filter(gene %in% dup_genes) %>% group_by(gene) %>% filter(std.logFC.rank == min(std.logFC.rank))
#keep_dups = top_1vsAll_marker_df %>% filter(gene %in% dup_genes) %>% group_by(gene) %>% filter(MeanRatio.rank == min(MeanRatio.rank))
top_1vsAll_marker_df = top_1vsAll_marker_df %>% filter(!gene %in% dup_genes)
top_1vsAll_marker_df = rbind(top_1vsAll_marker_df, keep_dups)
top_1vsAll_marker_df = top_1vsAll_marker_df %>% arrange(cellType.target)

top_1vsAll_marker_df %>% group_by(cellType.target) %>% summarise(n = n())





# #MHb1 top _markers
# genes_to_sum = top_1vsAll_marker_df |> filter(cellType.target == 'MHb.1') |> pull(gene)
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# genes_to_sum = genes_to_sum[!is.na(genes_to_sum )]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'MHb1 top markers', alpha_value = 0.1)

# #MHb2 top _markers
# genes_to_sum = top_1vsAll_marker_df |> filter(cellType.target == 'MHb.2') |> pull(gene)
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# genes_to_sum = genes_to_sum[!is.na(genes_to_sum )]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'MHb2 top markers', alpha_value = 0.1)

# #MHb1.2 top _markers
# genes_to_sum = top_1vsAll_marker_df |> filter(cellType.target == 'MHb.1.2') |> pull(gene)
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# genes_to_sum = genes_to_sum[!is.na(genes_to_sum )]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'MHb1.2 top markers', alpha_value = 0.1)

# #MHb3 top _markers
# genes_to_sum = top_1vsAll_marker_df |> filter(cellType.target == 'MHb.3') |> pull(gene)
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# genes_to_sum = genes_to_sum[!is.na(genes_to_sum )]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'MHb3 top markers', alpha_value = 0.1)


# genes_to_sum = c('COL25A1','GALR1','CBLN2','RFTN1','CHRM2','CALN1', 'PRKD1')
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'LHb2.7 top markers', alpha_value = 0.1)

# genes_to_sum = c('MMRN1','GRIK4','SEMA3D','EBF1','HTR4', 'KCNH8', 'TENM1', 'SLC35F3')
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'LHb1.3.4 top markers', alpha_value = 0.1)

# genes_to_sum = c('DAB1','NOVA1','GABRG3','ARPP21', 'HS3ST4', 'SEMA5A', 'GABRB1', 'CACNA2D1', 'SLIT1')
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'LHb4 top markers', alpha_value = 0.1)

# genes_to_sum = c('SLC32A1','PNOC','SIX3','PAX7')
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'Inhib4.2 top markers', alpha_value = 0.1)

# genes_to_sum = c('TFAP2B','NXPH1','NRXN3','GAD2','GRM8', 'PAX3', 'LHX1', 'EMX2', 'ADRA1A', 'GATA3')
# genes_to_sum = rowData(spe)$gene_id[match(genes_to_sum, rowData(spe)$gene_name)]
# plot_marker_signature_spatial(spe, sample_id = all_samples[3], genes_to_sum, px_per_plot, 
#                                plot_title = 'Inhib4.1 top markers', alpha_value = 0.1)






plot_marker_enrichment(spe, sample_id = all_samples[1], marker_set = top_1vsAll_marker_df, px_per_plot = 600, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[2], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[3], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[4], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[5], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[6], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[7], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[8], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[9], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[10], marker_set = top_1vsAll_marker_df, px_per_plot, flip = FALSE)







session_info()

  