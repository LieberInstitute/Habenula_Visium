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


plot_combo = function(spe, combo_name,region_name, target_celltypes, alpha_value = 1, flip = FALSE, is_pdf = FALSE) {
    stopifnot(!is.null(spe[[combo_name]]))

    spe[[combo_name]] = factor(
        spe[[combo_name]], levels = names(cluster_combos[[combo_name]])
    )

    all_samples = unique(spe$sample_id)

    p_list = list()
    for (this_sample_id in all_samples) {
        #   Run twice to overcome a bug with different behavior on the first
        #   plot
        for (i in seq_len(2)) {
            p_list[[this_sample_id]] = vis_clus_improved(
                    spe, sampleid = this_sample_id, clustervar = combo_name,
                    colors = cluster_combos[[combo_name]], target_celltypes, alpha_value, flip

                )
        }
    }
    p = plot_grid(plotlist = p_list, nrow = 1)
    
    if (is_pdf) {
        pdf(
            file.path(plot_dir, sprintf('%s_%s.pdf', combo_name, region_name)),
            width = 12, height = 6
        )
    } else {
        png(
            file.path(plot_dir, sprintf('%s_%s.png', combo_name, region_name)),
            width = px_per_plot * length(all_samples) / 1, height = px_per_plot * 1
        )
    }
    
    print(p)
    dev.off()
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

my_colors_mid = c(Excit.Thal = "#4d55b7",
  LHb.4 = "#00607A",
  Inhib.Thal = "#9a9fe7",
  Astrocyte = "#532222",
  MHb.1.2 = "#92007C",
  LHb.1 = "#008092",
  OPC = "#829454",
  Oligo =  "#384a08",
  Microglia = "#141b02",
  LHb.2.7 = "#6C9FA9",
  Endo = "#d95f02",
  LHb.1.3.4 = "#A8B8BC",
  MHb.1 = "#A86A9A",
  MHb.2 = "#BCA6B6",
  MHb.3 = "#56204eff",
  LHb.1.3 = "#C6C6C6",
  Inhib_LHb_4.1 = "#8B0000",
  Inhib_LHb_4.2 = "#DC143C",
  Ependymal = "#f5a105ff"
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
        Inhib_LHb_4.2 = my_colors_mid[["Inhib_LHb_4.2"]],
        Other = other_color
    ),
    MHb.2 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Inhib_LHb_4.2 = my_colors_mid[["Inhib_LHb_4.2"]],
        Other = other_color
    ),
    LHb.2.7 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Inhib_LHb_4.2 = my_colors_mid[["Inhib_LHb_4.2"]],
        Other = other_color
    ),
    LHb.4 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Inhib_LHb_4.2 = my_colors_mid[["Inhib_LHb_4.2"]],
        Other = other_color
    ),
    Excit_LHb = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Inhib_LHb_4.2 = my_colors_mid[["Inhib_LHb_4.2"]],
        Other = other_color
    ),
    Inhib_LHb_4.2 = c(
        MHb.1 = my_colors_mid[["MHb.1"]],
        MHb.2 = my_colors_mid[["MHb.2"]],
        LHb.2.7 = my_colors_mid[["LHb.2.7"]],
        LHb.4 = my_colors_mid[["LHb.4"]],
        Excit_LHb = "#0587f9",
        Inhib_LHb_4.2 = my_colors_mid[["Inhib_LHb_4.2"]],
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
    )
)



spe$MHb.1 = case_when(
    spe$banksy == 8 ~ "MHb.1",
    spe$banksy == 15 ~ "MHb.2",
    spe$banksy == 19 ~ "LHb.2.7",
    spe$banksy %in% c(7, 25) ~ "LHb.4",
    spe$banksy == 23 ~ "Excit_LHb",
    spe$banksy %in% c(5, 9) ~ "Inhib_LHb_4.2",
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

genes = c('GFAP', 'AQP4')
genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
plot_marker_spatial(spe, sample_id = all_samples[1], genes, px_per_plot, plot_title = 'astrocyte_markers_cells', alpha_value = .1, flip = FALSE)

genes = c('FOXJ1', 'PIFO')
genes = rowData(spe)$gene_id[match(genes, rowData(spe)$gene_name)]
plot_marker_spatial(spe, sample_id = all_samples[1], genes, px_per_plot, plot_title = 'ependymal_markers_cells', alpha_value = .1, flip = FALSE)

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





session_info()

