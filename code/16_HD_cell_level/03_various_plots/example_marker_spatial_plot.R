#Example script to visualize an aggregate marker set in visium spatial data

library(here)
library(tidyverse)
library(spatialLIBD)
library(sessioninfo)


#Paths for SPE object and cluster annotations
cluster_path = here(
    'processed-data', 'PATH_TO_CLUSTERS'
)
spe_path = here(
    'processed-data', 'PATH_TO_SPE_OBJECT'
)
anno_path = here(
    'processed-data', 'PATH_TO_CLUSTER_ANNOTATIONS'
)

plot_dir = here(
    'plots', 'PATH_TO_SAVE_PLOTS'
)


#Path to markers, these are computed from the single-nucleus data
deconvo_marker_path = file.path(dirname(here()), 'PATH_TO_MARKERS')
marker_stats_1vAll = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_1vAll.rds'))


#Viz function
#Takes in the SPE object, subsets to a particular sample
#I usually use CPM normalization, but the actual marker enrichment computations are derived from log1p(cpm), 
# should be very similar to which ever logcounts normalization you're likely using

#The actual data being plotted is the z-score of the marker enrichment score, in lines 59-61
#Assuming your spe object is segmented nuclei, for each nuclei you first compute the average log1(cpm) for each marker set.
#So for 10 marker sets, each nuclei will have 10 values. This is done with the MetaMarkers::score_cells() function
#Then each nuclei gets an enrichment value for each averaged marker set, essentially avg_marker_1 / avg(all marker sets), and so on
#This is the MetaMarkers::compute_marker_enrichment() function
#You could spatially plot these enrichment values, but I've found better success plotting the z-scored transformation of the enrichments

###########IMPORTANT
#This is currently set to cap the plot color scale at a max of 5, which was fine for the habenula data, but check out what it looks like for your dataset
#Lines 87-88
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

    #Fix up column names
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

    #Plot enrichment values for each marker set separately
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


#Load in SPE object and add in cluster annotations if needed.

spe = readRDS(spe_path)

anno_df = read_csv(anno_path, show_col_types = FALSE)
spe$banksy = tibble(key = spe$key) |>
    left_join(read_csv(cluster_path, show_col_types = FALSE), by = 'key') |>
    pull(banksy)
stopifnot(!any(is.na(spe$banksy)))




#Define your top marker gene sets
#I usually start with 50 and go from there
#This is using deconvoBuddies 1vsAll markers, but you can use which ever marker DE approach. 
top_1vsAll_marker_df = marker_stats_1vAll %>% group_by(cellType.target) %>% filter(std.logFC.rank <= 50)

#I've also found that doing the below to handle duplicate markers across cell-types tends to clean up signal.
#Though, this may end up filtering out the majority of markers for a cell-type if it is not distinct and overlaps with a bunch of others
#Usually not a problem if you're only filtering out a handful per cell-type, but the enrichments will be less robust when the sample size is more variable,
#like using gene sets that range in size from 10 to 50. So I like to check the duplicate filtering after

#For each duplicate, assign it to the cell-type with the better (minimum) rank
dup_genes = top_1vsAll_marker_df$gene[which(duplicated(top_1vsAll_marker_df$gene))]
keep_dups = top_1vsAll_marker_df %>% filter(gene %in% dup_genes) %>% group_by(gene) %>% filter(std.logFC.rank == min(std.logFC.rank))
top_1vsAll_marker_df = top_1vsAll_marker_df %>% filter(!gene %in% dup_genes)
top_1vsAll_marker_df = rbind(top_1vsAll_marker_df, keep_dups)
top_1vsAll_marker_df = top_1vsAll_marker_df %>% arrange(cellType.target)
#Check duplicate filtering
top_1vsAll_marker_df %>% group_by(cellType.target) %>% summarise(n = n())


#And plot
all_samples = unique(spe$sample_id)

plot_marker_enrichment(spe, sample_id = all_samples[1], marker_set = top_1vsAll_marker_df, px_per_plot = 600, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[2], marker_set = top_1vsAll_marker_df, px_per_plot = 600, flip = FALSE)
plot_marker_enrichment(spe, sample_id = all_samples[3], marker_set = top_1vsAll_marker_df, px_per_plot = 600, flip = FALSE)
#......
#......
#and so on
#......
#......
#......
