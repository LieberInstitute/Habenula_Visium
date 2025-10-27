# Make the plot for the 4 metrics for registration figure
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# Broad plot
plot_path = here(
    'plots', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'summary_score_4-metric_lineplot_broad.png'
)
summary_score_broad_batch = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'cleany',
        'heatmap_score_snRNAseq_broad.csv'
    ),
    show_col_types = FALSE
)

summary_score_broad_batch <- summary_score_broad_batch %>%
    pivot_longer(
        cols = c(
            unique_cell_type_count, 
            unique_mhb_lhb_counts, 
            one_to_one_cell_type_count, 
            multiple_mhb_lhb_counts
        ),
        names_to = "type",
        values_to = "value" 
    ) %>%
    select(original_cluster_count, type, value) %>%
    mutate(type = recode(type, 
                         unique_cell_type_count = "Unique Cell Type Count", 
                         unique_mhb_lhb_counts = "Unique MHB LHB Counts", 
                         one_to_one_cell_type_count = "One to One Cell Type Count", 
                         multiple_mhb_lhb_counts = "Multiple MHB LHB Counts")) %>%
    mutate(type2="Cleany")

summary_score_broad_libsize = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'normalized',
        'heatmap_score_snRNAseq_broad.csv'
    ),
    show_col_types = FALSE
)

summary_score_broad_libsize <- summary_score_broad_libsize %>%
    pivot_longer(
        cols = c(
            unique_cell_type_count, 
            unique_mhb_lhb_counts, 
            one_to_one_cell_type_count, 
            multiple_mhb_lhb_counts
        ),
        names_to = "type",
        values_to = "value" 
    ) %>%
    select(original_cluster_count, type, value) %>%
    mutate(type = recode(type, 
                         unique_cell_type_count = "Unique Cell Type Count", 
                         unique_mhb_lhb_counts = "Unique MHB LHB Counts", 
                         one_to_one_cell_type_count = "One to One Cell Type Count", 
                         multiple_mhb_lhb_counts = "Multiple MHB LHB Counts")) %>%
    mutate(type2="Normalized")

summary_score_broad <- rbind(summary_score_broad_batch, summary_score_broad_libsize)

summary_score_broad <- summary_score_broad %>% mutate(type=factor(type, 
                         levels = c("Unique Cell Type Count", 
                                    "One to One Cell Type Count",
                                    "Multiple MHB LHB Counts",
                                    "Unique MHB LHB Counts")))
# Make the barplot figure for broad

p1 = ggplot(summary_score_broad, aes(x = original_cluster_count, y = value, color = type2)) +
  geom_line() +
  geom_point(size=1) +
  labs(x = "K",
       y = "N of cell type",
       color = "Score Type",
       title = "With broad RNAseq reference") +
    scale_color_manual(values = c("Cleany" = "#1f77b4", "Normalized" = "#ff7f0e")) +
    facet_wrap(~ type , scales = "free") +
    theme_classic()

png(plot_path, width = 10, height = 10, units = "in", res = 300)
print(p1)
dev.off()

 # Fine plot
 plot_path = here(
    'plots', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'summary_score_4-metric_lineplot_fine.png'
)

summary_score_fine_batch = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'cleany',
        'heatmap_score_snRNAseq_fine.csv'
    ),
    show_col_types = FALSE
)

summary_score_fine_batch <- summary_score_fine_batch %>%
    pivot_longer(
        cols = c(
            unique_cell_type_count, 
            unique_mhb_lhb_counts, 
            one_to_one_cell_type_count, 
            multiple_mhb_lhb_counts
        ),
        names_to = "type",
        values_to = "value" 
    ) %>%
    select(original_cluster_count, type, value) %>%
    mutate(type = recode(type, 
                         unique_cell_type_count = "Unique Cell Type Count", 
                         unique_mhb_lhb_counts = "Unique MHB LHB Counts", 
                         one_to_one_cell_type_count = "One to One Cell Type Count", 
                         multiple_mhb_lhb_counts = "Multiple MHB LHB Counts")) %>%
    mutate(type2="Cleany")

summary_score_fine_libsize = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'normalized',
        'heatmap_score_snRNAseq_fine.csv'
    ),
    show_col_types = FALSE
)

summary_score_fine_libsize <- summary_score_fine_libsize %>%
    pivot_longer(
        cols = c(
            unique_cell_type_count, 
            unique_mhb_lhb_counts, 
            one_to_one_cell_type_count, 
            multiple_mhb_lhb_counts
        ),
        names_to = "type",
        values_to = "value" 
    ) %>%
    select(original_cluster_count, type, value) %>%
    mutate(type = recode(type, 
                         unique_cell_type_count = "Unique Cell Type Count", 
                         unique_mhb_lhb_counts = "Unique MHB LHB Counts", 
                         one_to_one_cell_type_count = "One to One Cell Type Count", 
                         multiple_mhb_lhb_counts = "Multiple MHB LHB Counts")) %>%
    mutate(type2="Normalized")

summary_score_fine <- rbind(summary_score_fine_batch, summary_score_fine_libsize)

summary_score_fine <- summary_score_fine %>% mutate(type=factor(type, 
                         levels = c("Unique Cell Type Count", 
                                    "One to One Cell Type Count",
                                    "Multiple MHB LHB Counts",
                                    "Unique MHB LHB Counts")))

# Make the barplot figure for fine
p2 = ggplot(summary_score_fine, aes(x = original_cluster_count, y = value, color = type2)) +
  geom_line() +
  geom_point(size=1) +
  labs(x = "K",
       y = "N of cell type",
       color = "Score Type",
       title = "With fine RNAseq reference") +
    scale_color_manual(values = c("Cleany" = "#1f77b4", "Normalized" = "#ff7f0e")) +
    facet_wrap(~ type , scales = "free") +
    theme_classic()

png(plot_path, width = 10, height = 10, units = "in", res = 300)
print(p2)
dev.off()
