# Make the plot for the 4 metrics for registration figure
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)

 # Broad plot
plot_path = here(
    'plots', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'summary_score_4-metric_barplot.png'
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
    mutate(original_cluster_count = as.factor(original_cluster_count), type2="Cleany") %>%
   group_by(type) %>%
   slice_max(order_by = value, n = 3, with_ties = FALSE) %>%
   ungroup()

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
    mutate(original_cluster_count = as.factor(original_cluster_count), type2="Normalized") %>%
    group_by(type) %>%
   slice_max(order_by = value, n = 3, with_ties = FALSE) %>%
   ungroup()

summary_score_broad <- rbind(summary_score_broad_batch, summary_score_broad_libsize)

 # Fine plot
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
    mutate(original_cluster_count = as.factor(original_cluster_count), type2="Cleany") %>%
   group_by(type) %>%
   slice_max(order_by = value, n = 3, with_ties = FALSE) %>%
   ungroup()

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
    mutate(original_cluster_count = as.factor(original_cluster_count), type2="Normalized")%>%
   group_by(type) %>%
   slice_max(order_by = value, n = 3, with_ties = FALSE) %>%
   ungroup()

summary_score_fine <- rbind(summary_score_fine_batch, summary_score_fine_libsize)

# Broad plot
summary_score_top10 <- summary_score_broad %>%
    mutate(type=factor(type, 
                         levels = c("Unique Cell Type Count", 
                                    "One to One Cell Type Count",
                                    "Multiple MHB LHB Counts",
                                    "Unique MHB LHB Counts"))) %>%
  group_by(type, type2) %>%
  mutate(bar_id = paste(type2, row_number(), sep = "_")) %>%
  ungroup()

# Make the barplot figure for broad
p1 <- ggplot(summary_score_top10, aes(x = bar_id, y = value, fill = type2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = original_cluster_count,y= 0), 
            vjust = 1.5, size = 3.5, color = "black") +
  facet_grid(. ~ type, scales = "free_x", space = "free_x") + 
  labs(title = "With broad RNAseq reference",
       x = NULL,
       y = "N of cell type",
       fill = "Score Type") +
  scale_fill_manual(values = c("Cleany" = "#1f77b4", "Normalized" = "#ff7f0e")) +
    theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    panel.spacing = unit(1.5, "lines")
  ) + scale_y_continuous(limits = c(0, 9),breaks = seq(0, 9, 2)) 

# Fine plot 

summary_score_top10 <- summary_score_fine %>%
    mutate(type=factor(type, 
                         levels = c("Unique Cell Type Count", 
                                    "One to One Cell Type Count",
                                    "Multiple MHB LHB Counts",
                                    "Unique MHB LHB Counts"))) %>%
  group_by(type, type2) %>%
  mutate(bar_id = paste(type2, row_number(), sep = "_")) %>%
  ungroup()

# Make the barplot figure for broad
p2 <- ggplot(summary_score_top10, aes(x = bar_id, y = value, fill = type2)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = original_cluster_count,y= 0), 
            vjust = 1.5, size = 3.5, color = "black") +
  facet_grid(. ~ type, scales = "free_x", space = "free_x") + 
  labs(title = "With fine RNAseq reference",
       x = NULL,
       y = "N of cell type",
       fill = "Score Type") +
  scale_fill_manual(values = c("Cleany" = "#1f77b4", "Normalized" = "#ff7f0e")) +
    theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    panel.spacing = unit(1.5, "lines")
  ) + scale_y_continuous(limits = c(0, 17),breaks = seq(0, 17, 2)) 

# Save the plot
# p = ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
png(plot_path, width = 9, height = 5, units = "in", res = 300)
print(p2)
dev.off()


