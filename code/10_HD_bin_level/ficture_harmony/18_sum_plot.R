# Make the plot for the summary metrics for registration figure
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(ggforce)
library(dplyr)

custom_theme <- theme(
panel.background = element_blank(),
axis.ticks = element_blank(),
axis.text = element_text(size = 12),
axis.title = element_text(size = 14, face = "bold"), axis.line = element_line(linewidth = 0.5) )

# Panel b: Summary score line plot
    # Broad plot
plot_path = here(
    'plots', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
    'summary_score_line_plot.png'
)
summary_score_broad_batch = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'cleany',
        'heatmap_score_snRNAseq_broad.csv'
    ),
    show_col_types = FALSE
)
summary_score_broad_libsize = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'normalized',
        'heatmap_score_snRNAseq_broad.csv'
    ),
    show_col_types = FALSE
)
summary_score_broad <- summary_score_broad_batch %>%
    dplyr::left_join(summary_score_broad_libsize, by = "original_cluster_count", 
                   suffix = c("_cleany", "_normalized")) %>%
    mutate(unique_cell_type_count_cleany=unique_cell_type_count_cleany/9,
           unique_mhb_lhb_counts_cleany=unique_mhb_lhb_counts_cleany/2,
           one_to_one_cell_type_count_cleany=one_to_one_cell_type_count_cleany/9,
           multiple_mhb_lhb_counts_cleany=multiple_mhb_lhb_counts_cleany/2,
           score_total_cleany=unique_cell_type_count_cleany+unique_mhb_lhb_counts_cleany+one_to_one_cell_type_count_cleany+multiple_mhb_lhb_counts_cleany) %>%
    mutate(unique_cell_type_count_normalized=unique_cell_type_count_normalized/9,
           unique_mhb_lhb_counts_normalized=unique_mhb_lhb_counts_normalized/2,
           one_to_one_cell_type_count_normalized=one_to_one_cell_type_count_normalized/9,
           multiple_mhb_lhb_counts_normalized=multiple_mhb_lhb_counts_normalized/2,
           score_total_normalized=unique_cell_type_count_normalized+unique_mhb_lhb_counts_normalized+one_to_one_cell_type_count_normalized+multiple_mhb_lhb_counts_normalized) %>%
    select(original_cluster_count, score_total_cleany, score_total_normalized)


# Read the fine summary score
summary_score_fine_batch = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'cleany',
        'heatmap_score_snRNAseq_fine.csv'
    ),
    show_col_types = FALSE
)
summary_score_fine_libsize = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'new_samples', 'ficture_harmony',
        'registration','sum_score', 'normalized',
        'heatmap_score_snRNAseq_fine.csv'
    ),
    show_col_types = FALSE
)
summary_score_fine <- summary_score_fine_batch %>%
    dplyr::left_join(summary_score_fine_libsize, by = "original_cluster_count", 
                   suffix = c("_cleany", "_normalized")) %>%
    mutate(unique_cell_type_count_cleany=unique_cell_type_count_cleany/17,
           unique_mhb_lhb_counts_cleany=unique_mhb_lhb_counts_cleany/10,
           one_to_one_cell_type_count_cleany=one_to_one_cell_type_count_cleany/17,
           multiple_mhb_lhb_counts_cleany=multiple_mhb_lhb_counts_cleany/10,
           score_total_cleany=unique_cell_type_count_cleany+unique_mhb_lhb_counts_cleany+one_to_one_cell_type_count_cleany+multiple_mhb_lhb_counts_cleany) %>%
    mutate(unique_cell_type_count_normalized=unique_cell_type_count_normalized/17,
           unique_mhb_lhb_counts_normalized=unique_mhb_lhb_counts_normalized/10,
           one_to_one_cell_type_count_normalized=one_to_one_cell_type_count_normalized/17,
           multiple_mhb_lhb_counts_normalized=multiple_mhb_lhb_counts_normalized/10,
           score_total_normalized=unique_cell_type_count_normalized+unique_mhb_lhb_counts_normalized+one_to_one_cell_type_count_normalized+multiple_mhb_lhb_counts_normalized) %>%
    select(original_cluster_count, score_total_cleany, score_total_normalized)

df_long <- summary_score_broad %>%
  pivot_longer(cols = c(score_total_cleany, score_total_normalized),
               names_to = "score_type",
               values_to = "score") %>% na.omit()

# Convert score_type to a factor with specific levels
df_long$score_type <- factor(df_long$score_type, 
                             levels = c("score_total_cleany", "score_total_normalized"),
                             labels = c("Cleany", "Normalized"))

head(df_long[order(-df_long$score),])
highlight_points <- df_long %>%
  filter(original_cluster_count %in% c(27, 29, 38)) %>%
  filter(score_type=="Cleany")

# Create the plot
p1 = ggplot(df_long, aes(x = original_cluster_count, y = score, color = score_type)) +
  geom_line() +
  geom_point(size=1) +
  labs(x = "K",
       y = "4-metric summary score",
       color = "Score Type",
       title = "With broad RNAseq reference") +
    scale_color_manual(values = c("Cleany" = "#1f77b4", "Normalized" = "#ff7f0e")) +
    scale_y_continuous(limits = c(0, 4))+
   geom_text(data = highlight_points,
            aes(label = original_cluster_count),
            vjust = -1, size = 3) +
    theme_classic()

# Panel c: Fine plot

df_long <- summary_score_fine %>%
    pivot_longer(cols = c(score_total_cleany, score_total_normalized),
               names_to = "score_type",
               values_to = "score") %>% na.omit()

# Convert score_type to a factor with specific levels
df_long$score_type <- factor(df_long$score_type, 
                             levels = c("score_total_cleany", "score_total_normalized"),
                             labels = c("Cleany", "Normalized"))

head(df_long[order(-df_long$score),])
highlight_points <- df_long %>%
  filter(original_cluster_count %in% c(29, 70, 40)) %>%
  filter(score_type=="Cleany")
# Create the plot
p2 = ggplot(df_long, aes(x = original_cluster_count, y = score, color = score_type)) +
    geom_line() +
    geom_point(size=1) +
     labs(x = "K",
       y = "4-metric summary score",
       color = "Score Type",
       title = "With fine RNAseq reference") +
    scale_color_manual(values = c("Cleany" = "#1f77b4", "Normalized" = "#ff7f0e")) +
    scale_y_continuous(limits = c(0, 4))+
    geom_text(data = highlight_points,
            aes(label = original_cluster_count),
            vjust = -1, size = 3) +
    theme_classic()

# Save the plots
# Combine the two plots
p = ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
png(plot_path, width = 7, height = 5, units = "in", res = 300)
print(p)
dev.off()






