# Make the paired test plot
# 4 matrics

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)


plot_path = here(
    'plots', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
    'summary_score_4-metric_paired_plot.png'
)

summary_score_fine_batch = read_csv(
    here(
        'processed-data', '10_HD_bin_level', 'probe_fix', 'ficture_harmony',
        'registration','sum_score', 'cleany_normalized',
        'heatmap_score_snRNAseq_fine_cor0.37.csv'
    ),
    show_col_types = FALSE
)

summary_score_fine <- summary_score_fine_batch %>%
    pivot_longer(
        cols = c(
            unique_cell_type_count_cleany, 
            unique_mhb_lhb_counts_cleany,
            one_to_one_cell_type_count_cleany,
            multiple_mhb_lhb_counts_cleany,
            unique_cell_type_count_normalized,
            unique_mhb_lhb_counts_normalized,
            one_to_one_cell_type_count_normalized,
            multiple_mhb_lhb_counts_normalized
        ),
        names_to = c("metric", "type"),
    names_pattern = "(.*)_(cleany|normalized)",
    values_to = "value"
                ) %>%
                pivot_wider(
    names_from = type,
    values_from = value) %>%
    mutate(type = recode(metric, 
                         unique_cell_type_count = "Unique Cell Type Count", 
                         unique_mhb_lhb_counts = "Unique MHB LHB Counts", 
                         one_to_one_cell_type_count = "One to One Cell Type Count", 
                         multiple_mhb_lhb_counts = "Multiple MHB LHB Counts")) %>%
    mutate(original_cluster_count = as.numeric(original_cluster_count)) %>%
    mutate(type = factor(type, 
                         levels = c("Unique Cell Type Count", 
                                    "One to One Cell Type Count",
                                    "Multiple MHB LHB Counts",
                                    "Unique MHB LHB Counts"))) %>%
    select(original_cluster_count, type, cleany, normalized)
    

# Make the barplot figure for fine
summary_score_fine<-as.data.frame(summary_score_fine)
p1 = ggpaired(
  summary_score_fine,
  "cleany",
  "normalized",
  x = NULL,
  y = NULL,
  id = NULL,
  fill = "condition",
  palette = c("#1f77b4", "#ff7f0e")
) + facet_wrap(~ type , scales = "free")


png(plot_path, width = 10, height = 10, units = "in", res = 300)
print(p1)
dev.off()


t.test_results <- summary_score_fine %>%
  group_by(type) %>%
  summarise(
    t_statistic = t.test(cleany, normalized, paired = TRUE)$statistic,
    p_value = t.test(cleany, normalized, paired = TRUE)$p.value
  )

#   type                       t_statistic p_value
#   <fct>                            <dbl>   <dbl>
# 1 Unique Cell Type Count            1.49  0.145 
# 2 One to One Cell Type Count        1.54  0.132 
# 3 Multiple MHB LHB Counts           1.99  0.0535
# 4 Unique MHB LHB Counts             1.28  0.208
