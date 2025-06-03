#   Read in Jaccard plots for all k and plot as a multi-page PDF

library(here)
library(sessioninfo)

plot_path = here(
    'plots', '10_HD_bin_level', 'probe_fix', 'cell_environment', 'jaccard.pdf'
)
in_paths = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'temp_jaccard', 'k_%d.rds'
)
k_values = c(seq(4, 40, 2), 70, 100)

plot_list = lapply(k_values, function(k) readRDS(sprintf(in_paths, k)))
pdf(plot_path, width = 9)
print(plot_list)
dev.off()

session_info()
