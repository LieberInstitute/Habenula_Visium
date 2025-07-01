library(here)
library(tidyverse)
library(crawdad)
library(sessioninfo)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
plot_path = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'crawdad', 'dot_plot_%s.pdf'
)
result_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad',
    '%s_results.csv'
)

sample_ids = readLines(sample_id_path)[1:3]

result_list = list()
for (sample_id in sample_ids) {
    result_list[[sample_id]] = sprintf(result_paths, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        mutate(sample_id = sample_id)

    result_list[[sample_id]]$Z_sig = correctZBonferroni(
        result_list[[sample_id]]
    )
}

do.call(rbind, result_list) |>
    #   First average Z-scores across permutations
    group_by(sample_id, neighbor, scale, reference, Z_sig) |>
    summarize(Z = mean(Z))
