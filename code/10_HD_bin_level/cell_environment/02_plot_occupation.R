library(here)
library(tidyverse)
library(segmented)
library(sessioninfo)

occupation_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'occupation', '%s.csv'
)
sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'cell_environment')

sample_ids = readLines(sample_id_path)

occupation_df_list = list()
for (sample_id in sample_ids) {
    occupation_df_list[[sample_id]] = sprintf(occupation_path, sample_id) |>
        read_csv(show_col_types = FALSE) |>
        mutate(sample_id = sample_id)
}
occupation_df = do.call(rbind, occupation_df_list)

#   Occupation fractions by sample
p = ggplot(
        occupation_df,
        aes(
            x = expansion_distance, y = occupation, color = sample_id,
            group = sample_id
        )
    ) +
    geom_point() +
    geom_line() +
    labs(
        x = 'Expansion Distance', y = 'Occupation Fraction',
        color = 'Sample ID'
    ) +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'occupation_individual.pdf'), width = 10)
print(p)
dev.off()

occupation_df = occupation_df |>
    group_by(expansion_distance) |>
    summarize(occupation = mean(occupation)) |>
    ungroup()

#   Mean occupation fractions
p = ggplot(occupation_df, aes(x = expansion_distance, y = occupation)) +
    geom_point() +
    geom_line() +
    labs(x = 'Expansion Distance', y = 'Occupation Fraction') +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'occupation_mean.pdf'))
print(p)
dev.off()

model = lm(occupation ~ expansion_distance, data = occupation_df)

message('Segmented model summary (occupation fraction vs. expansion distance):')
model |>
    segmented(
        seg.Z = ~expansion_distance, psi = list(expansion_distance = 6)
    ) |>
    summary()

session_info()
