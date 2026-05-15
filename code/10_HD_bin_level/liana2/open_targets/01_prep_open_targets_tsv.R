library(here)
library(tidyverse)

in_paths = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana2',
    'open_targets', '%s.tsv'
)
traits = c('MDD', 'substance')
risk_cols = c(
    'gwasCredibleSets', 'geneBurden', 'eva', 'genomicsEngland',
    'gene2Phenotype', 'uniprotLiterature', 'uniprotVariants',
    'orphanet', 'clingen'
)
plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'liana2',
    'open_targets'
)
cutoff_val = 0.1 # recommended for trait association

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

risk_df_list = list()
for (this_trait in traits) {
    risk_df_list[[this_trait]] = read_tsv(
            sprintf(in_paths, this_trait), show_col_types = FALSE
        ) |>
        select(symbol, all_of(risk_cols)) |>
        mutate(
            across(
                all_of(risk_cols), ~ as.numeric(if_else(. == 'No data', '0', .))
            ),
            max_risk = pmax(!!!syms(risk_cols)),
            trait = this_trait
        ) |>
        select(trait, symbol, max_risk)
}
risk_df = bind_rows(risk_df_list)

p = ggplot(risk_df, aes(x = max_risk)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(
        xintercept = cutoff_val, color = "red", linetype = "dashed", linewidth = 1
    ) +
    facet_wrap(~trait, nrow = 2, scales = 'free_y') +
    labs(
        title = "Distribution of Maximum Risk Scores by Trait",
        x = "Maximum Risk Score",
        y = "Number of Genes"
    )
pdf(file.path(plot_dir, 'max_risk_histogram.pdf'))
print(p)
dev.off()
