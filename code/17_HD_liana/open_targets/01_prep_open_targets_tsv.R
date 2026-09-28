library(here)
library(sessioninfo)
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
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana2',
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
    ) +
    theme_bw(base_size = 15)
pdf(file.path(plot_dir, 'max_risk_histogram.pdf'))
print(p)
dev.off()

risk_df |>
    filter(max_risk >= cutoff_val) |>
    dplyr::rename(gene = symbol) |>
    select(trait, gene) |>
    write_csv(file.path(out_path, 'risk_genes.csv'))

###--------------------------------------------
### Export OpenTargets risk genes in Python-compatible format
###--------------------------------------------

# Map short trait names in the R script to disease names used in the Python code
trait_name_map <- c(
    MDD = "MDD",
    substance = "Substance_dependence"
)

# Generate files for each trait
for (this_trait in traits) {

    disease_name <- trait_name_map[[this_trait]]

    this_risk_df <- risk_df |>
        filter(trait == this_trait, max_risk >= cutoff_val) |>
        transmute(
            gene = symbol,
            genetic_association = max_risk
        ) |>
        distinct(gene, .keep_all = TRUE) |>
        arrange(desc(genetic_association))

    # 1) Similar to Python output:
    #    Substance_dependence_risk_genes_01thr.csv
    #    First column is gene, so Python can read it with index_col=0
    write_csv(
        this_risk_df,
        file.path(
            out_path,
            paste0(disease_name, "_risk_genes_01thr.csv")
        )
    )

    # 2) Similar to Python output:
    #    Substance_dependence_risk_genes_01thr_hot1.csv
    #    Required by downstream Python code:
    #    tgts = pd.read_csv(...)
    #    tgts['gene'] = tgts['genes']
    this_hot1_df <- this_risk_df |>
        transmute(
            genes = gene,
            !!disease_name := 1
        )

    write_csv(
        this_hot1_df,
        file.path(
            out_path,
            paste0(disease_name, "_risk_genes_01thr_hot1.csv")
        )
    )
}

session_info()
