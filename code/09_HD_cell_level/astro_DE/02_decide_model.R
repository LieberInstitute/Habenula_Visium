library(tidyverse)
library(here)
library(sessioninfo)
library(SpatialExperiment)
library(variancePartition)
library(pheatmap)
library(qs2)

sce_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'astro_sce_pb.qs2'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'no_secondary', 'astro_DE',
    'decide_model'
)
out_dir = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'astro_DE'
)
cont_covariates = c('sum_umi', 'ncells', 'expr_chrM_ratio')

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

vp_wrapper = function(sce, pd, this_formula_vp, plot_path, out_path) {
    vp_obj = fitExtractVarPartModel(logcounts(sce), this_formula_vp, pd) |>
        sortCols()

    p = plotVarPart(vp_obj,) +
        theme_bw(base_size = 15) +
        theme(
            axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90),
            legend.position = 'none'
        )
    pdf(plot_path, width = 6, height = 8)
    print(p)
    dev.off()

    #   Also write the data to CSV, since it took a while to compute
    apply(vp_obj, 2, summary) |>
        as.data.frame() |>
        rownames_to_column("metric") |>
        pivot_longer(!metric) |>
        write_csv(out_path)
}

################################################################################
#   Load data and prepare covariates and models
################################################################################

sce = qs_read(sce_path)
pd = as.data.frame(colData(sce))

#   Will use centered and scaled versions of continuous covariates in the model
for (this_covariate in cont_covariates) {
    pd[[paste(this_covariate, 'scaled', sep = '_')]] = as.numeric(
        scale(pd[[this_covariate]])
    )
}
cont_covariates_scaled = paste(cont_covariates, 'scaled', sep = '_')

#   Appropriate formulas for canCorPairs and fitExtractVarPartModel,
#   respectively
this_formula = as.formula(
    paste(
        '~ astro_label + sample_id + donor +',
        paste0(cont_covariates_scaled, collapse = ' + ')
    )
)

this_formula_vp = list()
this_formula_vp[['before']] = as.formula(
    paste(
        '~ (1 | astro_label) + (1 | sample_id) + (1 | donor) +',
        paste0(cont_covariates_scaled, collapse = ' + ')
    )
)
this_formula_vp[['after']] = as.formula(
    '~ (1 | astro_label) + ncells_scaled + expr_chrM_ratio_scaled'
)

################################################################################
#   Canonical correlation analysis
################################################################################

p = canCorPairs(this_formula, pd) |>
    pheatmap(
        color = hcl.colors(50, "YlOrRd", rev = TRUE),
        fontsize = 18,
        border_color = "black",
        angle_col = 90,
        display_numbers = TRUE
    )
pdf(file.path(plot_dir, 'CCA_heatmap.pdf'))
print(p)
dev.off()

################################################################################
#   Variance-explained plots
################################################################################

#   'before' is trying out all potential covariates; 'after' is the proposed
#   final model for DE
for (model_time in c('before', 'after')) {
    vp_wrapper(
        sce, pd, this_formula_vp[[model_time]],
        file.path(plot_dir, sprintf('variance_explained_%s.pdf', model_time)),
        file.path(
            out_dir, sprintf('variance_explained_summary_%s.csv', model_time)
        )
    )
}

session_info()
