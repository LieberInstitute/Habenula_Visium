library("here")
library("SummarizedExperiment")
library("sessioninfo")

## Create output directories
dir_plots <- here("plots", "05_snRNA-seq_model_stats", "BayesSpace")
# dir_rdata <- here("processed-data", "05_snRNA-seq_model_stats")
# dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Load Habenula snRNA-seq data from the
## https://github.com/LieberInstitute/Habenula_Pilot
## project.

## If we were inside that project
## https://github.com/LieberInstitute/Habenula_Pilot/blob/85bfdf36505c829523358a98328a3e2b22bf2042/code/99_paper_figs/12_longer_Heatmap.R#L17-L18
##

habenula_pilot_path <- "~/Dropbox/Code/Habenula_Pilot"
if (!file.exists(habenula_pilot_path)) {
    habenula_pilot_path <-
        "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula"
}
input_sce_path <-
    file.path(
        habenula_pilot_path,
        "processed-data",
        "04_snRNA-seq",
        "sce_objects",
        "sce_final.Rdata"
    )
load(input_sce_path, verbose = TRUE)

sce_pseudo <-
    registration_pseudobulk(
        sce,
        var_registration = var_registration,
        var_sample_id = var_sample_id,
        min_ncells = min_ncells,
        pseudobulk_rds_file = pseudobulk_rds_file
    )

registration_mod <-
    registration_model(sce_pseudo, covars = covars)

block_cor <-
    registration_block_cor(sce_pseudo, registration_model = registration_mod)

results_enrichment <-
    registration_stats_enrichment(
        sce_pseudo,
        block_cor = block_cor,
        covars = covars,
        gene_ensembl = gene_ensembl,
        gene_name = gene_name
    )
