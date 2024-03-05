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

