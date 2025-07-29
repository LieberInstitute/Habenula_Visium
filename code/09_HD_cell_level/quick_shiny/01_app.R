library(spatialLIBD)
library(HDF5Array)
library(markdown)
library(here)

spe_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'spe_norm_filtered'
)
docs_dir = here('code', '09_HD_cell_level', 'quick_shiny', 'www')
sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Load the object and subset to the samples with good H&E images
spe <- loadHDF5SummarizedExperiment(spe_dir)
spe = spe[, spe$sample_id %in% readLines(sample_id_path)[1:3]]

## Deploy the website
run_app(
    spe,
    sce_layer = NULL,
    modeling_results = NULL,
    sig_genes = NULL,
    title = "habenula_atlas_HD",
    spe_discrete_vars = c(
        "ManualAnnotation",
        "labels_joint_source"
    ),
    spe_continuous_vars = c(
        "sum_umi",
        "sum_gene",
        "expr_chrM",
        "expr_chrM_ratio"
    ),
    default_cluster = "labels_joint_source",
    docs_path = docs_dir,
    is_stitched = TRUE
)
