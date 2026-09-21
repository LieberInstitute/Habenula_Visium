# Habenula_Visium

Welcome to the Hb_multiome project! Here you will find all code used to analyze the data generated as part of the manuscript [Manuscript name].  

[Project Website](https://research.libd.org/Hb_multiome/)

There are 2 GitHub repos associated with this study.
They are: 
1. [Hb_Multiome](https://github.com/LieberInstitute/Hb_multiome)
2. [Habenula_Visium](https://github.com/LieberInstitute/Habenula_Visium)


# Internal

JHPCE location: `/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium`

<br>

### Public available resources

<br><br>

| File name | Description     |
|:------------|:----------------|
| https://libd.shinyapps.io/Habenula_Visium_raw/  | Shiny app version with 12 Habenula samples, featuring both in-tissue and out-of-tissue spots with quality controls. |
|  https://libd.shinyapps.io/Habenula_Visium_Sp09/  | Shiny app version with 12 Habenula samples, featuring only in-tissue spots after standard quality control, manual annotations and Layer-level data for pseudobulk BayesSpace k=9  |


<br>

### JHPCE paths

<br><br>

| Description                              | JHPCE Directory                                                                 |
|:-----------------------------------------|:---------------------------------------------------------------------------|
| snRNAseq enrichment t-stat objects | <span style="color:red">`~/processed-data/05_snRNA-seq_model_stats/`</span> |
| Multiome enrichment t-stat objects | <span style="color:red">`~/processed-data/05_snRNA-seq_model_stats/` (Symbolic link: Hb_multiome/processed-data/08_spatial_registration_vs_multiome_snRNA-seq/ </span> |
| Visium Bayes-Space modeling-results | <span style="color:red">`~/processed-data/05_brain_area_differential_expression/modeling_results_BS/`</span> |

<br>


| Description                              | JHPCE Directory                                                                 |
|:-----------------------------------------|:---------------------------------------------------------------------------|
| Slurm to build references | <span style="color:red">`~/code/run_build_rna_references.sh`</span> |
| Slurm to compute Spatial-Registration        | <span style="color:red">`~/code/run_spatial_registration.sh`</span>  |



