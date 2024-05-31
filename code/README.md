
# Habenula_Visium

<br>

### Directory: *processed-data/*

##### SPE object descriptions

<br>

| Sub-directory           | File name                 | Description     |
|:------------            |:------------              |:----------------|
| 02_build_spe            | spe_raw.rds               | Basic *SpatialExperiment* (spe) object with donor information. It includes the n=5 samples, featuring both in-tissue and out-of-tissue spots, before quality control.   |
| 02_build_spe            | spe.rds                   | *spe* with spots in-tissue, and removed spots with not umi counts or not gene counts.   |
| 02_build_spe            | spe_qc.rds                | *spe* with spots in-tissue, and added scran variables before applying quality control.    |
| 04_harmony_BayesSpace   | spe_qcED_spatialLIBD.rds  | *spe* with in-tissue spots, QCed and excluding those manually annotated as low quality or having tissue artifacts.    |
| 04_harmony_BayesSpace   | spe_qcED_spatialLIBD_log.rds          | previous *spe* with log normalized counts.    |
| 04_harmony_BayesSpace   | spe_qcED_spatialLIBD_log_GLM-PCA.rds  | previous *spe* computed GLM-PCA.    |
| 04_harmony_BayesSpace   | assays.h5                             | previous *spe* with HDF5-backed object to control memory later.   |

<!--| 02_build_spe  | spe_qc_low_lib_edge.rds | *spe* with only spots in-tissue, and removed low-library spots at edge. This is the version used to deploy thespatialLIBD app: https://libd.shinyapps.io/Habenula_Visium/ |
| 04_harmony_BayesSpace  | spe_qc_low_spatialLIBD.rds | *spe* with only in-tissue spots, excluding those manually annotated as low quality or having artifacts. |-->

<br><br>

##### Other objects

<br>

| Sub-directory | File name       | Description     |
|:------------  |:------------    |:----------------|
| 04_harmony_BayesSpace  | top.hvgs.Rdata | Set of highly variable genes, based on variance modelling statistics from modelGeneVar. This object contains subsets of the top HVGs at FDR equal to 0.05 and 0.01, and proportion equal to 10%, 20% and 50%  |



<br><br>


Notes: 

1. code/03_spatialLIBD_app_prefiltered/ Directory contains the scripts to deploy the Habenula_Visium_raw/ shiny app.

2. code/03_spatialLIBD_app/ Directory contains the scripts to deploy the Habenula_Visium/ shiny app. 


Additional notes:

1. *spatialLIBD* version 1.15.4 is available on bioc-devel, not on bioc-release.
Aka via https://bioconductor.org/packages/devel/data/experiment/html/spatialLIBD.html + the installation instructions listed there. Be sure to have the correct credential access to install the devel version. This link has help resources to setup the PAT access https://usethis.r-lib.org/articles/git-credentials.html

2. The new *spatialLIBD* new version app start R (version “4.4”). Thus, if you are working in the *JHPCE* you we need to update your *bashrc* profile for JHPCE to load the proper R module. Ex. $ vim ~/.bashrc; Replace conda_R/4.3.x with conda_R/4.4

<br><br>


