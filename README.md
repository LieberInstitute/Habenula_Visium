
# Habenula_Visium

# Internal

JHPCE location: `/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium`

<br>

### Processed-data Directory

<br><br>

Script directory: *code/02_build_spe/*
<br>

| File name | Description     |
|:------------|:----------------|
| spe_raw.rds  | Basic *SpatialExperiment* (spe) object with donor information        |
| spe.rds  | Pre filtered *spe* object with spots in-tissue and removed spots with no counts       |
| spe_qc.rds |  Pre filtered *spe* object with the scran additional QC metrics       |
| spe_qc_low_lib_edge.rds | *spe* object removed low-library spots at edge. This will be used with the *spatialLIBD* app  |
 <!--| spe_qc_low_lib_edge_HighM.rds  | Filtered *spe.rds* object removed *low library at edge* and *high chrM percent* by spot       | -->


<br><br>

### code/03_spatialLIBD_app/ Directory

Notes: *spatialLIBD* version 1.15.4 is available on bioc-devel, not on bioc-release.
Aka via https://bioconductor.org/packages/devel/data/experiment/html/spatialLIBD.html + the installation instructions listed there.
<br>


1. The new *spatialLIBD* new version app start R (version “4.4”). Thus, if you are working in the *JHPCE* you we need to update your *bashrc* profile for JHPCE to load the proper R module. Ex. $ vim ~/.bashrc; Replace conda_R/4.3.x with conda_R/4.4

2. Working directory need to be on the current directory. 

3. Before running the *app.R* script, you need to create a symbolic link in the current directory pointing to the *spe.R* object to wrap in the *spatialLIBD* app.


<br><br>


