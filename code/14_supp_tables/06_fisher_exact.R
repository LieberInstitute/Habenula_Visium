#   Fisher's exact test: checking enrichment of our markers among reference DEGs
#   (or the reverse; we actually only use the p-value here, which should be
#   symmetric)
#
#   Markers include HD cell types, HD extracellular, and Multiome RNA
#   DEGs include fentanyl rat (habenula and amygdala separately) and habenula
#       pilot
#   Gene universe is the intersection of genes input to the DE and those after
#       basic QC in our data

library(here)
library(tidyverse)
library(sessioninfo)
library(SpatialExperiment)
library(Seurat)
library(qs2)

################################################################################
#   Paths
################################################################################

#-------------------------------------------------------------------------------
#   Markers
#-------------------------------------------------------------------------------

hd_cell_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA', 'gene_sets',
    'fine.tsv'
)
hd_extra_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'gene_sets', 'fine.tsv'
)
multiome_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/RNA/gene_sets/fine.tsv'

#-------------------------------------------------------------------------------
#   Cell-type mappings
#-------------------------------------------------------------------------------

multiome_map_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/cell_type_map.csv'
hd_map_path = here('raw-data', 'sample_info', 'hd_cell_type_map.csv')

#-------------------------------------------------------------------------------
#   Main SpatialExperiment-like objects (to later form universe genes)
#-------------------------------------------------------------------------------

hd_cell_spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary',
    'spe_norm_filtered_split.rds'
)
hd_extra_spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.qs2'
)
multiome_seur_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/16_shiny_app/02_prep_cell_level/atlas_seur_minimal.qs2'

#-------------------------------------------------------------------------------
#   DEGs
#-------------------------------------------------------------------------------

# From some detective work we can see that the full set of input genes is
# present at this path, unfiltered by significance. Check these snippets:
#    https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/4284ebac39e53fe23aea2a23812dc71441c14e10/code/05_DEA/01_Modeling.R#L81-L89
#    https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/4284ebac39e53fe23aea2a23812dc71441c14e10/code/05_DEA/01_Modeling.R#L103-L104
# Same logic was used for the amygdala as habenula. Proof here:
#    https://github.com/LieberInstitute/fentanyl_rat_hb_amy/blob/4284ebac39e53fe23aea2a23812dc71441c14e10/code/05_DEA/01_Modeling.R#L141-L142
fent_hb_deg_path = '/dcs04/lieber/marmaypag/fentanylRat_LIBD4270/fentanyl_rat_hb_amy/processed-data/05_DEA/results_Substance_all_vars_habenula.Rdata'
fent_amyg_deg_path = '/dcs04/lieber/marmaypag/fentanylRat_LIBD4270/fentanyl_rat_hb_amy/processed-data/05_DEA/results_Substance_all_vars_amygdala.Rdata'

pilot_hb_deg_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/10_DEA/04_DEA/DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
