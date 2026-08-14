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
    'processed-data', '12_apps_and_sharing', '01_prep_objects', 'spe_shiny.qs2'
)
hd_extra_spe_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'extracellular', 'spe_norm_filtered.qs2'
)
multiome_seur_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/11_link_prep/02_rebuild_atac_assay/cell_level_seur.qs2'

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

################################################################################
#   Functions
################################################################################

prep_marker_df = function(cell_path, all_genes, map_df) {
    hd_marker_df = read_tsv(cell_path, show_col_types = FALSE) |>
        left_join(map_df, by = c('set_id' = 'old_cell_type_dot')) |>
        dplyr::rename(cluster = new_cell_type) |>
        select(cluster, gene_id) |>
        mutate(is_marker = TRUE)

    gene_df = tibble(
            cluster = rep(
                unique(hd_marker_df$cluster), each = length(all_genes)
            ),
            gene_id = rep(
                all_genes, times = length(unique(hd_marker_df$cluster))
            )
        ) |>
        left_join(hd_marker_df, by = c('cluster', 'gene_id')) |>
        mutate(is_marker = coalesce(is_marker, FALSE))

    return(gene_df)
}

################################################################################
#   Main
################################################################################

#-------------------------------------------------------------------------------
#   For each dataset, prep a tibble for Fisher's exact test
#-------------------------------------------------------------------------------

#   For each dataset we're essentially forming a tibble with columns
#   cluster, gene_id, (is_marker or is_deg), then later testing all combinations
#   of reference and query clusters using the intersection of gene_id

hd_map_df = read_csv(hd_map_path, show_col_types = FALSE) |>
    mutate(old_cell_type_dot = str_replace(old_cell_type, '/', '.'))

#...............................................................................
#   HD cell types
#...............................................................................

spe_hd_cell = qs_read(hd_cell_spe_path)
hd_gene_df = prep_marker_df(
    hd_cell_path, rowData(spe_hd_cell)$gene_id, hd_map_df
)
rm(spe_hd_cell); gc()

#...............................................................................
#   HD extracellular cell types
#...............................................................................

spe_hd_extra = qs_read(hd_extra_spe_path)
spe_hd_extra$cell_type = tibble(
        old_cell_type = as.character(spe_hd_extra$cell_type)
    ) |>
    left_join(hd_map_df, by = c('old_cell_type')) |>
    pull(new_cell_type)

hd_extra_gene_df = prep_marker_df(
    hd_extra_path, rowData(spe_hd_extra)$gene_id, hd_map_df
)
rm(spe_hd_extra); gc()

#...............................................................................
#   Multiome
#...............................................................................

multiome_map_df = read_csv(multiome_map_path, show_col_types = FALSE) |>
    dplyr::rename(old_cell_type_dot = old_cell_type)

seur = qs_read(multiome_seur_path)
#   TODO: Ensembl IDs are not in the object (just symbols. Convert)
multiome_gene_df = prep_marker_df(multiome_path, genes_here, multiome_map_df)
rm(seur); gc()
