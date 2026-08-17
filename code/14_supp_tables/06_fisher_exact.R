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
library(Signac)
library(qs2)
library(rtracklayer)
library(orthogene)
library(AnnotationDbi)
library(org.Hs.eg.db)

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

#-------------------------------------------------------------------------------
#   Other
#-------------------------------------------------------------------------------

gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'

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

#   For Fisher's exact test, we want ENSEMBL IDs from the GTF. If we naively
#   join with the GTF by gene symbol, 9 genes fail to join. This is due to an
#   upstream issue (that really should've been fixed earlier) where because
#   non-unique gene symbols are used for the rownames, make.unique()
#   automatically added ".1" and similar suffixes.
#   Returns a named character vector of Ensembl IDs, one per row of seur[["RNA"]].
#   The RNA assay uses gene symbols as rownames; in a few cases make.unique() has
#   appended ".1" etc. to duplicate symbols, creating names not present in the GTF.
#   Legitimate gene names that end in .[0-9]+ (e.g. "AL627309.1") ARE directly in
#   the GTF and are left as-is. For true make.unique() artifacts the suffix is
#   stripped to recover the base symbol, and the nth occurrence of that symbol in
#   the rownames is mapped to the nth Ensembl ID for that symbol in the GTF.
grab_multiome_genes = function(seur, gtf_path) {
    gtf = import(gtf_path) |>
        as.data.frame() |>
        as_tibble() |>
        filter(type == "gene") |>
        select(gene_name, gene_id) |>
        distinct()

    gtf_symbols = unique(gtf$gene_name)

    #   Named list: gene symbol -> character vector of Ensembl IDs in GTF order
    ensembls_by_symbol = gtf |>
        group_by(gene_name) |>
        summarise(gene_ids = list(gene_id), .groups = "drop") |>
        deframe()

    original_rownames    = rownames(seur[["RNA"]])
    times_symbol_seen    = integer(0)
    resolved_ensembl_ids = character(length(original_rownames))

    for (i in seq_along(original_rownames)) {
        rn_i = original_rownames[i]
        #   Legitimate names like "AL627309.1" are directly in the GTF.
        #   make.unique() artifacts like "TBCE.1" are not; strip the suffix.
        base_symbol = if (rn_i %in% gtf_symbols) rn_i else sub("\\.[0-9]+$", "", rn_i)
        stopifnot(base_symbol %in% gtf_symbols)

        n_prior = times_symbol_seen[base_symbol]
        if (is.na(n_prior)) n_prior = 0L
        ens = ensembls_by_symbol[[base_symbol]]
        stopifnot(n_prior + 1L <= length(ens))

        resolved_ensembl_ids[i]        = ens[n_prior + 1L]
        times_symbol_seen[base_symbol] = n_prior + 1L
    }

    names(resolved_ensembl_ids) = original_rownames
    resolved_ensembl_ids
}

#   Given rat Ensembl IDs, find orthologs where possible, then convert from
#   human gene symbol to Ensembl ID. This is fragile in a couple ways, but
#   should be sufficient for our purposes
rat_to_human = function(gene_df) {
    mapped <- convert_orthologs(
        gene_df        = data.frame(row.names = unique(gene_df$gene_id)),
        gene_input     = "rownames",
        gene_output    = "dict",
        input_species  = "rat",
        output_species = "human",
        non121_strategy = "keep_popular",
        method         = "gprofiler"
    )
    
    map_df = tibble(rat_id = names(mapped), human_symbol = unname(mapped))
    human_ens <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = map_df$human_symbol,
        keytype = "SYMBOL",
        column = "ENSEMBL",
        multiVals = "first"
    )
    map_df = tibble(
            human_symbol = names(human_ens), human_id = unname(human_ens)
        ) |>
        distinct(human_symbol, .keep_all = TRUE) |>
        left_join(map_df, by = "human_symbol") |>
        dplyr::select(rat_id, human_id) |>
        filter(!is.na(human_id))

    gene_df = gene_df |>
        inner_join(map_df, by = c("gene_id" = "rat_id")) |>
        dplyr::select(cluster, human_id, is_deg) |>
        dplyr::rename(gene_id = human_id)

    return(gene_df)
}

rat_to_human2 = function(rat_ids) {
    library(biomaRt)

    rat <- useEnsembl(
        biomart = "genes",
        dataset = "rnorvegicus_gene_ensembl"
    )

    human <- useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl"
    )

    orthologs <- getLDS(
        attributes  = "ensembl_gene_id",
        filters     = "ensembl_gene_id",
        values      = rat_ids,
        mart        = rat,

        attributesL = c("ensembl_gene_id", "external_gene_name"),
        martL       = human
    )

    return(orthologs)
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

#   Unexpressed genes shouldn't be part of the universe for Fisher's exact test
keep_genes <- rownames(seur[['RNA']])[
    rowSums(GetAssayData(seur, assay = "RNA", layer = "counts") > 0) > 0
]
seur[["RNA"]] <- subset(seur[["RNA"]], features = keep_genes)

multiome_ensembl_ids = grab_multiome_genes(seur, gtf_path)
multiome_gene_df = prep_marker_df(
    multiome_path, multiome_ensembl_ids, multiome_map_df
)
rm(seur); gc()

#...............................................................................
#   Fentanyl habenula DEGs
#...............................................................................

fent_hb_gene_df = get(load(fent_hb_deg_path))[[1]] |>
    as_tibble() |>
    dplyr::rename(gene_id = ensemblID) |>
    mutate(is_deg = adj.P.Val < 0.05, cluster = 'global') |>
    dplyr::select(cluster, gene_id, is_deg) |>
    rat_to_human()

#...............................................................................
#   Fentanyl amygdala DEGs
#...............................................................................

fent_amyg_gene_df = get(load(fent_amyg_deg_path))[[1]] |>
    as_tibble() |>
    dplyr::rename(gene_id = ensemblID) |>
    mutate(is_deg = adj.P.Val < 0.05, cluster = 'global') |>
    dplyr::select(cluster, gene_id, is_deg) |>
    rat_to_human()
