#!/usr/bin/env Rscript


# =============================================================================
# Pair-specific GO enrichment
#
# FINAL WORKFLOW
#
# For EACH exact directional source -> target pair:
#
#   1. Rename cell types using latest hd_cell_type_map.csv
#
#   2. For every exact:
#          source + target + ligand + receptor
#      require detection in >= 4 samples
#
#   3. For each LR pair, calculate source-target enrichment:
#
#          enrichment =
#              interaction score in this source-target
#              /
#              mean interaction score of the same LR
#              across all reproducible source-target combinations
#
#   4. Require enrichment >= 3
#
#   5. Within each exact source-target pair:
#          rank ALL LR pairs passing enrichment >= 3 by interaction score
#
#   6. KEEP ALL enriched LR pairs (NO Top-N LR restriction)
#
#   7. From ALL enriched LR pairs:
#          extract ligand + receptor genes
#          follow LR rank
#          remove duplicated genes
#          KEEP ALL UNIQUE genes (NO Top-N gene restriction)
#
#   8. GO enrichment:
#          BP / MF / CC
#          custom universe = all genes present in cellular_annotated.h5ad
#          H5AD Ensembl IDs are stripped of version suffixes before mapping
#
#   9. GO is completely pair-specific.
#      Comparison labels are used ONLY for plotting.
#
#   10. GO IDs are retained internally and in CSV,
#       but removed from y-axis labels in the figure.
# =============================================================================


suppressPackageStartupMessages({

    library(data.table)

    library(clusterProfiler)

    library(org.Hs.eg.db)

    library(ggplot2)

    library(here)

    library(zellkonverter)
})


# =============================================================================
# 1. Paths
# =============================================================================

INPUT_FILE <- paste0(

    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",

    "Habenula_Visium/processed-data/10_HD_bin_level/",

    "no_secondary/liana2/table/",

    "significant_interactions_across_donors_5000.0.csv"
)


CELL_TYPE_MAP_FILE <- paste0(

    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",

    "Habenula_Visium/raw-data/sample_info/",

    "hd_cell_type_map.csv"
)


OUTPUT_DIR <- paste0(

    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",

    "Habenula_Visium/processed-data/10_HD_bin_level/",

    "no_secondary/liana2/GO/",

    "pair_specific_enrichment3_allLR_allGenes_h5adUniverse_3"
)


dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# =============================================================================
# 2. Analysis settings
# =============================================================================

MIN_SAMPLES <- 4L

ENRICHMENT_CUTOFF <- 3


FDR_CUTOFF <- 0.05


# -----------------------------------------------------------------------------
# GO background universe
# -----------------------------------------------------------------------------
# Use ALL genes present in the AnnData feature dimension (var) of the same
# cellular_annotated.h5ad dataset as the GO universe (~17,000 genes).
#
# H5AD var_names are Ensembl gene IDs. The code removes any Ensembl version
# suffix (e.g., ENSG00000123456.7 -> ENSG00000123456), converts ENSEMBL ->
# SYMBOL using org.Hs.eg.db, and uses the resulting unique SYMBOLs as universe.

H5AD_FILE <- here(
    "processed-data",
    "10_HD_bin_level",
    "new_samples2",
    "liana",
    "input_habenula",
    "cellular_annotated.h5ad"
)


# plotting only
TOP_N_GO_ACROSS_ONTOLOGIES_PER_PAIR <- 3L


ONTOLOGY_ORDER <- c(
    "BP",
    "MF",
    "CC"
)


# =============================================================================
# 3. Selected exact directional pairs
# =============================================================================

PAIR_INFO <- data.table(

    cell_pair = c(

        # ---------------------------------------------------------------------
        # Astrocyte -> Habenula
        # ---------------------------------------------------------------------

        "Astrocyte -> MHb_A",

        "Astrocyte -> MHb_B",

        "Astrocyte -> Excit_LHb",

        "Astrocyte -> LHb_A",

        "Astrocyte -> LHb_C",


        # ---------------------------------------------------------------------
        # MHb <-> Ependymal / Subependymal
        # ---------------------------------------------------------------------

        "MHb_A -> Ependymal",

        "MHb_A -> Subependymal",

        "MHb_B -> Ependymal",

        "MHb_B -> Subependymal",


        "Ependymal -> MHb_A",

        "Ependymal -> MHb_B",

        "Subependymal -> MHb_A",

        "Subependymal -> MHb_B",


        # ---------------------------------------------------------------------
        # MHb <-> LHb
        # ---------------------------------------------------------------------

        "MHb_A -> Excit_LHb",

        "MHb_A -> LHb_A",

        "MHb_B -> Excit_LHb",

        "MHb_B -> LHb_A",


        "Excit_LHb -> MHb_A",

        "Excit_LHb -> MHb_B",

        "LHb_A -> MHb_A",

        "LHb_A -> MHb_B",


        # ---------------------------------------------------------------------
        # Oligo -> Habenula
        # ---------------------------------------------------------------------

        "Oligo -> MHb_A",

        "Oligo -> MHb_B",

        "Oligo -> Excit_LHb",

        "Oligo -> LHb_A",

        "Oligo -> LHb_C"
    ),


    comparison = c(

        rep(
            "Astrocyte -> Habenula",
            5
        ),

        rep(
            "Ependymal/Subependymal <-> MHb",
            8
        ),

        rep(
            "MHb <-> LHb",
            8
        ),

        rep(
            "Oligo -> Habenula",
            5
        )
    )
)


PAIR_INFO[
    ,
    pair_order := seq_len(.N)
]


stopifnot(
    nrow(PAIR_INFO) == 26L
)


stopifnot(
    !anyDuplicated(
        PAIR_INFO$cell_pair
    )
)


COMPARISON_ORDER <- c(

    "MHb <-> LHb",

    "Ependymal/Subependymal <-> MHb",

    "Astrocyte -> Habenula",

    "Oligo -> Habenula"
)


# =============================================================================
# 4. Read LIANA results
# =============================================================================

dt <- fread(
    INPUT_FILE
)


required_columns <- c(

    "source",

    "ligand",

    "receptor",

    "target",

    "mean",

    "pval",

    "donor_id"
)


missing_columns <- setdiff(
    required_columns,
    names(dt)
)


if (length(missing_columns) > 0L) {

    stop(

        "Missing required columns: ",

        paste(
            missing_columns,
            collapse = ", "
        )
    )
}


dt[
    ,
    `:=`(

        source =
            trimws(
                as.character(source)
            ),

        ligand =
            trimws(
                as.character(ligand)
            ),

        receptor =
            trimws(
                as.character(receptor)
            ),

        target =
            trimws(
                as.character(target)
            ),

        donor_id =
            trimws(
                as.character(donor_id)
            ),

        mean =
            suppressWarnings(
                as.numeric(mean)
            ),

        pval =
            suppressWarnings(
                as.numeric(pval)
            )
    )
]


dt <- dt[

    !is.na(source) &
    source != "" &

    !is.na(target) &
    target != "" &

    !is.na(ligand) &
    ligand != "" &

    !is.na(receptor) &
    receptor != "" &

    !is.na(donor_id) &
    donor_id != "" &

    !is.na(mean) &
    is.finite(mean)
]


cat("\n")
cat("============================================================\n")
cat("INPUT SUMMARY\n")
cat("============================================================\n")


cat(
    "Rows:",
    nrow(dt),
    "\n"
)


cat(
    "Samples:",
    uniqueN(dt$donor_id),
    "\n"
)


cat(
    "Original source cell types:",
    uniqueN(dt$source),
    "\n"
)


cat(
    "Original target cell types:",
    uniqueN(dt$target),
    "\n"
)


# =============================================================================
# 5. Read latest cell-type mapping
# =============================================================================

cell_type_map <- fread(
    CELL_TYPE_MAP_FILE
)


find_first_column <- function(
    column_names,
    candidates
) {

    hit <- candidates[
        candidates %in%
        column_names
    ]

    if (length(hit) == 0L) {

        return(
            NA_character_
        )
    }

    hit[1]
}


old_col <- find_first_column(

    names(cell_type_map),

    c(
        "old_cell_type",
        "old",
        "old_name",
        "cell_type_old",
        "original_cell_type"
    )
)


new_col <- find_first_column(

    names(cell_type_map),

    c(
        "new_cell_type",
        "new",
        "new_name",
        "cell_type_new",
        "mapped_cell_type"
    )
)


if (
    is.na(old_col) ||
    is.na(new_col)
) {

    cat("\nColumns in mapping file:\n")

    print(
        names(cell_type_map)
    )


    stop(
        paste0(
            "\nCould not identify old/new cell-type columns.\n",
            "Expected something like:\n",
            "  old_cell_type\n",
            "  new_cell_type\n"
        )
    )
}


setnames(
    cell_type_map,
    old_col,
    "old_cell_type"
)


setnames(
    cell_type_map,
    new_col,
    "new_cell_type"
)


cell_type_map[
    ,
    `:=`(

        old_cell_type =
            trimws(
                as.character(old_cell_type)
            ),

        new_cell_type =
            trimws(
                as.character(new_cell_type)
            )
    )
]


cell_type_map <- cell_type_map[

    !is.na(old_cell_type) &
    old_cell_type != "" &

    !is.na(new_cell_type) &
    new_cell_type != ""
]


if (
    anyDuplicated(
        cell_type_map$old_cell_type
    )
) {

    cat("\nDuplicated mapping values:\n")

    print(

        cell_type_map[

            duplicated(old_cell_type) |

            duplicated(
                old_cell_type,
                fromLast = TRUE
            )
        ]
    )


    stop(
        "Duplicated old cell types in mapping file."
    )
}


cell_type_lookup <- setNames(

    cell_type_map$new_cell_type,

    cell_type_map$old_cell_type
)


# =============================================================================
# 6. Apply mapping to source and target
# =============================================================================

dt[
    ,
    `:=`(

        source_original =
            source,

        target_original =
            target
    )
]


dt[
    ,
    source_new :=
        unname(
            cell_type_lookup[
                source_original
            ]
        )
]


dt[
    ,
    target_new :=
        unname(
            cell_type_lookup[
                target_original
            ]
        )
]


dt[
    is.na(source_new) |
    source_new == "",
    source_new :=
        source_original
]


dt[
    is.na(target_new) |
    target_new == "",
    target_new :=
        target_original
]


dt[
    ,
    `:=`(

        source =
            source_new,

        target =
            target_new
    )
]


# =============================================================================
# 7. Mapping audit
# =============================================================================

mapping_audit <- unique(

    rbindlist(

        list(

            dt[
                ,
                .(
                    old_cell_type =
                        source_original,

                    new_cell_type =
                        source
                )
            ],

            dt[
                ,
                .(
                    old_cell_type =
                        target_original,

                    new_cell_type =
                        target
                )
            ]
        )
    )
)


setorder(
    mapping_audit,
    old_cell_type
)


cat("\n")
cat("============================================================\n")
cat("CELL TYPE MAPPING\n")
cat("============================================================\n")


print(
    mapping_audit
)


fwrite(

    mapping_audit,

    file.path(
        OUTPUT_DIR,
        "cell_type_mapping_audit.csv"
    )
)


# =============================================================================
# 8. Create source -> target label
#
# DO NOT restrict to selected 26 pairs yet.
# =============================================================================

dt[
    ,
    cell_pair :=
        paste(
            source,
            "->",
            target
        )
]


cat("\n")
cat(
    "Mapped source-target combinations:",
    uniqueN(dt$cell_pair),
    "\n"
)


# =============================================================================
# 9. Collapse duplicate rows within donor
# =============================================================================

donor_lr <- dt[
    ,
    .(
        donor_score =
            mean(
                mean,
                na.rm = TRUE
            )
    ),
    by = .(
        source,
        target,
        cell_pair,
        ligand,
        receptor,
        donor_id
    )
]


# =============================================================================
# 10. Collapse across donors
# =============================================================================

lr_all <- donor_lr[
    ,
    .(

        n_samples =
            uniqueN(
                donor_id
            ),

        interaction_score =
            mean(
                donor_score,
                na.rm = TRUE
            ),

        median_interaction_score =
            median(
                donor_score,
                na.rm = TRUE
            )
    ),
    by = .(
        source,
        target,
        cell_pair,
        ligand,
        receptor
    )
]


cat("\n")
cat("============================================================\n")
cat("LR REPRODUCIBILITY BEFORE FILTER\n")
cat("============================================================\n")


print(
    table(
        lr_all$n_samples
    )
)


cat(
    "Total source-target-specific LR interactions:",
    nrow(lr_all),
    "\n"
)


# =============================================================================
# 11. Require >=4 samples
# =============================================================================

lr_reproducible <- lr_all[
    n_samples >=
        MIN_SAMPLES
]


cat(
    "LR interactions after >= ",
    MIN_SAMPLES,
    " samples: ",
    nrow(lr_reproducible),
    "\n",
    sep = ""
)


if (
    nrow(lr_reproducible) == 0L
) {

    stop(
        "No LR interactions remain after sample reproducibility filtering."
    )
}


# =============================================================================
# 12. Compute LR-specific background score
# =============================================================================

lr_background <- lr_reproducible[
    ,
    .(

        mean_score_all_combos =
            mean(
                interaction_score,
                na.rm = TRUE
            ),

        n_source_target_combos =
            uniqueN(
                cell_pair
            )
    ),
    by = .(
        ligand,
        receptor
    )
]


lr_reproducible <- merge(

    lr_reproducible,

    lr_background,

    by = c(
        "ligand",
        "receptor"
    ),

    all.x = TRUE,

    sort = FALSE
)


# =============================================================================
# 13. Calculate source-target enrichment
# =============================================================================

lr_reproducible[
    ,
    enrichment :=
        interaction_score /
        mean_score_all_combos
]


lr_reproducible[
    ,
    lr_pair :=
        paste(
            ligand,
            receptor,
            sep = " -> "
        )
]


lr_reproducible <- lr_reproducible[

    !is.na(enrichment) &

    is.finite(enrichment) &

    mean_score_all_combos > 0
]


cat("\n")
cat("============================================================\n")
cat("LR ENRICHMENT SUMMARY\n")
cat("============================================================\n")


print(
    summary(
        lr_reproducible$enrichment
    )
)


cat(
    "LR/source-target combinations with enrichment >= ",
    ENRICHMENT_CUTOFF,
    ": ",
    sum(
        lr_reproducible$enrichment >=
            ENRICHMENT_CUTOFF
    ),
    "\n",
    sep = ""
)


fwrite(

    lr_reproducible,

    file.path(
        OUTPUT_DIR,
        "all_reproducible_LR_with_source_target_enrichment.csv"
    )
)


# =============================================================================
# 14. NOW restrict to selected 26 pairs
# =============================================================================

available_pairs <- unique(
    lr_reproducible$cell_pair
)


missing_pairs <- setdiff(

    PAIR_INFO$cell_pair,

    available_pairs
)


if (
    length(missing_pairs) > 0L
) {

    cat("\n")
    cat("Requested pairs missing after >=4-sample filter:\n")

    cat(
        paste(
            missing_pairs,
            collapse = "\n"
        ),
        "\n"
    )
}


selected_lr <- lr_reproducible[
    cell_pair %chin%
        PAIR_INFO$cell_pair
]


selected_lr[
    ,
    comparison :=
        PAIR_INFO$comparison[
            match(
                cell_pair,
                PAIR_INFO$cell_pair
            )
        ]
]


selected_lr[
    ,
    pair_order :=
        PAIR_INFO$pair_order[
            match(
                cell_pair,
                PAIR_INFO$cell_pair
            )
        ]
]


# =============================================================================
# 15. Require >=3-fold enrichment
# =============================================================================

selected_enriched_lr <- selected_lr[

    enrichment >=
        ENRICHMENT_CUTOFF
]


cat("\n")
cat("============================================================\n")
cat(">=3-FOLD ENRICHED LR PAIRS\n")
cat("============================================================\n")


enrichment_summary <- selected_enriched_lr[
    ,
    .(

        n_enriched_LR =
            .N,

        max_enrichment =
            max(
                enrichment,
                na.rm = TRUE
            ),

        median_enrichment =
            median(
                enrichment,
                na.rm = TRUE
            )
    ),
    by = .(
        comparison,
        pair_order,
        cell_pair
    )
]


setorder(
    enrichment_summary,
    pair_order
)


print(
    enrichment_summary
)


fwrite(

    selected_enriched_lr,

    file.path(
        OUTPUT_DIR,
        "selected_pairs_enrichment3_all_LR.csv"
    )
)


# =============================================================================
# 16. Rank ALL enriched LR pairs by interaction score
#
# IMPORTANT:
#   There is NO Top-N LR restriction.
#   Every LR pair with enrichment >= ENRICHMENT_CUTOFF is retained.
# =============================================================================

setorder(
    selected_enriched_lr,
    pair_order,
    -interaction_score,
    -enrichment,
    ligand,
    receptor
)


selected_enriched_lr[
    ,
    LR_rank :=
        seq_len(.N),
    by = cell_pair
]


all_enriched_lr <- copy(
    selected_enriched_lr
)


setorder(
    all_enriched_lr,
    pair_order,
    LR_rank
)


cat("\n")
cat("============================================================\n")
cat("ALL ENRICHED LR PAIRS (NO TOP-N LR FILTER)\n")
cat("============================================================\n")


all_enriched_lr_summary <- all_enriched_lr[
    ,
    .(

        n_LR =
            .N,

        n_unique_ligands =
            uniqueN(
                ligand
            ),

        n_unique_receptors =
            uniqueN(
                receptor
            ),

        n_unique_genes_available =
            length(
                union(
                    ligand,
                    receptor
                )
            ),

        strongest_score =
            max(
                interaction_score,
                na.rm = TRUE
            ),

        weakest_score =
            min(
                interaction_score,
                na.rm = TRUE
            ),

        min_enrichment =
            min(
                enrichment,
                na.rm = TRUE
            ),

        max_enrichment =
            max(
                enrichment,
                na.rm = TRUE
            )
    ),
    by = .(
        comparison,
        pair_order,
        cell_pair
    )
]


setorder(
    all_enriched_lr_summary,
    pair_order
)


print(
    all_enriched_lr_summary
)


fwrite(
    all_enriched_lr,
    file.path(
        OUTPUT_DIR,
        "pair_specific_enrichment3_all_LR.csv"
    )
)


fwrite(
    all_enriched_lr_summary,
    file.path(
        OUTPUT_DIR,
        "pair_specific_enrichment3_all_LR_summary.csv"
    )
)


# =============================================================================
# 17. Read GO universe directly from cellular_annotated.h5ad
#
# IMPORTANT:
#   H5AD var_names are Ensembl IDs.
#   Remove Ensembl version suffix first:
#       ENSG00000123456.7 -> ENSG00000123456
#   Then convert ENSEMBL -> SYMBOL using org.Hs.eg.db.
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("READING GO UNIVERSE FROM H5AD\n")
cat("============================================================\n")


if (
    !file.exists(
        H5AD_FILE
    )
) {

    stop(
        "H5AD file does not exist:\n",
        H5AD_FILE
    )
}


cat(
    "H5AD file:\n",
    H5AD_FILE,
    "\n"
)


# -----------------------------------------------------------------------------
# Read H5AD.
# use_hdf5 = TRUE keeps assay data HDF5-backed where possible.
# We only need the feature names here, not the expression matrix itself.
# -----------------------------------------------------------------------------

sce <- zellkonverter::readH5AD(
    H5AD_FILE,
    use_hdf5 = TRUE
)


cat(
    "H5AD dimensions (genes x observations): ",
    nrow(sce),
    " x ",
    ncol(sce),
    "\n",
    sep = ""
)


if (
    nrow(sce) == 0L
) {

    stop(
        "No genes/features were found in the H5AD file."
    )
}


# -----------------------------------------------------------------------------
# Extract Ensembl IDs from AnnData var_names / rownames(sce)
# -----------------------------------------------------------------------------

background_ensembl_raw <- rownames(
    sce
)


if (
    is.null(
        background_ensembl_raw
    )
) {

    stop(
        "H5AD rownames/var_names are NULL; cannot construct GO universe."
    )
}


background_ensembl_raw <- trimws(
    as.character(
        background_ensembl_raw
    )
)


background_ensembl_raw <- background_ensembl_raw[
    !is.na(background_ensembl_raw) &
        background_ensembl_raw != ""
]


cat(
    "Unique H5AD feature IDs before removing Ensembl versions: ",
    uniqueN(background_ensembl_raw),
    "\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# CRITICAL: remove Ensembl version suffix
#
# Examples:
#   ENSG00000187634.13 -> ENSG00000187634
#   ENSG00000188976.2  -> ENSG00000188976
# -----------------------------------------------------------------------------

background_ensembl <- sub(
    "\\.[0-9]+$",
    "",
    background_ensembl_raw
)


background_ensembl <- unique(
    background_ensembl
)


background_ensembl <- background_ensembl[
    !is.na(background_ensembl) &
        background_ensembl != ""
]


cat(
    "Unique Ensembl IDs after removing version suffix: ",
    length(background_ensembl),
    "\n",
    sep = ""
)


cat("First 20 cleaned Ensembl IDs:\n")
print(
    head(
        background_ensembl,
        20L
    )
)


# -----------------------------------------------------------------------------
# Sanity check: H5AD feature names should now look like Ensembl gene IDs.
# -----------------------------------------------------------------------------

n_ensembl_like <- sum(
    grepl(
        "^ENSG[0-9]+$",
        background_ensembl
    )
)


cat(
    "Features matching ^ENSG[0-9]+$: ",
    n_ensembl_like,
    " / ",
    length(background_ensembl),
    "\n",
    sep = ""
)


if (
    n_ensembl_like == 0L
) {

    stop(
        paste0(
            "No H5AD feature names look like Ensembl gene IDs after ",
            "removing version suffixes. Please inspect rownames(sce)."
        )
    )
}


# -----------------------------------------------------------------------------
# Convert ENSEMBL -> SYMBOL using org.Hs.eg.db
# -----------------------------------------------------------------------------

valid_ensembl_keys <- AnnotationDbi::keys(
    org.Hs.eg.db,
    keytype = "ENSEMBL"
)


ensembl_present_in_orgdb <- intersect(
    background_ensembl,
    valid_ensembl_keys
)


cat(
    "H5AD Ensembl IDs recognized by org.Hs.eg.db: ",
    length(ensembl_present_in_orgdb),
    " / ",
    length(background_ensembl),
    "\n",
    sep = ""
)


if (
    length(ensembl_present_in_orgdb) == 0L
) {

    stop(
        "No cleaned H5AD Ensembl IDs are recognized by org.Hs.eg.db."
    )
}


gene_mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ensembl_present_in_orgdb,
    keytype = "ENSEMBL",
    columns = c(
        "ENSEMBL",
        "SYMBOL"
    )
)


gene_mapping <- as.data.table(
    gene_mapping
)


gene_mapping <- gene_mapping[
    !is.na(ENSEMBL) &
        ENSEMBL != "" &
        !is.na(SYMBOL) &
        SYMBOL != ""
]


# One Ensembl ID can occasionally have more than one annotation row.
# Keep unique ENSEMBL-SYMBOL mappings, then use unique SYMBOLs as the universe.
gene_mapping <- unique(
    gene_mapping,
    by = c(
        "ENSEMBL",
        "SYMBOL"
    )
)


background_genes <- sort(
    unique(
        gene_mapping$SYMBOL
    )
)


mapped_ensembl <- unique(
    gene_mapping$ENSEMBL
)


unmapped_ensembl <- setdiff(
    background_ensembl,
    mapped_ensembl
)


cat("\n")
cat("H5AD GO UNIVERSE MAPPING SUMMARY\n")
cat("------------------------------------------------------------\n")

cat(
    "Raw H5AD feature IDs:                 ",
    uniqueN(background_ensembl_raw),
    "\n"
)

cat(
    "Cleaned unique Ensembl IDs:           ",
    length(background_ensembl),
    "\n"
)

cat(
    "Ensembl IDs mapped to >=1 SYMBOL:     ",
    length(mapped_ensembl),
    "\n"
)

cat(
    "Unique SYMBOLs in final GO universe:  ",
    length(background_genes),
    "\n"
)

cat(
    "Unmapped cleaned Ensembl IDs:         ",
    length(unmapped_ensembl),
    "\n"
)


# We expect roughly 17,000 genes after mapping for this dataset.
if (
    length(background_genes) < 10000L
) {

    warning(
        paste0(
            "Only ",
            length(background_genes),
            " unique gene SYMBOLs were generated from ",
            length(background_ensembl),
            " cleaned H5AD Ensembl IDs. Please inspect the mapping audit."
        ),
        call. = FALSE
    )
}


if (
    length(background_genes) == 0L
) {

    stop(
        "No H5AD Ensembl genes could be mapped to gene SYMBOLs."
    )
}


# -----------------------------------------------------------------------------
# Save full raw -> cleaned Ensembl audit
# -----------------------------------------------------------------------------

ensembl_cleaning_audit <- unique(
    data.table(
        H5AD_feature_raw = background_ensembl_raw,
        ENSEMBL = sub(
            "\\.[0-9]+$",
            "",
            background_ensembl_raw
        )
    )
)


fwrite(
    ensembl_cleaning_audit,
    file.path(
        OUTPUT_DIR,
        "GO_h5ad_Ensembl_version_cleaning_audit.csv"
    )
)


# -----------------------------------------------------------------------------
# Save ENSEMBL -> SYMBOL mapping
# -----------------------------------------------------------------------------

fwrite(
    gene_mapping,
    file.path(
        OUTPUT_DIR,
        "GO_h5ad_ENSEMBL_to_SYMBOL_mapping.csv"
    )
)


# -----------------------------------------------------------------------------
# Save exact SYMBOL universe used by enrichGO
# -----------------------------------------------------------------------------

fwrite(
    data.table(
        gene = background_genes
    ),
    file.path(
        OUTPUT_DIR,
        "GO_h5ad_gene_universe.csv"
    )
)


# -----------------------------------------------------------------------------
# Save cleaned Ensembl IDs that could not be mapped to SYMBOL
# -----------------------------------------------------------------------------

fwrite(
    data.table(
        ENSEMBL = unmapped_ensembl
    ),
    file.path(
        OUTPUT_DIR,
        "GO_h5ad_unmapped_ENSEMBL.csv"
    )
)


cat("============================================================\n")


# Free the large SingleCellExperiment object before GO analysis
rm(
    sce
)

gc()


# =============================================================================
# 18. Function:
#
# ALL enriched LR pairs -> ALL UNIQUE ligand/receptor genes
#
# IMPORTANT:
# No top-N gene restriction is applied.
# =============================================================================

get_all_unique_genes <- function(
    lr_dt
) {


    x <- copy(
        lr_dt
    )


    setorder(
        x,
        LR_rank
    )


    # NO Top-N LR filtering here.
    # Every enriched LR pair supplied in lr_dt is retained.

    if (
        nrow(x) == 0L
    ) {

        return(
            data.table()
        )
    }


    gene_rows <- vector(
        "list",
        nrow(x)
    )


    for (
        i in seq_len(
            nrow(x)
        )
    ) {


        gene_rows[[i]] <- data.table(

            gene = c(
                x$ligand[i],
                x$receptor[i]
            ),

            role = c(
                "Ligand",
                "Receptor"
            ),

            LR_rank =
                x$LR_rank[i],

            ligand =
                x$ligand[i],

            receptor =
                x$receptor[i],

            interaction_score =
                x$interaction_score[i],

            enrichment =
                x$enrichment[i],

            n_samples =
                x$n_samples[i]
        )
    }


    gene_dt <- rbindlist(
        gene_rows,
        use.names = TRUE,
        fill = TRUE
    )


    gene_dt <- gene_dt[
        !is.na(gene) &
            gene != ""
    ]


    # Follow LR rank and retain only the first occurrence
    # if a gene occurs in multiple LR pairs.
    gene_dt <- gene_dt[
        !duplicated(
            gene
        )
    ]


    # NO head(..., 30) here.
    # Keep every unique ligand/receptor gene represented by ALL enriched LR pairs.

    gene_dt[
        ,
        gene_rank :=
            seq_len(.N)
    ]


    setcolorder(

        gene_dt,

        c(
            "gene_rank",
            "gene",
            "role",
            "LR_rank",
            "ligand",
            "receptor",
            "interaction_score",
            "enrichment",
            "n_samples"
        )
    )


    gene_dt
}


# =============================================================================
# 19. Generate ALL unique genes from ALL enriched LR pairs
#     for each exact cell pair
# =============================================================================

go_gene_list <- list()

gene_counter <- 0L


for (
    pair_i in
    PAIR_INFO$cell_pair
) {


    lr_i <- all_enriched_lr[
        cell_pair ==
            pair_i
    ]


    if (
        nrow(lr_i) == 0L
    ) {

        warning(
            "No enriched LR pairs for: ",
            pair_i,
            call. = FALSE
        )

        next
    }


    genes_i <- get_all_unique_genes(
        lr_dt =
            lr_i
    )


    if (
        nrow(genes_i) == 0L
    ) {

        next
    }


    gene_counter <-
        gene_counter + 1L


    genes_i[
        ,
        `:=`(

            comparison =
                PAIR_INFO[
                    cell_pair ==
                        pair_i,
                    comparison
                ],

            pair_order =
                PAIR_INFO[
                    cell_pair ==
                        pair_i,
                    pair_order
                ],

            cell_pair =
                pair_i
        )
    ]


    go_gene_list[[gene_counter]] <-
        genes_i
}


if (
    length(go_gene_list) == 0L
) {

    stop(
        "No GO gene sets generated."
    )
}


go_gene_dt <- rbindlist(

    go_gene_list,

    use.names = TRUE,

    fill = TRUE
)


setorder(
    go_gene_dt,
    pair_order,
    gene_rank
)


# =============================================================================
# 20. GO input gene summary
# =============================================================================

go_gene_summary <- go_gene_dt[
    ,
    .(

        n_GO_genes_raw =
            uniqueN(
                gene
            ),

        n_GO_genes_in_background =
            uniqueN(
                gene[
                    gene %chin%
                        background_genes
                ]
            ),

        n_ligands =
            sum(
                role ==
                    "Ligand"
            ),

        n_receptors =
            sum(
                role ==
                    "Receptor"
            ),

        max_LR_rank_contributing =
            max(
                LR_rank,
                na.rm = TRUE
            )
    ),
    by = .(
        comparison,
        pair_order,
        cell_pair
    )
]


setorder(
    go_gene_summary,
    pair_order
)


go_gene_summary <- merge(

    go_gene_summary,

    all_enriched_lr_summary[
        ,
        .(
            pair_order,
            cell_pair,
            n_enriched_LR_pairs =
                n_LR
        )
    ],

    by = c(
        "pair_order",
        "cell_pair"
    ),

    all.x = TRUE
)


setorder(
    go_gene_summary,
    pair_order
)


cat("\n")
cat("============================================================\n")
cat("GO INPUT: ALL ENRICHED LR -> ALL UNIQUE LIGAND/RECEPTOR GENES\n")
cat("============================================================\n")


print(
    go_gene_summary
)


fwrite(

    go_gene_dt,

    file.path(
        OUTPUT_DIR,
        "pair_specific_enrichment3_allLR_allUniqueGenes.csv"
    )
)


fwrite(

    go_gene_summary,

    file.path(
        OUTPUT_DIR,
        "pair_specific_enrichment3_allLR_allUniqueGenes_summary.csv"
    )
)


# =============================================================================
# 21. Run GO enrichment
#
# Foreground:
#     all unique ligand/receptor genes from ALL enriched LR pairs
#
# Background:
#     all genes present in cellular_annotated.h5ad (~17,000)
# =============================================================================

go_result_list <- list()

go_counter <- 0L


for (
    pair_i in
    PAIR_INFO$cell_pair
) {


    genes_i_raw <- go_gene_dt[
        cell_pair ==
            pair_i,
        unique(
            gene
        )
    ]


    if (
        length(genes_i_raw) == 0L
    ) {

        next
    }


    # Foreground genes must belong to the H5AD gene universe.
    genes_i <- intersect(
        genes_i_raw,
        background_genes
    )


    genes_not_in_background <- setdiff(
        genes_i_raw,
        background_genes
    )


    if (
        length(genes_not_in_background) > 0L
    ) {

        warning(
            pair_i,
            ": ",
            length(genes_not_in_background),
            " LR genes are not present in the H5AD gene universe: ",
            paste(
                genes_not_in_background,
                collapse = ", "
            ),
            call. = FALSE
        )
    }


    if (
        length(genes_i) == 0L
    ) {

        warning(
            "No GO input genes remain after intersecting with background for: ",
            pair_i,
            call. = FALSE
        )

        next
    }


    comparison_i <- PAIR_INFO[
        cell_pair ==
            pair_i,
        comparison
    ]


    pair_order_i <- PAIR_INFO[
        cell_pair ==
            pair_i,
        pair_order
    ]


    n_lr_i <- nrow(
        all_enriched_lr[
            cell_pair ==
                pair_i
        ]
    )


    cat("\n")
    cat("------------------------------------------------------------\n")
    cat("GO:", pair_i, "\n")
    cat("Enriched LR pairs:", n_lr_i, "\n")

    cat(
        "Unique LR genes before universe intersection:",
        length(genes_i_raw),
        "\n"
    )

    cat(
        "Unique LR genes used for GO:",
        length(genes_i),
        "\n"
    )

    cat(
        "Expressed-gene universe:",
        length(background_genes),
        "\n"
    )

    cat("------------------------------------------------------------\n")


    for (
        ontology_i in
        ONTOLOGY_ORDER
    ) {


        go_obj <- enrichGO(

            gene =
                genes_i,

            # Custom universe = all genes expressed in this dataset.
            universe =
                background_genes,

            OrgDb =
                org.Hs.eg.db,

            keyType =
                "SYMBOL",

            ont =
                ontology_i,

            pAdjustMethod =
                "BH",

            pvalueCutoff =
                1,

            qvalueCutoff =
                1,

            readable =
                TRUE
        )


        go_df <- as.data.table(

            as.data.frame(
                go_obj
            )
        )


        if (
            nrow(go_df) == 0L
        ) {

            next
        }


        go_counter <-
            go_counter + 1L


        go_df[
            ,
            `:=`(

                comparison =
                    comparison_i,

                pair_order =
                    pair_order_i,

                cell_pair =
                    pair_i,

                ontology =
                    ontology_i,

                n_enriched_LR =
                    n_lr_i,

                n_input_genes_raw =
                    length(
                        genes_i_raw
                    ),

                n_input_genes =
                    length(
                        genes_i
                    ),

                n_background_genes =
                    length(
                        background_genes
                    )
            )
        ]


        go_result_list[[go_counter]] <-
            go_df
    }
}


if (
    length(go_result_list) == 0L
) {

    stop(
        "No GO enrichment results returned."
    )
}


# =============================================================================
# 22. Combine GO results
# =============================================================================

go_all <- rbindlist(

    go_result_list,

    use.names = TRUE,

    fill = TRUE
)


setorder(
    go_all,
    pair_order,
    ontology,
    p.adjust
)


fwrite(

    go_all,

    file.path(
        OUTPUT_DIR,
        "pair_specific_enrichment3_allLR_allGenes_h5adUniverse_GO_all.csv"
    )
)


# =============================================================================
# 23. Significant GO results
# =============================================================================

go_sig <- go_all[

    !is.na(p.adjust) &

    p.adjust <
        FDR_CUTOFF
]


fwrite(

    go_sig,

    file.path(
        OUTPUT_DIR,
        "pair_specific_enrichment3_allLR_allGenes_h5adUniverse_GO_FDR05.csv"
    )
)


# =============================================================================
# 23b. Supplementary table:
#      significant GO terms and constituent ligand/receptor genes
#
# One row = one:
#   source + target + GO term + constituent LR gene
#
# Final publication-ready columns:
#   source
#   target
#   ontology
#   term_id
#   term_description
#   gene
#   role
#   GeneRatio
#   BgRatio
#   Count
#   pvalue
#   FDR
#
# Only GO terms with FDR < 0.05 are included.
# =============================================================================

GO_GENE_TABLE_FILE <- file.path(
    OUTPUT_DIR,
    "GO_terms_and_constituent_LR_genes_FDR05.csv"
)


if (
    nrow(go_sig) > 0L
) {


    # -------------------------------------------------------------------------
    # Expand clusterProfiler geneID:
    #
    #     APOE/LRP1/APP
    #
    # becomes three rows:
    #
    #     APOE
    #     LRP1
    #     APP
    # -------------------------------------------------------------------------

    go_sig_genes <- go_sig[
        ,
        .(
            gene =
                unlist(
                    strsplit(
                        as.character(geneID),
                        "/",
                        fixed = TRUE
                    ),
                    use.names = FALSE
                )
        ),
        by = .(
            comparison,
            pair_order,
            cell_pair,
            ontology,
            ID,
            Description,
            GeneRatio,
            BgRatio,
            pvalue,
            p.adjust,
            Count
        )
    ]


    go_sig_genes[
        ,
        gene :=
            trimws(
                as.character(gene)
            )
    ]


    go_sig_genes <- go_sig_genes[
        !is.na(gene) &
            gene != ""
    ]


    # -------------------------------------------------------------------------
    # Determine whether each constituent gene appears as a ligand, receptor,
    # or both among the enriched LR pairs for that exact cell pair.
    #
    # Using all_enriched_lr is more complete than the deduplicated GO-input
    # gene table because the same gene can occur in multiple LR pairs and can
    # potentially serve in both roles.
    # -------------------------------------------------------------------------

    ligand_role_lookup <- unique(
        all_enriched_lr[
            !is.na(ligand) &
                ligand != "",
            .(
                cell_pair,
                gene = ligand,
                role = "Ligand"
            )
        ]
    )


    receptor_role_lookup <- unique(
        all_enriched_lr[
            !is.na(receptor) &
                receptor != "",
            .(
                cell_pair,
                gene = receptor,
                role = "Receptor"
            )
        ]
    )


    gene_role_lookup <- rbindlist(
        list(
            ligand_role_lookup,
            receptor_role_lookup
        ),
        use.names = TRUE,
        fill = TRUE
    )


    gene_role_lookup <- gene_role_lookup[
        ,
        .(
            role =
                paste(
                    sort(
                        unique(role)
                    ),
                    collapse = "/"
                )
        ),
        by = .(
            cell_pair,
            gene
        )
    ]


    # -------------------------------------------------------------------------
    # Add source and target cell types explicitly.
    # -------------------------------------------------------------------------

    cell_pair_lookup <- unique(
        all_enriched_lr[
            ,
            .(
                cell_pair,
                source,
                target
            )
        ]
    )


    go_sig_genes <- merge(
        go_sig_genes,
        gene_role_lookup,
        by = c(
            "cell_pair",
            "gene"
        ),
        all.x = TRUE,
        sort = FALSE
    )


    go_sig_genes <- merge(
        go_sig_genes,
        cell_pair_lookup,
        by = "cell_pair",
        all.x = TRUE,
        sort = FALSE
    )


    # -------------------------------------------------------------------------
    # Sort using pair_order internally, then retain only the final 12 columns.
    # -------------------------------------------------------------------------

    setorder(
        go_sig_genes,
        pair_order,
        ontology,
        p.adjust,
        ID,
        gene
    )


    go_sig_genes_supp <- go_sig_genes[
        ,
        .(
            source,
            target,
            ontology,
            term_id = ID,
            term_description = Description,
            gene,
            role,
            GeneRatio,
            BgRatio,
            Count,
            pvalue,
            FDR = p.adjust
        )
    ]


    fwrite(
        go_sig_genes_supp,
        GO_GENE_TABLE_FILE
    )


    cat(
        "Supplementary GO term-gene rows:",
        nrow(go_sig_genes_supp),
        "\n"
    )


    cat(
        "Supplementary GO term-gene table:\n",
        GO_GENE_TABLE_FILE,
        "\n"
    )


} else {


    warning(
        "No FDR-significant GO terms; constituent-gene supplementary table is empty.",
        call. = FALSE
    )


    empty_go_gene_table <- data.table(
        source = character(),
        target = character(),
        ontology = character(),
        term_id = character(),
        term_description = character(),
        gene = character(),
        role = character(),
        GeneRatio = character(),
        BgRatio = character(),
        Count = integer(),
        pvalue = numeric(),
        FDR = numeric()
    )


    fwrite(
        empty_go_gene_table,
        GO_GENE_TABLE_FILE
    )
}


cat("\n")
cat("============================================================\n")
cat("GO SUMMARY\n")
cat("============================================================\n")


cat(
    "Total GO result rows:",
    nrow(go_all),
    "\n"
)


cat(
    "FDR < 0.05:",
    nrow(go_sig),
    "\n"
)


# =============================================================================
# 24. Convert GeneRatio
# =============================================================================

ratio_to_numeric <- function(x) {


    vapply(

        strsplit(
            as.character(x),
            "/",
            fixed = TRUE
        ),

        function(z) {


            if (
                length(z) != 2L
            ) {

                return(
                    NA_real_
                )
            }


            as.numeric(z[1]) /
                as.numeric(z[2])
        },

        numeric(1)
    )
}


go_all[
    ,
    GeneRatio_numeric :=
        ratio_to_numeric(
            GeneRatio
        )
]


go_all[
    ,
    neglog10_FDR :=
        -log10(
            pmax(
                p.adjust,
                .Machine$double.xmin
            )
        )
]


# =============================================================================
# 25. Comparison dotplot
# =============================================================================

make_comparison_dotplot <- function(
    comparison_name
) {


    cat("\n")
    cat(
        "Creating plot:",
        comparison_name,
        "\n"
    )


    comparison_pairs <- PAIR_INFO[
        comparison ==
            comparison_name
    ][
        order(
            pair_order
        ),
        cell_pair
    ]


    comparison_pairs <- comparison_pairs[
        comparison_pairs %in%
            go_all$cell_pair
    ]


    if (
        length(comparison_pairs) == 0L
    ) {

        return(
            invisible(NULL)
        )
    }


    comparison_sig <- go_all[

        comparison ==
            comparison_name &

        !is.na(p.adjust) &

        p.adjust <
            FDR_CUTOFF
    ]


    if (
        nrow(comparison_sig) == 0L
    ) {

        warning(
            "No significant GO terms for: ",
            comparison_name,
            call. = FALSE
        )

        return(
            invisible(NULL)
        )
    }


    # -------------------------------------------------------------------------
    # Each exact cell pair contributes its top 3 significant GO terms ACROSS
    # all three ontologies combined (BP + MF + CC), NOT top 3 within each
    # ontology separately.
    #
    # Therefore, one cell pair contributes at most 3 displayed GO terms total.
    # Those 3 terms can come from any mixture of BP / MF / CC depending on FDR.
    # For example, a pair could contribute 2 BP + 1 CC, or 3 BP, etc.
    #
    # Selection is based primarily on FDR (p.adjust), with raw p-value used
    # only as a deterministic tie-breaker. If the same GO term is selected
    # by multiple cell pairs, it appears only once in the final term list,
    # while the dotplot still shows that term across all cell pairs in the
    # comparison.
    # -------------------------------------------------------------------------

    selected_terms <- comparison_sig[

        order(
            p.adjust,
            pvalue
        ),

        head(
            .SD,
            TOP_N_GO_ACROSS_ONTOLOGIES_PER_PAIR
        ),

        by = .(
            cell_pair
        )
    ][
        ,
        .(
            best_FDR =
                min(
                    p.adjust,
                    na.rm = TRUE
                )
        ),
        by = .(
            ontology,
            ID,
            Description
        )
    ]


    selected_terms[
        ,
        ontology_order :=
            match(
                ontology,
                ONTOLOGY_ORDER
            )
    ]


    setorder(
        selected_terms,
        ontology_order,
        best_FDR,
        Description
    )


    selected_terms[
        ,
        term_key :=
            seq_len(.N)
    ]


    # Keep GO ID internally so different GO terms remain uniquely identified
    selected_terms[
        ,
        term_label :=
            paste0(
                Description,
                " [",
                ID,
                "]"
            )
    ]


    # -------------------------------------------------------------------------
    # Complete cell-pair × GO-term grid
    # -------------------------------------------------------------------------

    grid_dt <- CJ(

        cell_pair =
            comparison_pairs,

        term_key =
            selected_terms$term_key,

        unique =
            TRUE
    )


    grid_dt <- merge(

        grid_dt,

        selected_terms[
            ,
            .(
                term_key,
                ontology,
                ID,
                Description,
                term_label
            )
        ],

        by =
            "term_key",

        all.x =
            TRUE
    )


    actual_dt <- go_all[
        comparison ==
            comparison_name,
        .(
            cell_pair,
            ontology,
            ID,
            p.adjust,
            GeneRatio_numeric
        )
    ]


    plot_dt <- merge(

        grid_dt,

        actual_dt,

        by = c(
            "cell_pair",
            "ontology",
            "ID"
        ),

        all.x =
            TRUE
    )


    plot_dt[
        ,
        neglog10_FDR :=
            fifelse(

                !is.na(p.adjust),

                -log10(
                    pmax(
                        p.adjust,
                        .Machine$double.xmin
                    )
                ),

                NA_real_
            )
    ]


    # -------------------------------------------------------------------------
    # GO term ordering
    # -------------------------------------------------------------------------

    term_order <- plot_dt[

        !is.na(neglog10_FDR),

        .(
            strongest_signal =
                max(
                    neglog10_FDR,
                    na.rm = TRUE
                )
        ),

        by = .(
            ontology,
            term_label
        )
    ]


    term_order[
        ,
        ontology_order :=
            match(
                ontology,
                ONTOLOGY_ORDER
            )
    ]


    setorder(
        term_order,
        ontology_order,
        -strongest_signal,
        term_label
    )


    plot_dt[
        ,
        term_label :=
            factor(

                term_label,

                levels =
                    rev(
                        term_order$term_label
                    )
            )
    ]


    plot_dt[
        ,
        cell_pair :=
            factor(

                cell_pair,

                levels =
                    comparison_pairs
            )
    ]


    plot_dt[
        ,
        ontology :=
            factor(

                ontology,

                levels =
                    ONTOLOGY_ORDER
            )
    ]


    # -------------------------------------------------------------------------
    # Safe output name
    # -------------------------------------------------------------------------

    safe_name <- gsub(
        "[^A-Za-z0-9]+",
        "_",
        comparison_name
    )


    # -------------------------------------------------------------------------
    # Save plotting data
    # -------------------------------------------------------------------------

    plot_save <- copy(
        plot_dt
    )


    plot_save[
        ,
        `:=`(

            cell_pair =
                as.character(
                    cell_pair
                ),

            ontology =
                as.character(
                    ontology
                ),

            term_label =
                as.character(
                    term_label
                )
        )
    ]


    fwrite(

        plot_save,

        file.path(

            OUTPUT_DIR,

            paste0(
                safe_name,
                "_GO_dotplot_data.csv"
            )
        )
    )


    # =============================================================================
    # Plot
    #
    # IMPORTANT:
    # GO IDs are removed ONLY from displayed y-axis labels.
    # Internal term_label still contains [GO:xxxxxxx].
    # =============================================================================

    p <- ggplot(

        plot_dt,

        aes(
            x =
                cell_pair,

            y =
                term_label
        )
    ) +


        geom_point(

            aes(

                size =
                    GeneRatio_numeric,

                color =
                    neglog10_FDR
            ),

            na.rm =
                TRUE
        ) +


        facet_grid(

            rows =
                vars(
                    ontology
                ),

            scales =
                "free_y",

            space =
                "free_y"
        ) +


        # ---------------------------------------------------------------------
        # Remove [GO:xxxxxxx] from displayed Y-axis labels
        # ---------------------------------------------------------------------

        scale_y_discrete(

            labels = function(x) {

                gsub(
                    "\\s*\\[GO:[^]]+\\]$",
                    "",
                    x
                )
            }
        ) +


        scale_color_gradient(

            low =
                "grey80",

            high =
                "#B40426",

            name =
                expression(
                    -log[10]("FDR")
                )
        ) +


        scale_size_continuous(

            name =
                "GeneRatio"
        ) +


        labs(

            # title =
            #     comparison_name,

            # subtitle =
            #     paste0(

            #         "Pair-specific GO; LR enrichment >= ",

            #         ENRICHMENT_CUTOFF,

            #         "; all enriched LR pairs; all unique ligand/receptor genes"
            #     ),

            x =
                NULL,

            y =
                NULL

            # caption =
            #     paste0(
            #         "LR interactions required in >= ",
            #         MIN_SAMPLES,
            #         " samples; GO universe = ",
            #         length(background_genes),
            #         " expressed genes"
            #     )
        ) +


        theme_bw(
            base_size =
                15
        ) +


        theme(

            panel.grid.major =
                element_line(
                    color =
                        "grey92",
                    linewidth =
                        0.3
                ),

            panel.grid.minor =
                element_blank(),

            axis.text.x =
                element_text(
                    angle =
                        45,
                    hjust =
                        1,
                    vjust =
                        1
                ),

            strip.background =
                element_rect(
                    fill =
                        "grey92",
                    color =
                        "grey60"
                ),

            strip.text =
                element_text(
                    face =
                        "bold"
                ),

            plot.title =
                element_text(
                    face =
                        "bold"
                ),

            plot.caption =
                element_text(
                    hjust =
                        0
                )
        )


    n_terms <- uniqueN(
        plot_dt$term_label
    )


    n_pairs <- length(
        comparison_pairs
    )


    # plot_width <- max(
    #     8,
    #     1.2 *
    #         n_pairs +
    #         5
    # )
    plot_width = 8


    # plot_height <- max(
    #     7,
    #     0.30 *
    #         n_terms +
    #         3
    # )
    plot_height = 7


    # -------------------------------------------------------------------------
    # PDF
    # -------------------------------------------------------------------------

    ggsave(

        filename =
            file.path(

                OUTPUT_DIR,

                paste0(
                    safe_name,
                    "_GO_dotplot.pdf"
                )
            ),

        plot =
            p,

        width =
            plot_width,

        height =
            plot_height,

        units =
            "in"
    )


    # -------------------------------------------------------------------------
    # PNG
    # -------------------------------------------------------------------------

    ggsave(

        filename =
            file.path(

                OUTPUT_DIR,

                paste0(
                    safe_name,
                    "_GO_dotplot.png"
                )
            ),

        plot =
            p,

        width =
            plot_width,

        height =
            plot_height,

        units =
            "in",

        dpi =
            300
    )
}


# =============================================================================
# 26. Generate four comparison figures
# =============================================================================

for (
    comparison_i in
    COMPARISON_ORDER
) {


    make_comparison_dotplot(
        comparison_i
    )
}


# =============================================================================
# 27. Final summary
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("FINAL GO ANALYSIS FINISHED\n")
cat("============================================================\n")


cat(
    "Minimum samples per LR:       ",
    MIN_SAMPLES,
    "\n"
)


cat(
    "Enrichment cutoff:            ",
    ENRICHMENT_CUTOFF,
    "x\n"
)


cat(
    "LR pairs per cell pair:       ALL with enrichment >= ",
    ENRICHMENT_CUTOFF,
    " (no Top-N LR filter)\n",
    sep = ""
)


cat(
    "GO foreground genes:          all unique ligand/receptor genes from ALL enriched LR pairs\n"
)


cat(
    "GO background genes:          ",
    length(background_genes),
    "\n"
)


cat(
    "Selected cell pairs analyzed: ",
    uniqueN(
        go_gene_dt$cell_pair
    ),
    " / 26\n"
)


cat(
    "FDR-significant GO rows:      ",
    nrow(
        go_sig
    ),
    "\n"
)


cat(
    "GO terms plotted per cell pair across BP/MF/CC: ",
    TOP_N_GO_ACROSS_ONTOLOGIES_PER_PAIR,
    " total (ontologies pooled for ranking)\n"
)


cat(
    "GO term-gene supp. table:     ",
    GO_GENE_TABLE_FILE,
    "\n"
)


cat(
    "Output directory:\n",
    OUTPUT_DIR,
    "\n"
)


cat("============================================================\n")