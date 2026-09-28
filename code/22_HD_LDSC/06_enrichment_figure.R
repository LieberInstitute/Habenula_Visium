#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(scales)
})


###############################################################################
# User settings
###############################################################################

input_file <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Habenula_Visium/processed-data/10_HD_bin_level/LDSC/",
    "04_sldsc_results/all_GWAS_DAR_sldsc_summary.tsv"
)

cell_type_map_file <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Hb_multiome/raw-data/cell_type_map.csv"
)

gwas_info_file <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Hb_multiome/processed-data/10_MAGMA/RNA/gwas_info.csv"
)

output_dir <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Habenula_Visium/processed-data/10_HD_bin_level/LDSC/",
    "04_sldsc_results/heatmap"
)


###############################################################################
# DAR sets to plot
###############################################################################

dar_sets_to_plot <- c(
    "open",
    "closed",
    "all"
)



###############################################################################
# Minimum common SNP count for display
###############################################################################

minimum_n_snps <- 100L


###############################################################################
# Plot only significant cells?
#
# FALSE:
#   color all ordinary GWAS cells by coefficient Z-score.
#
# TRUE:
#   only cells showing evidence for positive enrichment after cell-type-specific
#   BH-FDR correction are colored.
###############################################################################
show_only_significant <- FALSE


###############################################################################
# Significance measure for ORDINARY GWAS rows
#
# For EACH DAR set separately (open / closed / all):
#   within EACH cell type / annotation, raw coefficient_p_one_sided values
#   across all displayed GWAS traits are adjusted using Benjamini-Hochberg FDR.
#
# Therefore there is NO global adjustment across cell types.
###############################################################################

significance_column <- "coefficient_FDR_one_sided_within_cell_type"


###############################################################################
# Heatmap color limit
#
# Ordinary rows:
#   coefficient Z-score
###############################################################################

color_limit <- 4


###############################################################################
# Output directory
###############################################################################

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


###############################################################################
# Helper: require columns
###############################################################################

require_columns <- function(
    data,
    required,
    object_name
) {

    missing_columns <- setdiff(
        required,
        names(data)
    )

    if (length(missing_columns) > 0L) {

        stop(
            object_name,
            " is missing required columns: ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    }
}


###############################################################################
# Helper: significance stars
###############################################################################

make_significance_stars <- function(
    p,
    require_positive = FALSE,
    z = NULL
) {

    if (require_positive) {

        if (is.null(z)) {
            stop(
                "z must be supplied when require_positive = TRUE."
            )
        }

        positive <- (
            !is.na(z) &
            z > 0
        )

    } else {

        positive <- rep(
            TRUE,
            length(p)
        )
    }


    fifelse(
        positive &
        !is.na(p) &
        p < 0.001,
        "***",

        fifelse(
            positive &
            !is.na(p) &
            p < 0.01,
            "**",

            fifelse(
                positive &
                !is.na(p) &
                p < 0.05,
                "*",
                ""
            )
        )
    )
}


###############################################################################
# Read S-LDSC summary
###############################################################################

if (!file.exists(input_file)) {

    stop(
        "Input file does not exist: ",
        input_file
    )
}


if (
    is.na(
        file.info(input_file)$size
    ) ||
    file.info(input_file)$size == 0
) {

    stop(
        "Input file is empty: ",
        input_file
    )
}


dt <- fread(
    input_file,
    sep = "\t",
    quote = "",
    na.strings = c(
        "",
        "NA",
        "nan"
    )
)


###############################################################################
# Required columns
###############################################################################

required_columns <- c(
    "trait",
    "dar_set",
    "annotation",
    "source_resolution",
    "n_snps",

    # Raw one-sided P value used for cell-type-specific BH-FDR:
    "coefficient_p_one_sided"
)


require_columns(
    dt,
    required_columns,
    "S-LDSC summary"
)


###############################################################################
# Identify coefficient Z-score
###############################################################################

if (
    "coefficient_z" %in%
    names(dt)
) {

    dt[
        ,
        coefficient_z_plot :=
            suppressWarnings(
                as.numeric(
                    coefficient_z
                )
            )
    ]

    z_score_source <-
        "coefficient_z"


} else if (
    "Coefficient_z-score" %in%
    names(dt)
) {

    dt[
        ,
        coefficient_z_plot :=
            suppressWarnings(
                as.numeric(
                    `Coefficient_z-score`
                )
            )
    ]

    z_score_source <-
        "Coefficient_z-score"


} else if (
    all(
        c(
            "Coefficient",
            "Coefficient_std_error"
        ) %in% names(dt)
    )
) {

    dt[
        ,
        coefficient_z_plot :=
            suppressWarnings(
                as.numeric(
                    Coefficient
                )
            ) /
            suppressWarnings(
                as.numeric(
                    Coefficient_std_error
                )
            )
    ]

    z_score_source <-
        "Coefficient / Coefficient_std_error"


} else {

    stop(
        paste0(
            "No coefficient Z-score could be identified. ",
            "Expected coefficient_z, Coefficient_z-score, ",
            "or Coefficient + Coefficient_std_error."
        )
    )
}


###############################################################################
# Convert original one-sided P value to numeric
###############################################################################

dt[
    ,
    coefficient_p_one_sided :=
        suppressWarnings(
            as.numeric(
                coefficient_p_one_sided
            )
        )
]


###############################################################################
# Clean variables
###############################################################################

dt[
    ,
    `:=`(

        trait =
            trimws(
                as.character(
                    trait
                )
            ),

        dar_set =
            trimws(
                as.character(
                    dar_set
                )
            ),

        annotation =
            trimws(
                as.character(
                    annotation
                )
            ),

        source_resolution =
            trimws(
                as.character(
                    source_resolution
                )
            ),

        n_snps =
            suppressWarnings(
                as.numeric(
                    n_snps
                )
            )
    )
]


###############################################################################
# Ordinary-row FDR is intentionally NOT calculated here.
#
# It is calculated inside the DAR-set loop AFTER the minimum_n_snps filter.
# Within each DAR set, BH-FDR is applied separately for each annotation/cell
# type across all displayed GWAS traits.
###############################################################################


###############################################################################
# Read GWAS name mapping
#
# Current S-LDSC "trait" may use either gwas_info.csv "nickname" or
# "manuscript_name".  The reference-figure grouping code below resolves both.
###############################################################################

if (!file.exists(gwas_info_file)) {

    stop(
        "GWAS info file does not exist: ",
        gwas_info_file
    )
}


gwas_info <- fread(
    gwas_info_file,
    sep = ",",
    quote = "\"",
    na.strings = c(
        "",
        "NA",
        "nan"
    ),
    encoding = "UTF-8"
)


require_columns(
    gwas_info,
    c(
        "nickname",
        "manuscript_name"
    ),
    "GWAS info file"
)


gwas_info[
    ,
    `:=`(
        nickname =
            trimws(
                as.character(
                    nickname
                )
            ),

        manuscript_name =
            trimws(
                as.character(
                    manuscript_name
                )
            )
    )
]


gwas_info <- gwas_info[
    !is.na(nickname) &
    nickname != "" &
    !is.na(manuscript_name) &
    manuscript_name != ""
]


if (
    anyDuplicated(
        gwas_info$nickname
    )
) {

    stop(
        "Duplicated nickname values were found in gwas_info.csv."
    )
}


if (
    anyDuplicated(
        gwas_info$manuscript_name
    )
) {

    stop(
        "Duplicated manuscript_name values were found in gwas_info.csv."
    )
}


gwas_name_lookup <- setNames(
    gwas_info$manuscript_name,
    gwas_info$nickname
)


###############################################################################
# Exact GWAS grouping/order from the reference figure
#
# The labels below are copied exactly from the reference figure.  Importantly,
# we do NOT assume that these labels are necessarily the same naming convention
# used in the current S-LDSC input.  Both the reference labels and the current
# input traits are resolved through gwas_info.csv using BOTH nickname and
# manuscript_name.
###############################################################################

gwas_group_order <- c(
    "P",
    "Psychiatric",
    "Substance Use"
)


# Exact left-to-right labels shown in the reference figure.
gwas_trait_order <- c(

    # Factor
    "p_factor_Grotzinger",

    # Psychiatric
    "compulsive_F1_Grotzinger",
    "externalizing_Linnér",
    "internalizing_F4_Grotzinger",
    "MDD_Howard",
    "neurodev_F3_Grotzinger",
    "panic_Forster",
    "SCZ_Trubetskoy",
    "SCZ/BPD_F2_Grotzinger",

    # Substance Use
    "AUD_Zhou",
    "CUD_Johnson",
    "CUD_Pasman",
    "OUD_Deak",
    "SUD_F5_Grotzinger",
    "SUD_Hatoum",
    "SUD_Polimanti"
)


reference_gwas <- data.table(
    reference_label = gwas_trait_order,
    gwas_group = c(
        "P",
        rep("Psychiatric", 8L),
        rep("Substance Use", 7L)
    ),
    reference_order = seq_along(gwas_trait_order)
)


###############################################################################
# Name lookups from gwas_info.csv
###############################################################################

# nickname -> manuscript_name
gwas_name_lookup <- setNames(
    gwas_info$manuscript_name,
    gwas_info$nickname
)


# manuscript_name -> nickname
gwas_nickname_from_manuscript_lookup <- setNames(
    gwas_info$nickname,
    gwas_info$manuscript_name
)


###############################################################################
# Helper: resolve any GWAS label to canonical nickname
#
# A label is accepted if it matches EITHER:
#   - gwas_info.csv nickname
#   - gwas_info.csv manuscript_name
###############################################################################

resolve_to_nickname <- function(x) {

    x <- trimws(
        as.character(x)
    )

    out <- rep(
        NA_character_,
        length(x)
    )

    is_nickname <- x %chin% gwas_info$nickname

    out[is_nickname] <- x[is_nickname]

    is_manuscript <- (
        !is_nickname &
        x %chin% gwas_info$manuscript_name
    )

    out[is_manuscript] <- unname(
        gwas_nickname_from_manuscript_lookup[
            x[is_manuscript]
        ]
    )

    out
}


###############################################################################
# Resolve the REFERENCE-FIGURE labels through gwas_info.csv
###############################################################################

reference_gwas[
    ,
    trait_nickname :=
        resolve_to_nickname(
            reference_label
        )
]


unresolved_reference_labels <- reference_gwas[
    is.na(trait_nickname) |
    trait_nickname == "",
    reference_label
]


if (length(unresolved_reference_labels) > 0L) {

    warning(
        paste0(
            "The following labels copied from the reference figure match neither ",
            "nickname nor manuscript_name in the CURRENT gwas_info.csv:\n",
            paste(
                unresolved_reference_labels,
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}


# A canonical GWAS must not map to more than one reference label.
reference_nickname_duplicates <- reference_gwas[
    !is.na(trait_nickname),
    .N,
    by = trait_nickname
][
    N > 1L
]


if (nrow(reference_nickname_duplicates) > 0L) {

    stop(
        paste0(
            "Multiple reference-figure labels resolve to the same canonical GWAS ",
            "nickname:\n",
            paste(
                capture.output(
                    print(reference_nickname_duplicates)
                ),
                collapse = "\n"
            )
        )
    )
}


###############################################################################
# Resolve CURRENT S-LDSC input traits through the same mapping
###############################################################################

# Preserve the exact value in all_GWAS_DAR_sldsc_summary.tsv.
dt[
    ,
    trait_input :=
        trimws(
            as.character(
                trait
            )
        )
]


dt[
    ,
    trait_nickname :=
        resolve_to_nickname(
            trait_input
        )
]


unresolved_input_traits <- sort(
    unique(
        dt[
            is.na(trait_nickname) |
            trait_nickname == "",
            trait_input
        ]
    )
)


if (length(unresolved_input_traits) > 0L) {

    warning(
        paste0(
            "The following current S-LDSC input trait names match neither nickname ",
            "nor manuscript_name in gwas_info.csv and therefore cannot be matched ",
            "to the reference figure:\n",
            paste(
                unresolved_input_traits,
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}


###############################################################################
# Attach exact reference-figure group and display label by canonical nickname
###############################################################################

reference_resolved <- reference_gwas[
    !is.na(trait_nickname) &
    trait_nickname != ""
]


reference_label_lookup <- setNames(
    reference_resolved$reference_label,
    reference_resolved$trait_nickname
)


reference_group_lookup <- setNames(
    reference_resolved$gwas_group,
    reference_resolved$trait_nickname
)


reference_order_lookup <- setNames(
    reference_resolved$reference_order,
    reference_resolved$trait_nickname
)


dt[
    ,
    `:=`(
        trait_plot =
            unname(
                reference_label_lookup[
                    trait_nickname
                ]
            ),

        gwas_group =
            unname(
                reference_group_lookup[
                    trait_nickname
                ]
            ),

        reference_order =
            suppressWarnings(
                as.integer(
                    unname(
                        reference_order_lookup[
                            trait_nickname
                        ]
                    )
                )
            )
    )
]


###############################################################################
# Diagnose reference GWAS that are genuinely absent from the CURRENT input
###############################################################################

available_canonical_nicknames <- unique(
    dt[
        !is.na(trait_nickname),
        trait_nickname
    ]
)


missing_reference_gwas <- reference_resolved[
    !trait_nickname %chin% available_canonical_nicknames,
    reference_label
]


if (length(missing_reference_gwas) > 0L) {

    warning(
        paste0(
            "The following GWAS traits from the reference figure are genuinely ",
            "absent from the current S-LDSC input AFTER resolving both naming ",
            "systems through gwas_info.csv:\n",
            paste(
                missing_reference_gwas,
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}


###############################################################################
# Keep exactly the GWAS represented in the reference figure
###############################################################################

dt <- dt[
    !is.na(trait_plot) &
    trait_plot != "" &
    !is.na(gwas_group) &
    gwas_group != ""
]


if (nrow(dt) == 0L) {

    stop(
        paste0(
            "None of the reference-figure GWAS could be matched to the current ",
            "S-LDSC input. The script checked BOTH nickname and manuscript_name ",
            "for BOTH the reference labels and the input trait column."
        )
    )
}


###############################################################################
# Manuscript-friendly trait name for output tables
###############################################################################

dt[
    ,
    trait :=
        unname(
            gwas_name_lookup[
                trait_nickname
            ]
        )
]


dt[
    is.na(trait) |
    trait == "",
    trait :=
        trait_nickname
]


dt[
    ,
    trait_mapped := TRUE
]


###############################################################################
# GWAS naming/grouping audit
###############################################################################

gwas_name_audit <- unique(
    dt[
        ,
        .(
            trait_input,
            trait_nickname,
            trait,
            trait_plot,
            gwas_group,
            reference_order,
            trait_mapped
        )
    ]
)


setorder(
    gwas_name_audit,
    reference_order,
    trait_input
)


output_gwas_name_map <- file.path(
    output_dir,
    "GWAS_name_mapping.tsv"
)


fwrite(
    gwas_name_audit,
    output_gwas_name_map,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Basic filtering
#
# Keep rows with a valid raw one-sided S-LDSC P value.
# BH-FDR is calculated later within each DAR set, separately for each cell type.
###############################################################################

dt <- dt[
    !is.na(trait) &
    trait != "" &

    !is.na(dar_set) &
    dar_set != "" &

    !is.na(annotation) &
    annotation != "" &

    !is.na(n_snps) &

    is.finite(
        coefficient_z_plot
    ) &

    !is.na(
        coefficient_p_one_sided
    ) &

    is.finite(
        coefficient_p_one_sided
    ) &

    coefficient_p_one_sided >= 0 &
    coefficient_p_one_sided <= 1
]


if (nrow(dt) == 0L) {

    stop(
        "No complete rows remained after input validation."
    )
}




###############################################################################
# Validate DAR sets
###############################################################################

expected_dar_sets <- c(
    "open",
    "closed",
    "all"
)


unexpected_dar_sets <- setdiff(
    unique(
        dt$dar_set
    ),
    expected_dar_sets
)


if (
    length(
        unexpected_dar_sets
    ) > 0L
) {

    stop(
        "Unexpected dar_set values: ",
        paste(
            unexpected_dar_sets,
            collapse = ", "
        )
    )
}


###############################################################################
# Check unique result per:
#
# trait × dar_set × annotation
###############################################################################

duplicate_results <- dt[
    ,
    .N,
    by = .(
        trait,
        dar_set,
        annotation
    )
][
    N > 1L
]


if (
    nrow(
        duplicate_results
    ) > 0L
) {

    stop(
        paste0(
            "Multiple rows exist for some ",
            "trait × dar_set × annotation combinations:\n",
            paste(
                capture.output(
                    print(
                        duplicate_results
                    )
                ),
                collapse = "\n"
            )
        )
    )
}


###############################################################################
# Initial information
###############################################################################

message("")
message(
    "Z-score source: ",
    z_score_source
)

message(
    "Significance column: ",
    significance_column
)


###############################################################################
# Read cell-type mapping
###############################################################################

if (
    !file.exists(
        cell_type_map_file
    )
) {

    stop(
        "Cell-type mapping file does not exist: ",
        cell_type_map_file
    )
}


if (
    is.na(
        file.info(
            cell_type_map_file
        )$size
    ) ||
    file.info(
        cell_type_map_file
    )$size == 0
) {

    stop(
        "Cell-type mapping file is empty: ",
        cell_type_map_file
    )
}


cell_type_map <- fread(
    cell_type_map_file,
    sep = ",",
    quote = "\"",
    na.strings = c(
        "",
        "NA",
        "nan"
    )
)


###############################################################################
# Validate mapping columns
###############################################################################

required_map_columns <- c(
    "old_cell_type",
    "new_cell_type"
)


require_columns(
    cell_type_map,
    required_map_columns,
    "Cell-type mapping file"
)


###############################################################################
# Clean mapping
###############################################################################

cell_type_map[
    ,
    `:=`(

        old_cell_type =
            trimws(
                as.character(
                    old_cell_type
                )
            ),

        new_cell_type =
            trimws(
                as.character(
                    new_cell_type
                )
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
    nrow(
        cell_type_map
    ) == 0L
) {

    stop(
        "No valid rows remained in the cell-type mapping file."
    )
}


###############################################################################
# Check old-name uniqueness
###############################################################################

duplicated_old_names <- cell_type_map[
    duplicated(
        old_cell_type
    ) |
    duplicated(
        old_cell_type,
        fromLast = TRUE
    )
]


if (
    nrow(
        duplicated_old_names
    ) > 0L
) {

    duplicate_summary <- duplicated_old_names[
        ,
        .(
            mapped_names =
                paste(
                    unique(
                        new_cell_type
                    ),
                    collapse = " | "
                )
        ),
        by = old_cell_type
    ]


    message("")
    message(
        "Duplicated old_cell_type entries:"
    )

    print(
        duplicate_summary
    )


    stop(
        "Each old_cell_type must appear only once in the mapping file."
    )
}


###############################################################################
# Apply corrected cell-type names
###############################################################################

name_lookup <- setNames(
    cell_type_map$new_cell_type,
    cell_type_map$old_cell_type
)


dt[
    ,
    annotation_original :=
        annotation
]


dt[
    ,
    annotation_mapped :=
        annotation_original %chin%
        cell_type_map$old_cell_type
]


dt[
    ,
    annotation_display :=
        unname(
            name_lookup[
                annotation_original
            ]
        )
]


dt[
    is.na(
        annotation_display
    ) |
    annotation_display == "",
    annotation_display :=
        annotation_original
]


###############################################################################
# Mapping audit
###############################################################################

mapping_audit <- unique(
    dt[
        ,
        .(
            dar_set,
            annotation_original,
            annotation_display,
            annotation_mapped
        )
    ]
)


mapping_audit[
    ,
    dar_set_order :=
        match(
            dar_set,
            expected_dar_sets
        )
]


setorder(
    mapping_audit,
    dar_set_order,
    -annotation_mapped,
    annotation_display,
    annotation_original
)


mapping_audit[
    ,
    dar_set_order := NULL
]


output_annotation_map <- file.path(
    output_dir,
    "DAR_sldsc_annotation_name_mapping.tsv"
)


fwrite(
    mapping_audit,
    output_annotation_map,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Report unmapped names
###############################################################################

unmapped_annotations <- sort(
    unique(
        dt[
            annotation_mapped == FALSE,
            annotation_original
        ]
    )
)


if (
    length(
        unmapped_annotations
    ) > 0L
) {

    warning(
        paste0(
            "The following annotations were not found in ",
            basename(
                cell_type_map_file
            ),
            " and retain their original names:\n",
            paste(
                unmapped_annotations,
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}


###############################################################################
# Biological x-axis order
###############################################################################

custom_annotation_display_order <- c(

    # Whole habenula regions -- standalone rows
    "MHb",
    "LHb",

    # Medial habenula subpopulations
    "MHb_A",
    "MHb_B",
    "MHb_C",
    "MHb_D",

    # Lateral habenula subpopulations
    "LHb_A",
    "LHb_B",
    "LHb_C",
    "GABA_LHb_C.1",

    # Thalamic neuronal populations
    "Thalamus",
    "Excit. Thal",
    "Inhib. Thal",

    # Non-neuronal populations
    "Astrocyte",
    "OPC",
    "Oligo",
    "Microglia",
    "Ependymal",
    "Endo"
)

# Draw a narrow white separator line below the whole-region rows (MHb and LHb)
# without creating an extra blank y-axis label row.
whole_hb_separator_linewidth <- 5.0

# Thin black rules above and below the white gap.
# This creates a clean grouped-panel look without boxing the entire heatmap.
whole_hb_separator_border_linewidth <- 0.4


###############################################################################
# Significance description for ordinary rows
###############################################################################

significance_description <- paste0(
    "* FDR<0.05, ",
    "** FDR<0.01, ",
    "*** FDR<0.001; ",
    "BH-FDR across GWAS traits separately within each cell type"
)


###############################################################################
# DAR set titles
###############################################################################

dar_set_titles <- c(
    open = "Open DARs",
    closed = "Closed DARs",
    all = "All DARs"
)


###############################################################################
# Store outputs across DAR sets
###############################################################################

all_plot_data <- list()



###############################################################################
# Generate one heatmap per DAR set
###############################################################################

for (
    dar_set_i in
    dar_sets_to_plot
) {


    ###########################################################################
    # Skip absent sets
    ###########################################################################

    if (
        !dar_set_i %in%
        dt$dar_set
    ) {

        warning(
            "Skipping absent DAR set: ",
            dar_set_i,
            call. = FALSE
        )

        next
    }


    message("")
    message("============================================================")
    message(
        "Creating heatmap: ",
        dar_set_i
    )
    message("============================================================")


    ###########################################################################
    # Output paths
    ###########################################################################

    dar_output_dir <- file.path(
        output_dir,
        dar_set_i
    )


    dir.create(
        dar_output_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )


    output_pdf <- file.path(
        dar_output_dir,
        paste0(
            "DAR_sldsc_coefficient_z_heatmap_",
            dar_set_i,
            ".pdf"
        )
    )

    output_png <- file.path(
        dar_output_dir,
        paste0(
            "DAR_sldsc_coefficient_z_heatmap_",
            dar_set_i,
            ".png"
        )
    )

    output_plot_data <- file.path(
        dar_output_dir,
        paste0(
            "DAR_sldsc_heatmap_plot_data_",
            dar_set_i,
            ".tsv"
        )
    )



    ###########################################################################
    # Select DAR set
    ###########################################################################

    plot_dt <- copy(
        dt[
            dar_set ==
            dar_set_i
        ]
    )


    ###########################################################################
    # Filter annotations by SNP count
    ###########################################################################

    plot_dt <- plot_dt[
        n_snps >=
        minimum_n_snps
    ]


    if (
        nrow(
            plot_dt
        ) == 0L
    ) {

        warning(
            "No rows remain for ",
            dar_set_i,
            " after minimum_n_snps = ",
            minimum_n_snps,
            ".",
            call. = FALSE
        )

        next
    }


    ###########################################################################
    # BH-FDR for ordinary GWAS heatmap cells, separately by cell type
    #
    # IMPORTANT:
    # The current plot_dt contains ONE DAR set only and has already passed the
    # minimum_n_snps display filter.
    #
    # For EACH annotation/cell type separately, p.adjust() corrects the raw
    # one-sided S-LDSC P values across all displayed GWAS traits for that
    # cell type. Thus each cell type has its own FDR family within this DAR set.
    ###########################################################################

    plot_dt[
        ,
        coefficient_FDR_one_sided_within_cell_type :=
            p.adjust(
                coefficient_p_one_sided,
                method = "BH"
            ),
        by = annotation
    ]


    plot_dt[
        ,
        significance_value :=
            coefficient_FDR_one_sided_within_cell_type
    ]


    fdr_family_summary <- plot_dt[
        ,
        .(
            n_GWAS_tests = .N
        ),
        by = annotation
    ]

    message(
        "[", dar_set_i, "] FDR correction performed separately for ",
        nrow(fdr_family_summary),
        " cell types."
    )

    ###########################################################################
    # Ensure unique ordinary heatmap cells
    ###########################################################################

    duplicate_plot_cells <- plot_dt[
        ,
        .N,
        by = .(
            trait,
            annotation
        )
    ][
        N > 1L
    ]


    if (
        nrow(
            duplicate_plot_cells
        ) > 0L
    ) {

        stop(
            paste0(
                "Duplicate heatmap cells in ",
                dar_set_i,
                ":\n",
                paste(
                    capture.output(
                        print(
                            duplicate_plot_cells
                        )
                    ),
                    collapse = "\n"
                )
            )
        )
    }


    ###########################################################################
    # Annotation metadata
    ###########################################################################

    annotation_info <- unique(
        plot_dt[
            ,
            .(
                annotation,
                annotation_original,
                annotation_display,
                annotation_mapped,
                source_resolution,
                n_snps
            )
        ]
    )


    ###########################################################################
    # Validate metadata consistency WITHIN current DAR set
    ###########################################################################

    annotation_metadata_check <- annotation_info[
        ,
        .(
            n_resolutions =
                uniqueN(
                    source_resolution
                ),

            n_snp_counts =
                uniqueN(
                    n_snps
                ),

            resolution_values =
                paste(
                    sort(
                        unique(
                            source_resolution
                        )
                    ),
                    collapse = " | "
                ),

            n_snp_values =
                paste(
                    sort(
                        unique(
                            n_snps
                        )
                    ),
                    collapse = " | "
                )
        ),
        by = annotation
    ][
        n_resolutions > 1L |
        n_snp_counts > 1L
    ]


    if (
        nrow(
            annotation_metadata_check
        ) > 0L
    ) {

        print(
            annotation_metadata_check
        )

        stop(
            "Inconsistent annotation metadata within DAR set ",
            dar_set_i,
            "."
        )
    }


    annotation_info <- unique(
        annotation_info,
        by = "annotation"
    )


    ###########################################################################
    # Biological annotation order
    ###########################################################################

    annotation_info[
        ,
        custom_order_rank :=
            match(
                annotation_display,
                custom_annotation_display_order
            )
    ]


    annotation_info[
        ,
        annotation_is_unlisted :=
            is.na(
                custom_order_rank
            )
    ]


    annotation_info[
        annotation_is_unlisted == TRUE,
        custom_order_rank :=
            length(
                custom_annotation_display_order
            ) +
            frank(
                annotation_display,
                ties.method = "dense"
            )
    ]


    annotation_order <- annotation_info[
        order(
            annotation_is_unlisted,
            custom_order_rank,
            annotation_display,
            annotation
        ),
        annotation
    ]


    ###########################################################################
    # Ordinary GWAS significance
    ###########################################################################

    plot_dt[
        ,
        positive_enrichment_significant :=
            coefficient_z_plot > 0 &
            significance_value < 0.05
    ]


    plot_dt[
        ,
        significance_label :=
            make_significance_stars(
                p =
                    significance_value,
                require_positive =
                    TRUE,
                z =
                    coefficient_z_plot
            )
    ]


    ###########################################################################
    # Ordinary GWAS color
    ###########################################################################

    plot_dt[
        ,
        coefficient_z_clipped :=
            pmax(
                pmin(
                    coefficient_z_plot,
                    color_limit
                ),
                -color_limit
            )
    ]


    if (
        show_only_significant
    ) {

        plot_dt[
            ,
            fill_value :=
                fifelse(
                    positive_enrichment_significant,
                    coefficient_z_clipped,
                    NA_real_
                )
        ]

    } else {

        plot_dt[
            ,
            fill_value :=
                coefficient_z_clipped
        ]
    }


    ###########################################################################
    # Mark ordinary rows
    ###########################################################################

    plot_dt[
        ,
        row_type := "GWAS"
    ]


    ###########################################################################
    # Exact GWAS trait order from the reference figure
    ###########################################################################

    regular_trait_order <- gwas_trait_order[
        gwas_trait_order %chin%
        unique(
            plot_dt$trait_plot
        )
    ]


    ###########################################################################
    # Use ordinary GWAS rows only
    ###########################################################################

    plot_dt_full <- copy(
        plot_dt
    )


    ###########################################################################
    # Factor ordering
    #
    # GWAS blocks and GWAS order are fixed to match the reference figure.
    # Annotations are reversed so the first biological annotation appears at top.
    ###########################################################################

    plot_dt_full[
        ,
        trait_plot :=
            factor(
                trait_plot,
                levels =
                    regular_trait_order
            )
    ]


    plot_dt_full[
        ,
        gwas_group :=
            factor(
                gwas_group,
                levels =
                    gwas_group_order
            )
    ]


    plot_dt_full[
        ,
        annotation :=
            factor(
                annotation,
                levels =
                    rev(
                        annotation_order
                    )
            )
    ]


    ###########################################################################
    # Annotation labels
    ###########################################################################

    annotation_label_info <- annotation_info[
        match(
            annotation_order,
            annotation
        )
    ]


    annotation_labels <- setNames(
        paste0(
            annotation_label_info$annotation_display,
            "\n(n=",
            format(
                annotation_label_info$n_snps,
                big.mark = ",",
                scientific = FALSE,
                trim = TRUE
            ),
            ")"
        ),
        annotation_label_info$annotation
    )


    ###########################################################################
    # Narrow separator position
    #
    # Place a single thin white line below LHb, so MHb/LHb are visually
    # separated from the finer subpopulations without adding an extra row.
    ###########################################################################

    separator_yintercept <- NULL

    lhb_y_position <- unique(
        as.numeric(
            plot_dt_full[
                annotation_display == "LHb",
                annotation
            ]
        )
    )

    if (length(lhb_y_position) == 1L) {
        separator_yintercept <- lhb_y_position - 0.5
    }


    ###########################################################################
    # Sort plotting data
    ###########################################################################

    setorder(
        plot_dt_full,
        gwas_group,
        trait_plot,
        annotation
    )


    ###########################################################################
    # Save plot data
    ###########################################################################

    fwrite(
        plot_dt_full,
        output_plot_data,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )


    ###########################################################################
    # Store combined plotting data
    #
    # Convert factors back to character before combining DAR sets.
    ###########################################################################

    plot_dt_store <- copy(
        plot_dt_full
    )


    plot_dt_store[
        ,
        `:=`(
            trait =
                as.character(
                    trait
                ),
            trait_plot =
                as.character(
                    trait_plot
                ),
            gwas_group =
                as.character(
                    gwas_group
                ),
            annotation =
                as.character(
                    annotation
                )
        )
    ]


    all_plot_data[[dar_set_i]] <- plot_dt_store


    ###########################################################################
    # Subtitle
    ###########################################################################

    subtitle_text <- paste0(
        "Color: coefficient Z-score; ",
        significance_description,
        ". Color scale clipped at ",
        color_limit,
        "; n_snps ≥ ",
        format(
            minimum_n_snps,
            big.mark = ","
        )
    )

    if (
        show_only_significant
    ) {

        subtitle_text <- paste0(
            subtitle_text,
            "; nonsignificant cells shown in gray"
        )
    }


    ###########################################################################
    # Caption
    ###########################################################################

    caption_text <- paste0(
        "Red = positive coefficient, blue = negative coefficient, ",
        "white = coefficient Z = 0. ",
        "Stars indicate positive enrichment after BH-FDR correction ",
        "across GWAS traits separately within each cell type ",
        "and within the current DAR set."
    )


    ###########################################################################
    # Plot title
    ###########################################################################

    plot_title <- paste0(
        "S-LDSC across GWAS traits: ",
        dar_set_titles[
            dar_set_i
        ]
    )


    ###########################################################################
    # Draw heatmap
    #
    # Layout follows the reference figure:
    #   X = GWAS traits
    #   top facet strips = P / Psychiatric / Substance Use
    #   Y = DAR annotations
    ###########################################################################

    p <- ggplot(
        plot_dt_full,
        aes(
            x =
                trait_plot,
            y =
                annotation,
            fill =
                fill_value
        )
    ) +

        geom_tile(
            color = "white",
            linewidth = 0.35
        )

    if (!is.null(separator_yintercept)) {

        # White gap
        p <- p +
            geom_hline(
                yintercept = separator_yintercept,
                color = "white",
                linewidth = whole_hb_separator_linewidth
            ) +

            # Thin black rule at the upper edge of the gap
            geom_hline(
                yintercept = separator_yintercept + 0.10,
                color = "black",
                linewidth = whole_hb_separator_border_linewidth
            ) +

            # Thin black rule at the lower edge of the gap
            geom_hline(
                yintercept = separator_yintercept - 0.10,
                color = "black",
                linewidth = whole_hb_separator_border_linewidth
            )
    }

    p <- p +

        geom_text(
            aes(
                label =
                    significance_label
            ),
            size = 3.4,
            fontface = "bold",
            na.rm = TRUE
        ) +

        #######################################################################
        # Exact three GWAS blocks from the reference figure
        #######################################################################

        facet_grid(
            cols = vars(gwas_group),
            scales = "free_x",
            space = "free_x",
            drop = TRUE
        ) +

        scale_x_discrete(
            drop = TRUE
        ) +

        scale_y_discrete(
            labels =
                annotation_labels,
            drop = FALSE
        ) +

        scale_fill_gradient2(
            low = "#3B4CC0",
            mid = "white",
            high = "#B40426",
            midpoint = 0,
            limits = c(
                -color_limit,
                color_limit
            ),
            breaks = seq(
                -color_limit,
                color_limit,
                by = 1
            ),
            oob = squish,
            na.value = "grey90",
            name = "Z-score"
        ) +

        labs(
            title =
                plot_title,

            subtitle =
                subtitle_text,

            x =
                "GWAS trait",

            y =
                paste0(
                    dar_set_titles[
                        dar_set_i
                    ],
                    " annotation and number of common SNPs"
                ),

            caption =
                caption_text
        ) +

        theme_bw(
            base_size = 11
        ) +

        theme(

            panel.grid =
                element_blank(),

            axis.text.x =
                element_text(
                    angle = 90,
                    hjust = 1,
                    vjust = 0.5,
                    size = 9
                ),

            axis.text.y =
                element_text(
                    size = 9
                ),

            axis.title =
                element_text(
                    face = "bold"
                ),

            strip.background.x =
                element_rect(
                    fill = "grey85",
                    color = "grey40",
                    linewidth = 0.7
                ),

            strip.text.x =
                element_text(
                    size = 11
                ),

            panel.spacing.x =
                grid::unit(
                    0.10,
                    "in"
                ),

            plot.title =
                element_text(
                    face = "bold",
                    size = 14
                ),

            plot.subtitle =
                element_text(
                    size = 9.5
                ),

            plot.caption =
                element_text(
                    hjust = 0,
                    size = 8
                ),

            legend.title =
                element_text(
                    size = 10
                )
        )


    ###########################################################################
    # Dimensions
    ###########################################################################

    n_regular_traits <- uniqueN(
        plot_dt$trait_nickname
    )


    n_plot_rows <- n_regular_traits

    n_annotations <- uniqueN(
        plot_dt$annotation
    )


    plot_width <- 
        0.4 *
        n_regular_traits +
        2



    plot_height <- 
        max(7,0.36 *
        n_annotations +
        1.5)


    ###########################################################################
    # Save PDF
    ###########################################################################

    ggsave(
        filename =
            output_pdf,
        plot =
            p,
        width =
            plot_width,
        height =
            plot_height,
        units =
            "in"
    )


    ###########################################################################
    # Save PNG
    ###########################################################################

    ggsave(
        filename =
            output_png,
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


    ###########################################################################
    # Per-DAR-set summary
    ###########################################################################

    message("")
    message(
        "DAR set:              ",
        dar_set_i
    )

    message(
        "Heatmap PDF:          ",
        output_pdf
    )

    message(
        "Heatmap PNG:          ",
        output_png
    )

    message(
        "Plot data:            ",
        output_plot_data
    )


    message(
        "GWAS traits:          ",
        n_regular_traits
    )

    message(
        "Annotations:          ",
        n_annotations
    )
}


###############################################################################
# Check output
###############################################################################

if (
    length(
        all_plot_data
    ) == 0L
) {

    stop(
        "No DAR-set heatmaps were generated."
    )
}


###############################################################################
# Combined heatmap plotting data
###############################################################################

combined_plot_data <- rbindlist(
    all_plot_data,
    use.names = TRUE,
    fill = TRUE
)


combined_plot_data[
    ,
    dar_set_order :=
        match(
            dar_set,
            expected_dar_sets
        )
]


setorder(
    combined_plot_data,
    dar_set_order,
    trait,
    annotation
)


combined_plot_data[
    ,
    dar_set_order := NULL
]


combined_plot_data_file <- file.path(
    output_dir,
    "DAR_sldsc_heatmap_plot_data_all_sets.tsv"
)


fwrite(
    combined_plot_data,
    combined_plot_data_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Publication-ready S-LDSC supplementary table
#
# IMPORTANT:
# This table uses EXACTLY the same analysis subset and significance definition
# as the heatmaps above:
#
#   1. Only GWAS traits represented in gwas_trait_order are retained.
#   2. Only annotations with n_snps >= minimum_n_snps are retained.
#   3. Cell-type names use annotation_display after mapping.
#   4. Within EACH DAR set separately (open / closed / all), BH-FDR is applied
#      separately within EACH cell type across the displayed GWAS traits.
#   5. Significant positive enrichment is defined as:
#
#          coefficient_z > 0
#          AND
#          coefficient_FDR_one_sided_within_cell_type < 0.05
#
# The full supplementary table contains only the manuscript-facing columns
# selected below.  Plotting-only and repeated analysis-metadata fields are not
# included.
###############################################################################


###############################################################################
# Required source columns for the supplementary table
###############################################################################

supp_table_source_columns <- c(
    "trait_plot",
    "dar_set",
    "annotation_display",
    "n_snps",
    "Prop._SNPs",
    "Prop._h2",
    "Prop._h2_std_error",
    "Enrichment",
    "Enrichment_std_error",
    "Coefficient",
    "Coefficient_std_error",
    "coefficient_z_plot",
    "coefficient_p_one_sided",
    "coefficient_FDR_one_sided_within_cell_type",
    "positive_enrichment_significant"
)


missing_supp_columns <- setdiff(
    supp_table_source_columns,
    names(combined_plot_data)
)


if (length(missing_supp_columns) > 0L) {

    stop(
        paste0(
            "Cannot create the final S-LDSC supplementary table because these ",
            "required columns are missing from combined_plot_data: ",
            paste(
                missing_supp_columns,
                collapse = ", "
            )
        )
    )
}


###############################################################################
# Create a working table
#
# Keep the significance indicator temporarily so that the FDR<0.05 table can
# be generated.  It will NOT be included in the final 14-column output.
###############################################################################

ldsc_supp_work <- copy(
    combined_plot_data[
        ,
        ..supp_table_source_columns
    ]
)


###############################################################################
# Rename columns to manuscript-facing names
###############################################################################

setnames(
    ldsc_supp_work,
    old = c(
        "trait_plot",
        "annotation_display",
        "coefficient_z_plot",
        "coefficient_FDR_one_sided_within_cell_type"
    ),
    new = c(
        "GWAS_label",
        "cell_type",
        "coefficient_z",
        "coefficient_FDR"
    )
)


###############################################################################
# Explicit final column set requested for the supplementary table
###############################################################################

supp_table_columns <- c(
    "GWAS_label",
    "dar_set",
    "cell_type",
    "n_snps",
    "Prop._SNPs",
    "Prop._h2",
    "Prop._h2_std_error",
    "Enrichment",
    "Enrichment_std_error",
    "Coefficient",
    "Coefficient_std_error",
    "coefficient_z",
    "coefficient_p_one_sided",
    "coefficient_FDR"
)


###############################################################################
# Check uniqueness before writing
#
# Exactly one row should exist per:
#   GWAS × DAR set × cell type
###############################################################################

duplicate_supp_rows <- ldsc_supp_work[
    ,
    .N,
    by = .(
        GWAS_label,
        dar_set,
        cell_type
    )
][
    N > 1L
]


if (nrow(duplicate_supp_rows) > 0L) {

    stop(
        paste0(
            "Duplicate rows detected in the supplementary S-LDSC table:\n",
            paste(
                capture.output(
                    print(
                        duplicate_supp_rows
                    )
                ),
                collapse = "\n"
            )
        )
    )
}


###############################################################################
# Add ordering variables temporarily
###############################################################################

ldsc_supp_work[
    ,
    dar_set_order :=
        match(
            dar_set,
            expected_dar_sets
        )
]


ldsc_supp_work[
    ,
    gwas_order :=
        match(
            GWAS_label,
            gwas_trait_order
        )
]


ldsc_supp_work[
    ,
    cell_type_order :=
        match(
            cell_type,
            custom_annotation_display_order
        )
]


# Put cell types not listed in custom_annotation_display_order at the end.
ldsc_supp_work[
    is.na(cell_type_order),
    cell_type_order :=
        length(
            custom_annotation_display_order
        ) +
        frank(
            cell_type,
            ties.method = "dense"
        )
]


setorder(
    ldsc_supp_work,
    dar_set_order,
    gwas_order,
    cell_type_order,
    cell_type
)


###############################################################################
# Save significant positive-enrichment rows before dropping the temporary
# significance field.
###############################################################################

ldsc_sig_work <- copy(
    ldsc_supp_work[
        positive_enrichment_significant == TRUE
    ]
)


###############################################################################
# Final 14-column supplementary table
###############################################################################

ldsc_supp_table <- copy(
    ldsc_supp_work[
        ,
        ..supp_table_columns
    ]
)


ldsc_sig_table <- copy(
    ldsc_sig_work[
        ,
        ..supp_table_columns
    ]
)


###############################################################################
# Output paths
###############################################################################

ldsc_supp_table_file <- file.path(
    output_dir,
    "DAR_sldsc_supplementary_table.tsv"
)


ldsc_sig_table_file <- file.path(
    output_dir,
    "DAR_sldsc_supplementary_table_FDR05.tsv"
)


###############################################################################
# Write complete supplementary table
###############################################################################

fwrite(
    ldsc_supp_table,
    ldsc_supp_table_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Write significant positive-enrichment results only
#
# Same significance definition as heatmap stars:
#   coefficient_z > 0 AND coefficient_FDR < 0.05
###############################################################################

fwrite(
    ldsc_sig_table,
    ldsc_sig_table_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Supplementary-table summary
###############################################################################

message("")
message("Supplementary S-LDSC table:")
message("  All results:          ", ldsc_supp_table_file)
message("  Significant results:  ", ldsc_sig_table_file)
message("  Rows in full table:    ", nrow(ldsc_supp_table))
message("  FDR < 0.05 positive:   ", nrow(ldsc_sig_table))
message("  Columns in final table: ", paste(supp_table_columns, collapse = ", "))


###############################################################################
# Final summary
###############################################################################

message("")
message("============================================================")
message("ALL S-LDSC HEATMAPS FINISHED")
message("============================================================")
message("")

message(
    "Root output directory: ",
    output_dir
)

message(
    "Annotation mapping:    ",
    output_annotation_map
)

message(
    "GWAS name mapping:     ",
    output_gwas_name_map
)

message(
    "Combined plot data:    ",
    combined_plot_data_file
)

message(
    "Supplementary table:   ",
    ldsc_supp_table_file
)

message(
    "Significant table:     ",
    ldsc_sig_table_file
)


message("")
message(
    "Heatmaps generated for: ",
    paste(
        names(
            all_plot_data
        ),
        collapse = ", "
    )
)



message("")
message("Done.")