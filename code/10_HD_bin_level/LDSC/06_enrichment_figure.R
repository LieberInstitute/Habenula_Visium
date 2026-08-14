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
# Fisher meta-analysis row label
###############################################################################

meta_trait_label <- "Fisher meta-analysis"


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
#   only cells showing evidence for positive enrichment are colored.
#
# Fisher meta row follows the same rule using Fisher meta-analytic FDR.
###############################################################################

show_only_significant <- FALSE


###############################################################################
# Significance measure for ORDINARY GWAS rows
#
# For EACH DAR set separately (open / closed / all):
#   all raw coefficient_p_one_sided values across ALL GWAS traits × ALL
#   displayed annotations are adjusted together using Benjamini-Hochberg FDR.
#
# There is NO adjustment by GWAS and NO adjustment by cell type.
#
# Fisher meta-analysis:
#   raw coefficient_p_one_sided -> Fisher P per annotation -> BH-FDR across
#   annotations within the CURRENT DAR set.
###############################################################################

significance_column <- "coefficient_FDR_one_sided_within_dar_set"


###############################################################################
# Heatmap color limit
#
# Ordinary rows:
#   coefficient Z-score
#
# Fisher row:
#   nonnegative normal-equivalent Z converted from Fisher FDR.
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

    # Fisher meta-analysis always uses this:
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
# It is calculated inside the DAR-set loop AFTER the minimum_n_snps filter,
# so each plotted heatmap (open / closed / all) has its own FDR family:
# all displayed GWAS × annotation S-LDSC P values in that DAR set.
###############################################################################


###############################################################################
# Read GWAS name mapping
#
# Input S-LDSC "trait" is expected to match gwas_info.csv "nickname".
# For plotting/output, replace it with "manuscript_name".
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


dt[
    ,
    trait_nickname :=
        trait
]


dt[
    ,
    trait_mapped :=
        trait_nickname %chin%
        gwas_info$nickname
]


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


unmapped_gwas_traits <- sort(
    unique(
        dt[
            trait_mapped == FALSE,
            trait_nickname
        ]
    )
)


if (
    length(
        unmapped_gwas_traits
    ) > 0L
) {

    warning(
        paste0(
            "The following GWAS traits were not found in gwas_info.csv ",
            "and retain their original names:\n",
            paste(
                unmapped_gwas_traits,
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}


gwas_name_audit <- unique(
    dt[
        ,
        .(
            trait_nickname,
            trait,
            trait_mapped
        )
    ]
)


setorder(
    gwas_name_audit,
    trait_nickname
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
# Ordinary-row BH-FDR is calculated later within each DAR set.
# Fisher meta-analysis also uses these ORIGINAL one-sided P values.
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
# Prevent collision with the artificial Fisher row
###############################################################################

if (
    meta_trait_label %in%
    dt$trait
) {

    stop(
        "A real GWAS trait already has the reserved name: ",
        meta_trait_label
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
    "Ordinary-row significance column: ",
    significance_column
)

message(
    "Fisher meta-analysis input: coefficient_p_one_sided; ",
    "display significance: BH-FDR across annotations"
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

    # Medial habenula
    "MHb",
    "MHb_A",
    "MHb_B",
    "MHb_C",
    "MHb_D",

    # Lateral habenula
    "LHb",
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


###############################################################################
# Significance description for ordinary rows
###############################################################################

significance_description <- paste0(
    "* FDR<0.05, ",
    "** FDR<0.01, ",
    "*** FDR<0.001; ",
    "BH-FDR across all GWAS × annotations within this DAR set"
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

all_fisher_results <- list()


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
            "_with_Fisher_meta.pdf"
        )
    )


    output_png <- file.path(
        dar_output_dir,
        paste0(
            "DAR_sldsc_coefficient_z_heatmap_",
            dar_set_i,
            "_with_Fisher_meta.png"
        )
    )


    output_plot_data <- file.path(
        dar_output_dir,
        paste0(
            "DAR_sldsc_heatmap_plot_data_",
            dar_set_i,
            "_with_Fisher_meta.tsv"
        )
    )


    output_fisher <- file.path(
        dar_output_dir,
        paste0(
            "DAR_sldsc_Fisher_meta_",
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
    # BH-FDR for ordinary GWAS heatmap cells
    #
    # IMPORTANT:
    # The current plot_dt contains ONE DAR set only and has already passed the
    # minimum_n_snps display filter. Therefore p.adjust() below corrects ALL
    # displayed raw S-LDSC P values together across:
    #
    #   all GWAS traits × all displayed annotations
    #
    # within THIS DAR set. There is deliberately no `by = trait` and no
    # `by = annotation`.
    ###########################################################################

    plot_dt[
        ,
        coefficient_FDR_one_sided_within_dar_set :=
            p.adjust(
                coefficient_p_one_sided,
                method = "BH"
            )
    ]


    plot_dt[
        ,
        significance_value :=
            coefficient_FDR_one_sided_within_dar_set
    ]


    message(
        "[", dar_set_i, "] ordinary S-LDSC tests in FDR family: ",
        nrow(plot_dt)
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
    # Order ordinary GWAS traits
    ###########################################################################

    regular_trait_order <- plot_dt[
        ,
        .(
            ordering_score =
                max(
                    abs(
                        coefficient_z_clipped
                    ),
                    na.rm = TRUE
                )
        ),
        by = trait
    ][
        order(
            -ordering_score,
            trait
        ),
        trait
    ]


    ###########################################################################
    # Fisher meta-analysis
    #
    # For each annotation:
    #
    # X^2 = -2 * sum(log(P_i))
    #
    # df = 2k
    #
    # where P_i is coefficient_p_one_sided for GWAS i.
    #
    # NOTE:
    # Fisher meta-analysis uses ORIGINAL one-sided P values,
    # never FDR-adjusted values.
    ###########################################################################

    fisher_meta <- plot_dt[
        !is.na(
            coefficient_p_one_sided
        ) &
        is.finite(
            coefficient_p_one_sided
        ) &
        coefficient_p_one_sided >= 0 &
        coefficient_p_one_sided <= 1,
        {

            ###################################################################
            # Extract P values
            ###################################################################

            p_values <- as.numeric(
                coefficient_p_one_sided
            )


            ###################################################################
            # Numerical protection against log(0)
            ###################################################################

            p_values_safe <- pmax(
                p_values,
                .Machine$double.xmin
            )


            ###################################################################
            # Fisher statistic
            ###################################################################

            fisher_statistic <-
                -2 *
                sum(
                    log(
                        p_values_safe
                    )
                )


            ###################################################################
            # Degrees of freedom
            ###################################################################

            fisher_n_traits <-
                length(
                    p_values_safe
                )


            fisher_df <-
                2 *
                fisher_n_traits


            ###################################################################
            # Fisher combined P
            ###################################################################

            fisher_meta_p <-
                pchisq(
                    fisher_statistic,
                    df =
                        fisher_df,
                    lower.tail =
                        FALSE
                )


            list(

                fisher_meta_p =
                    fisher_meta_p,

                fisher_statistic =
                    fisher_statistic,

                fisher_df =
                    fisher_df,

                fisher_n_traits =
                    fisher_n_traits
            )
        },
        by = .(
            dar_set,
            annotation,
            annotation_original,
            annotation_display,
            annotation_mapped,
            source_resolution,
            n_snps
        )
    ]


    ###########################################################################
    # BH-FDR for Fisher meta-analysis P values
    #
    # Fisher itself is calculated from ORIGINAL one-sided P values.
    # Then its combined P values are adjusted across annotations within
    # the current DAR set.
    ###########################################################################

    fisher_meta[
        ,
        fisher_meta_FDR := {

            fdr_values <- rep(
                NA_real_,
                .N
            )

            valid_meta_p <- (
                !is.na(fisher_meta_p) &
                is.finite(fisher_meta_p) &
                fisher_meta_p >= 0 &
                fisher_meta_p <= 1
            )

            fdr_values[valid_meta_p] <- p.adjust(
                fisher_meta_p[valid_meta_p],
                method = "BH"
            )

            fdr_values
        }
    ]


    ###########################################################################
    # Convert Fisher FDR to a nonnegative normal-equivalent visualization Z
    ###########################################################################

    fisher_meta[
        ,
        fisher_FDR_for_z :=
            pmin(
                pmax(
                    fisher_meta_FDR,
                    .Machine$double.xmin
                ),
                1 -
                .Machine$double.eps
            )
    ]


    fisher_meta[
        ,
        fisher_z_equivalent :=
            pmax(
                qnorm(
                    fisher_FDR_for_z,
                    lower.tail = FALSE
                ),
                0
            )
    ]


    fisher_meta[
        ,
        fisher_FDR_for_z := NULL
    ]


    ###########################################################################
    # Make sure every displayed annotation received a Fisher result
    ###########################################################################

    missing_fisher_annotations <- setdiff(
        annotation_order,
        fisher_meta$annotation
    )


    if (
        length(
            missing_fisher_annotations
        ) > 0L
    ) {

        warning(
            paste0(
                "[",
                dar_set_i,
                "] Fisher meta-analysis could not be calculated for:\n",
                paste(
                    missing_fisher_annotations,
                    collapse = "\n"
                )
            ),
            call. = FALSE
        )
    }


    ###########################################################################
    # Number of ordinary GWAS traits contributing to each annotation
    ###########################################################################

    expected_trait_counts <- plot_dt[
        ,
        .(
            available_GWAS_traits =
                uniqueN(
                    trait
                )
        ),
        by = annotation
    ]


    fisher_meta <- merge(
        fisher_meta,
        expected_trait_counts,
        by = "annotation",
        all.x = TRUE,
        sort = FALSE
    )


    fisher_meta[
        ,
        all_available_traits_used :=
            fisher_n_traits ==
            available_GWAS_traits
    ]


    ###########################################################################
    # Warn if Fisher uses fewer traits due to missing P values
    ###########################################################################

    incomplete_fisher <- fisher_meta[
        all_available_traits_used ==
        FALSE
    ]


    if (
        nrow(
            incomplete_fisher
        ) > 0L
    ) {

        warning(
            paste0(
                "[",
                dar_set_i,
                "] Some Fisher meta-analyses use fewer GWAS traits ",
                "because coefficient_p_one_sided was missing:\n",
                paste(
                    capture.output(
                        print(
                            incomplete_fisher[
                                ,
                                .(
                                    annotation,
                                    fisher_n_traits,
                                    available_GWAS_traits
                                )
                            ]
                        )
                    ),
                    collapse = "\n"
                )
            ),
            call. = FALSE
        )
    }


    ###########################################################################
    # Save exact Fisher meta-analysis table
    ###########################################################################

    setorder(
        fisher_meta,
        fisher_meta_FDR,
        fisher_meta_p,
        annotation
    )


    fwrite(
        fisher_meta,
        output_fisher,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )


    all_fisher_results[[dar_set_i]] <- copy(
    fisher_meta
)


    ###########################################################################
    # Convert Fisher results into heatmap rows
    ###########################################################################

    fisher_plot <- copy(
        fisher_meta
    )


    fisher_plot[
        ,
        `:=`(

            trait =
                meta_trait_label,

            coefficient_z_plot =
                fisher_z_equivalent,

            significance_value =
                fisher_meta_FDR,

            positive_enrichment_significant =
                fisher_meta_FDR < 0.05,

            significance_label =
                make_significance_stars(
                    p =
                        fisher_meta_FDR,
                    require_positive =
                        FALSE
                ),

            row_type =
                "Fisher_meta"
        )
    ]


    ###########################################################################
    # Clip Fisher Z-equivalent
    ###########################################################################

    fisher_plot[
        ,
        coefficient_z_clipped :=
            pmin(
                fisher_z_equivalent,
                color_limit
            )
    ]


    ###########################################################################
    # Fisher fill
    ###########################################################################

    if (
        show_only_significant
    ) {

        fisher_plot[
            ,
            fill_value :=
                fifelse(
                    positive_enrichment_significant,
                    coefficient_z_clipped,
                    NA_real_
                )
        ]

    } else {

        fisher_plot[
            ,
            fill_value :=
                coefficient_z_clipped
        ]
    }


    ###########################################################################
    # Combine ordinary GWAS rows + Fisher row
    ###########################################################################

    plot_dt_full <- rbindlist(
        list(
            plot_dt,
            fisher_plot
        ),
        use.names = TRUE,
        fill = TRUE
    )


    ###########################################################################
    # Factor ordering
    #
    # ggplot discrete Y places the FIRST factor level at the bottom.
    #
    # Therefore:
    #
    #   Fisher meta-analysis = first level = bottom row
    #
    # Ordinary traits remain ordered by strongest absolute signal.
    ###########################################################################

    y_levels <- c(
        meta_trait_label,
        rev(
            regular_trait_order
        )
    )


    plot_dt_full[
        ,
        trait :=
            factor(
                trait,
                levels =
                    y_levels
            )
    ]


    plot_dt_full[
        ,
        annotation :=
            factor(
                annotation,
                levels =
                    annotation_order
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
    # Sort plotting data
    ###########################################################################

    setorder(
        plot_dt_full,
        trait,
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
        "GWAS rows: coefficient Z-score; ",
        significance_description,
        ". Fisher row: Fisher combination of original one-sided P values, ",
        "followed by BH-FDR across annotations; color and stars use Fisher FDR. ",
        "Color scale clipped at ",
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
        "Ordinary GWAS rows: red = positive coefficient, ",
        "blue = negative coefficient, white = coefficient Z = 0; ",
        "stars use BH-FDR within each DAR set across all GWAS × annotations. ",
        "The Fisher row combines original coefficient_p_one_sided values ",
        "across GWAS traits, then applies BH-FDR across annotations within ",
        "each DAR set. Because Fisher's method has no effect direction, ",
        "its FDR-derived visualization score is nonnegative."
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
    ###########################################################################

    p <- ggplot(
        plot_dt_full,
        aes(
            x =
                annotation,
            y =
                trait,
            fill =
                fill_value
        )
    ) +

        geom_tile(
            color = "white",
            linewidth = 0.35
        ) +

        geom_text(
            aes(
                label =
                    significance_label
            ),
            size = 3.4,
            fontface = "bold",
            na.rm = TRUE
        ) +

        scale_x_discrete(
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
                paste0(
                    dar_set_titles[
                        dar_set_i
                    ],
                    " annotation and number of common SNPs"
                ),

            y =
                "GWAS trait",

            caption =
                caption_text
        ) +

        coord_fixed(
            ratio = 0.8
        ) +

        theme_bw(
            base_size = 11
        ) +

        theme(

            panel.grid =
                element_blank(),

            axis.text.x =
                element_text(
                    angle = 45,
                    hjust = 1,
                    vjust = 1,
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
        plot_dt$trait
    )


    n_plot_rows <-
        n_regular_traits + 1L


    n_annotations <- uniqueN(
        plot_dt$annotation
    )


    plot_width <- max(
        10,
        0.75 *
        n_annotations +
        4
    )


    plot_height <- max(
        6,
        0.36 *
        n_plot_rows +
        2.5
    )


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
        "Fisher results:       ",
        output_fisher
    )

    message(
        "GWAS traits:          ",
        n_regular_traits
    )

    message(
        "Annotations:          ",
        n_annotations
    )

    message(
        "Fisher FDR < 0.05:    ",
        fisher_meta[
            fisher_meta_FDR <
            0.05,
            .N
        ]
    )

    message(
        "Fisher FDR < 0.01:    ",
        fisher_meta[
            fisher_meta_FDR <
            0.01,
            .N
        ]
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
    "DAR_sldsc_heatmap_plot_data_all_sets_with_Fisher_meta.tsv"
)


fwrite(
    combined_plot_data,
    combined_plot_data_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Combined Fisher meta-analysis table
###############################################################################

combined_fisher_results <- rbindlist(
    all_fisher_results,
    use.names = TRUE,
    fill = TRUE
)


combined_fisher_results[
    ,
    dar_set_order :=
        match(
            dar_set,
            expected_dar_sets
        )
]


setorder(
    combined_fisher_results,
    dar_set_order,
    fisher_meta_FDR,
    fisher_meta_p,
    annotation
)


combined_fisher_results[
    ,
    dar_set_order := NULL
]


combined_fisher_file <- file.path(
    output_dir,
    "DAR_sldsc_Fisher_meta_all_sets.tsv"
)


fwrite(
    combined_fisher_results,
    combined_fisher_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


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
    "Combined Fisher table: ",
    combined_fisher_file
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


###############################################################################
# Fisher summary
###############################################################################

message("")
message(
    "Fisher meta-analysis summary:"
)


print(
    combined_fisher_results[
        ,
        .(
            n_annotations =
                .N,

            fisher_FDR_lt_0_05 =
                sum(
                    fisher_meta_FDR <
                    0.05,
                    na.rm = TRUE
                ),

            fisher_FDR_lt_0_01 =
                sum(
                    fisher_meta_FDR <
                    0.01,
                    na.rm = TRUE
                ),

            fisher_FDR_lt_0_001 =
                sum(
                    fisher_meta_FDR <
                    0.001,
                    na.rm = TRUE
                )
        ),
        by = dar_set
    ][
        match(
            dar_set,
            expected_dar_sets
        )
    ]
)


message("")
message("Done.")