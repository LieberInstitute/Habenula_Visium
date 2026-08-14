#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})


###############################################################################
# Paths
###############################################################################

project_dir <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Habenula_Visium/processed-data/10_HD_bin_level/LDSC"
)


###############################################################################
# S-LDSC task file
#
# Expected columns:
#
#   trait
#   dar_set
#   annotation
#   sumstats
#   samp_prev
#   pop_prev
#   result_prefix
###############################################################################

task_file <- file.path(
    project_dir,
    "03_GWAS",
    "sldsc_tasks.tsv"
)


###############################################################################
# Combined DAR manifest
#
# Expected to contain all three DAR sets:
#
#   open
#   closed
#   all
###############################################################################

dar_manifest_file <- file.path(
    project_dir,
    "01_DAR_beds_hg19",
    "DAR_hg19_manifest_all_sets.tsv"
)


###############################################################################
# LD-score directory
#
# Expected structure:
#
# 02_DAR_ldscores_hg19/
# ├── open/
# │   ├── Astrocyte/
# │   └── ...
# ├── closed/
# │   ├── Astrocyte/
# │   └── ...
# └── all/
#     ├── Astrocyte/
#     └── ...
###############################################################################

ldscore_dir <- file.path(
    project_dir,
    "02_DAR_ldscores_hg19"
)


###############################################################################
# S-LDSC results
###############################################################################

result_dir <- file.path(
    project_dir,
    "04_sldsc_results"
)


###############################################################################
# Gathered output
###############################################################################

output_file <- file.path(
    result_dir,
    "all_GWAS_DAR_sldsc_summary.tsv"
)


###############################################################################
# Task-level gathering status
###############################################################################

status_file <- file.path(
    result_dir,
    "all_GWAS_DAR_sldsc_gather_status.tsv"
)


dir.create(
    result_dir,
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
# Helper: normalize column names
###############################################################################

normalize_column_name <- function(x) {

    y <- tolower(
        trimws(x)
    )

    y <- gsub(
        "[^a-z0-9]+",
        "_",
        y
    )

    y <- gsub(
        "_+",
        "_",
        y
    )

    y <- gsub(
        "^_|_$",
        "",
        y
    )

    y
}


###############################################################################
# Helper: identify a column in LDSC output
###############################################################################

find_unique_column <- function(
    column_names,
    exact_names = character(),
    regex = NULL,
    required = TRUE
) {

    normalized_columns <- normalize_column_name(
        column_names
    )

    ###########################################################################
    # Try exact names first
    ###########################################################################

    if (length(exact_names) > 0L) {

        normalized_exact <- normalize_column_name(
            exact_names
        )

        exact_index <- match(
            normalized_exact,
            normalized_columns,
            nomatch = 0L
        )

        exact_index <- exact_index[
            exact_index > 0L
        ]

        if (length(exact_index) > 0L) {

            return(
                column_names[
                    exact_index[1]
                ]
            )
        }
    }


    ###########################################################################
    # Then try regex
    ###########################################################################

    if (!is.null(regex)) {

        matched <- grep(
            regex,
            column_names,
            ignore.case = TRUE,
            value = TRUE
        )

        if (length(matched) == 1L) {

            return(
                matched
            )
        }

        if (length(matched) > 1L) {

            stop(
                "Multiple columns matched ",
                regex,
                ": ",
                paste(
                    matched,
                    collapse = ", "
                )
            )
        }
    }


    ###########################################################################
    # Not found
    ###########################################################################

    if (required) {

        stop(
            "Could not identify the requested column from: ",
            paste(
                column_names,
                collapse = ", "
            )
        )
    }

    NA_character_
}


###############################################################################
# Helper: read total M / M_5_50 counts
#
# Each annotation is now uniquely defined by:
#
#   dar_set + annotation
#
# Example:
#
#   open + Astrocyte
#   closed + Astrocyte
#   all + Astrocyte
###############################################################################

read_ldscore_m_total <- function(
    dar_set,
    annotation,
    suffix
) {

    prefix <- file.path(
        ldscore_dir,
        dar_set,
        annotation,
        annotation
    )

    files <- paste0(
        prefix,
        ".",
        seq_len(22L),
        suffix
    )


    ###########################################################################
    # Check files
    ###########################################################################

    missing <- files[
        !file.exists(files) |
        is.na(file.info(files)$size) |
        file.info(files)$size == 0
    ]

    if (length(missing) > 0L) {

        warning(
            "Missing SNP-count files for ",
            dar_set,
            " / ",
            annotation,
            ": ",
            paste(
                missing,
                collapse = ", "
            )
        )

        return(
            NA_real_
        )
    }


    ###########################################################################
    # Read each chromosome
    ###########################################################################

    values <- vapply(
        files,
        function(file) {

            x <- scan(
                file,
                what = numeric(),
                quiet = TRUE
            )

            if (length(x) == 0L) {

                return(
                    NA_real_
                )
            }

            sum(x)
        },
        numeric(1)
    )


    ###########################################################################
    # Check completeness
    ###########################################################################

    if (anyNA(values)) {

        warning(
            "Incomplete SNP count for ",
            dar_set,
            " / ",
            annotation,
            " using ",
            suffix,
            "."
        )

        return(
            NA_real_
        )
    }


    ###########################################################################
    # Sum chromosomes 1-22
    ###########################################################################

    sum(values)
}


###############################################################################
# Check input files
###############################################################################

for (file in c(
    task_file,
    dar_manifest_file
)) {

    if (!file.exists(file)) {

        stop(
            "Missing input file: ",
            file
        )
    }

    if (
        is.na(file.info(file)$size) ||
        file.info(file)$size == 0
    ) {

        stop(
            "Input file is empty: ",
            file
        )
    }
}


###############################################################################
# Read inputs
###############################################################################

tasks <- fread(
    task_file,
    sep = "\t",
    quote = "",
    na.strings = c(
        "",
        "NA"
    )
)


dar_manifest <- fread(
    dar_manifest_file,
    sep = "\t",
    quote = "",
    na.strings = c(
        "",
        "NA"
    )
)


###############################################################################
# Validate task-file columns
###############################################################################

require_columns(
    tasks,
    c(
        "trait",
        "dar_set",
        "annotation",
        "result_prefix"
    ),
    "S-LDSC task file"
)


###############################################################################
# Validate DAR-manifest columns
###############################################################################

require_columns(
    dar_manifest,
    c(
        "dar_set",
        "annotation"
    ),
    "DAR manifest"
)


###############################################################################
# Check task file
###############################################################################

if (nrow(tasks) == 0L) {

    stop(
        "The S-LDSC task file contains no tasks."
    )
}


###############################################################################
# Standardize DAR-set variables
###############################################################################

tasks[
    ,
    dar_set := trimws(
        as.character(
            dar_set
        )
    )
]


dar_manifest[
    ,
    dar_set := trimws(
        as.character(
            dar_set
        )
    )
]


###############################################################################
# Validate DAR-set values
###############################################################################

expected_dar_sets <- c(
    "open",
    "closed",
    "all"
)


unexpected_task_sets <- setdiff(
    unique(
        tasks$dar_set
    ),
    expected_dar_sets
)


if (length(unexpected_task_sets) > 0L) {

    stop(
        "Unexpected dar_set values in task file: ",
        paste(
            unexpected_task_sets,
            collapse = ", "
        )
    )
}


unexpected_manifest_sets <- setdiff(
    unique(
        dar_manifest$dar_set
    ),
    expected_dar_sets
)


if (length(unexpected_manifest_sets) > 0L) {

    stop(
        "Unexpected dar_set values in DAR manifest: ",
        paste(
            unexpected_manifest_sets,
            collapse = ", "
        )
    )
}


###############################################################################
# Check duplicate S-LDSC tasks
#
# Each:
#
# trait × dar_set × annotation
#
# should occur exactly once.
###############################################################################

duplicate_tasks <- tasks[
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


if (nrow(duplicate_tasks) > 0L) {

    stop(
        "Duplicate trait × dar_set × annotation tasks detected:\n",
        paste(
            capture.output(
                print(
                    duplicate_tasks
                )
            ),
            collapse = "\n"
        )
    )
}


###############################################################################
# Count SNPs in each custom annotation
#
# n_snps:
#   Common reference SNPs from .l2.M_5_50 files.
#
# n_snps_all:
#   Reference SNPs from .l2.M files.
#
# Important:
#
# Astrocyte-open,
# Astrocyte-closed,
# Astrocyte-all
#
# are THREE different custom annotations here.
###############################################################################

annotations_for_counts <- unique(
    tasks[
        ,
        .(
            dar_set,
            annotation
        )
    ]
)


###############################################################################
# Order DAR sets
###############################################################################

annotations_for_counts[
    ,
    dar_set_order := match(
        dar_set,
        expected_dar_sets
    )
]


setorder(
    annotations_for_counts,
    dar_set_order,
    annotation
)


annotations_for_counts[
    ,
    dar_set_order := NULL
]


###############################################################################
# Read SNP counts
###############################################################################

snp_counts <- rbindlist(
    lapply(
        seq_len(
            nrow(
                annotations_for_counts
            )
        ),
        function(i) {

            dar_set_i <- as.character(
                annotations_for_counts$dar_set[i]
            )

            annotation_i <- as.character(
                annotations_for_counts$annotation[i]
            )


            message(
                "Reading SNP counts: ",
                dar_set_i,
                " / ",
                annotation_i
            )


            data.table(
                dar_set = dar_set_i,
                annotation = annotation_i,

                n_snps = read_ldscore_m_total(
                    dar_set = dar_set_i,
                    annotation = annotation_i,
                    suffix = ".l2.M_5_50"
                ),

                n_snps_all = read_ldscore_m_total(
                    dar_set = dar_set_i,
                    annotation = annotation_i,
                    suffix = ".l2.M"
                )
            )
        }
    ),
    use.names = TRUE,
    fill = TRUE
)


###############################################################################
# Check SNP-count uniqueness
###############################################################################

duplicate_snp_counts <- snp_counts[
    ,
    .N,
    by = .(
        dar_set,
        annotation
    )
][
    N > 1L
]


if (nrow(duplicate_snp_counts) > 0L) {

    stop(
        "Duplicate SNP-count rows detected for DAR set × annotation."
    )
}


###############################################################################
# Gather one custom annotation row per S-LDSC result
###############################################################################

result_list <- vector(
    mode = "list",
    length = nrow(tasks)
)


status_list <- vector(
    mode = "list",
    length = nrow(tasks)
)


###############################################################################
# Loop through task file
###############################################################################

for (i in seq_len(
    nrow(tasks)
)) {

    ###########################################################################
    # Task identifiers
    ###########################################################################

    trait_i <- as.character(
        tasks$trait[i]
    )

    dar_set_i <- as.character(
        tasks$dar_set[i]
    )

    annotation_i <- as.character(
        tasks$annotation[i]
    )

    result_prefix_i <- as.character(
        tasks$result_prefix[i]
    )


    ###########################################################################
    # Result file
    ###########################################################################

    result_file <- paste0(
        result_prefix_i,
        ".results"
    )


    status_i <- "success"

    message_i <- NA_character_


    message(
        "[",
        i,
        "/",
        nrow(tasks),
        "] ",
        trait_i,
        " / ",
        dar_set_i,
        " / ",
        annotation_i
    )


    ###########################################################################
    # Check result file
    ###########################################################################

    if (
        !file.exists(result_file) ||
        is.na(file.info(result_file)$size) ||
        file.info(result_file)$size == 0
    ) {

        status_i <- "missing_result"

        message_i <- result_file


        status_list[[i]] <- data.table(
            trait = trait_i,
            dar_set = dar_set_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }


    ###########################################################################
    # Read result file
    ###########################################################################

    result <- tryCatch(
        fread(
            result_file,
            check.names = FALSE,
            na.strings = c(
                "",
                "NA",
                "nan"
            )
        ),
        error = function(e) {
            e
        }
    )


    ###########################################################################
    # Read failure
    ###########################################################################

    if (inherits(
        result,
        "error"
    )) {

        status_i <- "read_error"

        message_i <- conditionMessage(
            result
        )


        status_list[[i]] <- data.table(
            trait = trait_i,
            dar_set = dar_set_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }


    ###########################################################################
    # Empty result table
    ###########################################################################

    if (nrow(result) == 0L) {

        status_i <- "empty_results_table"

        message_i <- paste0(
            "The .results file contained no rows."
        )


        status_list[[i]] <- data.table(
            trait = trait_i,
            dar_set = dar_set_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }


    ###########################################################################
    # Identify Category column
    ###########################################################################

    category_col <- tryCatch(
        find_unique_column(
            names(result),
            exact_names = "Category",
            regex = "^Category$"
        ),
        error = function(e) {
            e
        }
    )


    ###########################################################################
    # Category column missing
    ###########################################################################

    if (inherits(
        category_col,
        "error"
    )) {

        status_i <- "missing_category_column"

        message_i <- conditionMessage(
            category_col
        )


        status_list[[i]] <- data.table(
            trait = trait_i,
            dar_set = dar_set_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }


    ###########################################################################
    # Select the single custom annotation row
    #
    # The S-LDSC command uses:
    #
    #   --ref-ld-chr baselineLD,custom_annotation,pooled_background
    #
    # Therefore LDSC labels:
    #
    #   *_0 = baseline-LD annotation set
    #   *_1 = custom cell-type DAR annotation (foreground)
    #   *_2 = pooled same-DAR-set background annotation
    #
    # The foreground row therefore still ends in "_1". Its coefficient is
    # conditional on both baselineLD and the pooled background.
    ###########################################################################

    result[
        ,
        category_original := as.character(
            get(
                category_col
            )
        )
    ]


    custom_row <- result[
        grepl(
            "_1$",
            category_original
        )
    ]


    ###########################################################################
    # No custom annotation row
    ###########################################################################

    if (nrow(custom_row) == 0L) {

        status_i <- "custom_annotation_not_found"


        message_i <- paste0(
            "No Category ending in _1 was found. ",
            "Available categories include: ",
            paste(
                head(
                    unique(
                        result$category_original
                    ),
                    20L
                ),
                collapse = ", "
            )
        )


        status_list[[i]] <- data.table(
            trait = trait_i,
            dar_set = dar_set_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }


    ###########################################################################
    # More than one custom annotation row
    ###########################################################################

    if (nrow(custom_row) > 1L) {

        status_i <- "multiple_custom_annotation_rows"


        message_i <- paste0(
            "Expected exactly one Category ending in _1, but found ",
            nrow(custom_row),
            ": ",
            paste(
                unique(
                    custom_row$category_original
                ),
                collapse = ", "
            )
        )


        status_list[[i]] <- data.table(
            trait = trait_i,
            dar_set = dar_set_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }


    ###########################################################################
    # Add task identifiers
    #
    # dar_set is essential because annotation names repeat between:
    #
    # open
    # closed
    # all
    ###########################################################################

    custom_row[
        ,
        `:=`(
            trait = trait_i,
            dar_set = dar_set_i,
            annotation = annotation_i,
            ldsc_category = category_original,
            result_file = result_file
        )
    ]


    ###########################################################################
    # Store successful result
    ###########################################################################

    result_list[[i]] <- custom_row


    ###########################################################################
    # Store status
    ###########################################################################

    status_list[[i]] <- data.table(
        trait = trait_i,
        dar_set = dar_set_i,
        annotation = annotation_i,
        result_file = result_file,
        status = status_i,
        message = message_i
    )
}


###############################################################################
# Write task-level gathering status
###############################################################################

status <- rbindlist(
    status_list,
    use.names = TRUE,
    fill = TRUE
)


###############################################################################
# Order status
###############################################################################

status[
    ,
    dar_set_order := match(
        dar_set,
        expected_dar_sets
    )
]


setorder(
    status,
    status,
    dar_set_order,
    trait,
    annotation
)


status[
    ,
    dar_set_order := NULL
]


###############################################################################
# Write status
###############################################################################

fwrite(
    status,
    status_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Combine successful results
###############################################################################

results <- rbindlist(
    result_list,
    use.names = TRUE,
    fill = TRUE
)


if (nrow(results) == 0L) {

    stop(
        "No S-LDSC results could be gathered. See: ",
        status_file
    )
}


###############################################################################
# Identify coefficient z-score
###############################################################################

z_col <- find_unique_column(
    names(results),
    exact_names = c(
        "Coefficient_z-score",
        "Coefficient_z_score",
        "Coefficient z-score",
        "Coefficient z score"
    ),
    regex = "^Coefficient.*[Zz].*[Ss]core$"
)


###############################################################################
# Convert coefficient z-score to numeric
###############################################################################

results[
    ,
    coefficient_z := suppressWarnings(
        as.numeric(
            get(
                z_col
            )
        )
    )
]


if (all(
    is.na(
        results$coefficient_z
    )
)) {

    stop(
        "The coefficient z-score column was found, ",
        "but all values are NA: ",
        z_col
    )
}


###############################################################################
# Coefficient P values
#
# One-sided:
#
#   H1: custom annotation coefficient > 0
#
# Two-sided:
#
#   H1: custom annotation coefficient != 0
###############################################################################

results[
    ,
    coefficient_p_one_sided :=
        pnorm(
            coefficient_z,
            lower.tail = FALSE
        )
]


results[
    ,
    coefficient_p_two_sided :=
        2 * pnorm(
            abs(
                coefficient_z
            ),
            lower.tail = FALSE
        )
]


###############################################################################
# FDR within each:
#
#   GWAS × DAR set
#
# This treats:
#
#   open DAR LDSC
#   closed DAR LDSC
#   all DAR LDSC
#
# as three separate families of cell-type tests.
###############################################################################

results[
    ,
    coefficient_FDR_one_sided_by_trait_dar_set :=
        p.adjust(
            coefficient_p_one_sided,
            method = "BH"
        ),
    by = .(
        trait,
        dar_set
    )
]


results[
    ,
    coefficient_FDR_two_sided_by_trait_dar_set :=
        p.adjust(
            coefficient_p_two_sided,
            method = "BH"
        ),
    by = .(
        trait,
        dar_set
    )
]


###############################################################################
# Global FDR across ALL:
#
# GWAS × DAR set × cell-type tests
###############################################################################

results[
    ,
    coefficient_FDR_one_sided_global :=
        p.adjust(
            coefficient_p_one_sided,
            method = "BH"
        )
]


results[
    ,
    coefficient_FDR_two_sided_global :=
        p.adjust(
            coefficient_p_two_sided,
            method = "BH"
        )
]


###############################################################################
# Add DAR metadata
#
# IMPORTANT:
#
# Metadata must now be joined using:
#
#   dar_set + annotation
#
# not annotation alone.
###############################################################################

desired_metadata_columns <- c(
    "dar_set",
    "annotation",
    "cell_type",
    "source_resolution",
    "input_peaks",
    "uniquely_lifted_peaks",
    "dropped_unmapped_peaks",
    "dropped_multimapped_peaks",
    "dropped_unmapped_or_multimapped",
    "dropped_non_autosomal_mappings",
    "hg19_autosomal_intervals_before_reduce",
    "hg19_autosomal_intervals",
    "hg19_total_bp",
    "bed_file"
)


dar_metadata_columns <- intersect(
    desired_metadata_columns,
    names(
        dar_manifest
    )
)


###############################################################################
# Ensure merge keys exist
###############################################################################

required_metadata_keys <- c(
    "dar_set",
    "annotation"
)


if (!all(
    required_metadata_keys %in%
        dar_metadata_columns
)) {

    stop(
        "The DAR manifest must contain both ",
        "dar_set and annotation columns."
    )
}


###############################################################################
# Extract unique metadata
###############################################################################

dar_metadata <- unique(
    dar_manifest[
        ,
        ..dar_metadata_columns
    ]
)


###############################################################################
# Check metadata uniqueness
#
# One row per:
#
# dar_set × annotation
###############################################################################

duplicate_metadata <- dar_metadata[
    ,
    .N,
    by = .(
        dar_set,
        annotation
    )
][
    N > 1L
]


if (nrow(
    duplicate_metadata
) > 0L) {

    stop(
        "DAR metadata contains multiple rows for these ",
        "DAR set × annotation combinations:\n",
        paste(
            capture.output(
                print(
                    duplicate_metadata
                )
            ),
            collapse = "\n"
        )
    )
}


###############################################################################
# Merge DAR metadata
###############################################################################

n_before_merge <- nrow(
    results
)


results <- merge(
    results,
    dar_metadata,
    by = c(
        "dar_set",
        "annotation"
    ),
    all.x = TRUE,
    sort = FALSE
)


###############################################################################
# Ensure merge did not duplicate rows
###############################################################################

if (nrow(results) != n_before_merge) {

    stop(
        "Row count changed after merging DAR metadata: ",
        n_before_merge,
        " -> ",
        nrow(results)
    )
}


###############################################################################
# Check for unmatched metadata
###############################################################################

if (
    "cell_type" %in% names(results) &&
    anyNA(results$cell_type)
) {

    missing_metadata <- unique(
        results[
            is.na(cell_type),
            .(
                dar_set,
                annotation
            )
        ]
    )

    warning(
        "Some S-LDSC results could not be matched to DAR metadata:\n",
        paste(
            capture.output(
                print(
                    missing_metadata
                )
            ),
            collapse = "\n"
        )
    )
}


###############################################################################
# Add annotation SNP counts
#
# Again, merge using:
#
# dar_set + annotation
###############################################################################

n_before_snp_merge <- nrow(
    results
)


results <- merge(
    results,
    snp_counts,
    by = c(
        "dar_set",
        "annotation"
    ),
    all.x = TRUE,
    sort = FALSE
)


###############################################################################
# Ensure SNP-count merge did not duplicate rows
###############################################################################

if (
    nrow(results) !=
    n_before_snp_merge
) {

    stop(
        "Row count changed after merging SNP counts: ",
        n_before_snp_merge,
        " -> ",
        nrow(results)
    )
}


###############################################################################
# Warn about missing SNP counts
###############################################################################

if (
    anyNA(results$n_snps) ||
    anyNA(results$n_snps_all)
) {

    missing_snp_counts <- unique(
        results[
            is.na(n_snps) |
            is.na(n_snps_all),
            .(
                dar_set,
                annotation
            )
        ]
    )

    warning(
        "Missing SNP-count information for:\n",
        paste(
            capture.output(
                print(
                    missing_snp_counts
                )
            ),
            collapse = "\n"
        )
    )
}


###############################################################################
# Check one result per:
#
# trait × dar_set × annotation
###############################################################################

duplicate_results <- results[
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


if (nrow(
    duplicate_results
) > 0L) {

    stop(
        "Multiple gathered rows exist for some ",
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
}


###############################################################################
# Arrange output
#
# Desired DAR-set order:
#
# open
# closed
# all
###############################################################################

results[
    ,
    dar_set_order := match(
        dar_set,
        expected_dar_sets
    )
]


setorder(
    results,
    dar_set_order,
    trait,
    coefficient_p_one_sided,
    annotation
)


results[
    ,
    dar_set_order := NULL
]


###############################################################################
# Priority columns
###############################################################################

priority_columns <- c(

    # Analysis identifiers
    "trait",
    "dar_set",
    "annotation",

    # DAR metadata
    "cell_type",
    "source_resolution",

    # SNP counts
    "n_snps",
    "n_snps_all",
    "Prop._SNPs",

    # Main S-LDSC tests
    "coefficient_z",

    "coefficient_p_one_sided",

    "coefficient_FDR_one_sided_by_trait_dar_set",

    "coefficient_FDR_one_sided_global",

    "coefficient_p_two_sided",

    "coefficient_FDR_two_sided_by_trait_dar_set",

    "coefficient_FDR_two_sided_global",

    # LDSC identifiers
    "ldsc_category",

    # Result file
    "result_file",

    # Additional DAR metadata
    "input_peaks",
    "uniquely_lifted_peaks",
    "hg19_autosomal_intervals",
    "hg19_total_bp",
    "bed_file"
)


###############################################################################
# Retain only columns that actually exist
###############################################################################

priority_columns <- intersect(
    priority_columns,
    names(results)
)


###############################################################################
# Reorder columns
###############################################################################

setcolorder(
    results,
    c(
        priority_columns,
        setdiff(
            names(results),
            priority_columns
        )
    )
)


###############################################################################
# Write gathered summary
###############################################################################

fwrite(
    results,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Summary statistics
###############################################################################

n_success <- status[
    status == "success",
    .N
]


n_failed <- status[
    status != "success",
    .N
]


###############################################################################
# Print final summary
###############################################################################

message("")
message("============================================================")
message("S-LDSC GATHER FINISHED")
message("============================================================")
message("")

message(
    "Summary output: ",
    output_file
)

message(
    "Gather status:  ",
    status_file
)

message("")

message(
    "Tasks expected:        ",
    nrow(tasks)
)

message(
    "Results gathered:      ",
    nrow(results)
)

message(
    "Successful tasks:      ",
    n_success
)

message(
    "Failed/missing:        ",
    n_failed
)

message(
    "Traits:                ",
    uniqueN(
        results$trait
    )
)

message(
    "DAR sets:              ",
    paste(
        expected_dar_sets[
            expected_dar_sets %in%
                unique(results$dar_set)
        ],
        collapse = ", "
    )
)

message(
    "Cell annotations:      ",
    uniqueN(
        results$annotation
    )
)

message(
    "DAR set × annotations: ",
    uniqueN(
        results[
            ,
            paste(
                dar_set,
                annotation,
                sep = "::"
            )
        ]
    )
)

message(
    "SNP-count fields: n_snps = M_5_50; n_snps_all = M"
)


###############################################################################
# Results by DAR set
###############################################################################

message("")
message("Results by DAR set:")


print(
    results[
        ,
        .(
            n_results = .N,
            n_traits = uniqueN(
                trait
            ),
            n_annotations = uniqueN(
                annotation
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


###############################################################################
# Significant results summary
#
# One-sided FDR within each:
#
# trait × dar_set
###############################################################################

message("")
message(
    "FDR < 0.05 results by DAR set ",
    "(one-sided, within trait × DAR set):"
)


print(
    results[
        !is.na(
            coefficient_FDR_one_sided_by_trait_dar_set
        ) &
        coefficient_FDR_one_sided_by_trait_dar_set < 0.05,
        .N,
        by = dar_set
    ][
        match(
            dar_set,
            expected_dar_sets
        )
    ]
)


###############################################################################
# Non-success task status
###############################################################################

if (n_failed > 0L) {

    message("")
    message(
        "Non-success task statuses:"
    )

    print(
        status[
            status != "success",
            .N,
            by = .(
                status,
                dar_set
            )
        ][
            order(
                -N
            )
        ]
    )
}


###############################################################################
# Done
###############################################################################

message("")
message("Done.")