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

task_file <- file.path(
    project_dir,
    "03_GWAS",
    "sldsc_tasks.tsv"
)

dar_manifest_file <- file.path(
    project_dir,
    "01_DAR_beds_hg19",
    "DAR_hg19_manifest.tsv"
)

ldscore_dir <- file.path(
    project_dir,
    "02_DAR_ldscores_hg19"
)

result_dir <- file.path(
    project_dir,
    "04_sldsc_results"
)

output_file <- file.path(
    result_dir,
    "all_GWAS_DAR_sldsc_summary.tsv"
)

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
# Helpers
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


find_unique_column <- function(
    column_names,
    exact_names = character(),
    regex = NULL,
    required = TRUE
) {

    normalized_columns <- normalize_column_name(
        column_names
    )

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

    if (!is.null(regex)) {

        matched <- grep(
            regex,
            column_names,
            ignore.case = TRUE,
            value = TRUE
        )

        if (length(matched) == 1L) {
            return(matched)
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


read_ldscore_m_total <- function(annotation, suffix) {

    prefix <- file.path(ldscore_dir, annotation, annotation)
    files <- paste0(prefix, ".", seq_len(22L), suffix)
    missing <- files[!file.exists(files)]

    if (length(missing) > 0L) {
        warning(
            "Missing SNP-count files for ", annotation, ": ",
            paste(missing, collapse = ", ")
        )
        return(NA_real_)
    }

    values <- vapply(
        files,
        function(file) {
            x <- scan(file, what = numeric(), quiet = TRUE)
            if (length(x) == 0L) {
                return(NA_real_)
            }
            sum(x)
        },
        numeric(1)
    )

    if (anyNA(values)) {
        warning(
            "Incomplete SNP count for ", annotation,
            " using ", suffix, "."
        )
        return(NA_real_)
    }

    sum(values)
}


###############################################################################
# Read inputs
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

    if (file.info(file)$size == 0) {
        stop(
            "Input file is empty: ",
            file
        )
    }
}

tasks <- fread(
    task_file,
    sep = "\t",
    quote = "",
    na.strings = c("", "NA")
)

dar_manifest <- fread(
    dar_manifest_file,
    sep = "\t",
    quote = "",
    na.strings = c("", "NA")
)

require_columns(
    tasks,
    c(
        "trait",
        "annotation",
        "result_prefix"
    ),
    "S-LDSC task file"
)

require_columns(
    dar_manifest,
    "annotation",
    "DAR manifest"
)

if (nrow(tasks) == 0L) {
    stop("The S-LDSC task file contains no tasks.")
}

###############################################################################
# Count SNPs in each custom annotation
#
# n_snps:
#   Common reference SNPs from .l2.M_5_50 files.
#
# n_snps_all:
#   Reference SNPs from .l2.M files.
###############################################################################

annotations_for_counts <- sort(
    unique(as.character(tasks$annotation))
)

snp_counts <- rbindlist(
    lapply(
        annotations_for_counts,
        function(annotation_i) {
            data.table(
                annotation = annotation_i,
                n_snps = read_ldscore_m_total(
                    annotation_i,
                    ".l2.M_5_50"
                ),
                n_snps_all = read_ldscore_m_total(
                    annotation_i,
                    ".l2.M"
                )
            )
        }
    ),
    use.names = TRUE,
    fill = TRUE
)

###############################################################################
# Gather one custom annotation row per result file
###############################################################################

result_list <- vector(
    mode = "list",
    length = nrow(tasks)
)

status_list <- vector(
    mode = "list",
    length = nrow(tasks)
)

for (i in seq_len(nrow(tasks))) {

    trait_i <- as.character(
        tasks$trait[i]
    )

    annotation_i <- as.character(
        tasks$annotation[i]
    )

    result_prefix_i <- as.character(
        tasks$result_prefix[i]
    )

    result_file <- paste0(
        result_prefix_i,
        ".results"
    )

    status_i <- "success"
    message_i <- NA_character_

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
            na.strings = c("", "NA", "nan")
        ),
        error = function(e) {
            e
        }
    )

    if (inherits(result, "error")) {

        status_i <- "read_error"
        message_i <- conditionMessage(result)

        status_list[[i]] <- data.table(
            trait = trait_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }

    if (nrow(result) == 0L) {

        status_i <- "empty_results_table"
        message_i <- "The .results file contained no rows."

        status_list[[i]] <- data.table(
            trait = trait_i,
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

    if (inherits(category_col, "error")) {

        status_i <- "missing_category_column"
        message_i <- conditionMessage(category_col)

        status_list[[i]] <- data.table(
            trait = trait_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }

    ###########################################################################
    # Select the custom annotation row
    #
    # The S-LDSC command uses:
    #
    #   --ref-ld-chr baselineLD,custom_annotation
    #
    # LDSC therefore labels:
    #   *_0 = baselineLD categories
    #   *_1 = the single custom DAR annotation
    #
    # In the current output, the custom row is named "L2_1".
    # The task file records whether that row represents Astrocyte, MHb, etc.
    ###########################################################################

    result[
        ,
        category_original := as.character(
            get(category_col)
        )
    ]

    custom_row <- result[
        grepl(
            "_1$",
            category_original
        )
    ]

    if (nrow(custom_row) == 0L) {

        status_i <- "custom_annotation_not_found"

        message_i <- paste0(
            "No Category ending in _1 was found. Available categories include: ",
            paste(
                head(
                    unique(result$category_original),
                    20L
                ),
                collapse = ", "
            )
        )

        status_list[[i]] <- data.table(
            trait = trait_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }

    if (nrow(custom_row) > 1L) {

        status_i <- "multiple_custom_annotation_rows"

        message_i <- paste0(
            "Expected exactly one Category ending in _1, but found ",
            nrow(custom_row),
            ": ",
            paste(
                unique(custom_row$category_original),
                collapse = ", "
            )
        )

        status_list[[i]] <- data.table(
            trait = trait_i,
            annotation = annotation_i,
            result_file = result_file,
            status = status_i,
            message = message_i
        )

        next
    }

    ###########################################################################
    # Add task identifiers
    ###########################################################################

    custom_row[
        ,
        `:=`(
            trait = trait_i,
            annotation = annotation_i,
            ldsc_category = category_original,
            result_file = result_file
        )
    ]

    result_list[[i]] <- custom_row

    status_list[[i]] <- data.table(
        trait = trait_i,
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

setorder(
    status,
    status,
    trait,
    annotation
)

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

results[
    ,
    coefficient_z := suppressWarnings(
        as.numeric(
            get(z_col)
        )
    )
]

if (all(is.na(results$coefficient_z))) {
    stop(
        "The coefficient z-score column was found, but all values are NA: ",
        z_col
    )
}

###############################################################################
# Coefficient P values
#
# One-sided:
#   tests for a positive annotation coefficient.
#
# Two-sided:
#   tests for a nonzero annotation coefficient in either direction.
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
            abs(coefficient_z),
            lower.tail = FALSE
        )
]

###############################################################################
# FDR within each GWAS
###############################################################################

results[
    ,
    coefficient_FDR_one_sided_by_trait :=
        p.adjust(
            coefficient_p_one_sided,
            method = "BH"
        ),
    by = trait
]

results[
    ,
    coefficient_FDR_two_sided_by_trait :=
        p.adjust(
            coefficient_p_two_sided,
            method = "BH"
        ),
    by = trait
]

###############################################################################
# Global FDR across all GWAS × DAR tests
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
###############################################################################

desired_metadata_columns <- c(
    "annotation",
    "cell_type",
    "source_resolution",
    "input_peaks",
    "uniquely_lifted_peaks",
    "hg19_autosomal_intervals",
    "hg19_total_bp"
)

dar_metadata_columns <- intersect(
    desired_metadata_columns,
    names(dar_manifest)
)

if (!"annotation" %in% dar_metadata_columns) {
    stop(
        "The DAR manifest does not contain an annotation column."
    )
}

dar_metadata <- unique(
    dar_manifest[
        ,
        ..dar_metadata_columns
    ]
)

duplicate_metadata <- dar_metadata[
    ,
    .N,
    by = annotation
][
    N > 1L
]

if (nrow(duplicate_metadata) > 0L) {

    stop(
        "DAR metadata contains multiple rows for these annotations: ",
        paste(
            duplicate_metadata$annotation,
            collapse = ", "
        )
    )
}

n_before_merge <- nrow(results)

results <- merge(
    results,
    dar_metadata,
    by = "annotation",
    all.x = TRUE,
    sort = FALSE
)

if (nrow(results) != n_before_merge) {

    stop(
        "Row count changed after merging DAR metadata: ",
        n_before_merge,
        " -> ",
        nrow(results)
    )
}

###############################################################################
# Add annotation SNP counts
###############################################################################

n_before_snp_merge <- nrow(results)

results <- merge(
    results,
    snp_counts,
    by = "annotation",
    all.x = TRUE,
    sort = FALSE
)

if (nrow(results) != n_before_snp_merge) {
    stop(
        "Row count changed after merging SNP counts: ",
        n_before_snp_merge,
        " -> ",
        nrow(results)
    )
}

###############################################################################
# Arrange output
###############################################################################

setorder(
    results,
    trait,
    coefficient_p_one_sided,
    annotation
)

priority_columns <- c(
    "trait",
    "annotation",
    "cell_type",
    "source_resolution",
    "n_snps",
    "n_snps_all",
    "Prop._SNPs",
    "coefficient_z",
    "coefficient_p_one_sided",
    "coefficient_FDR_one_sided_by_trait",
    "coefficient_FDR_one_sided_global",
    "coefficient_p_two_sided",
    "coefficient_FDR_two_sided_by_trait",
    "coefficient_FDR_two_sided_global",
    "ldsc_category",
    "result_file"
)

priority_columns <- intersect(
    priority_columns,
    names(results)
)

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
# Write output
###############################################################################

fwrite(
    results,
    output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

###############################################################################
# Summary
###############################################################################

n_success <- status[
    status == "success",
    .N
]

n_failed <- status[
    status != "success",
    .N
]

message("")
message("Summary output: ", output_file)
message("Gather status:  ", status_file)
message("")
message("Tasks expected:  ", nrow(tasks))
message("Results gathered:", nrow(results))
message("Successful tasks:", n_success)
message("Failed/missing:  ", n_failed)
message("Traits:          ", uniqueN(results$trait))
message("Annotations:     ", uniqueN(results$annotation))
message("SNP-count fields: n_snps = M_5_50; n_snps_all = M")

if (n_failed > 0L) {

    message("")
    message("Non-success task statuses:")

    print(
        status[
            status != "success",
            .N,
            by = status
        ][
            order(-N)
        ]
    )
}