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

inventory_file <- file.path(
    project_dir,
    "03_GWAS",
    "GWAS_inventory.tsv"
)

manifest_file <- file.path(
    project_dir,
    "03_GWAS",
    "GWAS_manifest.tsv"
)

###############################################################################
# Fixed sample-size overrides
#
# Used only when the original GWAS does not contain an N column.
###############################################################################

fixed_N_overrides <- c(
    MDD2019 = 807553
)

###############################################################################
# Read inventory
###############################################################################

if (!file.exists(inventory_file)) {
    stop(
        "GWAS inventory does not exist: ",
        inventory_file
    )
}

inventory <- fread(
    inventory_file,
    sep = "\t",
    quote = "",
    na.strings = c("", "NA")
)

required_inventory_columns <- c(
    "trait",
    "sumstats_file",
    "file_name",
    "file_size_bytes",
    "header"
)

missing_inventory_columns <- setdiff(
    required_inventory_columns,
    names(inventory)
)

if (length(missing_inventory_columns) > 0L) {
    stop(
        "Inventory is missing columns: ",
        paste(
            missing_inventory_columns,
            collapse = ", "
        )
    )
}

if (nrow(inventory) == 0L) {
    stop("The GWAS inventory is empty.")
}

###############################################################################
# Require one selected input file per trait
###############################################################################

duplicate_traits <- inventory[
    ,
    .N,
    by = trait
][
    N > 1L
]

if (nrow(duplicate_traits) > 0L) {
    stop(
        "More than one GWAS file was found for these traits:\n",
        paste(
            capture.output(
                print(duplicate_traits)
            ),
            collapse = "\n"
        )
    )
}

setorder(
    inventory,
    trait
)

###############################################################################
# Helper functions
###############################################################################

normalize_name <- function(x) {

    y <- trimws(
        as.character(x)
    )

    # Remove VCF header prefixes, e.g. #CHROM -> CHROM.
    y <- sub(
        "^#+",
        "",
        y
    )

    y <- toupper(y)

    # Normalize punctuation:
    # P-value -> P_VALUE
    # Effect.N -> EFFECT_N
    y <- gsub(
        "[^A-Z0-9]+",
        "_",
        y
    )

    y <- gsub(
        "_+",
        "_",
        y
    )

    y <- gsub(
        "^_+|_+$",
        "",
        y
    )

    y
}


parse_header_line <- function(header_line) {

    header_line <- sub(
        "^\ufeff",
        "",
        header_line
    )

    header_line <- sub(
        "\r$",
        "",
        header_line
    )

    if (grepl(
        "\t",
        header_line,
        fixed = TRUE
    )) {

        columns <- strsplit(
            header_line,
            "\t",
            fixed = TRUE
        )[[1]]

        delimiter <- "tab"

    } else {

        columns <- strsplit(
            trimws(header_line),
            "[[:space:]]+"
        )[[1]]

        delimiter <- "whitespace"
    }

    columns <- trimws(columns)

    if (length(columns) > 0L) {
        columns[1] <- sub(
            "^#+",
            "",
            columns[1]
        )
    }

    list(
        columns = columns,
        delimiter = delimiter
    )
}


read_actual_header <- function(
    file_path,
    max_lines = 5000L
) {

    if (!file.exists(file_path)) {
        stop(
            "GWAS file does not exist: ",
            file_path
        )
    }

    compressed <- grepl(
        "\\.(gz|bgz)$",
        file_path,
        ignore.case = TRUE
    )

    connection <- if (compressed) {

        gzfile(
            file_path,
            open = "rt"
        )

    } else {

        base::file(
            file_path,
            open = "rt"
        )
    }

    on.exit(
        close(connection),
        add = TRUE
    )

    lines <- readLines(
        connection,
        n = max_lines,
        warn = FALSE
    )

    if (length(lines) == 0L) {
        stop(
            "Could not read any lines from: ",
            file_path
        )
    }

    lines <- sub(
        "^\ufeff",
        "",
        lines
    )

    lines <- sub(
        "\r$",
        "",
        lines
    )

    # Skip:
    #   empty lines
    #   PGC VCF metadata lines beginning with ##
    candidate_indices <- which(
        nzchar(trimws(lines)) &
        !grepl("^##", lines)
    )

    if (length(candidate_indices) == 0L) {
        stop(
            "Could not identify the actual header in: ",
            file_path
        )
    }

    header_index <- candidate_indices[1]
    header_line <- lines[header_index]

    parsed <- parse_header_line(
        header_line
    )

    list(
        header_line_number = header_index,
        metadata_lines = header_index - 1L,
        header_line = header_line,
        columns = parsed$columns,
        delimiter = parsed$delimiter
    )
}


pick_column <- function(
    raw_columns,
    aliases
) {

    normalized_columns <- normalize_name(
        raw_columns
    )

    normalized_aliases <- normalize_name(
        aliases
    )

    index <- match(
        normalized_aliases,
        normalized_columns,
        nomatch = 0L
    )

    index <- index[
        index > 0L
    ]

    if (length(index) == 0L) {
        return(NA_character_)
    }

    raw_columns[
        index[1]
    ]
}


column_exists <- function(
    column_name,
    raw_columns
) {

    if (
        is.na(column_name) ||
        column_name == ""
    ) {
        return(FALSE)
    }

    normalize_name(column_name) %in%
        normalize_name(raw_columns)
}


collapse_or_na <- function(x) {

    x <- unique(
        x[
            !is.na(x) &
            x != ""
        ]
    )

    if (length(x) == 0L) {
        return(NA_character_)
    }

    paste(
        x,
        collapse = ";"
    )
}

###############################################################################
# Column aliases
###############################################################################

snp_aliases <- c(
    "SNP",
    "SNP_ID",
    "SNPID",
    "RSID",
    "RS_ID",
    "RS",
    "ID",
    "MARKERNAME",
    "MARKER_NAME",
    "VARIANT_ID"
)

chr_aliases <- c(
    "CHR",
    "CHROM",
    "CHROMOSOME",
    "CHROMSOME",
    "CHROSOME"
)

pos_aliases <- c(
    "BP",
    "POS",
    "POSITION"
)

# Ordinary P-value fields.
p_aliases <- c(
    "P",
    "PVALUE",
    "P_VALUE",
    "PVAL",
    "P_VAL",
    "P-VALUE",
    "GC_PVALUE"
)

# Some PGC sumstats VCF files store -log10(P) as LP.
neglog10_p_aliases <- c(
    "LP",
    "LOG10P",
    "NEG_LOG10_P",
    "MINUS_LOG10_P"
)

# Preference order for signed statistics.
signed_precedence <- c(
    "Z",
    "ZSCORE",
    "Z_SCORE",
    "ZSTAT",
    "Z_STAT",
    "BETA",
    "LOGOR",
    "LOG_OR",
    "LOG_ODDS",
    "OR",
    "ES",
    "EFFECT",
    "EFFECT_SIZE",
    "EFFECTS"
)

a1_aliases <- c(
    "A1",
    "ALLELE1",
    "ALLELE_1",
    "EFFECT_ALLELE",
    "EA",
    "ALT",
    "INC_ALLELE"
)

a2_aliases <- c(
    "A2",
    "ALLELE2",
    "ALLELE_2",
    "OTHER_ALLELE",
    "NON_EFFECT_ALLELE",
    "NEA",
    "REF",
    "DEC_ALLELE"
)

n_aliases <- c(
    "N",
    "TOTAL_N",
    "TOTALN",
    "SAMPLESIZE",
    "SAMPLE_SIZE",
    "EFFECTIVE_N",
    "EFFECTIVEN",
    "N_EFF",
    "NEFF",
    "SS",
    "OBS_CT",
    "WEIGHT"
)

n_case_aliases <- c(
    "NCASE",
    "N_CASE",
    "N_CASES",
    "NCAS",
    "N_CAS",
    "TOTAL_NCASE",
    "TOTAL_N_CASE",
    "CASES_N"
)

n_control_aliases <- c(
    "NCONTROL",
    "N_CONTROL",
    "N_CONTROLS",
    "NCON",
    "N_CON",
    "TOTAL_NCONTROL",
    "TOTAL_N_CONTROL",
    "CONTROLS_N"
)

info_aliases <- c(
    "INFO",
    "INFO_SCORE",
    "RSQ",
    "R2",
    "IMPUTATION_INFO"
)

frq_aliases <- c(
    "FREQ",
    "FREQUENCY",
    "FREQ1",
    "EAF",
    "FRQ",
    "MAF",
    "AF",
    "A1FREQ"
)

###############################################################################
# Exact overrides for known files
###############################################################################

column_overrides <- list(

    AUD = list(
        snp_col = "SNP_ID",
        p_col = "PValue",
        n_col = "SampleSize"
    ),

    ext_cannabis = list(
        snp_col = "SNP",
        p_col = "PVAL",
        n_col = "N"
    ),

    OUD = list(
        snp_col = "SNP_ID",
        p_col = "PValue",
        n_col = "Effective_N"
    ),

    MDD2019 = list(
        snp_col = "MarkerName",
        p_col = "P",
        n_col = NA_character_,
        signed_col = "LogOR"
    ),

    SUD2020 = list(
        snp_col = "rsID",
        p_col = "P-value",
        n_col = "Total_N",
        signed_col = "Zscore"
    )
)

###############################################################################
# Build manifest
###############################################################################

manifest_list <- vector(
    mode = "list",
    length = nrow(inventory)
)

for (i in seq_len(nrow(inventory))) {

    trait_i <- inventory$trait[i]
    sumstats_file_i <- inventory$sumstats_file[i]
    file_name_i <- inventory$file_name[i]

    message(
        "[",
        i,
        "/",
        nrow(inventory),
        "] Reading actual header: ",
        trait_i
    )

    header_info <- read_actual_header(
        sumstats_file_i
    )

    raw_columns <- header_info$columns

    if (length(raw_columns) < 2L) {
        stop(
            "Fewer than two header columns were detected for ",
            trait_i,
            ": ",
            paste(
                raw_columns,
                collapse = "|"
            )
        )
    }

    ###########################################################################
    # Detect variant identifier and genomic coordinates
    ###########################################################################

    snp_col <- pick_column(
        raw_columns,
        snp_aliases
    )

    chr_col <- pick_column(
        raw_columns,
        chr_aliases
    )

    pos_col <- pick_column(
        raw_columns,
        pos_aliases
    )

    ###########################################################################
    # Detect P value
    ###########################################################################

    p_col <- pick_column(
        raw_columns,
        p_aliases
    )

    p_transform <- "identity"

    if (is.na(p_col)) {

        p_col <- pick_column(
            raw_columns,
            neglog10_p_aliases
        )

        if (!is.na(p_col)) {
            p_transform <- "neglog10"
        }
    }

    ###########################################################################
    # Detect signed statistic
    ###########################################################################

    signed_col <- pick_column(
        raw_columns,
        signed_precedence
    )

    signed_null <- NA_real_

    if (!is.na(signed_col)) {

        signed_null <- if (
            normalize_name(signed_col) == "OR"
        ) {
            1
        } else {
            0
        }
    }

    ###########################################################################
    # Detect alleles
    ###########################################################################

    a1_col <- pick_column(
        raw_columns,
        a1_aliases
    )

    a2_col <- pick_column(
        raw_columns,
        a2_aliases
    )

    ###########################################################################
    # Detect sample-size columns
    ###########################################################################

    n_col <- pick_column(
        raw_columns,
        n_aliases
    )

    n_case_col <- pick_column(
        raw_columns,
        n_case_aliases
    )

    n_control_col <- pick_column(
        raw_columns,
        n_control_aliases
    )

    ###########################################################################
    # Optional fields
    ###########################################################################

    info_col <- pick_column(
        raw_columns,
        info_aliases
    )

    frq_col <- pick_column(
        raw_columns,
        frq_aliases
    )

    ###########################################################################
    # Apply known overrides
    ###########################################################################

    override_i <- column_overrides[[trait_i]]

    if (!is.null(override_i)) {

        if ("snp_col" %in% names(override_i)) {
            snp_col <- override_i$snp_col
        }

        if ("p_col" %in% names(override_i)) {
            p_col <- override_i$p_col
            p_transform <- "identity"
        }

        if ("n_col" %in% names(override_i)) {
            n_col <- override_i$n_col
        }

        if ("signed_col" %in% names(override_i)) {

            signed_col <- override_i$signed_col

            signed_null <- if (
                !is.na(signed_col) &&
                normalize_name(signed_col) == "OR"
            ) {
                1
            } else {
                0
            }
        }
    }

    ###########################################################################
    # Validate selected columns
    ###########################################################################

    selected_columns <- list(
        snp_col = snp_col,
        chr_col = chr_col,
        pos_col = pos_col,
        p_col = p_col,
        signed_col = signed_col,
        a1_col = a1_col,
        a2_col = a2_col,
        n_col = n_col,
        n_case_col = n_case_col,
        n_control_col = n_control_col,
        info_col = info_col,
        frq_col = frq_col
    )

    for (column_label in names(selected_columns)) {

        value <- selected_columns[[column_label]]

        if (
            !is.na(value) &&
            !column_exists(
                value,
                raw_columns
            )
        ) {

            warning(
                trait_i,
                ": selected ",
                column_label,
                " was not found in actual header: ",
                value
            )

            selected_columns[[column_label]] <- NA_character_
        }
    }

    snp_col <- selected_columns$snp_col
    chr_col <- selected_columns$chr_col
    pos_col <- selected_columns$pos_col
    p_col <- selected_columns$p_col

    signed_col <- selected_columns$signed_col
    a1_col <- selected_columns$a1_col
    a2_col <- selected_columns$a2_col

    n_col <- selected_columns$n_col
    n_case_col <- selected_columns$n_case_col
    n_control_col <- selected_columns$n_control_col

    info_col <- selected_columns$info_col
    frq_col <- selected_columns$frq_col

    if (is.na(signed_col)) {
        signed_null <- NA_real_
    }

    ###########################################################################
    # Fixed sample size
    ###########################################################################

    N_fixed <- NA_real_

    if (trait_i %in% names(fixed_N_overrides)) {
        N_fixed <- as.numeric(
            fixed_N_overrides[[trait_i]]
        )
    }

    ###########################################################################
    # File format
    ###########################################################################

    file_format <- if (
        grepl(
            "\\.vcf\\.tsv\\.gz$",
            sumstats_file_i,
            ignore.case = TRUE
        )
    ) {

        "pgc_sumstats_vcf"

    } else if (
        grepl(
            "\\.(gz|bgz)$",
            sumstats_file_i,
            ignore.case = TRUE
        )
    ) {

        "compressed_text"

    } else {

        "plain_text"
    }

    ###########################################################################
    # All variant identifier fields are treated as SNP fields
    ###########################################################################

    id_mode <- "snp"

    ###########################################################################
    # Statistical mode
    #
    # Prefer signed statistics when available.
    # Otherwise derive unsigned Z from P.
    ###########################################################################

    if (!is.na(signed_col)) {

        stat_mode <- "signed"

        use_signed_stat <- 1L
        use_unsigned_from_p <- 0L

        signed_sumstats_arg <- paste0(
            signed_col,
            ",",
            signed_null
        )

        a1_inc <- 0L

    } else {

        stat_mode <- "unsigned_from_p"

        use_signed_stat <- 0L
        use_unsigned_from_p <- 1L

        signed_sumstats_arg <- NA_character_

        a1_inc <- 1L
    }

    # Alleles are not required for partitioned heritability.
    no_alleles <- 1L

    ###########################################################################
    # Sample-size mode
    ###########################################################################

    if (!is.na(n_col)) {

        sample_size_mode <- "N_column"

    } else if (
        !is.na(n_case_col) &&
        !is.na(n_control_col)
    ) {

        sample_size_mode <- "case_control_columns"

    } else if (!is.na(N_fixed)) {

        sample_size_mode <- "fixed_N"

    } else {

        sample_size_mode <- "missing"
    }

    ###########################################################################
    # Inclusion checks
    ###########################################################################

    reasons <- character()

    if (is.na(snp_col)) {
        reasons <- c(
            reasons,
            "missing_SNP_column"
        )
    }

    if (is.na(p_col)) {
        reasons <- c(
            reasons,
            "missing_P_or_LP_column"
        )
    }

    has_sample_size <- (
        !is.na(n_col) ||
        (
            !is.na(n_case_col) &&
            !is.na(n_control_col)
        ) ||
        !is.na(N_fixed)
    )

    if (!has_sample_size) {
        reasons <- c(
            reasons,
            "sample_size_required"
        )
    }

    include <- as.integer(
        length(reasons) == 0L
    )

    ###########################################################################
    # Informational preprocessing notes
    ###########################################################################

    notes <- character()

    if (header_info$metadata_lines > 0L) {

        notes <- c(
            notes,
            paste0(
                "skip_",
                header_info$metadata_lines,
                "_metadata_lines"
            )
        )
    }

    if (p_transform == "neglog10") {

        notes <- c(
            notes,
            "convert_neglog10P_to_P"
        )
    }

    if (stat_mode == "unsigned_from_p") {

        notes <- c(
            notes,
            "derive_unsigned_Z_from_P"
        )
    }

    ###########################################################################
    # Save manifest row
    ###########################################################################

    manifest_list[[i]] <- data.table(
        trait = trait_i,
        sumstats_file = sumstats_file_i,
        file_name = file_name_i,

        file_format = file_format,
        delimiter = header_info$delimiter,
        header_line_number = header_info$header_line_number,
        metadata_lines = header_info$metadata_lines,

        id_mode = id_mode,
        snp_col = snp_col,
        chr_col = chr_col,
        pos_col = pos_col,

        p_col = p_col,
        p_transform = p_transform,

        stat_mode = stat_mode,
        signed_col = signed_col,
        signed_null = signed_null,
        signed_sumstats_arg = signed_sumstats_arg,

        use_signed_stat = use_signed_stat,
        use_unsigned_from_p = use_unsigned_from_p,

        a1_col = a1_col,
        a2_col = a2_col,
        no_alleles = no_alleles,
        a1_inc = a1_inc,

        n_col = n_col,
        n_case_col = n_case_col,
        n_control_col = n_control_col,
        N_fixed = N_fixed,
        sample_size_mode = sample_size_mode,

        info_col = info_col,
        frq_col = frq_col,

        samp_prev = NA_real_,
        pop_prev = NA_real_,

        include = include,
        review_reason = collapse_or_na(
            reasons
        ),
        preprocessing_notes = collapse_or_na(
            notes
        ),

        actual_header = paste(
            raw_columns,
            collapse = "|"
        )
    )
}

###############################################################################
# Combine and write manifest
###############################################################################

manifest <- rbindlist(
    manifest_list,
    use.names = TRUE,
    fill = TRUE
)

setorder(
    manifest,
    trait
)

fwrite(
    manifest,
    manifest_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

###############################################################################
# Final summary
###############################################################################

message("")
message("Manifest written to: ", manifest_file)
message("Traits detected: ", nrow(manifest))
message("Ready: ", sum(manifest$include == 1L))
message("Need review: ", sum(manifest$include == 0L))

message("")
message("Manifest summary:")

print(
    manifest[
        ,
        .(
            trait,
            snp_col,
            p_col,
            p_transform,
            signed_col,
            stat_mode,
            n_col,
            N_fixed,
            sample_size_mode,
            include,
            review_reason,
            preprocessing_notes
        )
    ]
)

if (any(manifest$include == 0L)) {

    message("")
    message("Traits requiring review:")

    print(
        manifest[
            include == 0L,
            .(
                trait,
                sumstats_file,
                actual_header,
                review_reason
            )
        ]
    )
}