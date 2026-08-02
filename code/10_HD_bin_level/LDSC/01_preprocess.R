#!/usr/bin/env Rscript

# Preprocess fine-resolution DARs plus selected mid-resolution DARs
# from hg38 to hg19 for LDSC.
#
# Input files:
#   peak, cell_type, DA_direction
#
# DA_direction is intentionally ignored. All DARs from the same cell type
# are combined into one annotation.

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(IRanges)
    library(rtracklayer)
})

###############################################################################
# Input and output paths
###############################################################################

fine_dar_file <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Hb_multiome/processed-data/15_DARs/03_gather/",
    "DARs_LDSC_fine.csv.gz"
)

mid_dar_file <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Hb_multiome/processed-data/15_DARs/03_gather/",
    "DARs_LDSC_mid.csv.gz"
)

chain_file <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Hb_multiome/processed-data/10_MAGMA/",
    "hg38ToHg19.over.chain"
)

out_dir <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Habenula_Visium/processed-data/10_HD_bin_level/LDSC/",
    "01_DAR_beds_hg19"
)

# Add only these cell types from the mid-resolution DAR file.
mid_cell_types <- c(
    "MHb",
    "LHb",
    "Thalamus"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

###############################################################################
# Helper functions
###############################################################################

safe_filename <- function(x) {

    y <- gsub(
        "[^A-Za-z0-9._-]+",
        "_",
        x
    )

    y <- gsub(
        "^_+|_+$",
        "",
        y
    )

    y
}


read_and_prepare_dars <- function(
    dar_file,
    expected_resolution,
    keep_cell_types = NULL
) {

    message("Reading DAR file: ", dar_file)

    if (!file.exists(dar_file)) {
        stop("DAR file does not exist: ", dar_file)
    }

    dar <- fread(dar_file)

    ###########################################################################
    # Required columns
    #
    # DA_direction may be present in the file, but it is not used for
    # filtering, grouping, annotation naming, or BED generation.
    ###########################################################################

    required_columns <- c(
        "peak",
        "cell_type"
    )

    missing_columns <- setdiff(
        required_columns,
        names(dar)
    )

    if (length(missing_columns) > 0L) {
        stop(
            "Missing columns in ",
            dar_file,
            ": ",
            paste(missing_columns, collapse = ", ")
        )
    }

    ###########################################################################
    # Clean variables
    ###########################################################################

    dar[, peak := trimws(as.character(peak))]
    dar[, cell_type := trimws(as.character(cell_type))]

    ###########################################################################
    # For the mid-resolution file, retain only selected cell types
    ###########################################################################

    if (!is.null(keep_cell_types)) {

        available_cell_types <- sort(
            unique(dar$cell_type)
        )

        missing_cell_types <- setdiff(
            keep_cell_types,
            available_cell_types
        )

        if (length(missing_cell_types) > 0L) {
            stop(
                "The following requested mid-resolution cell types ",
                "were not found: ",
                paste(missing_cell_types, collapse = ", "),
                "\nAvailable mid-resolution cell types are: ",
                paste(available_cell_types, collapse = ", ")
            )
        }

        dar <- dar[
            cell_type %in% keep_cell_types
        ]
    }

    ###########################################################################
    # Parse chr-start-end coordinates
    ###########################################################################

    coords <- tstrsplit(
        dar$peak,
        "-",
        fixed = TRUE
    )

    if (length(coords) != 3L) {
        stop(
            "Could not parse peak coordinates as chr-start-end in: ",
            dar_file
        )
    }

    dar[, `:=`(
        chr = coords[[1]],
        start = suppressWarnings(
            as.integer(coords[[2]])
        ),
        end = suppressWarnings(
            as.integer(coords[[3]])
        ),
        source_resolution = expected_resolution,
        source_file = dar_file
    )]

    ###########################################################################
    # Basic validity filtering only
    #
    # No filtering is performed using:
    #   resolution
    #   p_val_adj
    #   avg_log2FC
    #   DA_direction
    ###########################################################################

    n_before_filter <- nrow(dar)

    dar <- dar[
        !is.na(chr) &
        chr != "" &
        !is.na(start) &
        !is.na(end) &
        !is.na(cell_type) &
        cell_type != "" &
        start >= 1L &
        start <= end
    ]

    message(
        "Invalid rows removed from ",
        expected_resolution,
        ": ",
        n_before_filter - nrow(dar)
    )

    ###########################################################################
    # Remove duplicate peaks within each cell type
    #
    # DA_direction is deliberately excluded from the duplicate definition.
    # If the same peak appears more than once for the same cell type, it is
    # retained once regardless of direction.
    ###########################################################################

    n_before_deduplication <- nrow(dar)

    dar <- unique(
        dar,
        by = c(
            "peak",
            "cell_type",
            "source_resolution"
        )
    )

    message(
        "Duplicate rows removed from ",
        expected_resolution,
        ": ",
        n_before_deduplication - nrow(dar)
    )

    message(
        "DAR rows retained from ",
        expected_resolution,
        ": ",
        nrow(dar)
    )

    message(
        "Cell types retained from ",
        expected_resolution,
        ": ",
        paste(
            sort(unique(dar$cell_type)),
            collapse = ", "
        )
    )

    dar
}

###############################################################################
# Read fine-resolution and selected mid-resolution DARs
###############################################################################

fine_dar <- read_and_prepare_dars(
    dar_file = fine_dar_file,
    expected_resolution = "fine"
)

mid_dar <- read_and_prepare_dars(
    dar_file = mid_dar_file,
    expected_resolution = "mid",
    keep_cell_types = mid_cell_types
)

if (nrow(fine_dar) == 0L) {
    stop("No fine-resolution DARs remain after processing.")
}

if (nrow(mid_dar) == 0L) {
    stop("No selected mid-resolution DARs remain after processing.")
}

###############################################################################
# Ensure fine and mid annotations do not have duplicate cell-type names
###############################################################################

overlapping_cell_types <- intersect(
    unique(fine_dar$cell_type),
    unique(mid_dar$cell_type)
)

if (length(overlapping_cell_types) > 0L) {
    stop(
        "These cell-type names occur in both fine and selected mid DARs: ",
        paste(overlapping_cell_types, collapse = ", "),
        "\nThe script requires each cell-type annotation to come from ",
        "only one resolution."
    )
}

###############################################################################
# Combine fine and selected mid DARs
###############################################################################

dar <- rbindlist(
    list(
        fine_dar,
        mid_dar
    ),
    use.names = TRUE,
    fill = TRUE
)

# One annotation per cell type.
# DA_direction is not included in the annotation name.
dar[, annotation := safe_filename(cell_type)]

cell_map <- unique(
    dar[
        ,
        .(
            cell_type,
            annotation,
            source_resolution,
            source_file
        )
    ]
)

setorder(
    cell_map,
    source_resolution,
    cell_type
)

###############################################################################
# Check annotation file names
###############################################################################

if (any(cell_map$annotation == "")) {
    stop("At least one cell type produced an empty annotation filename.")
}

duplicated_annotation_names <- cell_map[
    duplicated(annotation) |
    duplicated(annotation, fromLast = TRUE)
]

if (nrow(duplicated_annotation_names) > 0L) {
    stop(
        "Multiple cell types produce the same sanitized annotation name:\n",
        paste(
            capture.output(
                print(duplicated_annotation_names)
            ),
            collapse = "\n"
        )
    )
}

message("Combined DAR rows: ", nrow(dar))
message("Total annotations: ", nrow(cell_map))

message(
    "Fine-resolution annotations: ",
    sum(cell_map$source_resolution == "fine")
)

message(
    "Mid-resolution annotations: ",
    sum(cell_map$source_resolution == "mid")
)

message("DAR counts by annotation:")

print(
    dar[
        ,
        .N,
        by = .(
            source_resolution,
            cell_type
        )
    ][
        order(
            source_resolution,
            cell_type
        )
    ]
)

###############################################################################
# Import hg38-to-hg19 chain file
###############################################################################

if (!file.exists(chain_file)) {
    stop("Chain file does not exist: ", chain_file)
}

message("Importing chain file: ", chain_file)

chain <- import.chain(chain_file)

autosomes <- paste0(
    "chr",
    1:22
)

###############################################################################
# Process each cell-type annotation
###############################################################################

manifest <- vector(
    mode = "list",
    length = nrow(cell_map)
)

for (i in seq_len(nrow(cell_map))) {

    ct_i <- cell_map$cell_type[i]
    annotation_i <- cell_map$annotation[i]
    resolution_i <- cell_map$source_resolution[i]
    source_file_i <- cell_map$source_file[i]

    message(
        "[",
        i,
        "/",
        nrow(cell_map),
        "] Processing: ",
        ct_i,
        " (",
        resolution_i,
        ")"
    )

    dat <- dar[
        cell_type == ct_i &
        source_resolution == resolution_i
    ]

    if (nrow(dat) == 0L) {
        warning(
            "No DAR rows found for: ",
            ct_i,
            " (",
            resolution_i,
            ")"
        )
        next
    }

    ###########################################################################
    # Convert hg38 coordinates to GRanges
    #
    # Assumption:
    # Peak coordinates are 1-based, closed genomic intervals.
    ###########################################################################

    gr38 <- GRanges(
        seqnames = dat$chr,
        ranges = IRanges(
            start = dat$start,
            end = dat$end
        )
    )

    names(gr38) <- dat$peak

    ###########################################################################
    # Lift over from hg38 to hg19
    ###########################################################################

    lifted_list <- liftOver(
        gr38,
        chain
    )

    n_mappings <- elementNROWS(
        lifted_list
    )

    n_unmapped <- sum(
        n_mappings == 0L
    )

    n_multimapped <- sum(
        n_mappings > 1L
    )

    # Retain only peaks with exactly one hg19 mapping.
    uniquely_mapped <- n_mappings == 1L

    if (!any(uniquely_mapped)) {
        warning(
            "No uniquely mapped hg19 intervals for: ",
            ct_i
        )
        next
    }

    gr19_all <- unlist(
        lifted_list[uniquely_mapped],
        use.names = FALSE
    )

    ###########################################################################
    # Retain autosomes
    ###########################################################################

    keep_autosome <- as.character(
        seqnames(gr19_all)
    ) %chin% autosomes

    gr19_autosomal <- gr19_all[
        keep_autosome
    ]

    n_non_autosomal_mappings <- sum(
        !keep_autosome
    )

    if (length(gr19_autosomal) == 0L) {
        warning(
            "No autosomal hg19 intervals retained for: ",
            ct_i
        )
        next
    }

    ###########################################################################
    # Merge overlapping or directly adjacent intervals
    ###########################################################################

    gr19_reduced <- reduce(
        gr19_autosomal,
        ignore.strand = TRUE
    )

    if (length(gr19_reduced) == 0L) {
        warning(
            "No intervals remain after reducing overlaps for: ",
            ct_i
        )
        next
    }

    ###########################################################################
    # Convert GRanges to BED
    #
    # GRanges:
    #   1-based, closed
    #
    # BED:
    #   0-based, half-open
    #
    # BED start = GRanges start - 1
    # BED end   = GRanges end
    ###########################################################################

    bed <- data.table(
        chr = as.character(
            seqnames(gr19_reduced)
        ),
        start = start(gr19_reduced) - 1L,
        end = end(gr19_reduced)
    )

    bed[, chr_number := as.integer(
        sub("^chr", "", chr)
    )]

    setorder(
        bed,
        chr_number,
        start,
        end
    )

    bed[, chr_number := NULL]

    ###########################################################################
    # BED QC
    ###########################################################################

    if (anyNA(bed$start) || anyNA(bed$end)) {
        stop(
            "Missing BED coordinates generated for: ",
            ct_i
        )
    }

    if (any(bed$start < 0L)) {
        stop(
            "Negative BED start coordinate generated for: ",
            ct_i
        )
    }

    if (any(bed$start >= bed$end)) {
        stop(
            "Invalid BED interval generated for: ",
            ct_i
        )
    }

    ###########################################################################
    # Write BED file
    ###########################################################################

    bed_path <- file.path(
        out_dir,
        paste0(
            annotation_i,
            ".bed"
        )
    )

    fwrite(
        bed,
        file = bed_path,
        sep = "\t",
        col.names = FALSE,
        quote = FALSE
    )

    ###########################################################################
    # Save manifest information
    ###########################################################################

    manifest[[i]] <- data.table(
        cell_type = ct_i,
        annotation = annotation_i,
        source_resolution = resolution_i,
        source_file = source_file_i,
        input_peaks = nrow(dat),
        uniquely_lifted_peaks = sum(uniquely_mapped),
        dropped_unmapped_peaks = n_unmapped,
        dropped_multimapped_peaks = n_multimapped,
        dropped_unmapped_or_multimapped = sum(!uniquely_mapped),
        dropped_non_autosomal_mappings = n_non_autosomal_mappings,
        hg19_autosomal_intervals_before_reduce = length(
            gr19_autosomal
        ),
        hg19_autosomal_intervals = length(
            gr19_reduced
        ),
        hg19_total_bp = sum(
            width(gr19_reduced)
        ),
        bed_file = bed_path
    )
}

###############################################################################
# Write manifest
###############################################################################

manifest <- rbindlist(
    manifest,
    use.names = TRUE,
    fill = TRUE
)

if (nrow(manifest) == 0L) {
    stop("No BED files were generated.")
}

setorder(
    manifest,
    source_resolution,
    cell_type
)

manifest_path <- file.path(
    out_dir,
    "DAR_hg19_manifest.tsv"
)

fwrite(
    manifest,
    file = manifest_path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

###############################################################################
# Final checks
###############################################################################

expected_mid_annotations <- safe_filename(
    mid_cell_types
)

missing_mid_annotations <- setdiff(
    expected_mid_annotations,
    manifest[
        source_resolution == "mid",
        annotation
    ]
)

if (length(missing_mid_annotations) > 0L) {
    stop(
        "The following selected mid-resolution BED files were not generated: ",
        paste(missing_mid_annotations, collapse = ", ")
    )
}

missing_annotations <- setdiff(
    cell_map$annotation,
    manifest$annotation
)

if (length(missing_annotations) > 0L) {
    warning(
        "The following annotations did not produce BED files: ",
        paste(missing_annotations, collapse = ", ")
    )
}

###############################################################################
# Final summary
###############################################################################

message("Finished.")
message("BED directory: ", out_dir)
message("Manifest: ", manifest_path)
message("BED files generated: ", nrow(manifest))

message(
    "Fine-resolution BED files: ",
    sum(manifest$source_resolution == "fine")
)

message(
    "Mid-resolution BED files: ",
    sum(manifest$source_resolution == "mid")
)

message(
    "Mid-resolution annotations included: ",
    paste(
        manifest[
            source_resolution == "mid",
            cell_type
        ],
        collapse = ", "
    )
)

message(
    "Total input peaks: ",
    sum(manifest$input_peaks)
)

message(
    "Total uniquely lifted peaks: ",
    sum(manifest$uniquely_lifted_peaks)
)

message(
    "Total hg19 autosomal intervals after reducing: ",
    sum(manifest$hg19_autosomal_intervals)
)