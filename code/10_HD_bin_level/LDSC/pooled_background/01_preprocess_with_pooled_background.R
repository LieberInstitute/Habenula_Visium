#!/usr/bin/env Rscript

# Preprocess fine-resolution DARs plus selected mid-resolution DARs
# from hg38 to hg19 for LDSC.
#
# Generate three separate DAR annotation sets:
#
#   1. open   = DA_direction == "Up"
#   2. closed = DA_direction == "Down"
#   3. all    = Up + Down combined
#
# Output structure:
#
# 01_DAR_beds_hg19/
# ├── open/
# │   ├── Astrocyte.bed
# │   ├── Microglia.bed
# │   └── ...
# ├── closed/
# │   ├── Astrocyte.bed
# │   ├── Microglia.bed
# │   └── ...
# └── all/
#     ├── Astrocyte.bed
#     ├── Microglia.bed
#     └── ...
#
# These three directories are intended for three separate LDSC runs.


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


###############################################################################
# Pooled DAR background annotation
#
# For each DAR set separately:
#   open   -> union of open DARs across all cell types
#   closed -> union of closed DARs across all cell types
#   all    -> union of all DARs across all cell types
#
# This BED is used only as a control/background annotation in S-LDSC.
# It is intentionally NOT added to the ordinary cell-type manifest.
###############################################################################

pooled_background_annotation <- "__DAR_BACKGROUND__"


###############################################################################
# Add only these cell types from the mid-resolution DAR file
###############################################################################

mid_cell_types <- c(
    "MHb",
    "LHb",
    "Thalamus"
)


###############################################################################
# Create root output directory
###############################################################################

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


###############################################################################
# Helper function: safe file names
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


###############################################################################
# Helper function: read and prepare DAR files
###############################################################################

read_and_prepare_dars <- function(
    dar_file,
    expected_resolution,
    keep_cell_types = NULL
) {

    message("")
    message("Reading DAR file: ", dar_file)

    if (!file.exists(dar_file)) {
        stop(
            "DAR file does not exist: ",
            dar_file
        )
    }


    ###########################################################################
    # Read file
    ###########################################################################

    dar <- fread(dar_file)


    ###########################################################################
    # Required columns
    ###########################################################################

    required_columns <- c(
        "peak",
        "cell_type",
        "DA_direction"
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
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    }


    ###########################################################################
    # Clean variables
    ###########################################################################

    dar[, peak := trimws(
        as.character(peak)
    )]

    dar[, cell_type := trimws(
        as.character(cell_type)
    )]

    dar[, DA_direction := trimws(
        as.character(DA_direction)
    )]


    ###########################################################################
    # Check DA_direction
    #
    # Up   -> open DAR
    # Down -> closed DAR
    ###########################################################################

    observed_directions <- sort(
        unique(
            dar[
                !is.na(DA_direction) &
                DA_direction != "",
                DA_direction
            ]
        )
    )

    message(
        "Observed DA_direction values: ",
        paste(
            observed_directions,
            collapse = ", "
        )
    )

    unexpected_directions <- setdiff(
        observed_directions,
        c(
            "Up",
            "Down"
        )
    )

    if (length(unexpected_directions) > 0L) {

        stop(
            "Unexpected DA_direction values in ",
            dar_file,
            ": ",
            paste(
                unexpected_directions,
                collapse = ", "
            ),
            "\nExpected values are: Up, Down"
        )
    }


    ###########################################################################
    # For the mid-resolution file, retain only selected cell types
    ###########################################################################

    if (!is.null(keep_cell_types)) {

        available_cell_types <- sort(
            unique(
                dar$cell_type
            )
        )

        missing_cell_types <- setdiff(
            keep_cell_types,
            available_cell_types
        )

        if (length(missing_cell_types) > 0L) {

            stop(
                "The following requested mid-resolution cell types ",
                "were not found: ",
                paste(
                    missing_cell_types,
                    collapse = ", "
                ),
                "\nAvailable mid-resolution cell types are: ",
                paste(
                    available_cell_types,
                    collapse = ", "
                )
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
            as.integer(
                coords[[2]]
            )
        ),

        end = suppressWarnings(
            as.integer(
                coords[[3]]
            )
        ),

        source_resolution = expected_resolution,

        source_file = dar_file

    )]


    ###########################################################################
    # Basic validity filtering
    ###########################################################################

    n_before_filter <- nrow(dar)

    dar <- dar[

        !is.na(chr) &
        chr != "" &

        !is.na(start) &
        !is.na(end) &

        !is.na(cell_type) &
        cell_type != "" &

        !is.na(DA_direction) &
        DA_direction %in% c(
            "Up",
            "Down"
        ) &

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
    # Remove duplicate peaks within:
    #
    # cell type × direction × resolution
    #
    # Direction is retained here because open and closed DARs will later
    # be analyzed separately.
    ###########################################################################

    n_before_deduplication <- nrow(dar)

    dar <- unique(
        dar,
        by = c(
            "peak",
            "cell_type",
            "DA_direction",
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
            sort(
                unique(
                    dar$cell_type
                )
            ),
            collapse = ", "
        )
    )


    ###########################################################################
    # Direction counts
    ###########################################################################

    message(
        "DAR counts by direction from ",
        expected_resolution,
        ":"
    )

    print(
        dar[
            ,
            .N,
            by = DA_direction
        ][
            order(
                DA_direction
            )
        ]
    )


    dar
}


###############################################################################
# Read fine-resolution DARs
###############################################################################

fine_dar <- read_and_prepare_dars(
    dar_file = fine_dar_file,
    expected_resolution = "fine"
)


###############################################################################
# Read selected mid-resolution DARs
###############################################################################

mid_dar <- read_and_prepare_dars(
    dar_file = mid_dar_file,
    expected_resolution = "mid",
    keep_cell_types = mid_cell_types
)


###############################################################################
# Make sure data remain after filtering
###############################################################################

if (nrow(fine_dar) == 0L) {

    stop(
        "No fine-resolution DARs remain after processing."
    )
}

if (nrow(mid_dar) == 0L) {

    stop(
        "No selected mid-resolution DARs remain after processing."
    )
}


###############################################################################
# Ensure fine and mid annotations do not have duplicate cell-type names
###############################################################################

overlapping_cell_types <- intersect(
    unique(
        fine_dar$cell_type
    ),
    unique(
        mid_dar$cell_type
    )
)

if (length(overlapping_cell_types) > 0L) {

    stop(
        "These cell-type names occur in both fine and selected mid DARs: ",
        paste(
            overlapping_cell_types,
            collapse = ", "
        ),
        "\nThe script requires each cell-type annotation to come from ",
        "only one resolution."
    )
}


###############################################################################
# Combine fine-resolution + selected mid-resolution DARs
###############################################################################

dar <- rbindlist(
    list(
        fine_dar,
        mid_dar
    ),
    use.names = TRUE,
    fill = TRUE
)


###############################################################################
# Create annotation name
#
# Annotation remains CELL TYPE ONLY.
#
# We do NOT create names such as:
#
#   Astrocyte_Up
#   Astrocyte_Down
#
# because open/closed/all are three separate LDSC runs.
###############################################################################

dar[, annotation := safe_filename(
    cell_type
)]


###############################################################################
# Check annotation names
###############################################################################

if (any(dar$annotation == "")) {

    stop(
        "At least one cell type produced an empty annotation filename."
    )
}


###############################################################################
# Master cell-type map
###############################################################################

cell_map_master <- unique(
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
    cell_map_master,
    source_resolution,
    cell_type
)


###############################################################################
# Check sanitized annotation names
###############################################################################

duplicated_annotation_names <- cell_map_master[
    duplicated(annotation) |
    duplicated(
        annotation,
        fromLast = TRUE
    )
]

if (nrow(duplicated_annotation_names) > 0L) {

    stop(
        "Multiple cell types produce the same sanitized annotation name:\n",
        paste(
            capture.output(
                print(
                    duplicated_annotation_names
                )
            ),
            collapse = "\n"
        )
    )
}


###############################################################################
# Basic combined-data summary
###############################################################################

message("")
message("============================================================")
message("Combined DAR data summary")
message("============================================================")

message(
    "Combined DAR rows: ",
    nrow(dar)
)

message(
    "Total cell-type annotations: ",
    nrow(cell_map_master)
)

message(
    "Fine-resolution annotations: ",
    sum(
        cell_map_master$source_resolution == "fine"
    )
)

message(
    "Mid-resolution annotations: ",
    sum(
        cell_map_master$source_resolution == "mid"
    )
)


message("")
message("DAR counts by cell type and direction:")

print(
    dar[
        ,
        .N,
        by = .(
            source_resolution,
            cell_type,
            DA_direction
        )
    ][
        order(
            source_resolution,
            cell_type,
            DA_direction
        )
    ]
)


###############################################################################
# Create three DAR sets
#
# open:
#   DA_direction == Up
#
# closed:
#   DA_direction == Down
#
# all:
#   union of Up and Down within each cell type
#
# For "all", direction is ignored during deduplication.
###############################################################################

dar_open <- dar[
    DA_direction == "Up"
]

dar_closed <- dar[
    DA_direction == "Down"
]

dar_all <- unique(
    dar,
    by = c(
        "peak",
        "cell_type",
        "source_resolution"
    )
)


dar_sets <- list(
    open = dar_open,
    closed = dar_closed,
    all = dar_all
)


###############################################################################
# Check DAR set sizes
###############################################################################

message("")
message("============================================================")
message("DAR set sizes")
message("============================================================")

message(
    "Open DAR rows:   ",
    nrow(
        dar_sets$open
    )
)

message(
    "Closed DAR rows: ",
    nrow(
        dar_sets$closed
    )
)

message(
    "All DAR rows:    ",
    nrow(
        dar_sets$all
    )
)


###############################################################################
# Print counts by DAR set and cell type
###############################################################################

for (set_name in names(dar_sets)) {

    message("")
    message(
        "Counts for DAR set: ",
        set_name
    )

    print(
        dar_sets[[set_name]][
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
}


###############################################################################
# Import hg38-to-hg19 chain file
###############################################################################

if (!file.exists(chain_file)) {

    stop(
        "Chain file does not exist: ",
        chain_file
    )
}

message("")
message(
    "Importing chain file: ",
    chain_file
)

chain <- import.chain(
    chain_file
)


###############################################################################
# Autosomes
###############################################################################

autosomes <- paste0(
    "chr",
    1:22
)


###############################################################################
# Store manifests from all three runs
###############################################################################

all_manifests <- list()


###############################################################################
# Process each DAR set separately
#
# open
# closed
# all
###############################################################################

for (dar_set_name in names(dar_sets)) {


    ###########################################################################
    # Start DAR set
    ###########################################################################

    message("")
    message("============================================================")
    message(
        "Processing DAR set: ",
        dar_set_name
    )
    message("============================================================")


    dar_current <- copy(
        dar_sets[[dar_set_name]]
    )


    if (nrow(dar_current) == 0L) {

        warning(
            "DAR set is empty: ",
            dar_set_name
        )

        next
    }


    ###########################################################################
    # Output directory for this DAR set
    ###########################################################################

    current_out_dir <- file.path(
        out_dir,
        dar_set_name
    )

    dir.create(
        current_out_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )


    ###########################################################################
    # Build annotation map for this DAR set
    #
    # Some cell types could theoretically have only Up or only Down DARs.
    # Therefore each run uses the annotations actually available in that set.
    ###########################################################################

    cell_map_current <- unique(
        dar_current[
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
        cell_map_current,
        source_resolution,
        cell_type
    )


    ###########################################################################
    # Check missing cell types relative to the master annotation list
    ###########################################################################

    missing_cell_types_this_set <- setdiff(
        cell_map_master$cell_type,
        cell_map_current$cell_type
    )

    if (length(missing_cell_types_this_set) > 0L) {

        warning(
            "The following cell types have no ",
            dar_set_name,
            " DARs and will not produce BED files: ",
            paste(
                missing_cell_types_this_set,
                collapse = ", "
            )
        )
    }


    ###########################################################################
    # Initialize manifest
    ###########################################################################

    manifest <- vector(
        mode = "list",
        length = nrow(
            cell_map_current
        )
    )


    ###########################################################################
    # Process each cell-type annotation
    ###########################################################################

    for (i in seq_len(
        nrow(
            cell_map_current
        )
    )) {


        #######################################################################
        # Annotation information
        #######################################################################

        ct_i <- cell_map_current$cell_type[i]

        annotation_i <- cell_map_current$annotation[i]

        resolution_i <- cell_map_current$source_resolution[i]

        source_file_i <- cell_map_current$source_file[i]


        message(
            "[",
            dar_set_name,
            "] [",
            i,
            "/",
            nrow(
                cell_map_current
            ),
            "] Processing: ",
            ct_i,
            " (",
            resolution_i,
            ")"
        )


        #######################################################################
        # Extract DARs for this cell type
        #######################################################################

        dat <- dar_current[
            cell_type == ct_i &
            source_resolution == resolution_i
        ]


        if (nrow(dat) == 0L) {

            warning(
                "No DAR rows found for: ",
                ct_i,
                " (",
                resolution_i,
                "), DAR set = ",
                dar_set_name
            )

            next
        }


        #######################################################################
        # Convert hg38 coordinates to GRanges
        #
        # Input peak coordinates are assumed:
        #
        #   1-based
        #   closed intervals
        #######################################################################

        gr38 <- GRanges(

            seqnames = dat$chr,

            ranges = IRanges(
                start = dat$start,
                end = dat$end
            )
        )


        names(gr38) <- dat$peak


        #######################################################################
        # LiftOver hg38 -> hg19
        #######################################################################

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


        #######################################################################
        # Retain peaks having exactly one hg19 mapping
        #######################################################################

        uniquely_mapped <- (
            n_mappings == 1L
        )


        if (!any(uniquely_mapped)) {

            warning(
                "No uniquely mapped hg19 intervals for: ",
                ct_i,
                ", DAR set = ",
                dar_set_name
            )

            next
        }


        gr19_all <- unlist(
            lifted_list[
                uniquely_mapped
            ],
            use.names = FALSE
        )


        #######################################################################
        # Retain autosomes only
        #######################################################################

        keep_autosome <- as.character(
            seqnames(
                gr19_all
            )
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
                ct_i,
                ", DAR set = ",
                dar_set_name
            )

            next
        }


        #######################################################################
        # Merge overlapping or directly adjacent intervals
        #
        # This is done separately for each:
        #
        # DAR set × cell type
        #
        # Therefore:
        #
        # open:
        #   reduce only Up peaks
        #
        # closed:
        #   reduce only Down peaks
        #
        # all:
        #   reduce union of Up + Down peaks
        #######################################################################

        gr19_reduced <- reduce(
            gr19_autosomal,
            ignore.strand = TRUE
        )


        if (length(gr19_reduced) == 0L) {

            warning(
                "No intervals remain after reducing overlaps for: ",
                ct_i,
                ", DAR set = ",
                dar_set_name
            )

            next
        }


        #######################################################################
        # Convert GRanges to BED
        #
        # GRanges:
        #
        #   1-based, closed
        #
        # BED:
        #
        #   0-based, half-open
        #
        # Therefore:
        #
        # BED start = GRanges start - 1
        # BED end   = GRanges end
        #######################################################################

        bed <- data.table(

            chr = as.character(
                seqnames(
                    gr19_reduced
                )
            ),

            start = start(
                gr19_reduced
            ) - 1L,

            end = end(
                gr19_reduced
            )
        )


        #######################################################################
        # Sort chromosomes numerically
        #######################################################################

        bed[, chr_number := as.integer(
            sub(
                "^chr",
                "",
                chr
            )
        )]


        setorder(
            bed,
            chr_number,
            start,
            end
        )


        bed[, chr_number := NULL]


        #######################################################################
        # BED QC
        #######################################################################

        if (
            anyNA(bed$start) ||
            anyNA(bed$end)
        ) {

            stop(
                "Missing BED coordinates generated for: ",
                ct_i,
                ", DAR set = ",
                dar_set_name
            )
        }


        if (any(
            bed$start < 0L
        )) {

            stop(
                "Negative BED start coordinate generated for: ",
                ct_i,
                ", DAR set = ",
                dar_set_name
            )
        }


        if (any(
            bed$start >= bed$end
        )) {

            stop(
                "Invalid BED interval generated for: ",
                ct_i,
                ", DAR set = ",
                dar_set_name
            )
        }


        #######################################################################
        # Write BED file
        #######################################################################

        bed_path <- file.path(
            current_out_dir,
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


        #######################################################################
        # Save manifest information
        #######################################################################

        manifest[[i]] <- data.table(

            dar_set = dar_set_name,

            cell_type = ct_i,

            annotation = annotation_i,

            source_resolution = resolution_i,

            source_file = source_file_i,

            input_peaks = nrow(
                dat
            ),

            uniquely_lifted_peaks = sum(
                uniquely_mapped
            ),

            dropped_unmapped_peaks = n_unmapped,

            dropped_multimapped_peaks = n_multimapped,

            dropped_unmapped_or_multimapped = sum(
                !uniquely_mapped
            ),

            dropped_non_autosomal_mappings =
                n_non_autosomal_mappings,

            hg19_autosomal_intervals_before_reduce =
                length(
                    gr19_autosomal
                ),

            hg19_autosomal_intervals =
                length(
                    gr19_reduced
                ),

            hg19_total_bp =
                sum(
                    width(
                        gr19_reduced
                    )
                ),

            bed_file = bed_path
        )
    }


    ###########################################################################
    # Combine manifest rows
    ###########################################################################

    manifest <- rbindlist(
        manifest,
        use.names = TRUE,
        fill = TRUE
    )


    if (nrow(manifest) == 0L) {

        stop(
            "No BED files were generated for DAR set: ",
            dar_set_name
        )
    }


    ###########################################################################
    # Create pooled same-DAR-set background BED
    #
    # This is the union of all cell-type BEDs within the CURRENT DAR set.
    # Example for dar_set_name == "open":
    #
    #   background = union(open DARs from every cell type)
    #
    # The background is generated from the already-lifted hg19 BED files, so
    # target and background annotations use exactly the same liftOver/QC rules.
    ###########################################################################

    if (pooled_background_annotation %chin% manifest$annotation) {
        stop(
            "Reserved pooled background annotation name is already used: ",
            pooled_background_annotation
        )
    }

    pooled_bed_parts <- lapply(
        manifest$bed_file,
        function(bed_file_i) {

            if (!file.exists(bed_file_i) ||
                is.na(file.info(bed_file_i)$size) ||
                file.info(bed_file_i)$size == 0) {
                stop(
                    "Missing/empty cell-type BED while building pooled background: ",
                    bed_file_i
                )
            }

            x <- fread(
                bed_file_i,
                header = FALSE,
                col.names = c("chr", "start", "end")
            )

            if (ncol(x) != 3L || nrow(x) == 0L) {
                stop(
                    "Invalid cell-type BED while building pooled background: ",
                    bed_file_i
                )
            }

            x
        }
    )

    pooled_bed_raw <- rbindlist(
        pooled_bed_parts,
        use.names = TRUE,
        fill = FALSE
    )

    pooled_gr <- GRanges(
        seqnames = pooled_bed_raw$chr,
        ranges = IRanges(
            start = as.integer(pooled_bed_raw$start) + 1L,
            end = as.integer(pooled_bed_raw$end)
        )
    )

    pooled_gr <- reduce(
        pooled_gr,
        ignore.strand = TRUE
    )

    pooled_bed <- data.table(
        chr = as.character(seqnames(pooled_gr)),
        start = start(pooled_gr) - 1L,
        end = end(pooled_gr)
    )

    pooled_bed[, chr_num := suppressWarnings(
        as.integer(sub("^chr", "", chr))
    )]

    if (anyNA(pooled_bed$chr_num)) {
        stop(
            "Non-autosomal chromosome found in pooled background for DAR set: ",
            dar_set_name
        )
    }

    setorder(
        pooled_bed,
        chr_num,
        start,
        end
    )

    pooled_bed[, chr_num := NULL]

    if (nrow(pooled_bed) == 0L || any(pooled_bed$start >= pooled_bed$end)) {
        stop(
            "Invalid pooled background generated for DAR set: ",
            dar_set_name
        )
    }

    pooled_background_bed <- file.path(
        current_out_dir,
        paste0(
            pooled_background_annotation,
            ".bed"
        )
    )

    fwrite(
        pooled_bed[, .(chr, start, end)],
        file = pooled_background_bed,
        sep = "\t",
        col.names = FALSE,
        quote = FALSE
    )

    message(
        "Pooled background BED: ",
        pooled_background_bed
    )

    message(
        "Pooled background intervals: ",
        nrow(pooled_bed),
        "; total bp: ",
        sum(pooled_bed$end - pooled_bed$start)
    )


    ###########################################################################
    # Sort manifest
    ###########################################################################

    setorder(
        manifest,
        source_resolution,
        cell_type
    )


    ###########################################################################
    # Write manifest for this DAR set
    ###########################################################################

    manifest_path <- file.path(
        current_out_dir,
        paste0(
            "DAR_hg19_manifest_",
            dar_set_name,
            ".tsv"
        )
    )


    fwrite(
        manifest,
        file = manifest_path,
        sep = "\t",
        quote = FALSE,
        na = "NA"
    )


    ###########################################################################
    # Save manifest for final combined summary
    ###########################################################################

    all_manifests[[
        dar_set_name
    ]] <- copy(
        manifest
    )


    ###########################################################################
    # Check selected mid-resolution annotations
    ###########################################################################

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

        warning(
            "The following selected mid-resolution annotations ",
            "were not generated for DAR set ",
            dar_set_name,
            ": ",
            paste(
                missing_mid_annotations,
                collapse = ", "
            )
        )
    }


    ###########################################################################
    # Final summary for current DAR set
    ###########################################################################

    message("")
    message(
        "Finished DAR set: ",
        dar_set_name
    )

    message(
        "BED directory: ",
        current_out_dir
    )

    message(
        "Manifest: ",
        manifest_path
    )

    message(
        "BED files generated: ",
        nrow(
            manifest
        )
    )

    message(
        "Fine-resolution BED files: ",
        sum(
            manifest$source_resolution == "fine"
        )
    )

    message(
        "Mid-resolution BED files: ",
        sum(
            manifest$source_resolution == "mid"
        )
    )

    message(
        "Total input peaks: ",
        sum(
            manifest$input_peaks
        )
    )

    message(
        "Total uniquely lifted peaks: ",
        sum(
            manifest$uniquely_lifted_peaks
        )
    )

    message(
        "Total hg19 autosomal intervals after reducing: ",
        sum(
            manifest$hg19_autosomal_intervals
        )
    )

    message(
        "Total hg19 bp covered: ",
        sum(
            manifest$hg19_total_bp
        )
    )
}


###############################################################################
# Combined manifest across open / closed / all
###############################################################################

combined_manifest <- rbindlist(
    all_manifests,
    use.names = TRUE,
    fill = TRUE
)


combined_manifest_path <- file.path(
    out_dir,
    "DAR_hg19_manifest_all_sets.tsv"
)


fwrite(
    combined_manifest,
    file = combined_manifest_path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Final comparison summary
###############################################################################

summary_table <- combined_manifest[
    ,
    .(
        n_annotations = .N,

        total_input_peaks = sum(
            input_peaks
        ),

        total_uniquely_lifted_peaks = sum(
            uniquely_lifted_peaks
        ),

        total_hg19_intervals = sum(
            hg19_autosomal_intervals
        ),

        total_hg19_bp = sum(
            hg19_total_bp
        )
    ),
    by = dar_set
]


summary_table[, dar_set_order := match(
    dar_set,
    c(
        "open",
        "closed",
        "all"
    )
)]

setorder(
    summary_table,
    dar_set_order
)

summary_table[, dar_set_order := NULL]


summary_path <- file.path(
    out_dir,
    "DAR_hg19_summary_by_set.tsv"
)


fwrite(
    summary_table,
    file = summary_path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)


###############################################################################
# Final output
###############################################################################

message("")
message("============================================================")
message("ALL DAR BED GENERATION FINISHED")
message("============================================================")

message("")
message(
    "Root BED directory: ",
    out_dir
)

message("")
message("Three separate LDSC annotation sets:")

message(
    "  OPEN:   ",
    file.path(
        out_dir,
        "open"
    )
)

message(
    "  CLOSED: ",
    file.path(
        out_dir,
        "closed"
    )
)

message(
    "  ALL:    ",
    file.path(
        out_dir,
        "all"
    )
)

message("")
message(
    "Combined manifest: ",
    combined_manifest_path
)

message(
    "Summary table: ",
    summary_path
)

message("")
message("Summary by DAR set:")

print(
    summary_table
)

message("")
message("Done.")