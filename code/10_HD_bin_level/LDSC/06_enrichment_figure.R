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
    "Habenula_Visium/raw-data/sample_info/",
    "hd_cell_type_map.csv"
)

output_dir <- paste0(
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/",
    "Habenula_Visium/processed-data/10_HD_bin_level/LDSC/",
    "04_sldsc_results/heatmap"
)

# Minimum number of common SNPs required for display.
# Change to 0 to show every annotation.
minimum_n_snps <- 100L

# FALSE: show all cells and mark nominally significant cells with stars.
# TRUE: only nominally significant cells receive color; others are gray.
show_only_significant <- FALSE

# Maximum absolute value shown on the color scale.
# Larger values are clipped only for visualization.
color_limit <- 6

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

output_pdf <- file.path(
    output_dir,
    "DAR_sldsc_enrichment_heatmap_nominalP.pdf"
)

output_png <- file.path(
    output_dir,
    "DAR_sldsc_enrichment_heatmap_nominalP.png"
)

output_plot_data <- file.path(
    output_dir,
    "DAR_sldsc_enrichment_heatmap_plot_data.tsv"
)

output_annotation_map <- file.path(
    output_dir,
    "DAR_sldsc_annotation_name_mapping.tsv"
)

###############################################################################
# Read and validate S-LDSC input
###############################################################################

if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
}

if (file.info(input_file)$size == 0) {
    stop("Input file is empty: ", input_file)
}

dt <- fread(
    input_file,
    sep = "\t",
    quote = "",
    na.strings = c("", "NA", "nan")
)

required_columns <- c(
    "trait",
    "annotation",
    "source_resolution",
    "n_snps",
    "Enrichment",
    "Enrichment_p"
)

missing_columns <- setdiff(
    required_columns,
    names(dt)
)

if (length(missing_columns) > 0L) {
    stop(
        "Missing required columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

dt[
    ,
    `:=`(
        trait = trimws(as.character(trait)),
        annotation = trimws(as.character(annotation)),
        source_resolution = trimws(as.character(source_resolution)),
        n_snps = suppressWarnings(as.numeric(n_snps)),
        Enrichment = suppressWarnings(as.numeric(Enrichment)),
        Enrichment_p = suppressWarnings(as.numeric(Enrichment_p))
    )
]

dt <- dt[
    !is.na(trait) &
        trait != "" &
        !is.na(annotation) &
        annotation != "" &
        !is.na(n_snps) &
        !is.na(Enrichment) &
        !is.na(Enrichment_p)
]

if (nrow(dt) == 0L) {
    stop("No complete rows remained after input validation.")
}

###############################################################################
# Read and validate cell-type mapping
###############################################################################

if (!file.exists(cell_type_map_file)) {
    stop(
        "Cell-type mapping file does not exist: ",
        cell_type_map_file
    )
}

if (file.info(cell_type_map_file)$size == 0) {
    stop(
        "Cell-type mapping file is empty: ",
        cell_type_map_file
    )
}

cell_type_map <- fread(
    cell_type_map_file,
    sep = ",",
    quote = "\"",
    na.strings = c("", "NA", "nan")
)

required_map_columns <- c(
    "old_cell_type",
    "new_cell_type",
    "color"
)

missing_map_columns <- setdiff(
    required_map_columns,
    names(cell_type_map)
)

if (length(missing_map_columns) > 0L) {
    stop(
        "Missing required columns in cell-type mapping file: ",
        paste(missing_map_columns, collapse = ", ")
    )
}

cell_type_map[
    ,
    `:=`(
        old_cell_type = trimws(as.character(old_cell_type)),
        new_cell_type = trimws(as.character(new_cell_type)),
        color = trimws(as.character(color))
    )
]

cell_type_map <- cell_type_map[
    !is.na(old_cell_type) &
        old_cell_type != "" &
        !is.na(new_cell_type) &
        new_cell_type != ""
]

if (nrow(cell_type_map) == 0L) {
    stop(
        "No valid rows remained in the cell-type mapping file."
    )
}

###############################################################################
# Validate mapping-file uniqueness
###############################################################################

duplicated_old_names <- cell_type_map[
    duplicated(old_cell_type) |
        duplicated(old_cell_type, fromLast = TRUE)
]

if (nrow(duplicated_old_names) > 0L) {

    duplicate_summary <- duplicated_old_names[
        ,
        .(
            mapped_names = paste(
                unique(new_cell_type),
                collapse = " | "
            ),
            mapped_colors = paste(
                unique(color),
                collapse = " | "
            )
        ),
        by = old_cell_type
    ]

    message("")
    message("Duplicated old_cell_type entries:")
    print(duplicate_summary)

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

color_lookup <- setNames(
    cell_type_map$color,
    cell_type_map$old_cell_type
)

# Preserve the original annotation name.
dt[
    ,
    annotation_original := annotation
]

# Determine whether each annotation occurs in the mapping file.
dt[
    ,
    annotation_mapped :=
        annotation_original %chin% cell_type_map$old_cell_type
]

# Obtain the corrected display name.
dt[
    ,
    annotation_display := unname(
        name_lookup[annotation_original]
    )
]

# Obtain the corresponding cell-type color.
dt[
    ,
    annotation_color := unname(
        color_lookup[annotation_original]
    )
]

# Retain the original name when an annotation is not found in the mapping file.
dt[
    is.na(annotation_display) |
        annotation_display == "",
    annotation_display := annotation_original
]

###############################################################################
# Report mapped and unmapped annotation names
###############################################################################

mapping_audit <- unique(
    dt[
        ,
        .(
            annotation_original,
            annotation_display,
            annotation_color,
            annotation_mapped
        )
    ]
)

setorder(
    mapping_audit,
    -annotation_mapped,
    annotation_display,
    annotation_original
)

fwrite(
    mapping_audit,
    output_annotation_map,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

mapped_annotations <- sort(
    unique(
        dt[
            annotation_mapped == TRUE,
            annotation_original
        ]
    )
)

unmapped_annotations <- sort(
    unique(
        dt[
            annotation_mapped == FALSE,
            annotation_original
        ]
    )
)

message("")
message(
    "Mapped annotations: ",
    length(mapped_annotations)
)

message(
    "Unmapped annotations: ",
    length(unmapped_annotations)
)

if (length(unmapped_annotations) > 0L) {
    warning(
        paste0(
            "The following annotations were not found in ",
            basename(cell_type_map_file),
            " and will retain their original names:\n",
            paste(
                unmapped_annotations,
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}

###############################################################################
# Check whether multiple old names share the same new display name
###############################################################################

duplicated_display_names <- mapping_audit[
    ,
    .(
        n_original_names = uniqueN(annotation_original),
        original_names = paste(
            sort(unique(annotation_original)),
            collapse = " | "
        )
    ),
    by = annotation_display
][
    n_original_names > 1L
]

if (nrow(duplicated_display_names) > 0L) {
    warning(
        paste0(
            "Some corrected display names correspond to multiple original ",
            "annotation names. The columns remain separate internally, but ",
            "their displayed x-axis labels may be identical:\n",
            paste(
                paste0(
                    duplicated_display_names$annotation_display,
                    " <- ",
                    duplicated_display_names$original_names
                ),
                collapse = "\n"
            )
        ),
        call. = FALSE
    )
}

###############################################################################
# Filter annotations by common-SNP count
###############################################################################

plot_dt <- dt[
    n_snps >= minimum_n_snps
]

if (nrow(plot_dt) == 0L) {
    stop(
        "No rows remain after applying minimum_n_snps = ",
        minimum_n_snps,
        "."
    )
}

###############################################################################
# Nominal significance and plotting metric
###############################################################################

# Enrichment_p is treated as the nominal P value.
plot_dt[
    ,
    nominal_significant := Enrichment_p < 0.05
]

plot_dt[
    ,
    significance_label := fifelse(
        Enrichment_p < 0.001,
        "***",
        fifelse(
            Enrichment_p < 0.01,
            "**",
            fifelse(
                Enrichment_p < 0.05,
                "*",
                ""
            )
        )
    )
]

# Positive values indicate Enrichment > 1.
# Negative values indicate Enrichment < 1.
# Magnitude is -log10(nominal P).
plot_dt[
    ,
    signed_log10_p :=
        sign(Enrichment - 1) *
        -log10(
            pmax(
                Enrichment_p,
                1e-300
            )
        )
]

plot_dt[
    ,
    signed_log10_p_clipped := pmax(
        pmin(
            signed_log10_p,
            color_limit
        ),
        -color_limit
    )
]

if (show_only_significant) {

    plot_dt[
        ,
        fill_value := fifelse(
            nominal_significant,
            signed_log10_p_clipped,
            NA_real_
        )
    ]

} else {

    plot_dt[
        ,
        fill_value := signed_log10_p_clipped
    ]
}

###############################################################################
# Order traits
###############################################################################

# Traits with stronger nominal signals are placed first.
trait_order <- plot_dt[
    ,
    .(
        ordering_score = max(
            abs(signed_log10_p_clipped),
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

###############################################################################
# Collect and validate annotation-level information
###############################################################################

annotation_info <- unique(
    plot_dt[
        ,
        .(
            annotation,
            annotation_original,
            annotation_display,
            annotation_color,
            annotation_mapped,
            source_resolution,
            n_snps
        )
    ]
)

# Check that each annotation has one resolution and one SNP count.
annotation_metadata_check <- annotation_info[
    ,
    .(
        n_resolutions = uniqueN(source_resolution),
        n_snp_counts = uniqueN(n_snps),
        resolution_values = paste(
            sort(unique(source_resolution)),
            collapse = " | "
        ),
        n_snp_values = paste(
            sort(unique(n_snps)),
            collapse = " | "
        )
    ),
    by = annotation
][
    n_resolutions > 1L |
        n_snp_counts > 1L
]

if (nrow(annotation_metadata_check) > 0L) {
    message("")
    message(
        "Annotations with inconsistent resolution or n_snps values:"
    )

    print(annotation_metadata_check)

    stop(
        paste0(
            "Each annotation must have a consistent source_resolution ",
            "and n_snps value across traits."
        )
    )
}

annotation_info <- unique(
    annotation_info,
    by = "annotation"
)

###############################################################################
# Order annotations
###############################################################################

resolution_order <- c(
    "broad",
    "mid",
    "fine"
)

annotation_info[
    ,
    resolution_rank := match(
        tolower(source_resolution),
        resolution_order
    )
]

annotation_info[
    is.na(resolution_rank),
    resolution_rank := length(resolution_order) + 1L
]

annotation_order <- annotation_info[
    order(
        resolution_rank,
        -n_snps,
        annotation_display,
        annotation
    ),
    annotation
]

###############################################################################
# Convert plotting variables to factors
###############################################################################

plot_dt[
    ,
    trait := factor(
        trait,
        levels = rev(trait_order)
    )
]

plot_dt[
    ,
    annotation := factor(
        annotation,
        levels = annotation_order
    )
]

###############################################################################
# Construct corrected annotation labels
###############################################################################

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

###############################################################################
# Write the exact data used in the heatmap
###############################################################################

setorder(
    plot_dt,
    trait,
    annotation
)

fwrite(
    plot_dt,
    output_plot_data,
    sep = "\t",
    quote = FALSE,
    na = "NA"
)

###############################################################################
# Draw heatmap
###############################################################################

subtitle_text <- paste0(
    "Nominal Enrichment P values; ",
    "* P<0.05, ** P<0.01, *** P<0.001; ",
    "annotations with n_snps ≥ ",
    format(
        minimum_n_snps,
        big.mark = ","
    )
)

if (show_only_significant) {
    subtitle_text <- paste0(
        subtitle_text,
        "; nonsignificant cells shown in gray"
    )
}

p <- ggplot(
    plot_dt,
    aes(
        x = annotation,
        y = trait,
        fill = fill_value
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.35
    ) +
    geom_text(
        aes(
            label = significance_label
        ),
        size = 3.4,
        fontface = "bold",
        na.rm = TRUE
    ) +
    scale_x_discrete(
        labels = annotation_labels,
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
        oob = squish,
        na.value = "grey90",
        name = paste0(
            "Direction ×\n-log10(P)\n",
            "(clipped at ±",
            color_limit,
            ")"
        )
    ) +
    labs(
        title = paste0(
            "S-LDSC enrichment across GWAS traits ",
            "and DAR annotations"
        ),
        subtitle = subtitle_text,
        x = "DAR annotation and number of common SNPs",
        y = "GWAS trait",
        caption = paste0(
            "Red: Enrichment > 1; blue: Enrichment < 1. ",
            "Negative enrichment estimates should be interpreted ",
            "as out-of-bounds estimates."
        )
    ) +
    coord_fixed(
        ratio = 0.8
    ) +
    theme_bw(
        base_size = 11
    ) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1,
            size = 9
        ),
        axis.text.y = element_text(
            size = 9
        ),
        axis.title = element_text(
            face = "bold"
        ),
        plot.title = element_text(
            face = "bold",
            size = 14
        ),
        plot.subtitle = element_text(
            size = 10
        ),
        plot.caption = element_text(
            hjust = 0,
            size = 8
        ),
        legend.title = element_text(
            size = 10
        )
    )

###############################################################################
# Save outputs
###############################################################################

n_traits <- uniqueN(plot_dt$trait)
n_annotations <- uniqueN(plot_dt$annotation)

plot_width <- max(
    10,
    0.75 * n_annotations + 4
)

plot_height <- max(
    6,
    0.36 * n_traits + 2.5
)

ggsave(
    filename = output_pdf,
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in"
)

ggsave(
    filename = output_png,
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300
)

###############################################################################
# Completion messages
###############################################################################

message("")
message("Heatmap PDF:       ", output_pdf)
message("Heatmap PNG:       ", output_png)
message("Plot data:         ", output_plot_data)
message("Annotation mapping:", output_annotation_map)
message("Rows plotted:      ", nrow(plot_dt))
message("Traits:            ", n_traits)
message("Annotations:       ", n_annotations)
message(
    "Nominal P<0.05 cells: ",
    plot_dt[
        nominal_significant == TRUE,
        .N
    ]
)