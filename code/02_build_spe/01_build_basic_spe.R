library("SpatialExperiment")
library("spatialLIBD")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")

## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Define some info for the samples
sample_info <- data.frame(
    sample_id = c(
        "V12D07-075_C1"
    )
)
sample_info$subject <- "Br3854"
sample_info$sample_path <-
    file.path(
        here::here("processed-data", "01_spaceranger"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## https://github.com/LieberInstitute/Visium_SPG_AD/blob/master/raw-data/Visium_SPG_AD_ITG_MasterExcelSummarySheet.xlsx
## TODO Update this info!
donor_info <- data.frame(
    subject = c("Br3854"),
    age = c(65.75),
    sex = c("F"),
    race = "EA/CAUC",
    pmi = c(31.5),
    diagnosis = c("Control"),
    rin = c(7)
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)


## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE,
    reference_gtf = here("raw-data", "genes.gtf") ## Not needed at JHPCE
)
Sys.time()
# 2024-02-22 14:43:49.823746 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2024-02-22 14:43:56.531453 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2024-02-22 14:43:56.687378 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2024-02-22 14:43:56.86425 rtracklayer::import: reading the reference GTF file
# 2024-02-22 14:44:28.499147 adding gene information to the SPE object
# 2024-02-22 14:44:28.523581 adding information used by spatialLIBD
# [1] "2024-02-22 14:44:28 EST"

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

# ## Read in cell counts and segmentation results
# segmentations_list <-
#     lapply(sample_info$sample_id, function(sampleid) {
#         file <-
#             here(
#                 "processed-data",
#                 "spaceranger",
#                 sampleid,
#                 "outs",
#                 "spatial",
#                 "tissue_spot_counts.csv"
#             )
#         if (!file.exists(file)) {
#             return(NULL)
#         }
#         x <- read.csv(file)
#         x$key <- paste0(x$barcode, "_", sampleid)
#         return(x)
#     })
# ## Merge them (once the these files are done, this could be replaced by an rbind)
# segmentations <-
#     Reduce(function(...) {
#         merge(..., all = TRUE)
#     }, segmentations_list[lengths(segmentations_list) > 0])
#
# ## Add the information
# segmentation_match <- match(spe$key, segmentations$key)
# segmentation_info <-
#     segmentations[segmentation_match, -which(
#         colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
#     )]
# colData(spe) <- cbind(colData(spe), segmentation_info)

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 8748
length(no_expr) / nrow(spe) * 100
# [1] 23.90099
spe <- spe[-no_expr, ]

## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw_wholegenome <- spe

saveRDS(spe_raw_wholegenome, file.path(dir_rdata, "spe_raw_wholegenome.rds"))

## Size in Gb
lobstr::obj_size(spe_raw_wholegenome)
# 1.651702

## Now drop the spots outside the tissue
spe <- spe_raw_wholegenome[, spe_raw_wholegenome$in_tissue]
dim(spe)
# [1] 27853 38287
## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}


lobstr::obj_size(spe)
# 1.534376

saveRDS(spe, file.path(dir_rdata, "spe.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
