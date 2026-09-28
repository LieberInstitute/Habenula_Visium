#!/bin/bash

#SBATCH -p katun
#SBATCH --job-name=02_DAR_ldscore
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/02_DAR_ldscores_hg19_%a.out
#SBATCH --error=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/02_DAR_ldscores_hg19_%a.out

set -eo pipefail


###############################################################################
# Clean Python environment
###############################################################################

unset PYTHONHOME
unset PYTHONPATH


###############################################################################
# Activate LDSC Python 2 environment
###############################################################################

# Conda deactivate/activate hooks may reference undefined variables.
set +u

source \
    /jhpce/shared/jhpce/core/conda/miniconda3-24.3.0/etc/profile.d/conda.sh

conda activate py2_env


###############################################################################
# Load bedtools
###############################################################################

module load bedtools/2.31.0

set -u


###############################################################################
# Environment information
###############################################################################

echo "Conda environment: ${CONDA_DEFAULT_ENV:-unknown}"
echo "Python executable: $(which python)"
python --version

echo "Bedtools executable: $(which bedtools)"
bedtools --version


###############################################################################
# LDSC resources
###############################################################################

LDSC_DIR="/users/cliu3/Thesis/ldsc_clean"

LDSC="${LDSC_DIR}/ldsc.py"


###############################################################################
# 1000 Genomes EUR PLINK reference
###############################################################################

REF_PREFIX="/users/cliu3/Thesis/PRset/1000G_EUR_Phase3_plink/1000G.EUR.QC"


###############################################################################
# Baseline-LD scores
###############################################################################

BASELINE_PREFIX="/users/cliu3/Thesis/LDscoredata/1000G_Phase3_baselineLD_ldscores/baselineLD."


###############################################################################
# Check required environment variables
###############################################################################

if [[ -z "${TASK_FILE:-}" || ! -s "${TASK_FILE}" ]]; then
    echo "ERROR: TASK_FILE is missing or empty."
    exit 1
fi


if [[ -z "${OUT_DIR:-}" ]]; then
    echo "ERROR: OUT_DIR is not defined."
    exit 1
fi


if [[ -z "${SNP_LIST_DIR:-}" ]]; then
    echo "ERROR: SNP_LIST_DIR is not defined."
    exit 1
fi


###############################################################################
# Check software
###############################################################################

if [[ ! -s "${LDSC}" ]]; then

    echo "ERROR: LDSC script is missing:"
    echo "${LDSC}"

    exit 1
fi


if ! command -v bedtools >/dev/null 2>&1; then

    echo "ERROR: bedtools is not available in PATH."

    exit 1
fi


if ! command -v gzip >/dev/null 2>&1; then

    echo "ERROR: gzip is not available in PATH."

    exit 1
fi


if ! command -v zcat >/dev/null 2>&1; then

    echo "ERROR: zcat is not available in PATH."

    exit 1
fi


###############################################################################
# Check LDSC Python dependencies
###############################################################################

python - <<'PY'

import numpy
import scipy
import pandas
import bitarray

print("Core LDSC Python dependencies loaded successfully.")

PY


###############################################################################
# Read this Slurm array task
###############################################################################

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then

    echo "ERROR: SLURM_ARRAY_TASK_ID is not defined."

    echo "This script must be submitted as a Slurm array job."

    exit 1
fi


task_line=$(
    sed -n "${SLURM_ARRAY_TASK_ID}p" "${TASK_FILE}"
)


if [[ -z "${task_line}" ]]; then

    echo "ERROR: Could not read task ${SLURM_ARRAY_TASK_ID}"

    exit 1
fi


###############################################################################
# Task table columns:
#
#   1 DAR_SET
#   2 ANNOTATION
#   3 BED_FILE
#   4 CHR
###############################################################################

IFS=$'\t' read -r \
    DAR_SET \
    ANNOTATION \
    BED_FILE \
    CHR \
    <<< "${task_line}"


###############################################################################
# Validate task
###############################################################################

if [[ \
    -z "${DAR_SET}" || \
    -z "${ANNOTATION}" || \
    -z "${BED_FILE}" || \
    -z "${CHR}" \
]]; then

    echo "ERROR: Invalid task line:"
    echo "${task_line}"

    exit 1
fi


###############################################################################
# Validate DAR set
###############################################################################

case "${DAR_SET}" in

    open|closed|all)
        ;;

    *)
        echo "ERROR: Invalid DAR_SET:"
        echo "${DAR_SET}"
        exit 1
        ;;

esac


###############################################################################
# Print task information
###############################################################################

echo ""
echo "============================================================"
echo "DAR LDSC task"
echo "============================================================"

echo "DAR set:     ${DAR_SET}"
echo "Annotation:  ${ANNOTATION}"
echo "BED file:    ${BED_FILE}"
echo "Chromosome:  ${CHR}"


###############################################################################
# Paths
###############################################################################

BFILE="${REF_PREFIX}.${CHR}"


###############################################################################
# Baseline SNP list
###############################################################################

SNPLIST="${SNP_LIST_DIR}/baselineLD.${CHR}.snplist"


###############################################################################
# Baseline LD-score file
###############################################################################

BASELINE_LDSCORE="${BASELINE_PREFIX}${CHR}.l2.ldscore.gz"


###############################################################################
# Output directory
#
# Resulting structure:
#
# 02_DAR_ldscores_hg19/
#
#   open/
#       Astrocyte/
#           Astrocyte.1.annot.gz
#           Astrocyte.1.l2.ldscore.gz
#           ...
#
#   closed/
#       Astrocyte/
#           ...
#
#   all/
#       Astrocyte/
#           ...
###############################################################################

DAR_SET_DIR="${OUT_DIR}/${DAR_SET}"

ANNOTATION_DIR="${DAR_SET_DIR}/${ANNOTATION}"

ANNOT_FILE="${ANNOTATION_DIR}/${ANNOTATION}.${CHR}.annot.gz"

OUT_PREFIX="${ANNOTATION_DIR}/${ANNOTATION}.${CHR}"


mkdir -p \
    "${DAR_SET_DIR}" \
    "${ANNOTATION_DIR}"


###############################################################################
# Validate PLINK reference
###############################################################################

for ext in bed bim fam; do

    if [[ ! -s "${BFILE}.${ext}" ]]; then

        echo "ERROR: Missing PLINK file:"
        echo "${BFILE}.${ext}"

        exit 1
    fi

done


###############################################################################
# Validate DAR BED file
###############################################################################

if [[ ! -s "${BED_FILE}" ]]; then

    echo "ERROR: Missing or empty DAR BED file:"
    echo "${BED_FILE}"

    exit 1
fi


###############################################################################
# Validate baseline SNP list
###############################################################################

if [[ ! -s "${SNPLIST}" ]]; then

    echo "ERROR: Missing baseline SNP list:"
    echo "${SNPLIST}"

    exit 1
fi


###############################################################################
# Validate baseline LD-score file
###############################################################################

if [[ ! -s "${BASELINE_LDSCORE}" ]]; then

    echo "ERROR: Missing baseline LD-score file:"
    echo "${BASELINE_LDSCORE}"

    exit 1
fi


###############################################################################
# Step 1
#
# BIM SNPs -> BED
#
# BIM:
#
#   column 1 = chromosome
#   column 2 = SNP ID
#   column 3 = genetic position
#   column 4 = base-pair position
#
# A SNP at 1-based position P becomes:
#
#   [P-1, P)
#
# in BED coordinates.
###############################################################################

JOB_TAG="${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"


TMP_SNP_BED="${ANNOTATION_DIR}/.${ANNOTATION}.${CHR}.${JOB_TAG}.snps.tmp.bed"

TMP_OVERLAP="${ANNOTATION_DIR}/.${ANNOTATION}.${CHR}.${JOB_TAG}.overlap.tmp.tsv"

TMP_ANNOT="${ANNOTATION_DIR}/.${ANNOTATION}.${CHR}.${JOB_TAG}.annot.tmp.gz"


###############################################################################
# Cleanup temporary files on exit
###############################################################################

cleanup() {

    rm -f \
        "${TMP_SNP_BED}" \
        "${TMP_OVERLAP}" \
        "${TMP_ANNOT}"
}

trap cleanup EXIT


###############################################################################
# Create SNP BED
###############################################################################

echo ""
echo "Creating SNP BED file from:"
echo "${BFILE}.bim"


awk '
BEGIN {
    OFS = "\t"
}

{
    chromosome = $1

    if (chromosome !~ /^chr/) {
        chromosome = "chr" chromosome
    }


    if ($4 < 1) {

        print \
            "Invalid BIM position at row " NR ": " $4 \
            > "/dev/stderr"

        exit 1
    }


    # BED is 0-based and half-open.
    #
    # SNP at position P:
    #
    #   start = P - 1
    #   end   = P

    print \
        chromosome, \
        $4 - 1, \
        $4, \
        NR
}

' "${BFILE}.bim" \
    > "${TMP_SNP_BED}"


###############################################################################
# SNP BED QC
###############################################################################

N_BIM=$(wc -l < "${BFILE}.bim")

N_SNP_BED=$(wc -l < "${TMP_SNP_BED}")


echo "BIM SNPs:     ${N_BIM}"
echo "SNP BED rows: ${N_SNP_BED}"


if [[ "${N_BIM}" -ne "${N_SNP_BED}" ]]; then

    echo "ERROR: SNP BED row count does not match BIM row count."

    exit 1
fi


###############################################################################
# Step 2
#
# Determine which SNPs overlap DAR intervals
#
# Output:
#
# 1 chromosome
# 2 BED start
# 3 BED end
# 4 original BIM row
# 5 number of overlapping DAR intervals
###############################################################################

echo ""
echo "Intersecting SNPs with DAR intervals."

echo "DAR set: ${DAR_SET}"

echo "BED:"
echo "${BED_FILE}"


bedtools intersect \
    -a "${TMP_SNP_BED}" \
    -b "${BED_FILE}" \
    -c \
    > "${TMP_OVERLAP}"


###############################################################################
# bedtools output QC
###############################################################################

N_OVERLAP_ROWS=$(wc -l < "${TMP_OVERLAP}")


if [[ "${N_BIM}" -ne "${N_OVERLAP_ROWS}" ]]; then

    echo "ERROR: bedtools output row count does not match BIM row count."

    exit 1
fi


###############################################################################
# Step 3
#
# Create one-column thin annotation
#
# Header:
#
#   annotation name
#
# Values:
#
#   1 = SNP overlaps >= 1 DAR interval
#   0 = SNP does not overlap a DAR interval
###############################################################################

{

    printf "%s\n" "${ANNOTATION}"


    awk '
    {
        if ($5 > 0) {
            print 1
        } else {
            print 0
        }
    }
    ' "${TMP_OVERLAP}"

} | gzip -c > "${TMP_ANNOT}"


###############################################################################
# Only move into final location after successful creation
###############################################################################

mv \
    "${TMP_ANNOT}" \
    "${ANNOT_FILE}"


echo ""
echo "Created annotation:"
echo "${ANNOT_FILE}"


###############################################################################
# Annotation QC
###############################################################################

N_ANNOT=$(
    zcat "${ANNOT_FILE}" |
        awk '
        END {
            print NR - 1
        }
        '
)


N_ANNOTATED=$(
    zcat "${ANNOT_FILE}" |
        awk '
        NR > 1 && $1 != 0 {
            n++
        }

        END {
            print n + 0
        }
        '
)


echo ""
echo "BIM SNPs:        ${N_BIM}"
echo "Annotation rows: ${N_ANNOT}"
echo "Annotated SNPs:  ${N_ANNOTATED}"


if [[ "${N_BIM}" -ne "${N_ANNOT}" ]]; then

    echo "ERROR: Annotation row count does not match BIM."

    exit 1
fi


if [[ "${N_ANNOTATED}" -eq 0 ]]; then

    echo "WARNING: No SNPs were annotated."

    echo "DAR set:     ${DAR_SET}"
    echo "Annotation:  ${ANNOTATION}"
    echo "Chromosome:  ${CHR}"

fi


###############################################################################
# Step 4
#
# Compute custom partitioned LD scores
###############################################################################

echo ""
echo "Computing LD scores."

echo "DAR set:    ${DAR_SET}"
echo "Annotation: ${ANNOTATION}"
echo "Chromosome: ${CHR}"


python "${LDSC}" \
    --l2 \
    --bfile "${BFILE}" \
    --ld-wind-cm 1 \
    --annot "${ANNOT_FILE}" \
    --thin-annot \
    --print-snps "${SNPLIST}" \
    --out "${OUT_PREFIX}"


###############################################################################
# Check expected LDSC outputs
###############################################################################

for file in \
    "${OUT_PREFIX}.l2.ldscore.gz" \
    "${OUT_PREFIX}.l2.M" \
    "${OUT_PREFIX}.l2.M_5_50"; do


    if [[ ! -s "${file}" ]]; then

        echo "ERROR: Expected output is missing:"
        echo "${file}"

        exit 1
    fi

done


###############################################################################
# Compare custom LD-score SNP count with baseline-LD
###############################################################################

N_CUSTOM_LDSCORE=$(
    zcat "${OUT_PREFIX}.l2.ldscore.gz" |
        awk '
        END {
            print NR - 1
        }
        '
)


N_BASELINE_LDSCORE=$(
    zcat "${BASELINE_LDSCORE}" |
        awk '
        END {
            print NR - 1
        }
        '
)


echo ""
echo "Custom LD-score SNPs:   ${N_CUSTOM_LDSCORE}"
echo "Baseline LD-score SNPs: ${N_BASELINE_LDSCORE}"


if [[ "${N_CUSTOM_LDSCORE}" -ne "${N_BASELINE_LDSCORE}" ]]; then

    echo "ERROR: Custom and baseline LD scores contain different SNP counts."

    exit 1
fi


###############################################################################
# Compare SNP IDs AND order with baseline-LD
###############################################################################

if ! diff -q \
    <(
        zcat "${OUT_PREFIX}.l2.ldscore.gz" |
            awk '
            NR > 1 {
                print $2
            }
            '
    ) \
    <(
        zcat "${BASELINE_LDSCORE}" |
            awk '
            NR > 1 {
                print $2
            }
            '
    ) >/dev/null; then


    echo "ERROR: Custom and baseline LD scores contain different SNP IDs or order."

    exit 1
fi


###############################################################################
# Finished
###############################################################################

echo ""
echo "============================================================"
echo "Finished successfully"
echo "============================================================"

echo "DAR set:"
echo "${DAR_SET}"

echo ""

echo "Annotation:"
echo "${ANNOTATION}"

echo ""

echo "Chromosome:"
echo "${CHR}"

echo ""

echo "Output prefix:"
echo "${OUT_PREFIX}"