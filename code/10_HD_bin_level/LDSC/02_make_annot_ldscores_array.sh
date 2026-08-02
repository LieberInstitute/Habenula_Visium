#!/bin/bash

#SBATCH --job-name=02_DAR_ldscore
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/02_DAR_ldscores_hg19_%a.out
#SBATCH --error=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/02_DAR_ldscores_hg19_%a.out

set -eo pipefail

unset PYTHONHOME
unset PYTHONPATH

###############################################################################
# Activate LDSC Python 2 environment
###############################################################################

# Conda deactivate hooks may reference undefined variables.
set +u

source \
    /jhpce/shared/jhpce/core/conda/miniconda3-24.3.0/etc/profile.d/conda.sh

conda activate py2_env

# Load command-line bedtools.
module load bedtools/2.31.0

set -u

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

REF_PREFIX="/users/cliu3/Thesis/PRset/1000G_EUR_Phase3_plink/1000G.EUR.QC"

BASELINE_PREFIX="/users/cliu3/Thesis/LDscoredata/1000G_Phase3_baselineLD_ldscores/baselineLD."

###############################################################################
# Check inputs and software
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

# Only check packages required for LDSC itself.
python - <<'PY'
import numpy
import scipy
import pandas
import bitarray

print("Core LDSC Python dependencies loaded successfully.")
PY

###############################################################################
# Read this array task
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

IFS=$'\t' read -r ANNOTATION BED_FILE CHR <<< "${task_line}"

if [[ -z "${ANNOTATION}" || -z "${BED_FILE}" || -z "${CHR}" ]]; then
    echo "ERROR: Invalid task line:"
    echo "${task_line}"
    exit 1
fi

echo "Annotation: ${ANNOTATION}"
echo "BED file:   ${BED_FILE}"
echo "Chromosome: ${CHR}"

###############################################################################
# Paths
###############################################################################

BFILE="${REF_PREFIX}.${CHR}"

SNPLIST="${SNP_LIST_DIR}/baselineLD.${CHR}.snplist"

BASELINE_LDSCORE="${BASELINE_PREFIX}${CHR}.l2.ldscore.gz"

ANNOTATION_DIR="${OUT_DIR}/${ANNOTATION}"

ANNOT_FILE="${ANNOTATION_DIR}/${ANNOTATION}.${CHR}.annot.gz"

OUT_PREFIX="${ANNOTATION_DIR}/${ANNOTATION}.${CHR}"

mkdir -p "${ANNOTATION_DIR}"

###############################################################################
# Validate input files
###############################################################################

for ext in bed bim fam; do
    if [[ ! -s "${BFILE}.${ext}" ]]; then
        echo "ERROR: Missing PLINK file:"
        echo "${BFILE}.${ext}"
        exit 1
    fi
done

if [[ ! -s "${BED_FILE}" ]]; then
    echo "ERROR: Missing or empty DAR BED file:"
    echo "${BED_FILE}"
    exit 1
fi

if [[ ! -s "${SNPLIST}" ]]; then
    echo "ERROR: Missing baseline SNP list:"
    echo "${SNPLIST}"
    exit 1
fi

if [[ ! -s "${BASELINE_LDSCORE}" ]]; then
    echo "ERROR: Missing baseline LD-score file:"
    echo "${BASELINE_LDSCORE}"
    exit 1
fi

###############################################################################
# Step 1: BED -> thin annotation
#
# The BIM file has:
#   column 1: chromosome
#   column 2: SNP ID
#   column 3: genetic position
#   column 4: base-pair position
#
# Each BIM SNP is converted into a 1-bp BED interval.
# The order of the BIM SNPs is retained.
###############################################################################

JOB_TAG="${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID}"

TMP_SNP_BED="${ANNOTATION_DIR}/.${ANNOTATION}.${CHR}.${JOB_TAG}.snps.tmp.bed"
TMP_OVERLAP="${ANNOTATION_DIR}/.${ANNOTATION}.${CHR}.${JOB_TAG}.overlap.tmp.tsv"
TMP_ANNOT="${ANNOTATION_DIR}/.${ANNOTATION}.${CHR}.${JOB_TAG}.annot.tmp.gz"

cleanup() {
    rm -f \
        "${TMP_SNP_BED}" \
        "${TMP_OVERLAP}" \
        "${TMP_ANNOT}"
}

trap cleanup EXIT

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
        print "Invalid BIM position at row " NR ": " $4 > "/dev/stderr"
        exit 1
    }

    # BED is 0-based and half-open.
    # A SNP at 1-based position P becomes [P-1, P).
    print chromosome, $4 - 1, $4, NR
}
' "${BFILE}.bim" > "${TMP_SNP_BED}"

N_BIM=$(wc -l < "${BFILE}.bim")
N_SNP_BED=$(wc -l < "${TMP_SNP_BED}")

echo "BIM SNPs:     ${N_BIM}"
echo "SNP BED rows: ${N_SNP_BED}"

if [[ "${N_BIM}" -ne "${N_SNP_BED}" ]]; then
    echo "ERROR: SNP BED row count does not match BIM row count."
    exit 1
fi

###############################################################################
# Count DAR overlaps for each SNP.
#
# Output columns:
#   1 chromosome
#   2 BED start
#   3 BED end
#   4 original BIM row number
#   5 number of overlapping DAR intervals
###############################################################################

echo "Intersecting SNPs with DAR intervals."

bedtools intersect \
    -a "${TMP_SNP_BED}" \
    -b "${BED_FILE}" \
    -c \
    > "${TMP_OVERLAP}"

N_OVERLAP_ROWS=$(wc -l < "${TMP_OVERLAP}")

if [[ "${N_BIM}" -ne "${N_OVERLAP_ROWS}" ]]; then
    echo "ERROR: bedtools output row count does not match BIM row count."
    exit 1
fi

###############################################################################
# Create one-column thin annotation.
#
# The header is the annotation name.
# Each subsequent row is:
#   1 if the SNP overlaps at least one DAR
#   0 otherwise
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

# Move into place only after successful creation.
mv "${TMP_ANNOT}" "${ANNOT_FILE}"

echo "Created annotation:"
echo "${ANNOT_FILE}"

###############################################################################
# Annotation QC
###############################################################################

N_ANNOT=$(
    zcat "${ANNOT_FILE}" |
        awk 'END {print NR - 1}'
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

echo "BIM SNPs:        ${N_BIM}"
echo "Annotation rows: ${N_ANNOT}"
echo "Annotated SNPs:  ${N_ANNOTATED}"

if [[ "${N_BIM}" -ne "${N_ANNOT}" ]]; then
    echo "ERROR: Annotation row count does not match BIM."
    exit 1
fi

if [[ "${N_ANNOTATED}" -eq 0 ]]; then
    echo "WARNING: No SNPs were annotated on chromosome ${CHR}."
fi

###############################################################################
# Step 2: Compute custom partitioned LD scores
###############################################################################

echo "Computing LD scores."

python "${LDSC}" \
    --l2 \
    --bfile "${BFILE}" \
    --ld-wind-cm 1 \
    --annot "${ANNOT_FILE}" \
    --thin-annot \
    --print-snps "${SNPLIST}" \
    --out "${OUT_PREFIX}"

###############################################################################
# Output QC
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

N_CUSTOM_LDSCORE=$(
    zcat "${OUT_PREFIX}.l2.ldscore.gz" |
        awk 'END {print NR - 1}'
)

N_BASELINE_LDSCORE=$(
    zcat "${BASELINE_LDSCORE}" |
        awk 'END {print NR - 1}'
)

echo "Custom LD-score SNPs:   ${N_CUSTOM_LDSCORE}"
echo "Baseline LD-score SNPs: ${N_BASELINE_LDSCORE}"

if [[ "${N_CUSTOM_LDSCORE}" -ne "${N_BASELINE_LDSCORE}" ]]; then
    echo "ERROR: Custom and baseline LD scores contain different SNP counts."
    exit 1
fi

if ! diff -q \
    <(
        zcat "${OUT_PREFIX}.l2.ldscore.gz" |
            awk 'NR > 1 {print $2}'
    ) \
    <(
        zcat "${BASELINE_LDSCORE}" |
            awk 'NR > 1 {print $2}'
    ) >/dev/null; then

    echo "ERROR: Custom and baseline LD scores contain different SNP IDs or order."
    exit 1
fi

echo "Finished successfully:"
echo "${OUT_PREFIX}"