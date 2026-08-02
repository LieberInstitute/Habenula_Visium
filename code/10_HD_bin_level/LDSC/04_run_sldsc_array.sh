#!/bin/bash

#SBATCH --job-name=04_DAR_sldsc
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/04_DAR_sldsc_%a.out
#SBATCH --error=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/04_DAR_sldsc_%a.out

set -euo pipefail

###############################################################################
# Activate LDSC Python 2 environment
###############################################################################

unset PYTHONHOME
unset PYTHONPATH

set +u

source \
    /jhpce/shared/jhpce/core/conda/miniconda3-24.3.0/etc/profile.d/conda.sh

conda activate py2_env

set -u

###############################################################################
# Resources
###############################################################################

LDSC="/users/cliu3/Thesis/ldsc_clean/ldsc.py"

BASELINE_PREFIX="/users/cliu3/Thesis/LDscoredata/1000G_Phase3_baselineLD_ldscores/baselineLD."

WEIGHTS_PREFIX="/users/cliu3/Thesis/LDscoredata/weights_hm3_no_hla/weights."

FRQ_PREFIX="/users/cliu3/Thesis/LDscoredata/1000G_Phase3_frq/1000G.EUR.QC."

ARRAY_ID="${SLURM_ARRAY_TASK_ID:-1}"
JOB_ID="${SLURM_JOB_ID:-manual}"

###############################################################################
# Check environment
###############################################################################

echo "SLURM job:        ${JOB_ID}"
echo "SLURM array task: ${ARRAY_ID}"
echo "Node:             $(hostname)"
echo "Conda environment:${CONDA_DEFAULT_ENV:-unknown}"
echo "Python:           $(which python)"

python --version

###############################################################################
# Check exported variables
###############################################################################

for var in \
    TASK_FILE \
    LDSCORE_DIR; do

    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: ${var} is not defined."
        exit 1
    fi

done

###############################################################################
# Check fixed resources
###############################################################################

for file in \
    "${TASK_FILE}" \
    "${LDSC}"; do

    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Missing or empty file:"
        echo "${file}"
        exit 1
    fi

done

if [[ ! -d "${LDSCORE_DIR}" ]]; then
    echo "ERROR: Custom LD-score directory does not exist:"
    echo "${LDSCORE_DIR}"
    exit 1
fi

###############################################################################
# Determine whether TASK_FILE has a header
#
# Expected columns:
#
# trait
# annotation
# sumstats
# samp_prev
# pop_prev
# result_prefix
###############################################################################

FIRST_FIELD=$(
    awk -F'\t' 'NR == 1 {print $1}' "${TASK_FILE}"
)

if [[ "${FIRST_FIELD}" == "trait" ]]; then

    TASK_LINE_NUMBER=$((ARRAY_ID + 1))
    echo "Task file contains a header."

else

    TASK_LINE_NUMBER="${ARRAY_ID}"
    echo "Task file does not contain a header."

fi

###############################################################################
# Read task
###############################################################################

task_line=$(
    sed -n "${TASK_LINE_NUMBER}p" "${TASK_FILE}"
)

if [[ -z "${task_line}" ]]; then
    echo "ERROR: Could not read S-LDSC task ${ARRAY_ID}."
    echo "Requested task-file line: ${TASK_LINE_NUMBER}"
    exit 1
fi

IFS=$'\t' read -r \
    TRAIT \
    ANNOTATION \
    SUMSTATS \
    SAMP_PREV \
    POP_PREV \
    RESULT_PREFIX \
    <<< "${task_line}"

###############################################################################
# Validate task fields
###############################################################################

is_present() {
    [[ -n "${1:-}" && "${1}" != "NA" ]]
}

for field_name in \
    TRAIT \
    ANNOTATION \
    SUMSTATS \
    RESULT_PREFIX; do

    if ! is_present "${!field_name:-}"; then
        echo "ERROR: Task field ${field_name} is missing."
        echo "Task line:"
        echo "${task_line}"
        exit 1
    fi

done

echo
echo "Trait:         ${TRAIT}"
echo "Annotation:    ${ANNOTATION}"
echo "Sumstats:      ${SUMSTATS}"
echo "Sample prev:   ${SAMP_PREV}"
echo "Population prev:${POP_PREV}"
echo "Result prefix: ${RESULT_PREFIX}"
echo

###############################################################################
# Check munged GWAS input
###############################################################################

if [[ ! -s "${SUMSTATS}" ]]; then
    echo "ERROR: Munged summary-statistics file is missing or empty:"
    echo "${SUMSTATS}"
    exit 1
fi

echo "Munged sumstats preview:"

gzip -cd "${SUMSTATS}" |
    head || true

###############################################################################
# Custom annotation prefix
###############################################################################

CUSTOM_PREFIX="${LDSCORE_DIR}/${ANNOTATION}/${ANNOTATION}."

echo
echo "Custom LD-score prefix:"
echo "${CUSTOM_PREFIX}"

###############################################################################
# Check chromosome files
###############################################################################

echo
echo "Checking LD-score resources..."

for chr in $(seq 1 22); do

    for suffix in \
        "l2.ldscore.gz" \
        "l2.M" \
        "l2.M_5_50"; do

        custom_file="${CUSTOM_PREFIX}${chr}.${suffix}"

        if [[ ! -s "${custom_file}" ]]; then
            echo "ERROR: Missing custom LD-score file:"
            echo "${custom_file}"
            exit 1
        fi

    done

    for file in \
        "${BASELINE_PREFIX}${chr}.l2.ldscore.gz" \
        "${BASELINE_PREFIX}${chr}.l2.M" \
        "${BASELINE_PREFIX}${chr}.l2.M_5_50" \
        "${WEIGHTS_PREFIX}${chr}.l2.ldscore.gz" \
        "${FRQ_PREFIX}${chr}.frq"; do

        if [[ ! -s "${file}" ]]; then
            echo "ERROR: Missing reference file:"
            echo "${file}"
            exit 1
        fi

    done

done

echo "All chromosome resources were found."

###############################################################################
# Create result directory
###############################################################################

RESULT_DIR=$(
    dirname "${RESULT_PREFIX}"
)

mkdir -p "${RESULT_DIR}"

###############################################################################
# Build LDSC command
###############################################################################

CMD=(
    python
    "${LDSC}"

    --h2
    "${SUMSTATS}"

    --ref-ld-chr
    "${BASELINE_PREFIX},${CUSTOM_PREFIX}"

    --w-ld-chr
    "${WEIGHTS_PREFIX}"

    --overlap-annot

    --frqfile-chr
    "${FRQ_PREFIX}"

    --print-coefficients

    --out
    "${RESULT_PREFIX}"
)

###############################################################################
# Optional liability-scale transformation
###############################################################################

if is_present "${SAMP_PREV}" ||
   is_present "${POP_PREV}"; then

    if ! is_present "${SAMP_PREV}" ||
       ! is_present "${POP_PREV}"; then

        echo "ERROR: samp_prev and pop_prev must be supplied together."
        echo "samp_prev: ${SAMP_PREV}"
        echo "pop_prev:  ${POP_PREV}"
        exit 1
    fi

    if ! awk \
        -v samp="${SAMP_PREV}" \
        -v pop="${POP_PREV}" \
        'BEGIN {
            if (
                samp > 0 && samp < 1 &&
                pop > 0 && pop < 1
            ) {
                exit 0
            } else {
                exit 1
            }
        }'; then

        echo "ERROR: Prevalence values must both be between 0 and 1."
        echo "samp_prev: ${SAMP_PREV}"
        echo "pop_prev:  ${POP_PREV}"
        exit 1
    fi

    CMD+=(
        --samp-prev
        "${SAMP_PREV}"

        --pop-prev
        "${POP_PREV}"
    )

    echo
    echo "Liability-scale transformation will be applied."

else

    echo
    echo "Observed-scale heritability will be reported."

fi

###############################################################################
# Run S-LDSC
###############################################################################

echo
echo "Running command:"

printf "%q " "${CMD[@]}"
printf "\n\n"

"${CMD[@]}"

###############################################################################
# Check outputs
###############################################################################

for file in \
    "${RESULT_PREFIX}.log" \
    "${RESULT_PREFIX}.results"; do

    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Missing expected result:"
        echo "${file}"
        exit 1
    fi

done

###############################################################################
# Result preview
###############################################################################

echo
echo "Results preview:"

head -10 "${RESULT_PREFIX}.results" || true

echo
echo "Key log information:"

grep -E \
    "Total Observed scale h2|Total Liability scale h2|Lambda GC|Mean Chi|Intercept|Ratio|WARNING" \
    "${RESULT_PREFIX}.log" || true

echo
echo "Finished successfully:"
echo "${RESULT_PREFIX}.results"