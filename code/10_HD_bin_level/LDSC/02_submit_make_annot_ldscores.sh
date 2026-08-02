#!/bin/bash

set -euo pipefail

PROJECT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC"

BED_DIR="${PROJECT_DIR}/01_DAR_beds_hg19"
MANIFEST="${BED_DIR}/DAR_hg19_manifest.tsv"

OUT_DIR="${PROJECT_DIR}/02_DAR_ldscores_hg19"
TASK_FILE="${OUT_DIR}/DAR_ldscore_tasks.tsv"
SNP_LIST_DIR="${OUT_DIR}/baselineLD_snplists"

BASELINE_PREFIX="/users/cliu3/Thesis/LDscoredata/1000G_Phase3_baselineLD_ldscores/baselineLD."

SCRIPT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/LDSC"

ARRAY_SCRIPT="${SCRIPT_DIR}/02_make_annot_ldscores_array.sh"

mkdir -p \
    "${OUT_DIR}" \
    "${OUT_DIR}/logs" \
    "${SNP_LIST_DIR}"

if [[ ! -s "${MANIFEST}" ]]; then
    echo "Manifest not found or empty:"
    echo "${MANIFEST}"
    exit 1
fi

if [[ ! -s "${ARRAY_SCRIPT}" ]]; then
    echo "Array script not found:"
    echo "${ARRAY_SCRIPT}"
    exit 1
fi

###############################################################################
# Extract the exact SNP set used in the existing baseline-LD files.
###############################################################################

for chr in $(seq 1 22); do

    baseline_ldscore="${BASELINE_PREFIX}${chr}.l2.ldscore.gz"
    snplist="${SNP_LIST_DIR}/baselineLD.${chr}.snplist"

    if [[ ! -s "${baseline_ldscore}" ]]; then
        echo "Missing baseline-LD file:"
        echo "${baseline_ldscore}"
        exit 1
    fi

    zcat "${baseline_ldscore}" |
        awk 'NR > 1 {print $2}' \
        > "${snplist}"

    if [[ ! -s "${snplist}" ]]; then
        echo "Failed to create SNP list:"
        echo "${snplist}"
        exit 1
    fi

    echo "chr${chr}: $(wc -l < "${snplist}") baseline SNPs"
done

###############################################################################
# Construct cell type × chromosome task table.
###############################################################################

: > "${TASK_FILE}"

while IFS=$'\t' read -r annotation bed_file; do

    if [[ -z "${annotation}" || -z "${bed_file}" ]]; then
        continue
    fi

    if [[ ! -s "${bed_file}" ]]; then
        echo "BED file does not exist or is empty:"
        echo "${bed_file}"
        exit 1
    fi

    for chr in $(seq 1 22); do
        printf "%s\t%s\t%s\n" \
            "${annotation}" \
            "${bed_file}" \
            "${chr}" \
            >> "${TASK_FILE}"
    done

done < <(
    tail -n +2 "${MANIFEST}" |
        cut -f2,14
)

N_TASKS=$(wc -l < "${TASK_FILE}")

if [[ "${N_TASKS}" -eq 0 ]]; then
    echo "No tasks were generated."
    exit 1
fi

echo "Task file: ${TASK_FILE}"
echo "Number of tasks: ${N_TASKS}"

###############################################################################
# Limit concurrency to 44 simultaneous jobs.
###############################################################################

sbatch \
    --array="1-${N_TASKS}%44" \
    --export=ALL,TASK_FILE="${TASK_FILE}",OUT_DIR="${OUT_DIR}",SNP_LIST_DIR="${SNP_LIST_DIR}" \
    "${ARRAY_SCRIPT}"