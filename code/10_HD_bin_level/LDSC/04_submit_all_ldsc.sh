#!/bin/bash

set -euo pipefail

###############################################################################
# Paths
###############################################################################

PROJECT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC"

GWAS_MANIFEST="${PROJECT_DIR}/03_GWAS/GWAS_manifest.tsv"

DAR_MANIFEST="${PROJECT_DIR}/01_DAR_beds_hg19/DAR_hg19_manifest.tsv"

MUNGED_DIR="${PROJECT_DIR}/03_GWAS/munged"

LDSCORE_DIR="${PROJECT_DIR}/02_DAR_ldscores_hg19"

RESULT_DIR="${PROJECT_DIR}/04_sldsc_results"

TASK_FILE="${PROJECT_DIR}/03_GWAS/sldsc_tasks.tsv"

ARRAY_SCRIPT="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/LDSC/04_run_sldsc_array.sh"

LOG_DIR="${PROJECT_DIR}/logs"

###############################################################################
# Create directories
###############################################################################

mkdir -p \
    "${RESULT_DIR}" \
    "${LOG_DIR}" \
    "$(dirname "${TASK_FILE}")"

###############################################################################
# Check required files and directories
###############################################################################

for file in \
    "${GWAS_MANIFEST}" \
    "${DAR_MANIFEST}" \
    "${ARRAY_SCRIPT}"; do

    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Missing or empty input:"
        echo "${file}"
        exit 1
    fi

done

for directory in \
    "${MUNGED_DIR}" \
    "${LDSCORE_DIR}"; do

    if [[ ! -d "${directory}" ]]; then
        echo "ERROR: Missing directory:"
        echo "${directory}"
        exit 1
    fi

done

###############################################################################
# Temporary files
###############################################################################

TRAIT_FILE=$(mktemp)
ANNOTATION_FILE=$(mktemp)
VALID_ANNOTATION_FILE=$(mktemp)
SKIPPED_TRAIT_FILE=$(mktemp)
SKIPPED_ANNOTATION_FILE=$(mktemp)

cleanup() {

    rm -f \
        "${TRAIT_FILE}" \
        "${ANNOTATION_FILE}" \
        "${VALID_ANNOTATION_FILE}" \
        "${SKIPPED_TRAIT_FILE}" \
        "${SKIPPED_ANNOTATION_FILE}"
}

trap cleanup EXIT

###############################################################################
# Validate GWAS manifest columns
###############################################################################

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }

    required[1] = "trait"
    required[2] = "include"
    required[3] = "samp_prev"
    required[4] = "pop_prev"

    for (j = 1; j <= 4; j++) {
        if (!(required[j] in col)) {
            print "ERROR: Missing GWAS manifest column: " required[j] > "/dev/stderr"
            exit 1
        }
    }

    exit 0
}
' "${GWAS_MANIFEST}"

###############################################################################
# Validate DAR manifest columns
###############################################################################

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }

    if (!("annotation" in col)) {
        print "ERROR: Missing DAR manifest column: annotation" > "/dev/stderr"
        exit 1
    }

    exit 0
}
' "${DAR_MANIFEST}"

###############################################################################
# Extract included GWAS traits
#
# Output:
# trait    samp_prev    pop_prev
###############################################################################

awk -F'\t' '
BEGIN {
    OFS = "\t"
}

NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }
    next
}

$col["include"] == 1 {
    samp = $col["samp_prev"]
    pop = $col["pop_prev"]

    if (samp == "") {
        samp = "NA"
    }

    if (pop == "") {
        pop = "NA"
    }

    print $col["trait"], samp, pop
}
' "${GWAS_MANIFEST}" |
LC_ALL=C sort -t $'\t' -k1,1 \
> "${TRAIT_FILE}"

N_INCLUDED_TRAITS=$(
    awk 'NF > 0 {n++} END {print n + 0}' "${TRAIT_FILE}"
)

if [[ "${N_INCLUDED_TRAITS}" -eq 0 ]]; then
    echo "ERROR: No GWAS traits have include=1."
    exit 1
fi

###############################################################################
# Extract unique DAR annotations
###############################################################################

awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }
    next
}

{
    annotation = $col["annotation"]

    if (annotation != "" && annotation != "NA") {
        print annotation
    }
}
' "${DAR_MANIFEST}" |
LC_ALL=C sort -u \
> "${ANNOTATION_FILE}"

N_ANNOTATIONS=$(
    awk 'NF > 0 {n++} END {print n + 0}' "${ANNOTATION_FILE}"
)

if [[ "${N_ANNOTATIONS}" -eq 0 ]]; then
    echo "ERROR: No annotations were found in the DAR manifest."
    exit 1
fi

###############################################################################
# Check complete chromosome LD-score files
###############################################################################

echo "Checking custom LD-score files..."

while IFS= read -r annotation; do

    [[ -n "${annotation}" ]] || continue

    custom_prefix="${LDSCORE_DIR}/${annotation}/${annotation}."

    annotation_complete=1
    missing_file=""

    for chr in $(seq 1 22); do

        for suffix in \
            "l2.ldscore.gz" \
            "l2.M" \
            "l2.M_5_50"; do

            file="${custom_prefix}${chr}.${suffix}"

            if [[ ! -s "${file}" ]]; then
                annotation_complete=0
                missing_file="${file}"
                break 2
            fi

        done

    done

    if [[ "${annotation_complete}" -eq 1 ]]; then

        printf "%s\n" "${annotation}" \
            >> "${VALID_ANNOTATION_FILE}"

    else

        printf "%s\t%s\n" \
            "${annotation}" \
            "${missing_file}" \
            >> "${SKIPPED_ANNOTATION_FILE}"
    fi

done < "${ANNOTATION_FILE}"

N_VALID_ANNOTATIONS=$(
    awk 'NF > 0 {n++} END {print n + 0}' \
        "${VALID_ANNOTATION_FILE}"
)

if [[ "${N_VALID_ANNOTATIONS}" -eq 0 ]]; then
    echo "ERROR: No annotation has a complete set of LD-score files."
    exit 1
fi

###############################################################################
# Initialize task file
#
# This task file contains a header.
###############################################################################

printf "trait\tannotation\tsumstats\tsamp_prev\tpop_prev\tresult_prefix\n" \
    > "${TASK_FILE}"

###############################################################################
# Construct GWAS × annotation Cartesian product
###############################################################################

N_READY_TRAITS=0

while IFS=$'\t' read -r \
    trait \
    samp_prev \
    pop_prev; do

    [[ -n "${trait}" ]] || continue

    sumstats="${MUNGED_DIR}/${trait}/${trait}.sumstats.gz"

    if [[ ! -s "${sumstats}" ]]; then

        printf "%s\t%s\n" \
            "${trait}" \
            "${sumstats}" \
            >> "${SKIPPED_TRAIT_FILE}"

        continue
    fi

    N_READY_TRAITS=$((N_READY_TRAITS + 1))

    mkdir -p "${RESULT_DIR}/${trait}"

    while IFS= read -r annotation; do

        [[ -n "${annotation}" ]] || continue

        result_prefix="${RESULT_DIR}/${trait}/${trait}__${annotation}"

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${trait}" \
            "${annotation}" \
            "${sumstats}" \
            "${samp_prev}" \
            "${pop_prev}" \
            "${result_prefix}" \
            >> "${TASK_FILE}"

    done < "${VALID_ANNOTATION_FILE}"

done < "${TRAIT_FILE}"

###############################################################################
# Count tasks
###############################################################################

N_TASKS=$(
    awk 'NR > 1 && NF > 0 {n++} END {print n + 0}' \
        "${TASK_FILE}"
)

if [[ "${N_TASKS}" -eq 0 ]]; then
    echo "ERROR: No S-LDSC tasks were created."
    exit 1
fi

EXPECTED_TASKS=$((N_READY_TRAITS * N_VALID_ANNOTATIONS))

if [[ "${N_TASKS}" -ne "${EXPECTED_TASKS}" ]]; then

    echo "ERROR: Unexpected number of S-LDSC tasks."
    echo "Ready traits:      ${N_READY_TRAITS}"
    echo "Valid annotations: ${N_VALID_ANNOTATIONS}"
    echo "Expected tasks:    ${EXPECTED_TASKS}"
    echo "Observed tasks:    ${N_TASKS}"

    exit 1
fi

###############################################################################
# Summary
###############################################################################

echo
echo "Included traits in manifest: ${N_INCLUDED_TRAITS}"
echo "Traits with munged files:     ${N_READY_TRAITS}"
echo "Annotations in manifest:      ${N_ANNOTATIONS}"
echo "Complete annotations:         ${N_VALID_ANNOTATIONS}"
echo "S-LDSC tasks:                 ${N_TASKS}"

echo
echo "Task file:"
echo "${TASK_FILE}"

echo
echo "Ready GWAS traits:"

awk -F'\t' '
NR > 1 {
    trait[$1] = 1
}

END {
    for (x in trait) {
        print x
    }
}
' "${TASK_FILE}" |
LC_ALL=C sort |
nl -ba

echo
echo "Complete annotations:"

nl -ba "${VALID_ANNOTATION_FILE}"

###############################################################################
# Report skipped GWAS traits
###############################################################################

echo
echo "Included GWAS traits without munged output:"

if [[ -s "${SKIPPED_TRAIT_FILE}" ]]; then
    column -t -s $'\t' "${SKIPPED_TRAIT_FILE}"
else
    echo "None"
fi

###############################################################################
# Report incomplete annotations
###############################################################################

echo
echo "Annotations with incomplete LD-score files:"

if [[ -s "${SKIPPED_ANNOTATION_FILE}" ]]; then
    column -t -s $'\t' "${SKIPPED_ANNOTATION_FILE}"
else
    echo "None"
fi

###############################################################################
# Preview task file
###############################################################################

echo
echo "Task preview:"

head -6 "${TASK_FILE}" |
column -t -s $'\t'

###############################################################################
# Submit SLURM array
###############################################################################

echo
echo "Submitting S-LDSC array: 1-${N_TASKS}%20"

JOB_OUTPUT=$(
    sbatch \
        --array="1-${N_TASKS}%20" \
        --export=ALL,TASK_FILE="${TASK_FILE}",LDSCORE_DIR="${LDSCORE_DIR}" \
        "${ARRAY_SCRIPT}"
)

echo "${JOB_OUTPUT}"

echo
echo "Submission finished."