#!/bin/bash

set -euo pipefail

###############################################################################
# Paths
###############################################################################

PROJECT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC"

MANIFEST="${PROJECT_DIR}/03_GWAS/GWAS_manifest.tsv"

OUT_DIR="${PROJECT_DIR}/03_GWAS/munged"

TASK_FILE="${PROJECT_DIR}/03_GWAS/GWAS_munge_tasks.tsv"

ARRAY_SCRIPT="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/LDSC/03c_run_munge_GWAS_array.sh"

LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p \
    "${OUT_DIR}" \
    "${LOG_DIR}" \
    "$(dirname "${TASK_FILE}")"

###############################################################################
# Locate HapMap3 SNP list
###############################################################################

HM3_SNPLIST="/users/cliu3/Thesis/LDscoredata/LDSCORE_w_hm3.snplist"

if [[ ! -s "${HM3_SNPLIST}" ]]; then

    echo "Default HapMap3 SNP list was not found:"
    echo "${HM3_SNPLIST}"
    echo
    echo "Searching under /users/cliu3/Thesis ..."

    mapfile -t HM3_MATCHES < <(
        find /users/cliu3/Thesis \
            -type f \
            -name "w_hm3.snplist" \
            2>/dev/null |
        sort
    )

    if [[ "${#HM3_MATCHES[@]}" -eq 1 ]]; then

        HM3_SNPLIST="${HM3_MATCHES[0]}"

        echo "Using:"
        echo "${HM3_SNPLIST}"

    else

        echo "ERROR: Could not uniquely locate w_hm3.snplist."
        echo "Number of matches: ${#HM3_MATCHES[@]}"

        if [[ "${#HM3_MATCHES[@]}" -gt 0 ]]; then
            printf "  %s\n" "${HM3_MATCHES[@]}"
        fi

        exit 1
    fi
fi

###############################################################################
# Validate inputs
###############################################################################

for file in \
    "${MANIFEST}" \
    "${ARRAY_SCRIPT}" \
    "${HM3_SNPLIST}"; do

    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Missing or empty input:"
        echo "${file}"
        exit 1
    fi
done

###############################################################################
# Build task table from current manifest
#
# The task file intentionally has no header because SLURM_ARRAY_TASK_ID starts
# from 1 and the array script reads one task per line.
###############################################################################

TMP_TASK_FILE="${TASK_FILE}.tmp"

awk -F'\t' '
BEGIN {
    OFS = "\t"
}

NR == 1 {

    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }

    required[1]  = "trait"
    required[2]  = "sumstats_file"
    required[3]  = "file_format"
    required[4]  = "delimiter"
    required[5]  = "header_line_number"
    required[6]  = "snp_col"
    required[7]  = "p_col"
    required[8]  = "p_transform"
    required[9]  = "stat_mode"
    required[10] = "signed_col"
    required[11] = "signed_null"
    required[12] = "signed_sumstats_arg"
    required[13] = "use_signed_stat"
    required[14] = "use_unsigned_from_p"
    required[15] = "a1_col"
    required[16] = "a2_col"
    required[17] = "no_alleles"
    required[18] = "a1_inc"
    required[19] = "n_col"
    required[20] = "n_case_col"
    required[21] = "n_control_col"
    required[22] = "N_fixed"
    required[23] = "sample_size_mode"
    required[24] = "info_col"
    required[25] = "frq_col"
    required[26] = "samp_prev"
    required[27] = "pop_prev"
    required[28] = "include"

    for (j = 1; j <= 28; j++) {

        if (!(required[j] in col)) {

            print \
                "ERROR: Missing manifest column: " required[j] \
                > "/dev/stderr"

            exit 1
        }
    }

    next
}

$col["include"] == 1 {

    print \
        $col["trait"], \
        $col["sumstats_file"], \
        $col["file_format"], \
        $col["delimiter"], \
        $col["header_line_number"], \
        $col["snp_col"], \
        $col["p_col"], \
        $col["p_transform"], \
        $col["stat_mode"], \
        $col["signed_col"], \
        $col["signed_null"], \
        $col["signed_sumstats_arg"], \
        $col["use_signed_stat"], \
        $col["use_unsigned_from_p"], \
        $col["a1_col"], \
        $col["a2_col"], \
        $col["no_alleles"], \
        $col["a1_inc"], \
        $col["n_col"], \
        $col["n_case_col"], \
        $col["n_control_col"], \
        $col["N_fixed"], \
        $col["sample_size_mode"], \
        $col["info_col"], \
        $col["frq_col"], \
        $col["samp_prev"], \
        $col["pop_prev"]
}
' "${MANIFEST}" > "${TMP_TASK_FILE}"

mv "${TMP_TASK_FILE}" "${TASK_FILE}"

###############################################################################
# Count and validate tasks
###############################################################################

N_TASKS=$(
    awk '
        NF > 0 {
            n++
        }

        END {
            print n + 0
        }
    ' "${TASK_FILE}"
)

if [[ "${N_TASKS}" -eq 0 ]]; then
    echo "ERROR: No traits have include=1 in the manifest."
    exit 1
fi

###############################################################################
# Confirm every selected GWAS file exists
###############################################################################

while IFS=$'\t' read -r \
    TRAIT \
    SUMSTATS_FILE \
    REST; do

    if [[ ! -s "${SUMSTATS_FILE}" ]]; then

        echo "ERROR: Selected GWAS file is missing or empty."
        echo "Trait: ${TRAIT}"
        echo "File:  ${SUMSTATS_FILE}"

        exit 1
    fi

done < "${TASK_FILE}"

###############################################################################
# Summary
###############################################################################

echo
echo "Munge tasks: ${N_TASKS}"
echo "Task file:   ${TASK_FILE}"
echo "Output dir:  ${OUT_DIR}"
echo "HM3 file:    ${HM3_SNPLIST}"
echo "Array script:${ARRAY_SCRIPT}"

echo
echo "Selected traits:"

cut -f1 "${TASK_FILE}" |
    nl -ba

echo
echo "Excluded traits:"

EXCLUDED=$(
    awk -F'\t' '
    NR == 1 {

        for (i = 1; i <= NF; i++) {
            col[$i] = i
        }

        next
    }

    $col["include"] != 1 {

        reason = $col["review_reason"]

        if (reason == "" || reason == "NA") {
            reason = "include_not_equal_to_1"
        }

        print $col["trait"] "\t" reason
    }
    ' "${MANIFEST}"
)

if [[ -n "${EXCLUDED}" ]]; then
    printf "%s\n" "${EXCLUDED}" |
        column -t -s $'\t'
else
    echo "None"
fi

###############################################################################
# Submit SLURM array
###############################################################################

echo
echo "Submitting array 1-${N_TASKS}%10 ..."

JOB_OUTPUT=$(
    sbatch \
        --array="1-${N_TASKS}%10" \
        --export=ALL,TASK_FILE="${TASK_FILE}",OUT_DIR="${OUT_DIR}",HM3_SNPLIST="${HM3_SNPLIST}" \
        "${ARRAY_SCRIPT}"
)

echo "${JOB_OUTPUT}"
echo
echo "Submission finished."