#!/bin/bash

set -euo pipefail


###############################################################################
# Paths
###############################################################################

PROJECT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC"


###############################################################################
# GWAS manifest
###############################################################################

GWAS_MANIFEST="${PROJECT_DIR}/03_GWAS/GWAS_manifest.tsv"


###############################################################################
# DAR manifest
#
# This is the combined manifest containing:
#
#   open
#   closed
#   all
###############################################################################

DAR_MANIFEST="${PROJECT_DIR}/01_DAR_beds_hg19/DAR_hg19_manifest_all_sets.tsv"


###############################################################################
# Munged GWAS
###############################################################################

MUNGED_DIR="${PROJECT_DIR}/03_GWAS/munged"


###############################################################################
# Custom LD-score directory
#
# Expected:
#
# 02_DAR_ldscores_hg19/
# ├── open/
# ├── closed/
# └── all/
###############################################################################

LDSCORE_DIR="${PROJECT_DIR}/02_DAR_ldscores_hg19"

# Pooled same-DAR-set background annotation
BACKGROUND_ANNOTATION="__DAR_BACKGROUND__"


###############################################################################
# S-LDSC result directory
#
# Will become:
#
# 04_sldsc_results/
# ├── open/
# ├── closed/
# └── all/
###############################################################################

RESULT_DIR="${PROJECT_DIR}/04_sldsc_results"


###############################################################################
# Task file
###############################################################################

TASK_FILE="${PROJECT_DIR}/03_GWAS/sldsc_tasks.tsv"


###############################################################################
# Array script
###############################################################################

ARRAY_SCRIPT="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/code/10_HD_bin_level/LDSC/pooled_background/04_run_sldsc_array_with_pooled_background.sh"

###############################################################################
# Log directory
###############################################################################

LOG_DIR="${PROJECT_DIR}/logs"


###############################################################################
# Create directories
###############################################################################

mkdir -p \
    "${RESULT_DIR}" \
    "${RESULT_DIR}/open" \
    "${RESULT_DIR}/closed" \
    "${RESULT_DIR}/all" \
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

DAR_ANNOTATION_FILE=$(mktemp)

VALID_DAR_ANNOTATION_FILE=$(mktemp)

SKIPPED_TRAIT_FILE=$(mktemp)

SKIPPED_ANNOTATION_FILE=$(mktemp)


cleanup() {

    rm -f \
        "${TRAIT_FILE}" \
        "${DAR_ANNOTATION_FILE}" \
        "${VALID_DAR_ANNOTATION_FILE}" \
        "${SKIPPED_TRAIT_FILE}" \
        "${SKIPPED_ANNOTATION_FILE}"
}


trap cleanup EXIT


###############################################################################
# Validate GWAS manifest columns
#
# Required:
#
# trait
# include
# samp_prev
# pop_prev
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

            print \
                "ERROR: Missing GWAS manifest column: " required[j] \
                > "/dev/stderr"

            exit 1
        }
    }


    exit 0
}

' "${GWAS_MANIFEST}"


###############################################################################
# Validate DAR manifest columns
#
# Required:
#
# dar_set
# annotation
###############################################################################

awk -F'\t' '

NR == 1 {

    for (i = 1; i <= NF; i++) {
        col[$i] = i
    }


    required[1] = "dar_set"
    required[2] = "annotation"


    for (j = 1; j <= 2; j++) {

        if (!(required[j] in col)) {

            print \
                "ERROR: Missing DAR manifest column: " required[j] \
                > "/dev/stderr"

            exit 1
        }
    }


    exit 0
}

' "${DAR_MANIFEST}"


###############################################################################
# Extract included GWAS traits
#
# Output:
#
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


    print \
        $col["trait"], \
        samp, \
        pop
}

' "${GWAS_MANIFEST}" |
    LC_ALL=C sort -t $'\t' -k1,1 \
    > "${TRAIT_FILE}"


###############################################################################
# Number of included GWAS
###############################################################################

N_INCLUDED_TRAITS=$(
    awk '
    NF > 0 {
        n++
    }

    END {
        print n + 0
    }
    ' "${TRAIT_FILE}"
)


if [[ "${N_INCLUDED_TRAITS}" -eq 0 ]]; then

    echo "ERROR: No GWAS traits have include=1."

    exit 1
fi

###############################################################################
# Extract unique DAR set × annotation combinations
#
# Output:
#
# open      Astrocyte
# open      Microglia
# ...
# closed    Astrocyte
# ...
# all       Astrocyte
# ...
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

{
    dar_set = $col["dar_set"]
    annotation = $col["annotation"]

    if (dar_set != "" && dar_set != "NA" && annotation != "" && annotation != "NA") {

        if (dar_set != "open" && dar_set != "closed" && dar_set != "all") {
            print "ERROR: Unexpected dar_set: " dar_set > "/dev/stderr"
            exit 1
        }

        print dar_set, annotation
    }
}
' "${DAR_MANIFEST}" |
    LC_ALL=C sort -u \
    > "${DAR_ANNOTATION_FILE}"

###############################################################################
# Number of DAR set × annotation combinations
###############################################################################

N_DAR_ANNOTATIONS=$(
    awk '
    NF > 0 {
        n++
    }

    END {
        print n + 0
    }
    ' "${DAR_ANNOTATION_FILE}"
)


if [[ "${N_DAR_ANNOTATIONS}" -eq 0 ]]; then

    echo "ERROR: No DAR set × annotation combinations were found."

    exit 1
fi


###############################################################################
# Check complete chromosome LD-score files
#
# Expected:
#
# LDSCORE_DIR/
#   open/
#       Astrocyte/
#           Astrocyte.1.l2.ldscore.gz
#           ...
#
#   closed/
#
#   all/
###############################################################################

echo
echo "Checking custom LD-score files..."


while IFS=$'\t' read -r \
    dar_set \
    annotation; do


    [[ -n "${dar_set}" ]] || continue

    [[ -n "${annotation}" ]] || continue


    ###########################################################################
    # Custom LD-score prefix
    ###########################################################################

    custom_prefix="${LDSCORE_DIR}/${dar_set}/${annotation}/${annotation}."


    annotation_complete=1

    missing_file=""


    ###########################################################################
    # Require all chromosomes
    ###########################################################################

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


    ###########################################################################
    # Record complete / incomplete annotation
    ###########################################################################

    if [[ "${annotation_complete}" -eq 1 ]]; then

        printf "%s\t%s\n" \
            "${dar_set}" \
            "${annotation}" \
            >> "${VALID_DAR_ANNOTATION_FILE}"


    else

        printf "%s\t%s\t%s\n" \
            "${dar_set}" \
            "${annotation}" \
            "${missing_file}" \
            >> "${SKIPPED_ANNOTATION_FILE}"

    fi


done < "${DAR_ANNOTATION_FILE}"


###############################################################################
# Check pooled background LD-score files
#
# Every foreground regression will condition on the background from the same
# DAR set, so fail early here instead of launching many jobs that would fail.
###############################################################################

echo
echo "Checking pooled background LD-score files..."

for dar_set in open closed all; do

    background_prefix="${LDSCORE_DIR}/${dar_set}/${BACKGROUND_ANNOTATION}/${BACKGROUND_ANNOTATION}."

    for chr in $(seq 1 22); do
        for suffix in \
            "l2.ldscore.gz" \
            "l2.M" \
            "l2.M_5_50"; do

            background_file="${background_prefix}${chr}.${suffix}"

            if [[ ! -s "${background_file}" ]]; then
                echo "ERROR: Missing pooled background LD-score file:"
                echo "${background_file}"
                exit 1
            fi
        done
    done

    echo "Background complete: ${dar_set}"

done


###############################################################################
# Number of valid DAR set × annotation combinations
###############################################################################

N_VALID_DAR_ANNOTATIONS=$(
    awk '
    NF > 0 {
        n++
    }

    END {
        print n + 0
    }
    ' "${VALID_DAR_ANNOTATION_FILE}"
)


if [[ "${N_VALID_DAR_ANNOTATIONS}" -eq 0 ]]; then

    echo "ERROR: No DAR annotation has a complete set of LD-score files."

    exit 1
fi


###############################################################################
# Initialize task file
#
# This task file contains a header.
#
# Columns:
#
# trait
# dar_set
# annotation
# sumstats
# samp_prev
# pop_prev
# result_prefix
###############################################################################

printf \
    "trait\tdar_set\tannotation\tsumstats\tsamp_prev\tpop_prev\tresult_prefix\n" \
    > "${TASK_FILE}"


###############################################################################
# Construct:
#
# GWAS × DAR set × annotation
#
# Cartesian product
###############################################################################

N_READY_TRAITS=0


while IFS=$'\t' read -r \
    trait \
    samp_prev \
    pop_prev; do


    [[ -n "${trait}" ]] || continue


    ###########################################################################
    # Munged GWAS file
    ###########################################################################

    sumstats="${MUNGED_DIR}/${trait}/${trait}.sumstats.gz"


    if [[ ! -s "${sumstats}" ]]; then

        printf "%s\t%s\n" \
            "${trait}" \
            "${sumstats}" \
            >> "${SKIPPED_TRAIT_FILE}"

        continue
    fi


    N_READY_TRAITS=$((N_READY_TRAITS + 1))


    ###########################################################################
    # Loop over all valid DAR set × annotation pairs
    ###########################################################################

    while IFS=$'\t' read -r \
        dar_set \
        annotation; do


        [[ -n "${dar_set}" ]] || continue

        [[ -n "${annotation}" ]] || continue


        #######################################################################
        # Result directory
        #
        # Example:
        #
        # 04_sldsc_results/
        #     open/
        #         AD/
        #             AD__Astrocyte.results
        #######################################################################

        mkdir -p \
            "${RESULT_DIR}/${dar_set}/${trait}"


        #######################################################################
        # Result prefix
        #######################################################################

        result_prefix="${RESULT_DIR}/${dar_set}/${trait}/${trait}__${annotation}"


        #######################################################################
        # Add task
        #######################################################################

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${trait}" \
            "${dar_set}" \
            "${annotation}" \
            "${sumstats}" \
            "${samp_prev}" \
            "${pop_prev}" \
            "${result_prefix}" \
            >> "${TASK_FILE}"


    done < "${VALID_DAR_ANNOTATION_FILE}"


done < "${TRAIT_FILE}"


###############################################################################
# Count tasks
###############################################################################

N_TASKS=$(
    awk '
    NR > 1 && NF > 0 {
        n++
    }

    END {
        print n + 0
    }
    ' "${TASK_FILE}"
)


if [[ "${N_TASKS}" -eq 0 ]]; then

    echo "ERROR: No S-LDSC tasks were created."

    exit 1
fi


###############################################################################
# Expected number of tasks
###############################################################################

EXPECTED_TASKS=$((N_READY_TRAITS * N_VALID_DAR_ANNOTATIONS))


if [[ "${N_TASKS}" -ne "${EXPECTED_TASKS}" ]]; then

    echo "ERROR: Unexpected number of S-LDSC tasks."

    echo "Ready traits:             ${N_READY_TRAITS}"
    echo "Valid DAR annotations:    ${N_VALID_DAR_ANNOTATIONS}"
    echo "Expected tasks:           ${EXPECTED_TASKS}"
    echo "Observed tasks:           ${N_TASKS}"

    exit 1
fi


###############################################################################
# Summary
###############################################################################

echo
echo "============================================================"
echo "S-LDSC TASK SUMMARY"
echo "============================================================"

echo
echo "Included traits in manifest: ${N_INCLUDED_TRAITS}"

echo "Traits with munged files:     ${N_READY_TRAITS}"

echo "DAR set × annotations:        ${N_DAR_ANNOTATIONS}"

echo "Complete DAR annotations:     ${N_VALID_DAR_ANNOTATIONS}"

echo "S-LDSC tasks:                 ${N_TASKS}"


###############################################################################
# Counts by DAR set
###############################################################################

echo
echo "Complete annotations by DAR set:"


awk -F'\t' '

{
    n[$1]++
}

END {

    print "open:   " n["open"] + 0

    print "closed: " n["closed"] + 0

    print "all:    " n["all"] + 0
}

' "${VALID_DAR_ANNOTATION_FILE}"


###############################################################################
# Task counts by DAR set
###############################################################################

echo
echo "S-LDSC tasks by DAR set:"


awk -F'\t' '

NR > 1 {
    n[$2]++
}

END {

    print "open:   " n["open"] + 0

    print "closed: " n["closed"] + 0

    print "all:    " n["all"] + 0
}

' "${TASK_FILE}"


###############################################################################
# Task file
###############################################################################

echo
echo "Task file:"

echo "${TASK_FILE}"


###############################################################################
# Ready GWAS traits
###############################################################################

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


###############################################################################
# Complete annotations
###############################################################################

echo
echo "Complete DAR set × annotation combinations:"


column -t -s $'\t' \
    "${VALID_DAR_ANNOTATION_FILE}"


###############################################################################
# Report skipped GWAS traits
###############################################################################

echo
echo "Included GWAS traits without munged output:"


if [[ -s "${SKIPPED_TRAIT_FILE}" ]]; then

    column -t -s $'\t' \
        "${SKIPPED_TRAIT_FILE}"

else

    echo "None"

fi


###############################################################################
# Report incomplete DAR annotations
###############################################################################

echo
echo "DAR annotations with incomplete LD-score files:"


if [[ -s "${SKIPPED_ANNOTATION_FILE}" ]]; then

    column -t -s $'\t' \
        "${SKIPPED_ANNOTATION_FILE}"

else

    echo "None"

fi


###############################################################################
# Preview task file
###############################################################################

echo
echo "Task preview:"


head -10 "${TASK_FILE}" |
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