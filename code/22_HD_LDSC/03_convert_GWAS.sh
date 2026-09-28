# check GWAS files for LDSC

#!/bin/bash

set -euo pipefail

###############################################################################
# Paths
###############################################################################

ROOT="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA"

PROJECT_DIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC"

OUT_DIR="${PROJECT_DIR}/03_GWAS"

OUT_FILE="${OUT_DIR}/GWAS_inventory.tsv"

mkdir -p "${OUT_DIR}"

###############################################################################
# Output header
###############################################################################

printf "trait\tsumstats_file\tfile_name\tfile_size_bytes\theader\n" \
    > "${OUT_FILE}"

###############################################################################
# Directories excluded from the GWAS inventory
#
# MDD is intentionally excluded. We are not using:
#   MDD/p_values.tsv
###############################################################################

exclude_directory() {

    local trait="$1"

    case "${trait}" in

        first_test|\
        gene_sets|\
        habenula_pilot_gwas|\
        logs|\
        MDD|\
        new_cluster|\
        OUD_unused|\
        RNA|\
        trios)

            return 0
            ;;

        *)

            return 1
            ;;
    esac
}

###############################################################################
# Traits whose files are specified explicitly
###############################################################################

is_explicit_trait() {

    local trait="$1"

    case "${trait}" in

        MDD2019|\
        panic|\
        SCZ|\
        SUD2020)

            return 0
            ;;

        *)

            return 1
            ;;
    esac
}

###############################################################################
# Read the first line of compressed or uncompressed files
###############################################################################

read_header() {

    local file="$1"

    case "${file}" in

        *.gz|*.bgz)

            gzip -cd "${file}" 2>/dev/null |
                head -n 1 || true
            ;;

        *.bz2)

            bzip2 -cd "${file}" 2>/dev/null |
                head -n 1 || true
            ;;

        *.xz)

            xz -cd "${file}" 2>/dev/null |
                head -n 1 || true
            ;;

        *)

            head -n 1 "${file}" 2>/dev/null || true
            ;;
    esac
}

###############################################################################
# Add one GWAS file to the inventory
###############################################################################

add_inventory_row() {

    local trait="$1"
    local file="$2"

    if [[ ! -s "${file}" ]]; then
        echo "ERROR: GWAS file does not exist or is empty:"
        echo "${file}"
        return 1
    fi

    local file_name
    local file_size
    local header
    local header_clean

    file_name=$(basename "${file}")

    file_size=$(
        stat -Lc "%s" "${file}"
    )

    header=$(
        read_header "${file}"
    )

    if [[ -z "${header}" ]]; then
        echo "ERROR: Could not read header from:"
        echo "${file}"
        return 1
    fi

    # Replace tabs with "|" so that the inventory itself remains valid TSV.
    header_clean=$(
        printf "%s\n" "${header}" |
            sed 's/\r$//' |
            tr '\t' '|' |
            sed 's/^ *//; s/ *$//'
    )

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "${trait}" \
        "${file}" \
        "${file_name}" \
        "${file_size}" \
        "${header_clean}" \
        >> "${OUT_FILE}"

    echo "Added: ${trait} -> ${file_name}"
}

###############################################################################
# Automatically find standard *.txt.gz GWAS files
###############################################################################

echo "Searching standard GWAS files."
echo

for trait_dir in "${ROOT}"/*; do

    [[ -d "${trait_dir}" ]] || continue

    trait=$(basename "${trait_dir}")

    if exclude_directory "${trait}"; then
        echo "Skipping excluded directory: ${trait}"
        continue
    fi

    # These traits are added using explicit paths below.
    if is_explicit_trait "${trait}"; then
        echo "Using explicit file for: ${trait}"
        continue
    fi

    echo "Scanning: ${trait}"

    mapfile -t candidates < <(
        find -L "${trait_dir}" \
            -mindepth 1 \
            -maxdepth 1 \
            -type f \
            -name "*.txt.gz" \
            2>/dev/null |
            sort
    )

    if [[ "${#candidates[@]}" -eq 0 ]]; then
        echo "WARNING: no .txt.gz file found for ${trait}"
        continue
    fi

    if [[ "${#candidates[@]}" -gt 1 ]]; then

        echo "WARNING: multiple .txt.gz files found for ${trait}:"

        printf "  %s\n" "${candidates[@]}"

        echo "All files will be added. Review the inventory before continuing."
    fi

    for file in "${candidates[@]}"; do
        add_inventory_row "${trait}" "${file}"
    done

done

###############################################################################
# Add explicitly specified GWAS files
###############################################################################

echo
echo "Adding explicitly specified GWAS files."
echo

SPECIAL_TRAITS=(
    "MDD2019"
    "panic"
    "SCZ"
    "SUD2020"
)

SPECIAL_FILES=(
    "${ROOT}/MDD2019/MDD2019.txt"
    "${ROOT}/panic/panic.vcf.tsv.gz"
    "${ROOT}/SCZ/SCZ.vcf.tsv.gz"
    "${ROOT}/SUD2020/SUD2020.tbl"
)

if [[ "${#SPECIAL_TRAITS[@]}" -ne "${#SPECIAL_FILES[@]}" ]]; then
    echo "ERROR: SPECIAL_TRAITS and SPECIAL_FILES have different lengths."
    exit 1
fi

for i in "${!SPECIAL_TRAITS[@]}"; do

    trait="${SPECIAL_TRAITS[$i]}"
    file="${SPECIAL_FILES[$i]}"

    add_inventory_row "${trait}" "${file}"

done

###############################################################################
# Sort inventory by trait and file name
###############################################################################

{
    head -n 1 "${OUT_FILE}"

    tail -n +2 "${OUT_FILE}" |
        sort \
            -t $'\t' \
            -k1,1 \
            -k3,3

} > "${OUT_FILE}.tmp"

mv "${OUT_FILE}.tmp" "${OUT_FILE}"

###############################################################################
# Check for duplicate trait entries
###############################################################################

duplicate_traits=$(
    awk -F'\t' '
        NR > 1 {
            count[$1]++
        }

        END {
            for (trait in count) {
                if (count[trait] > 1) {
                    print trait
                }
            }
        }
    ' "${OUT_FILE}"
)

if [[ -n "${duplicate_traits}" ]]; then

    echo
    echo "WARNING: multiple input files were found for these traits:"
    echo "${duplicate_traits}"
    echo
    echo "Review GWAS_inventory.tsv before continuing."

fi

###############################################################################
# Summary
###############################################################################

echo
echo "Inventory written to:"
echo "${OUT_FILE}"

echo
echo "Files found:"

awk -F'\t' '
    NR > 1 {
        print $1 "\t" $3
    }
' "${OUT_FILE}" |
    column -t -s $'\t'

echo
echo "Candidate counts by trait:"

awk -F'\t' '
    NR > 1 {
        count[$1]++
    }

    END {
        for (trait in count) {
            print trait, count[trait]
        }
    }
' "${OUT_FILE}" |
    sort

N_FILES=$(
    awk 'END {print NR - 1}' "${OUT_FILE}"
)

N_TRAITS=$(
    awk -F'\t' '
        NR > 1 {
            trait[$1] = 1
        }

        END {
            print length(trait)
        }
    ' "${OUT_FILE}"
)

echo
echo "Total GWAS files: ${N_FILES}"
echo "Total GWAS traits: ${N_TRAITS}"

###############################################################################
# Confirm that MDD is not included
###############################################################################

if awk -F'\t' '
    NR > 1 && $1 == "MDD" {
        found = 1
    }

    END {
        exit(found ? 0 : 1)
    }
' "${OUT_FILE}"; then

    echo "ERROR: MDD was unexpectedly included in the inventory."
    exit 1

else

    echo "Confirmed: MDD is excluded."

fi

###############################################################################
# Report directories without selected files
###############################################################################

echo
echo "Included trait directories without a selected GWAS file:"

missing_count=0

for trait_dir in "${ROOT}"/*; do

    [[ -d "${trait_dir}" ]] || continue

    trait=$(basename "${trait_dir}")

    if exclude_directory "${trait}"; then
        continue
    fi

    if ! awk -F'\t' -v trait="${trait}" '
        NR > 1 && $1 == trait {
            found = 1
        }

        END {
            exit(found ? 0 : 1)
        }
    ' "${OUT_FILE}"; then

        echo "${trait}"
        missing_count=$((missing_count + 1))
    fi

done

if [[ "${missing_count}" -eq 0 ]]; then
    echo "None"
fi

echo
echo "Finished."