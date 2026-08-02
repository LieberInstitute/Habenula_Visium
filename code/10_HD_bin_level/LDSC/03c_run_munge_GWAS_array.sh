#!/bin/bash

#SBATCH --job-name=03_munge_GWAS
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/03_munge_GWAS_%a.out
#SBATCH --error=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LDSC/logs/03_munge_GWAS_%a.out

set -eo pipefail

###############################################################################
# Activate Python 2 LDSC environment
###############################################################################

unset PYTHONHOME
unset PYTHONPATH

set +u

source \
    /jhpce/shared/jhpce/core/conda/miniconda3-24.3.0/etc/profile.d/conda.sh

conda activate py2_env

set -u

###############################################################################
# Paths
###############################################################################

LDSC_DIR="/users/cliu3/Thesis/ldsc_clean"
MUNGE="${LDSC_DIR}/munge_sumstats.py"

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
# Check exported variables and files
###############################################################################

for var in TASK_FILE OUT_DIR HM3_SNPLIST; do

    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: ${var} is not defined."
        exit 1
    fi

done

for file in \
    "${TASK_FILE}" \
    "${HM3_SNPLIST}" \
    "${MUNGE}"; do

    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Missing or empty file:"
        echo "${file}"
        exit 1
    fi

done

###############################################################################
# Read one task
#
# Columns:
#  1 trait
#  2 sumstats_file
#  3 file_format
#  4 delimiter
#  5 header_line_number
#  6 snp_col
#  7 p_col
#  8 p_transform
#  9 stat_mode
# 10 signed_col
# 11 signed_null
# 12 signed_sumstats_arg
# 13 use_signed_stat
# 14 use_unsigned_from_p
# 15 a1_col
# 16 a2_col
# 17 no_alleles
# 18 a1_inc
# 19 n_col
# 20 n_case_col
# 21 n_control_col
# 22 N_fixed
# 23 sample_size_mode
# 24 info_col
# 25 frq_col
# 26 samp_prev
# 27 pop_prev
###############################################################################

task_line=$(
    sed -n "${ARRAY_ID}p" "${TASK_FILE}"
)

if [[ -z "${task_line}" ]]; then
    echo "ERROR: Could not read task ${ARRAY_ID}."
    exit 1
fi

IFS=$'\t' read -r \
    TRAIT \
    SUMSTATS \
    FILE_FORMAT \
    DELIMITER \
    HEADER_LINE_NUMBER \
    SNP_COL \
    P_COL \
    P_TRANSFORM \
    STAT_MODE \
    SIGNED_COL \
    SIGNED_NULL \
    SIGNED_SUMSTATS_ARG \
    USE_SIGNED_STAT \
    USE_UNSIGNED_FROM_P \
    A1_COL \
    A2_COL \
    NO_ALLELES \
    A1_INC \
    N_COL \
    N_CASE_COL \
    N_CONTROL_COL \
    N_FIXED \
    SAMPLE_SIZE_MODE \
    INFO_COL \
    FRQ_COL \
    SAMP_PREV \
    POP_PREV \
    <<< "${task_line}"

###############################################################################
# Helpers
###############################################################################

is_present() {
    [[ -n "${1:-}" && "${1}" != "NA" ]]
}

stream_sumstats() {

    case "${SUMSTATS}" in

        *.gz|*.bgz)
            gzip -cd "${SUMSTATS}"
            ;;

        *.bz2)
            bzip2 -cd "${SUMSTATS}"
            ;;

        *.xz)
            xz -cd "${SUMSTATS}"
            ;;

        *)
            cat "${SUMSTATS}"
            ;;
    esac
}

###############################################################################
# Display task
###############################################################################

echo
echo "Trait:               ${TRAIT}"
echo "Summary statistics:  ${SUMSTATS}"
echo "File format:         ${FILE_FORMAT}"
echo "Delimiter:           ${DELIMITER}"
echo "Header line:         ${HEADER_LINE_NUMBER}"
echo "SNP column:          ${SNP_COL}"
echo "P column:            ${P_COL}"
echo "P transformation:    ${P_TRANSFORM}"
echo "Statistic mode:      ${STAT_MODE}"
echo "Signed column:       ${SIGNED_COL}"
echo "Signed null:         ${SIGNED_NULL}"
echo "Sample-size mode:    ${SAMPLE_SIZE_MODE}"
echo "N column:            ${N_COL}"
echo "N case column:       ${N_CASE_COL}"
echo "N control column:    ${N_CONTROL_COL}"
echo "Fixed N:             ${N_FIXED}"
echo "No alleles:          ${NO_ALLELES}"
echo "A1 increment:        ${A1_INC}"
echo

###############################################################################
# Validate task
###############################################################################

if [[ ! -s "${SUMSTATS}" ]]; then
    echo "ERROR: Summary-statistics file is missing or empty:"
    echo "${SUMSTATS}"
    exit 1
fi

if ! [[ "${HEADER_LINE_NUMBER}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Invalid header line number: ${HEADER_LINE_NUMBER}"
    exit 1
fi

if [[ "${HEADER_LINE_NUMBER}" -lt 1 ]]; then
    echo "ERROR: Header line number must be at least 1."
    exit 1
fi

if ! is_present "${SNP_COL}"; then
    echo "ERROR: No SNP column was specified for ${TRAIT}."
    exit 1
fi

if ! is_present "${P_COL}"; then
    echo "ERROR: No P-value column was specified for ${TRAIT}."
    exit 1
fi

case "${P_TRANSFORM}" in

    identity|neglog10)
        ;;

    *)
        echo "ERROR: Unsupported P-value transformation:"
        echo "${P_TRANSFORM}"
        exit 1
        ;;
esac

case "${STAT_MODE}" in

    signed)

        if ! is_present "${SIGNED_COL}"; then
            echo "ERROR: Signed mode selected without signed_col."
            exit 1
        fi

        if ! is_present "${SIGNED_NULL}"; then
            echo "ERROR: Signed mode selected without signed_null."
            exit 1
        fi
        ;;

    unsigned_from_p)
        ;;

    *)
        echo "ERROR: Unsupported statistic mode:"
        echo "${STAT_MODE}"
        exit 1
        ;;
esac

case "${SAMPLE_SIZE_MODE}" in

    N_column)

        if ! is_present "${N_COL}"; then
            echo "ERROR: N_column mode selected without n_col."
            exit 1
        fi
        ;;

    case_control_columns)

        if ! is_present "${N_CASE_COL}" ||
           ! is_present "${N_CONTROL_COL}"; then

            echo "ERROR: case_control_columns requires both N columns."
            exit 1
        fi
        ;;

    fixed_N)

        if ! is_present "${N_FIXED}"; then
            echo "ERROR: fixed_N mode selected without N_fixed."
            exit 1
        fi
        ;;

    *)
        echo "ERROR: Unsupported sample-size mode:"
        echo "${SAMPLE_SIZE_MODE}"
        exit 1
        ;;
esac

###############################################################################
# Output paths
###############################################################################

TRAIT_DIR="${OUT_DIR}/${TRAIT}"
OUT_PREFIX="${TRAIT_DIR}/${TRAIT}"

STANDARD_FILE="${OUT_PREFIX}.standardized.tsv.gz"

JOB_TAG="${JOB_ID}_${ARRAY_ID}"

TMP_STANDARD="${TRAIT_DIR}/.${TRAIT}.${JOB_TAG}.standardized.tmp.tsv.gz"

mkdir -p "${TRAIT_DIR}"

cleanup() {
    rm -f "${TMP_STANDARD}"
}

trap cleanup EXIT

###############################################################################
# Standardize input
#
# Output:
#   unsigned: SNP P N
#   signed:   SNP P N SIGNED
#
# The first column of w_hm3.snplist is used as the SNP whitelist.
###############################################################################

echo "Standardizing summary statistics..."

stream_sumstats |
tail -n +"${HEADER_LINE_NUMBER}" |
awk \
    -v hm3_file="${HM3_SNPLIST}" \
    -v snp_name="${SNP_COL}" \
    -v p_name="${P_COL}" \
    -v p_transform="${P_TRANSFORM}" \
    -v stat_mode="${STAT_MODE}" \
    -v signed_name="${SIGNED_COL}" \
    -v sample_mode="${SAMPLE_SIZE_MODE}" \
    -v n_name="${N_COL}" \
    -v n_case_name="${N_CASE_COL}" \
    -v n_control_name="${N_CONTROL_COL}" \
    -v fixed_n="${N_FIXED}" \
    '
BEGIN {
    FS = "[[:space:]]+"
    OFS = "\t"
    log10_value = log(10)
}

function normalize_name(x, y) {
    y = toupper(x)
    sub(/^#+/, "", y)
    gsub(/[^A-Z0-9]+/, "_", y)
    gsub(/_+/, "_", y)
    sub(/^_+/, "", y)
    sub(/_+$/, "", y)
    return y
}

function clean_number(x, y) {
    y = x
    gsub(/,/, "", y)
    return y
}

function is_number(x, y) {
    y = clean_number(x)
    return (y ~ /^[-+]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eE][-+]?[0-9]+)?$/)
}

FILENAME == hm3_file {
    if ($1 != "" && normalize_name($1) != "SNP") {
        hm3[$1] = 1
    }
    next
}

FILENAME == "-" && FNR == 1 {

    target_snp = normalize_name(snp_name)
    target_p = normalize_name(p_name)
    target_signed = normalize_name(signed_name)
    target_n = normalize_name(n_name)
    target_n_case = normalize_name(n_case_name)
    target_n_control = normalize_name(n_control_name)

    for (i = 1; i <= NF; i++) {

        header_i = normalize_name($i)

        if (header_i == target_snp) {
            snp_i = i
        }

        if (header_i == target_p) {
            p_i = i
        }

        if (stat_mode == "signed" && header_i == target_signed) {
            signed_i = i
        }

        if (sample_mode == "N_column" && header_i == target_n) {
            n_i = i
        }

        if (sample_mode == "case_control_columns" &&
            header_i == target_n_case) {

            n_case_i = i
        }

        if (sample_mode == "case_control_columns" &&
            header_i == target_n_control) {

            n_control_i = i
        }
    }

    if (!snp_i) {
        print "ERROR: SNP column not found: " snp_name > "/dev/stderr"
        exit 2
    }

    if (!p_i) {
        print "ERROR: P-value column not found: " p_name > "/dev/stderr"
        exit 2
    }

    if (stat_mode == "signed" && !signed_i) {
        print "ERROR: Signed column not found: " signed_name > "/dev/stderr"
        exit 2
    }

    if (sample_mode == "N_column" && !n_i) {
        print "ERROR: N column not found: " n_name > "/dev/stderr"
        exit 2
    }

    if (sample_mode == "case_control_columns" &&
        (!n_case_i || !n_control_i)) {

        print "ERROR: Case/control N columns not found." > "/dev/stderr"
        exit 2
    }

    if (stat_mode == "signed") {
        print "SNP", "P", "N", "SIGNED"
    } else {
        print "SNP", "P", "N"
    }

    next
}

FILENAME == "-" {

    snp = $snp_i

    if (snp == "" || snp == "." || snp == "NA") {
        next
    }

    # rsIDs are retained; nonmatching identifiers are naturally removed.
    if (!(snp in hm3)) {
        next
    }

    ###########################################################################
    # P value
    ###########################################################################

    p_raw = clean_number($p_i)

    if (!is_number(p_raw)) {
        next
    }

    p_raw += 0

    if (p_transform == "neglog10") {

        if (p_raw < 0) {
            next
        }

        p = exp(-p_raw * log10_value)

        if (p == 0) {
            p = 1e-300
        }

    } else {

        p = p_raw

        if (p == 0) {
            p = 1e-300
        }
    }

    if (p <= 0 || p > 1) {
        next
    }

    ###########################################################################
    # Sample size
    ###########################################################################

    if (sample_mode == "N_column") {

        n_raw = clean_number($n_i)

        if (!is_number(n_raw)) {
            next
        }

        n_value = n_raw + 0

    } else if (sample_mode == "case_control_columns") {

        n_case_raw = clean_number($n_case_i)
        n_control_raw = clean_number($n_control_i)

        if (!is_number(n_case_raw) || !is_number(n_control_raw)) {
            next
        }

        n_case = n_case_raw + 0
        n_control = n_control_raw + 0

        if (n_case <= 0 || n_control <= 0) {
            next
        }

        n_value = 4 * n_case * n_control / (n_case + n_control)

    } else {

        n_raw = clean_number(fixed_n)

        if (!is_number(n_raw)) {
            next
        }

        n_value = n_raw + 0
    }

    if (n_value <= 0) {
        next
    }

    ###########################################################################
    # Signed statistic
    ###########################################################################

    if (stat_mode == "signed") {

        signed_raw = clean_number($signed_i)

        if (!is_number(signed_raw)) {
            next
        }

        signed_value = signed_raw + 0
    }

    ###########################################################################
    # Output
    ###########################################################################

    if (stat_mode == "signed") {

        printf "%s\t%.17g\t%.17g\t%.17g\n", \
            snp, p, n_value, signed_value

    } else {

        printf "%s\t%.17g\t%.17g\n", \
            snp, p, n_value
    }
}
' \
    "${HM3_SNPLIST}" \
    - |
gzip -c > "${TMP_STANDARD}"

###############################################################################
# Standardized-file QC
###############################################################################

if [[ ! -s "${TMP_STANDARD}" ]]; then
    echo "ERROR: Standardized output was not created."
    exit 1
fi

N_STANDARD=$(
    gzip -cd "${TMP_STANDARD}" |
    awk 'NR > 1 {n++} END {print n + 0}'
)

echo
echo "HM3 SNPs retained: ${N_STANDARD}"

if [[ "${N_STANDARD}" -eq 0 ]]; then
    echo "ERROR: No SNPs were retained for ${TRAIT}."
    exit 1
fi

if [[ "${N_STANDARD}" -lt 200000 ]]; then
    echo "WARNING: Fewer than 200,000 SNPs were retained."
fi

mv \
    "${TMP_STANDARD}" \
    "${STANDARD_FILE}"

echo
echo "Standardized file:"
echo "${STANDARD_FILE}"

echo
echo "Standardized preview:"

gzip -cd "${STANDARD_FILE}" |
head || true

###############################################################################
# Run munge_sumstats.py
###############################################################################

CMD=(
    python
    "${MUNGE}"

    --sumstats
    "${STANDARD_FILE}"

    --out
    "${OUT_PREFIX}"

    --snp
    SNP

    --p
    P

    --N-col
    N

    --no-alleles

    --chunksize
    500000
)

if [[ "${STAT_MODE}" == "signed" ]]; then

    CMD+=(
        --signed-sumstats
        "SIGNED,${SIGNED_NULL}"
    )

else

    CMD+=(
        --a1-inc
    )
fi

echo
echo "Running command:"

printf "%q " "${CMD[@]}"
printf "\n\n"

"${CMD[@]}"

###############################################################################
# Output QC
###############################################################################

for file in \
    "${OUT_PREFIX}.sumstats.gz" \
    "${OUT_PREFIX}.log"; do

    if [[ ! -s "${file}" ]]; then
        echo "ERROR: Expected output is missing:"
        echo "${file}"
        exit 1
    fi

done

echo
echo "Munged output preview:"

gzip -cd "${OUT_PREFIX}.sumstats.gz" |
head || true

echo
echo "LDSC munging summary:"

grep -E \
    "Read .* SNPs|SNPs remain|Mean chi|Lambda GC|Max chi|WARNING" \
    "${OUT_PREFIX}.log" || true

echo
echo "Finished successfully:"
echo "${OUT_PREFIX}.sumstats.gz"