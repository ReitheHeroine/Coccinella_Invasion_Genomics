#!/usr/bin/env bash
: <<'HEADER'
title: 'diversity.sh'
project: 'BIOL624 Final Project: Selection Detection in Lady Beetles'
author: 'Reina Hastings'
contact: 'reinahastings13@gmail.com'
date created: 11/20/2025
last modified: 11/20/2025
purpose: 'Compute per-chromosome Tajima’s D and nucleotide diversity (π) with vcftools, globally or per population.'
inputs:
  - Directory of filtered per-chromosome VCF files (e.g., VARIANTS_BY_CHR/FILTERED_CHR_VCFS/)
  - Optional population file (one sample ID per line) for per-population stats
  - Window size in bp (e.g., 10000)
outputs:
  - <outdir>/tajima/<chrom>[_<poplabel>].Tajima.D
  - <outdir>/pi/<chrom>[_<poplabel>].windowed.pi
  - <outdir>/pi_sites/<chrom>[_<poplabel>].sites.pi
required tools:
  - vcftools
example usage:
    chmod +x diversity.sh

    # Global (all individuals)
    ./diversity.sh -i VARIANTS_BY_CHR/FILTERED_CHR_VCFS -w 10000 -o results

    # Per-population (e.g., Europe only)
    ./diversity.sh -i VARIANTS_BY_CHR/FILTERED_CHR_VCFS \
                         -w 10000 \
                         -o results \
                         -p metadata/pop_Europe.txt
notes:
  - Uses non-overlapping windows of fixed physical size (bp) for windowed stats.
  - When a population file is provided, vcftools is run with --keep <popfile>.
HEADER

set -euo pipefail


# argument parsing
usage() {
    cat <<EOF
Usage: diversity.sh -i <vcf_dir> [-w window_bp] [-o outdir] [-p popfile] [-g vcf_glob]

Required:
  -i  Directory containing per-chromosome VCFs (e.g., VARIANTS_BY_CHR/FILTERED_CHR_VCFS)

Optional:
  -w  Window size in bp for Tajima's D and windowed pi (default: 10000)
  -o  Output base directory (default: results)
  -p  Population file (one sample ID per line); if provided, stats are per population
  -g  Glob pattern for VCFs inside -i (default: '*.vcf.gz')
  -h  Show this help message

Examples:
  # Global (all individuals)
  diversity.sh -i VARIANTS_BY_CHR/FILTERED_CHR_VCFS -w 10000 -o ../results

  # Europe-only stats
  diversity.sh -i VARIANTS_BY_CHR/FILTERED_CHR_VCFS -w 10000 -o ../results -p metadata/pop_Europe.txt
EOF
}

VCF_DIR=""
WIN=10000
OUTDIR="results"
POP_FILE=""
VCF_GLOB="*.vcf.gz"

while getopts ":i:w:o:p:g:h" opt; do
    case "${opt}" in
        i) VCF_DIR="${OPTARG}" ;;
        w) WIN="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;
        p) POP_FILE="${OPTARG}" ;;
        g) VCF_GLOB="${OPTARG}" ;;
        h) usage; exit 0 ;;
        \?) echo "ERROR: Invalid option -${OPTARG}" >&2; usage; exit 1 ;;
        :)  echo "ERROR: Option -${OPTARG} requires an argument." >&2; usage; exit 1 ;;
    esac
done

if [[ -z "${VCF_DIR}" ]]; then
    echo "ERROR: -i <vcf_dir> is required." >&2
    usage
    exit 1
fi

if [[ ! -d "${VCF_DIR}" ]]; then
    echo "ERROR: VCF directory not found: ${VCF_DIR}" >&2
    exit 1
fi

# Check vcftools is available
if ! command -v vcftools >/dev/null 2>&1; then
    echo "ERROR: vcftools not found in PATH." >&2
    exit 1
fi

# setup output dirs & logging
TAJ_DIR="${OUTDIR}/tajima"
PI_DIR="${OUTDIR}/pi"
PISITE_DIR="${OUTDIR}/pi_sites"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "${TAJ_DIR}" "${PI_DIR}" "${PISITE_DIR}" "${LOG_DIR}"

# create a timestamped run log
STAMP=$(date +"%Y%m%d_%H%M%S")
RUN_LOG="${LOG_DIR}/diversity_${STAMP}.log"

echo "Running diversity.sh" | tee -a "${RUN_LOG}"
echo "  VCF_DIR:   ${VCF_DIR}" | tee -a "${RUN_LOG}"
echo "  WINDOW:    ${WIN} bp" | tee -a "${RUN_LOG}"
echo "  OUTDIR:    ${OUTDIR}" | tee -a "${RUN_LOG}"
echo "  POP_FILE:  ${POP_FILE:-<none>}" | tee -a "${RUN_LOG}"
echo "  VCF_GLOB:  ${VCF_GLOB}" | tee -a "${RUN_LOG}"

# population handling
SUFFIX=""

if [[ -n "${POP_FILE}" ]]; then
    if [[ ! -f "${POP_FILE}" ]]; then
        echo "ERROR: Population file not found: ${POP_FILE}" | tee -a "${RUN_LOG}" >&2
        exit 1
    fi
    POP_LABEL=$(basename "${POP_FILE}")
    POP_LABEL="${POP_LABEL%.*}"
    SUFFIX="_${POP_LABEL}"
    echo "Per-population mode: using --keep ${POP_FILE} (label: ${POP_LABEL})" | tee -a "${RUN_LOG}"
else
    echo "Global mode: using all individuals in each VCF." | tee -a "${RUN_LOG}"
fi

# main loop over VCFs
shopt -s nullglob
VCF_LIST=("${VCF_DIR}"/${VCF_GLOB})

if [[ ${#VCF_LIST[@]} -eq 0 ]]; then
    echo "ERROR: No VCFs matching pattern '${VCF_GLOB}' found in ${VCF_DIR}" | tee -a "${RUN_LOG}" >&2
    exit 1
fi

for VCF in "${VCF_LIST[@]}"; do
    BASENAME=$(basename "${VCF}")
    CHR="${BASENAME%.vcf.gz}"
    CHR="${CHR%.vcf}"

    echo "Processing ${BASENAME} (chrom label: ${CHR}) ..." | tee -a "${RUN_LOG}"

    if [[ -n "${POP_FILE}" ]]; then
        echo "  vcftools --gzvcf ${VCF} --keep ${POP_FILE} --TajimaD ${WIN} --out ${TAJ_DIR}/${CHR}${SUFFIX}" >> "${RUN_LOG}"
        vcftools \
            --gzvcf "${VCF}" \
            --keep "${POP_FILE}" \
            --TajimaD "${WIN}" \
            --out "${TAJ_DIR}/${CHR}${SUFFIX}" \
            >> "${RUN_LOG}" 2>&1

        echo "  vcftools --gzvcf ${VCF} --keep ${POP_FILE} --window-pi ${WIN} --out ${PI_DIR}/${CHR}${SUFFIX}" >> "${RUN_LOG}"
        vcftools \
            --gzvcf "${VCF}" \
            --keep "${POP_FILE}" \
            --window-pi "${WIN}" \
            --out "${PI_DIR}/${CHR}${SUFFIX}" \
            >> "${RUN_LOG}" 2>&1

        echo "  vcftools --gzvcf ${VCF} --keep ${POP_FILE} --site-pi --out ${PISITE_DIR}/${CHR}${SUFFIX}" >> "${RUN_LOG}"
        vcftools \
            --gzvcf "${VCF}" \
            --keep "${POP_FILE}" \
            --site-pi \
            --out "${PISITE_DIR}/${CHR}${SUFFIX}" \
            >> "${RUN_LOG}" 2>&1

    else
        echo "  vcftools --gzvcf ${VCF} --TajimaD ${WIN} --out ${TAJ_DIR}/${CHR}" >> "${RUN_LOG}"
        vcftools \
            --gzvcf "${VCF}" \
            --TajimaD "${WIN}" \
            --out "${TAJ_DIR}/${CHR}" \
            >> "${RUN_LOG}" 2>&1

        echo "  vcftools --gzvcf ${VCF} --window-pi ${WIN} --out ${PI_DIR}/${CHR}" >> "${RUN_LOG}"
        vcftools \
            --gzvcf "${VCF}" \
            --window-pi "${WIN}" \
            --out "${PI_DIR}/${CHR}" \
            >> "${RUN_LOG}" 2>&1

        echo "  vcftools --gzvcf ${VCF} --site-pi --out ${PISITE_DIR}/${CHR}" >> "${RUN_LOG}"
        vcftools \
            --gzvcf "${VCF}" \
            --site-pi \
            --out "${PISITE_DIR}/${CHR}" \
            >> "${RUN_LOG}" 2>&1
    fi

    echo "  Finished ${CHR}${SUFFIX}" | tee -a "${RUN_LOG}"
done

echo "All done. See:" | tee -a "${RUN_LOG}"
echo "  Tajima's D : ${TAJ_DIR}" | tee -a "${RUN_LOG}"
echo "  Windowed pi: ${PI_DIR}" | tee -a "${RUN_LOG}"
echo "  Site pi    : ${PISITE_DIR}" | tee -a "${RUN_LOG}"
echo "  Log file   : ${RUN_LOG}" | tee -a "${RUN_LOG}"