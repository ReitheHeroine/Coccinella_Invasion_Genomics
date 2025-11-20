#!/usr/bin/env bash

# title: diversity.sh
# project: BIOL624 Final Project — Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-11-17
# last modified: 2025-11-20
#
# purpose:
#   Compute Tajima’s D, windowed nucleotide diversity (π), and site-level π
#   from per-chromosome VCFs using vcftools. Supports global (all samples) or
#   per-population analyses using pop files generated from create_metadata.py
#   (e.g., pop_CHI.txt, pop_WEU.txt, pop_EEU.txt, pop_USA.txt).
#
# inputs:
#   - Directory containing per-chromosome VCFs (*.vcf.gz)
#   - Optional population file (one sample ID per line)
#   - Optional window size (default: 5000 bp)
#
# outputs:
#   <outdir>/tajima/<chrom>[_<poplabel>].Tajima.D
#   <outdir>/pi/<chrom>[_<poplabel>].windowed.pi
#   <outdir>/pi_sites/<chrom>[_<poplabel>].sites.pi
#   <outdir>/logs/diversity_<timestamp>.log
#
# usage example:
#   # Global (whole population)
#   ./diversity.sh \
#       -i ../data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS \
#       -o ../results/wholepop_results \
#       -w 5000
#
#   # CHI-only diversity
#   ./diversity.sh \
#       -i ../data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS \
#       -p ../metadata/pop_CHI.txt \
#       -o ../results/CHI_results \
#       -w 5000

set -euo pipefail

# usage() message
usage() {
cat <<EOF
Usage: diversity.sh -i <vcf_dir> [-w window_bp] [-o outdir] [-p popfile] [-g vcf_glob]

Required:
  -i   Directory containing per-chromosome VCFs (*.vcf.gz)

Optional:
  -w   Window size in bp for Tajima's D & windowed pi (default: 5000)
  -o   Output directory (default: results/)
  -p   Population file (pop_CHI.txt, pop_WEU.txt, etc.)
  -g   Glob pattern for VCFs (default: "*.vcf.gz")
  -h   Show this help message

Examples:
  # Global diversity
  ./diversity.sh -i ../data/FILTERED_CHR_VCFS -o ../results/wholepop_results

  # Per-continent diversity
  ./diversity.sh -i ../data/FILTERED_CHR_VCFS -p ../metadata/pop_WEU.txt -o ../results/WEU_results
EOF
}

# argument parsing
VCF_DIR=""
WIN=5000
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
        \?) echo "ERROR: Unknown option -${OPTARG}" >&2; usage; exit 1 ;;
        :)  echo "ERROR: Option -${OPTARG} requires an argument." >&2; usage; exit 1 ;;
    esac
done

# sanity checks
if [[ -z "${VCF_DIR}" ]]; then
    echo "ERROR: -i <vcf_dir> is required." >&2
    usage
    exit 1
fi

if [[ ! -d "${VCF_DIR}" ]]; then
    echo "ERROR: VCF directory not found: ${VCF_DIR}" >&2
    exit 1
fi

if ! command -v vcftools >/dev/null 2>&1; then
    echo "ERROR: vcftools not found in PATH." >&2
    exit 1
fi

# setup output directories
TAJ_DIR="${OUTDIR}/tajima"
PI_DIR="${OUTDIR}/pi"
PISITE_DIR="${OUTDIR}/pi_sites"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "${TAJ_DIR}" "${PI_DIR}" "${PISITE_DIR}" "${LOG_DIR}"

STAMP=$(date +"%Y%m%d_%H%M%S")
RUN_LOG="${LOG_DIR}/diversity_${STAMP}.log"

# intro log
echo "Running diversity.sh" | tee -a "${RUN_LOG}"
echo "  VCF_DIR:   ${VCF_DIR}" | tee -a "${RUN_LOG}"
echo "  WINDOW:    ${WIN}" | tee -a "${RUN_LOG}"
echo "  OUTDIR:    ${OUTDIR}" | tee -a "${RUN_LOG}"
echo "  POP_FILE:  ${POP_FILE:-<none>}" | tee -a "${RUN_LOG}"
echo "  VCF_GLOB:  ${VCF_GLOB}" | tee -a "${RUN_LOG}"

# population label handling
SUFFIX=""

if [[ -n "${POP_FILE}" ]]; then
    if [[ ! -f "${POP_FILE}" ]]; then
        echo "ERROR: Population file not found: ${POP_FILE}" | tee -a "${RUN_LOG}" >&2
        exit 1
    fi

    # example: pop_CHI.txt → CHI
    POP_LABEL=$(basename "${POP_FILE}")
    POP_LABEL="${POP_LABEL%.*}"           # drop .txt
    POP_LABEL="${POP_LABEL#pop_}"         # drop pop_
    SUFFIX="_${POP_LABEL}"

    echo "Per-population mode: ${POP_LABEL}" | tee -a "${RUN_LOG}"
else
    echo "Global mode: whole population" | tee -a "${RUN_LOG}"
    POP_LABEL="wholepop"
fi

# main loop over chromosomes
shopt -s nullglob
VCF_LIST=("${VCF_DIR}"/${VCF_GLOB})

if [[ ${#VCF_LIST[@]} -eq 0 ]]; then
    echo "ERROR: No VCFs matching pattern '${VCF_GLOB}' found." | tee -a "${RUN_LOG}"
    exit 1
fi

for VCF in "${VCF_LIST[@]}"; do
    BASENAME=$(basename "${VCF}")
    CHR="${BASENAME%.vcf.gz}"
    CHR="${CHR%.vcf}"

    echo "Processing ${CHR} ..." | tee -a "${RUN_LOG}"

    # TajimaD
    vcftools \
        --gzvcf "${VCF}" \
        ${POP_FILE:+--keep "${POP_FILE}"} \
        --TajimaD "${WIN}" \
        --out "${TAJ_DIR}/${CHR}${SUFFIX}" \
        >> "${RUN_LOG}" 2>&1

    # windowed pi
    vcftools \
        --gzvcf "${VCF}" \
        ${POP_FILE:+--keep "${POP_FILE}"} \
        --window-pi "${WIN}" \
        --out "${PI_DIR}/${CHR}${SUFFIX}" \
        >> "${RUN_LOG}" 2>&1

    # site-pi
    vcftools \
        --gzvcf "${VCF}" \
        ${POP_FILE:+--keep "${POP_FILE}"} \
        --site-pi \
        --out "${PISITE_DIR}/${CHR}${SUFFIX}" \
        >> "${RUN_LOG}" 2>&1

    echo "  Finished ${CHR}${SUFFIX}" | tee -a "${RUN_LOG}"
done

# completion message
echo "All done. Output written to:" | tee -a "${RUN_LOG}"
echo "  TajimaD:   ${TAJ_DIR}" | tee -a "${RUN_LOG}"
echo "  window-pi: ${PI_DIR}" | tee -a "${RUN_LOG}"
echo "  site-pi:   ${PISITE_DIR}" | tee -a "${RUN_LOG}"
echo "  Log file:  ${RUN_LOG}" | tee -a "${RUN_LOG}"