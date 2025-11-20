#!/usr/bin/env bash

# title: fst.sh
# project: BIOL624 Final Project — Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2025-11-17
# last modified: 2025-11-20
#
# purpose:
#   Compute windowed Weir & Cockerham FST between two populations using vcftools.
#   Designed for continent-level comparisons (CHI, WEU, EEU, USA).
#
# inputs:
#   - Directory of per-chromosome VCFs (e.g., ../data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS/)
#   - Population file for group 1 (pop1)
#   - Population file for group 2 (pop2)
#   - Optional window size (default 5000 bp)
#
# outputs:
#   - <outdir>/<chrom>_FST.windowed.weir.fst
#   - <outdir>/logs/fst_<timestamp>.log
#
# usage example:
#   ./fst.sh \
#       -i ../data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS \
#       -p1 ../metadata/pop_CHI.txt \
#       -p2 ../metadata/pop_WEU.txt \
#       -w 5000 \
#       -o ../results/fst_CHI_vs_WEU_results

set -euo pipefail

# usage
usage() {
    cat <<EOF
Usage: fst.sh -i <vcf_dir> -p1 <pop1_file> -p2 <pop2_file> [-w window_bp] [-o outdir]

Required:
  -i   Directory containing per-chromosome VCFs
  -p1  Pop file for group 1
  -p2  Pop file for group 2

Optional:
  -w   Window size (default 5000 bp)
  -o   Output directory (auto-generated if not provided)

EOF
}

# parse arguments
WINDOW=5000
OUTDIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--vcf-dir) VCF_DIR="$2"; shift 2;;
        -p1|--pop1)   POP1="$2"; shift 2;;
        -p2|--pop2)   POP2="$2"; shift 2;;
        -w|--window)  WINDOW="$2"; shift 2;;
        -o|--outdir)  OUTDIR="$2"; shift 2;;
        -h|--help)    usage; exit 0;;
        *) echo "Unknown option: $1"; usage; exit 1;;
    esac
done

# validate required arguments
if [[ -z "${VCF_DIR:-}" || -z "${POP1:-}" || -z "${POP2:-}" ]]; then
    echo "ERROR: missing required arguments." >&2
    usage
    exit 1
fi

if [[ ! -d "$VCF_DIR" ]]; then
    echo "ERROR: VCF directory not found: $VCF_DIR" >&2
    exit 1
fi

if [[ ! -f "$POP1" || ! -f "$POP2" ]]; then
    echo "ERROR: One or both population files not found." >&2
    exit 1
fi

# auto-generate OUTDIR if not provided
if [[ -z "$OUTDIR" ]]; then
    P1_LABEL=$(basename "$POP1" .txt)
    P2_LABEL=$(basename "$POP2" .txt)
    OUTDIR="../results/fst_${P1_LABEL}_vs_${P2_LABEL}_results"
fi

mkdir -p "$OUTDIR/logs"

# logging
STAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="$OUTDIR/logs/fst_${STAMP}.log"

echo "Running fst.sh" | tee "$LOGFILE"
echo "  VCF_DIR : $VCF_DIR"   | tee -a "$LOGFILE"
echo "  POP1    : $POP1"      | tee -a "$LOGFILE"
echo "  POP2    : $POP2"      | tee -a "$LOGFILE"
echo "  OUTDIR  : $OUTDIR"    | tee -a "$LOGFILE"
echo "  WINDOW  : $WINDOW bp" | tee -a "$LOGFILE"

# run FST for each chromosome
shopt -s nullglob
VCF_FILES=("$VCF_DIR"/*.vcf.gz)
shopt -u nullglob

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .vcf.gz files found in $VCF_DIR" | tee -a "$LOGFILE"
    exit 1
fi

for VCF in "${VCF_FILES[@]}"; do
    CHR=$(basename "$VCF" .vcf.gz)

    echo "Processing $VCF (CHR: $CHR)" | tee -a "$LOGFILE"

    vcftools \
        --gzvcf "$VCF" \
        --weir-fst-pop "$POP1" \
        --weir-fst-pop "$POP2" \
        --fst-window-size "$WINDOW" \
        --fst-window-step "$WINDOW" \
        --out "$OUTDIR/${CHR}_FST" \
        2>&1 | tee -a "$LOGFILE"
done

echo "Done!" | tee -a "$LOGFILE"