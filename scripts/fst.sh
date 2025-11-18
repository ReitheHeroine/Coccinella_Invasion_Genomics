#!/usr/bin/env bash

# title: 'fst.sh'
# project: 'BIOL624 Final Project: Selection Detection in Lady Beetles'
# author: 'Reina Hastings'
# purpose: 'Compute windowed Weir & Cockerham FST between two populations.'
# usage:
#   ./fst.sh -i ../data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS \
#            -p1 ../metadata/pop_Europe.txt \
#            -p2 ../metadata/pop_Asia.txt \
#            -w 5000 \
#            -o ../results/fst_Eur_vs_Asia

set -euo pipefail

# defaults
WINDOW=5000

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--vcf-dir)
            VCF_DIR="$2"; shift 2;;
        -p1|--pop1)
            POP1="$2"; shift 2;;
        -p2|--pop2)
            POP2="$2"; shift 2;;
        -w|--window)
            WINDOW="$2"; shift 2;;
        -o|--outdir)
            OUTDIR="$2"; shift 2;;
        *)
            echo "Unknown option: $1"; exit 1;;
    esac
done

mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/logs"

LOGFILE="$OUTDIR/logs/fst_$(date +%Y%m%d_%H%M%S).log"

echo "FST analysis" | tee "$LOGFILE"
echo "VCF_DIR: $VCF_DIR" | tee -a "$LOGFILE"
echo "POP1:    $POP1"     | tee -a "$LOGFILE"
echo "POP2:    $POP2"     | tee -a "$LOGFILE"
echo "WINDOW:  $WINDOW"   | tee -a "$LOGFILE"
echo "OUTDIR:  $OUTDIR"   | tee -a "$LOGFILE"

# make sure VCF_DIR exists
if [[ ! -d "$VCF_DIR" ]]; then
    echo "Error: VCF_DIR '$VCF_DIR' does not exist." | tee -a "$LOGFILE"
    exit 1
fi

# expand VCF list
shopt -s nullglob
VCF_FILES=("$VCF_DIR"/*.vcf.gz)
shopt -u nullglob

if [[ ${#VCF_FILES[@]} -eq 0 ]]; then
    echo "Error: no .vcf.gz files found in '$VCF_DIR'." | tee -a "$LOGFILE"
    exit 1
fi

for VCF in "${VCF_FILES[@]}"; do
    CHR=$(basename "$VCF" .filtered.vcf.gz)
    echo "Processing $VCF ($CHR)" | tee -a "$LOGFILE"

    vcftools \
        --gzvcf "$VCF" \
        --weir-fst-pop "$POP1" \
        --weir-fst-pop "$POP2" \
        --fst-window-size "$WINDOW" \
        --fst-window-step "$WINDOW" \
        --out "$OUTDIR/${CHR}_FST" \
        2>&1 | tee -a "$LOGFILE"
done

echo "Done." | tee -a "$LOGFILE"