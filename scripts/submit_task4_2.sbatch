#!/usr/bin/env bash
# ============================================================================
# title:         submit_task4_2.sbatch
# project:       BIOL624 Population Genomics of Coccinella septempunctata
# author:        Reina Hastings (reinahastings13@gmail.com)
# date created:  2026-04-14
# last modified: 2026-04-14 (added submit-time email, 30-min monitor interval,
#                TIME_LIMIT mail-type; dropped cpus-per-task 16->8 and
#                walltime 2h->3h to dispatch faster on a busy debug partition)
# purpose:       SLURM submission script for Task 4.2 (gene-to-GO mapping) on
#                the anthill cluster. Runs scripts/gene_to_go_mapping.py under
#                the ladybug conda env, streaming NCBI gene2go (Route A) and
#                performing the two-step BLAST orthology transfer
#                Csep -> Haxy -> Tcas (Route B). Launches a background monitor
#                that emails a progress summary every 30 minutes. Also sends
#                a submit-time confirmation email via sendmail (Slurm has no
#                SUBMIT mail-type).
# inputs:        data/go_annotations/gene2go (9.9 GB, NCBI)
#                data/go_annotations/tribolium_gene2go.tsv
#                data/{csep,haxy,tcas}/.../protein.faa + genomic.gff
#                results/enrichment/foreground_genes.tsv
#                results/enrichment/background_genes.tsv
# outputs:       results/enrichment/gene2go_foreground.tsv
#                results/enrichment/gene2go_background.tsv
#                results/enrichment/gene_go_cache.tsv (+ meta)
#                results/enrichment/go_mapping_summary.txt
# usage:         sbatch scripts/submit_task4_2.sbatch
# cluster notes: anthill SLURM quirk - --mem must be 1 (M) or omitted.
#                Only partition available is 'debug'. Compute nodes have
#                >=16 cores and >=128 GB RAM; requesting 16 CPUs is safe.
# ============================================================================

#SBATCH --job-name=task4_2_gene2go
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=1
#SBATCH --time=03:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=rhastings5624@sdsu.edu
#SBATCH --output=results/enrichment/logs/task4_2.%j.out
#SBATCH --error=results/enrichment/logs/task4_2.%j.err

# --- Environment Setup ---
# Clean SLURM env, then bring in conda and activate ladybug (which has
# blastp 2.17 and biopython 1.87).
module purge
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate ladybug

set -uo pipefail

# --- Locate Project Dir ---
# Assumes sbatch is run from the project root on anthill, e.g.
# /usr/scratchF/rhastings/BIOL624_Project_Fall_2025. SLURM sets
# SLURM_SUBMIT_DIR to whichever directory sbatch was invoked from.
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_DIR}"

mkdir -p results/enrichment/logs

# --- Submit-Time Confirmation Email ---
# Slurm has no SUBMIT mail-type; send an explicit notice at the start of the
# job allocation so the user gets a record that submission succeeded. BEGIN
# still fires separately when the job is dispatched.
JOBID="${SLURM_JOB_ID:-local}"
MAIL_TO="rhastings5624@sdsu.edu"
{
    printf 'To: %s\n' "${MAIL_TO}"
    printf 'Subject: [Task 4.2] Job %s submitted\n' "${JOBID}"
    printf 'Content-Type: text/plain; charset=UTF-8\n\n'
    printf 'Task 4.2 job %s has been accepted by Slurm on %s.\n' "${JOBID}" "$(hostname)"
    printf 'Submit time: %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    printf 'Partition:   %s\n' "${SLURM_JOB_PARTITION:-debug}"
    printf 'Progress emails every 30 min once running.\n'
} | /usr/sbin/sendmail -t || true

# --- Launch Background Monitor (30-min email updates) ---
LOG_PATH="${PROJECT_DIR}/results/enrichment/logs/task4_2.${JOBID}.out"

bash scripts/monitor_task4_2.sh \
    --log "${LOG_PATH}" \
    --to "${MAIL_TO}" \
    --jobid "${JOBID}" \
    --interval 1800 &
MONITOR_PID=$!

cleanup() {
    if kill -0 "${MONITOR_PID}" 2>/dev/null; then
        kill "${MONITOR_PID}" 2>/dev/null || true
    fi
}
trap cleanup EXIT

# --- Run Task 4.2 ---
echo "=== Task 4.2 start: $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname) ==="
echo "Job ID:  ${JOBID}"
echo "Threads: ${SLURM_CPUS_PER_TASK:-unset}"
echo ''

python scripts/gene_to_go_mapping.py \
    --enrichment-dir results/enrichment \
    --data-dir data \
    --threads "${SLURM_CPUS_PER_TASK:-8}"

RC=$?
echo ''
echo "=== Task 4.2 end: $(date '+%Y-%m-%d %H:%M:%S %Z') (rc=${RC}) ==="
exit "${RC}"
