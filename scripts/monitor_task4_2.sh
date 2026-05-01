#!/bin/bash
# ============================================================================
# title:         monitor_task4_2.sh
# project:       BIOL624 Population Genomics of Coccinella septempunctata
# author:        Reina Hastings (reinahastings13@gmail.com)
# date created:  2026-04-14
# last modified: 2026-04-14
# purpose:       Background progress monitor for the Task 4.2 gene_to_go_mapping
#                job. Emails a progress summary at a fixed interval by tailing
#                the job's stdout log and reporting elapsed time, recent
#                progress lines, and any warnings/errors. Scheduler-agnostic:
#                works with SLURM (anthill) or any scheduler that writes
#                stdout to a file.
# inputs:        --log       Path to the job stdout log (.out)
#                --to        Recipient email address
#                --jobid     Job ID (for subject line; SLURM_JOB_ID on anthill)
#                --interval  Seconds between updates (default 1500 = 25 min)
# outputs:       Email messages via /usr/sbin/sendmail.
# usage:         monitor_task4_2.sh --log logs/x.out --to you@example.com \
#                    --jobid 12345 --interval 1500 &
# ============================================================================
set -uo pipefail

# --- Argument Parsing ---
LOG=''
TO=''
JOBID=''
INTERVAL=1500

while [[ $# -gt 0 ]]; do
    case "$1" in
        --log)      LOG="$2";      shift 2 ;;
        --to)       TO="$2";       shift 2 ;;
        --jobid)    JOBID="$2";    shift 2 ;;
        --interval) INTERVAL="$2"; shift 2 ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

if [[ -z "$LOG" || -z "$TO" || -z "$JOBID" ]]; then
    echo "Usage: $0 --log PATH --to EMAIL --jobid JOBID [--interval SEC]"
    exit 1
fi

SENDMAIL=/usr/sbin/sendmail
HOSTNAME=$(hostname)
START_TIME=$(date +%s)

# --- Email Helpers ---
send_update() {
    local subject="$1"
    local body_file="$2"
    {
        echo "To: ${TO}"
        echo "From: task4.2-monitor@${HOSTNAME}"
        echo "Subject: ${subject}"
        echo "Content-Type: text/plain; charset=UTF-8"
        echo ''
        cat "${body_file}"
    } | "${SENDMAIL}" -t
}

# --- Summary Builder ---
# Pulls elapsed time, current pipeline step marker (e.g. [3/7]), last log
# lines, and any warnings/errors. The script's main() prints numbered
# [N/7] headers, so we surface the most recent one to show current stage.
build_summary() {
    local out="$1"
    local now_ts elapsed elapsed_human
    now_ts=$(date +%s)
    elapsed=$(( now_ts - START_TIME ))
    elapsed_human=$(printf '%dh %02dm' $(( elapsed / 3600 )) $(( (elapsed % 3600) / 60 )))

    {
        echo "Task 4.2 gene_to_go_mapping progress update"
        echo "Job ID:    ${JOBID}"
        echo "Host:      ${HOSTNAME}"
        echo "Elapsed:   ${elapsed_human}"
        echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S %Z')"
        echo ''
        echo '--- Current stage ---'
        if [[ -f "$LOG" ]]; then
            grep -E '^\[[0-9]+/7\]' "$LOG" 2>/dev/null | tail -1 || echo '(no stage marker yet)'
        else
            echo "(log file not yet created: $LOG)"
        fi
        echo ''
        echo '--- Last 30 log lines ---'
        if [[ -f "$LOG" ]]; then
            tail -30 "$LOG"
        fi
        echo ''
        echo '--- Warnings / errors ---'
        if [[ -f "$LOG" ]]; then
            grep -E 'ERROR|WARN|Traceback|Exception' "$LOG" 2>/dev/null | tail -10 \
                || echo '(none)'
        fi
    } > "$out"
}

# --- Main Loop ---
TMPFILE=$(mktemp)
trap "rm -f $TMPFILE" EXIT

while true; do
    sleep "$INTERVAL"
    build_summary "$TMPFILE"
    send_update "Task 4.2 ${JOBID} progress ($(date '+%H:%M'))" "$TMPFILE"
done