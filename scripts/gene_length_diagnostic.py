#!/usr/bin/env python3
# =============================================================================
# Title:        gene_length_diagnostic.py
# Project:      Population genomics of Coccinella septempunctata (BIOL624)
# Author:       Reina Hastings <reinahastings13@gmail.com>
# Date created: 2026-04-29
# Last modified: 2026-04-29
#
# Purpose:
#   Diagnostic to reconcile gene-count and total-bp discrepancy between the
#   strict feature=="gene" computation and the NCBI annotation report figure
#   (233,535,789 bp / 25,480 genes). Tabulates feature types in column 3,
#   pools all top-level (no-Parent) loci, computes both naive and merged
#   non-redundant totals, and breaks counts down by seqid to detect anything
#   on alt/unplaced contigs.
#
# Inputs:
#   --gff   Path to the NCBI RefSeq GFF (uncompressed or .gz).
#   --outdir  Directory for diagnostic TSVs.
#
# Outputs:
#   <outdir>/feature_type_table.tsv    one row per column-3 feature type
#   <outdir>/toplevel_feature_table.tsv  same, restricted to no-Parent rows
#   <outdir>/seqid_gene_counts.tsv     gene/top-level counts per seqid
#   <outdir>/merged_totals.tsv         naive vs merged-interval totals
#
# Usage example:
#   python scripts/gene_length_diagnostic.py \
#     --gff data/coccinella_septempunctata_/ncbi_dataset/data/GCF_907165205.1/genomic.gff \
#     --outdir results/annotation_qc
# =============================================================================

from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import IO


# --- Argument parsing --------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Diagnostic for GFF gene tabulation.")
    p.add_argument("--gff", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


# --- Helpers -----------------------------------------------------------------

def open_gff(path: Path) -> IO[str]:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_attributes(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in attr_field.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def merge_intervals(intervals: list[tuple[int, int]]) -> int:
    """Sum of merged (non-redundant) bp covered by a list of [start,end] intervals.

    Intervals are 1-based inclusive. Two intervals are merged when they share
    at least one base (s <= cur_end). Adjacent intervals (s == cur_end + 1)
    are NOT merged but contribute the same total either way.
    """
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1
    return total


# --- Main --------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Per-feature-type aggregates (all rows, including child features)
    ft_count: Counter[str] = Counter()
    ft_total_bp: dict[str, int] = defaultdict(int)

    # Top-level (no-Parent) aggregates
    tl_count: Counter[str] = Counter()
    tl_total_bp: dict[str, int] = defaultdict(int)
    tl_intervals_per_seq: dict[str, list[tuple[int, int]]] = defaultdict(list)

    # Strict gene aggregates (feature == "gene")
    gene_intervals_per_seq: dict[str, list[tuple[int, int]]] = defaultdict(list)
    gene_total_bp_naive = 0

    # Per-seqid counts
    seq_gene_count: Counter[str] = Counter()
    seq_toplevel_count: Counter[str] = Counter()

    with open_gff(args.gff) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid = fields[0]
            ftype = fields[2]
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            length = end - start + 1
            attrs = parse_attributes(fields[8])
            has_parent = "Parent" in attrs

            ft_count[ftype] += 1
            ft_total_bp[ftype] += length

            if not has_parent:
                tl_count[ftype] += 1
                tl_total_bp[ftype] += length
                tl_intervals_per_seq[seqid].append((start, end))
                seq_toplevel_count[seqid] += 1

            if ftype == "gene":
                gene_intervals_per_seq[seqid].append((start, end))
                gene_total_bp_naive += length
                seq_gene_count[seqid] += 1

    # --- TSV 1: all column-3 feature types --------------------------------
    with open(args.outdir / "feature_type_table.tsv", "w") as out:
        out.write("feature_type\tn_rows\ttotal_bp_naive\n")
        for ft, n in sorted(ft_count.items(), key=lambda kv: -kv[1]):
            out.write(f"{ft}\t{n}\t{ft_total_bp[ft]}\n")

    # --- TSV 2: top-level (no-Parent) feature types -----------------------
    with open(args.outdir / "toplevel_feature_table.tsv", "w") as out:
        out.write("feature_type\tn_rows\ttotal_bp_naive\n")
        for ft, n in sorted(tl_count.items(), key=lambda kv: -kv[1]):
            out.write(f"{ft}\t{n}\t{tl_total_bp[ft]}\n")

    # --- TSV 3: per-seqid counts ------------------------------------------
    all_seqids = set(seq_gene_count) | set(seq_toplevel_count)
    with open(args.outdir / "seqid_gene_counts.tsv", "w") as out:
        out.write("seqid\tplacement\tn_gene_features\tn_toplevel_features\n")
        for seqid in sorted(all_seqids):
            placement = "chromosome" if seqid.startswith("NC_") else "unplaced"
            out.write(
                f"{seqid}\t{placement}\t"
                f"{seq_gene_count.get(seqid, 0)}\t"
                f"{seq_toplevel_count.get(seqid, 0)}\n"
            )

    # --- Merged-interval totals -------------------------------------------
    gene_merged_bp = sum(merge_intervals(ivs) for ivs in gene_intervals_per_seq.values())
    tl_merged_bp = sum(merge_intervals(ivs) for ivs in tl_intervals_per_seq.values())
    tl_total_naive = sum(tl_total_bp.values())
    tl_total_count = sum(tl_count.values())

    ncbi_total_bp = 233_535_789
    ncbi_n_genes = 25_480

    with open(args.outdir / "merged_totals.tsv", "w") as out:
        out.write("subset\tn_features\ttotal_bp_naive\ttotal_bp_merged\tdelta_vs_ncbi_233535789\n")
        out.write(
            f"strict_gene\t{sum(seq_gene_count.values())}\t"
            f"{gene_total_bp_naive}\t{gene_merged_bp}\t"
            f"{gene_merged_bp - ncbi_total_bp}\n"
        )
        out.write(
            f"all_toplevel\t{tl_total_count}\t"
            f"{tl_total_naive}\t{tl_merged_bp}\t"
            f"{tl_merged_bp - ncbi_total_bp}\n"
        )
        out.write(f"ncbi_report\t{ncbi_n_genes}\t{ncbi_total_bp}\tNA\t0\n")

    # --- Console summary --------------------------------------------------
    print("=== Column-3 feature-type table ===")
    print(f"{'feature_type':<25} {'n_rows':>10} {'total_bp':>14}")
    for ft, n in sorted(ft_count.items(), key=lambda kv: -kv[1]):
        print(f"{ft:<25} {n:>10} {ft_total_bp[ft]:>14}")

    print("\n=== Top-level (no Parent=) feature-type table ===")
    print(f"{'feature_type':<25} {'n_rows':>10} {'total_bp':>14}")
    for ft, n in sorted(tl_count.items(), key=lambda kv: -kv[1]):
        print(f"{ft:<25} {n:>10} {tl_total_bp[ft]:>14}")

    print("\n=== Merged-interval reconciliation ===")
    print(f"strict feature=='gene':   n={sum(seq_gene_count.values()):>6}  "
          f"naive_bp={gene_total_bp_naive:>12}  merged_bp={gene_merged_bp:>12}")
    print(f"all top-level loci:       n={tl_total_count:>6}  "
          f"naive_bp={tl_total_naive:>12}  merged_bp={tl_merged_bp:>12}")
    print(f"NCBI annotation report:   n={ncbi_n_genes:>6}  "
          f"total_bp={ncbi_total_bp:>12}")

    print("\n=== Per-seqid distribution (top 25 by gene count) ===")
    print(f"{'seqid':<20} {'placement':<12} {'n_gene':>8} {'n_toplevel':>11}")
    for seqid, _ in seq_gene_count.most_common(25):
        placement = "chromosome" if seqid.startswith("NC_") else "unplaced"
        print(f"{seqid:<20} {placement:<12} "
              f"{seq_gene_count[seqid]:>8} {seq_toplevel_count.get(seqid, 0):>11}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
