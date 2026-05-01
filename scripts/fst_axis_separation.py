#!/usr/bin/env python3
################################################################################
# title: fst_axis_separation.py
# project: BIOL624 Final Project -- Population genomics of Coccinella septempunctata
# author: Reina Hastings (reinahastings13@gmail.com)
# date created: 2026-04-20
# last modified: 2026-04-23 (Task 4.7 Step 7e: pure-SNP-only convention for
#                            max_fst_among_invasive and max_fst_chi; mixed-axis
#                            SNPs contribute NaN to both primary axes. True
#                            per-axis FSTs for mixed SNPs are looked up from
#                            source .weir.fst files and reported in two new
#                            aux columns (max_fst_among_inv_mixed_lookup,
#                            max_fst_chi_mixed_lookup). Fixes the systematic
#                            inflation identified in Task 4.6 (e.g., LOC123319966
#                            0.874 -> 0.320, LOC123318563 0.930 -> 0.314).)
#
# purpose:
#   Task 4.4: Separate per-SNP FST values by selection axis (among-invasive vs
#   CHI-vs-invasive) and aggregate to per-gene maxima across foreground genes.
#   The master candidate table stores a single percentile_fst_max per SNP that
#   mixes axes; reporting the wrong axis inflates effect sizes for POST_INVASION
#   genes (among-invasive genome-wide mean ~0.003-0.017 vs CHI ~0.06). This
#   script classifies each 2+ method SNP by axis using percentile_comparisons
#   and joins to the canonical SNP-to-gene mapping produced by Task 4.1
#   (snp_to_gene_mapping.tsv) to emit per-gene max FST on each axis. Using the
#   Task 4.1 mapping (rather than recomputing flanking) guarantees strand-aware
#   5 kb up / 2 kb down windows and keeps the candidate table coherent with the
#   foreground_genes.tsv definition.
#
#   Mixed-axis convention (per Task 4.7 Step 7e):
#     - max_fst_among_invasive / max_fst_chi use PURE axis SNPs only. A mixed
#       SNP (outlier in both among-inv and CHI comparisons) contributes NaN to
#       both primary maxima, because percentile_fst_max for that SNP is the
#       max across all its outlier comparisons and cannot be attributed to a
#       single axis honestly.
#     - For transparency, mixed-axis SNP positions are looked up in the source
#       .weir.fst files to get true per-axis FSTs. Per-gene max over these
#       looked-up values is reported in the aux columns
#       max_fst_among_inv_mixed_lookup / max_fst_chi_mixed_lookup. These are
#       informational (for manuscript narrative) -- they are NOT folded into
#       the primary max columns.
#
# inputs:
#   --master    results/cross_method_concordance/master_candidate_table.tsv
#               (filter n_methods >= 2 for the foreground SNP set)
#   --foreground results/enrichment/foreground_genes.tsv
#   --snp-to-gene results/enrichment/snp_to_gene_mapping.tsv
#               (canonical Task 4.1 SNP-to-gene mapping; strand-aware flanking)
#   --fst-dir   results/fst
#               (contains results/fst/<comp>/fst_genomewide/<chrom>.filtered.weir.fst
#               for the 6 pairwise comparisons; used for mixed-axis SNP lookup)
#
# outputs (all in results/enrichment/):
#   candidate_gene_fst_by_axis.tsv  One row per foreground gene.
#                                   Columns: gene_id, chrom, gene_start,
#                                   gene_end, consensus_pattern, n_outlier_snps,
#                                   n_2plus_method_snps, max_n_methods,
#                                   n_snps_among_inv, n_snps_chi, n_snps_mixed,
#                                   max_fst_among_invasive, max_fst_chi,
#                                   max_fst_among_inv_mixed_lookup,
#                                   max_fst_chi_mixed_lookup,
#                                   detecting_methods, location_type
#   fst_axis_separation_summary.txt Text log of axis classification counts and
#                                   per-gene coverage statistics.
#
# usage:
#   python scripts/fst_axis_separation.py \
#     --master results/cross_method_concordance/master_candidate_table.tsv \
#     --foreground results/enrichment/foreground_genes.tsv \
#     --snp-to-gene results/enrichment/snp_to_gene_mapping.tsv \
#     --fst-dir results/fst \
#     --out-dir results/enrichment
################################################################################

import argparse
import sys
from pathlib import Path

import pandas as pd

# --- Axis definitions ---
AMONG_INVASIVE_COMPS = {"EEU_vs_WEU", "EEU_vs_USA", "USA_vs_WEU"}
CHI_COMPS = {"CHI_vs_EEU", "CHI_vs_WEU", "CHI_vs_USA"}
ALL_COMPS = AMONG_INVASIVE_COMPS | CHI_COMPS


def lookup_mixed_axis_fst(
    mixed_positions: list[tuple[str, int]],
    fst_dir: Path,
) -> dict[tuple[str, int], dict[str, float]]:
    """Look up per-comparison FST for a small set of (chrom, pos) positions.

    For each chromosome represented in `mixed_positions`, load the
    corresponding .weir.fst file for each of the 6 comparisons and pull out
    the rows matching the target positions. Returns:
        {(chrom, pos): {comparison_name: fst_value, ...}}

    Missing positions or comparisons are simply absent from the inner dict
    (caller falls back to NaN). Only the chromosomes we need are read; this
    is the targeted-lookup approach specified in the Task 4.7 handoff
    (avoids scanning the full 800K-row files across all 6 comparisons).
    """
    if not mixed_positions:
        return {}

    # Organize target positions by chromosome so we read each chrom file once
    targets_by_chrom: dict[str, set[int]] = {}
    for chrom, pos in mixed_positions:
        targets_by_chrom.setdefault(chrom, set()).add(int(pos))

    out: dict[tuple[str, int], dict[str, float]] = {
        (c, p): {} for c, p in mixed_positions
    }

    for comp in sorted(ALL_COMPS):
        for chrom, target_posns in targets_by_chrom.items():
            f = fst_dir / comp / "fst_genomewide" / f"{chrom}.filtered.weir.fst"
            if not f.exists():
                print(f"[warn] missing .weir.fst: {f}", file=sys.stderr)
                continue
            df = pd.read_csv(f, sep="\t")
            hits = df[df["POS"].isin(target_posns)]
            for _, row in hits.iterrows():
                key = (chrom, int(row["POS"]))
                fst = row["WEIR_AND_COCKERHAM_FST"]
                try:
                    out[key][comp] = float(fst)
                except (TypeError, ValueError):
                    # vcftools sometimes writes '-nan' or similar
                    out[key][comp] = float("nan")

    return out


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description="Task 4.4: Separate per-SNP FST by axis and aggregate to genes.",
    )
    p.add_argument("--master", required=True, type=Path,
                   help="master_candidate_table.tsv (cross-method concordance output)")
    p.add_argument("--foreground", required=True, type=Path,
                   help="foreground_genes.tsv (Task 4.1 output)")
    p.add_argument("--snp-to-gene", required=True, type=Path,
                   help="snp_to_gene_mapping.tsv (canonical Task 4.1 mapping)")
    p.add_argument("--fst-dir", required=True, type=Path,
                   help="results/fst (contains <comp>/fst_genomewide/<chrom>.filtered.weir.fst)")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Directory for candidate_gene_fst_by_axis.tsv and log")
    return p.parse_args()


def classify_axis(comparisons: str) -> str:
    """Classify a SNP's detecting comparisons into an axis category.

    Returns one of: among_invasive, chi_vs_invasive, mixed, no_percentile.
    """
    if not isinstance(comparisons, str) or comparisons.strip() in ("", "NA"):
        return "no_percentile"
    comps = {c.strip() for c in comparisons.split(";") if c.strip()}
    has_among = bool(comps & AMONG_INVASIVE_COMPS)
    has_chi = bool(comps & CHI_COMPS)
    if has_among and has_chi:
        return "mixed"
    if has_among:
        return "among_invasive"
    if has_chi:
        return "chi_vs_invasive"
    return "no_percentile"


def assign_axis_fst(row: pd.Series) -> tuple[float, float]:
    """Return (fst_among_inv, fst_chi) for a SNP given its axis classification.

    Pure-SNP-only convention (Task 4.7 Step 7e): mixed-axis SNPs contribute
    NaN to BOTH primary axes. Their true per-axis FSTs are looked up from
    .weir.fst files and reported in aux columns downstream, not folded into
    these primary maxima.
    """
    fst = row["percentile_fst_max"]
    axis = row["axis"]
    if pd.isna(fst):
        return (float("nan"), float("nan"))
    if axis == "among_invasive":
        return (fst, float("nan"))
    if axis == "chi_vs_invasive":
        return (float("nan"), fst)
    # mixed or no_percentile: contribute nothing to either primary axis
    return (float("nan"), float("nan"))


def derive_state_and_block(row: pd.Series) -> tuple[str | None, str | None]:
    """Classify a gene by its axis state (A/B/C/D) and assign selection_block.

    See handoff §3.1 Task 4.7 Step 7f for state and block definitions:
        A -> AMONG_INVASIVE      (pure among-inv SNPs only)
        B -> NATIVE_VS_INVASIVE  (pure CHI SNPs only)
        C -> LAYERED_SELECTION   (pure SNPs on both axes)
        D -> POPULATION_SPECIFIC (only mixed-axis SNPs)
    A gene is assigned to a block iff location_type == "genic" (the single
    uniform filter; threshold-based filtering replaced with magnitude-based
    ranking via max_FST).
    """
    has_ami     = pd.notna(row["max_fst_among_invasive"])
    has_chi     = pd.notna(row["max_fst_chi"])
    has_ami_mix = pd.notna(row["max_fst_among_inv_mixed_lookup"])
    has_chi_mix = pd.notna(row["max_fst_chi_mixed_lookup"])

    if has_ami and not has_chi:
        state = "A"
    elif has_chi and not has_ami:
        state = "B"
    elif has_ami and has_chi:
        state = "C"
    elif has_ami_mix or has_chi_mix:
        state = "D"
    else:
        state = None

    is_genic = row.get("location_type") == "genic"
    if state is None or not is_genic:
        block = None
    else:
        block = {
            "A": "AMONG_INVASIVE",
            "B": "NATIVE_VS_INVASIVE",
            "C": "LAYERED_SELECTION",
            "D": "POPULATION_SPECIFIC",
        }[state]
    return state, block


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load inputs ---
    print(f"[load] master: {args.master}", file=sys.stderr)
    master = pd.read_csv(args.master, sep="\t", low_memory=False)
    print(f"[load] foreground: {args.foreground}", file=sys.stderr)
    fg = pd.read_csv(args.foreground, sep="\t")
    print(f"[load] snp_to_gene: {args.snp_to_gene}", file=sys.stderr)
    s2g = pd.read_csv(args.snp_to_gene, sep="\t")

    n_master_total = len(master)
    master = master[master["n_methods"] >= 2].copy()
    n_2plus = len(master)
    print(f"[filter] master: {n_master_total} -> {n_2plus} SNPs (n_methods >= 2)",
          file=sys.stderr)

    # --- Step 1: classify each 2+ method SNP by axis ---
    master["axis"] = master["percentile_comparisons"].apply(classify_axis)
    axis_counts = master["axis"].value_counts().to_dict()
    print(f"[axis] classification counts: {axis_counts}", file=sys.stderr)

    # --- Step 2: compute per-SNP axis-specific FST (pure-SNP-only) ---
    fst_pairs = master.apply(assign_axis_fst, axis=1, result_type="expand")
    fst_pairs.columns = ["fst_among_inv", "fst_chi"]
    master = pd.concat([master, fst_pairs], axis=1)

    n_mixed = int((master["axis"] == "mixed").sum())
    if n_mixed:
        print(f"[info] mixed-axis SNPs: n={n_mixed} (contribute NaN to primary "
              f"axes; true per-axis FSTs looked up below for aux columns)",
              file=sys.stderr)

    # --- Step 2b: lookup true per-axis FST for mixed SNPs from .weir.fst ---
    mixed = master[master["axis"] == "mixed"]
    mixed_positions = list(zip(mixed["chrom"], mixed["pos"]))
    print(f"[lookup] reading .weir.fst for {len(mixed_positions)} mixed-axis "
          f"SNPs on {len(set(mixed['chrom']))} chromosome(s)...",
          file=sys.stderr)
    mixed_lookup = lookup_mixed_axis_fst(mixed_positions, args.fst_dir)

    # Derive per-SNP lookup columns: max across among-inv comps, max across CHI comps.
    def _max_over(comps: set[str], by_comp: dict[str, float]) -> float:
        vals = [by_comp[c] for c in comps if c in by_comp and not pd.isna(by_comp[c])]
        return max(vals) if vals else float("nan")

    mixed_lookup_rows = []
    for (chrom, pos), by_comp in mixed_lookup.items():
        mixed_lookup_rows.append({
            "snp_id": f"{chrom}:{pos}",
            "fst_among_inv_mixed_lookup": _max_over(AMONG_INVASIVE_COMPS, by_comp),
            "fst_chi_mixed_lookup":       _max_over(CHI_COMPS, by_comp),
        })
    mixed_lookup_df = pd.DataFrame(
        mixed_lookup_rows,
        columns=["snp_id", "fst_among_inv_mixed_lookup", "fst_chi_mixed_lookup"],
    )
    master = master.merge(mixed_lookup_df, on="snp_id", how="left")

    # Sanity report: confirm the known cases from Task 4.6 QC report
    for sanity_snp in ("NC_058197.1:17290234", "NC_058196.1:13830930"):
        row = master[master["snp_id"] == sanity_snp]
        if len(row) == 1 and row["axis"].iat[0] == "mixed":
            print(f"[sanity] {sanity_snp}: among_inv_true="
                  f"{row['fst_among_inv_mixed_lookup'].iat[0]:.4f}, "
                  f"chi_true={row['fst_chi_mixed_lookup'].iat[0]:.4f}",
                  file=sys.stderr)

    # --- Step 3: per-gene aggregation via canonical snp_to_gene mapping ---
    # snp_to_gene_mapping.tsv carries one row per (SNP, gene) pair using Task 4.1's
    # strand-aware 5 kb upstream / 2 kb downstream flanking. Joining on snp_id then
    # grouping by gene_id is equivalent to recomputing the window join, without risk
    # of divergence from the foreground_genes.tsv definition.
    snp_cols = ["snp_id", "n_methods", "axis", "fst_among_inv", "fst_chi",
                "fst_among_inv_mixed_lookup", "fst_chi_mixed_lookup"]
    joined = s2g.merge(master[snp_cols], on="snp_id", how="inner")
    expected_pairs = len(s2g)
    print(f"[join] snp_to_gene pairs: {expected_pairs}; joined to master: {len(joined)}",
          file=sys.stderr)
    if len(joined) != expected_pairs:
        print(f"[warn] {expected_pairs - len(joined)} SNP-gene pairs failed to join to "
              f"master (unexpected: snp_to_gene_mapping is restricted to 2+ method SNPs)",
              file=sys.stderr)

    def agg_gene(g: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "n_2plus_method_snps": g["snp_id"].nunique(),
            "max_n_methods": int(g["n_methods"].max()),
            "n_snps_among_inv": int((g["axis"] == "among_invasive").sum()),
            "n_snps_chi": int((g["axis"] == "chi_vs_invasive").sum()),
            "n_snps_mixed": int((g["axis"] == "mixed").sum()),
            "max_fst_among_invasive": g["fst_among_inv"].max(skipna=True),
            "max_fst_chi": g["fst_chi"].max(skipna=True),
            "max_fst_among_inv_mixed_lookup": g["fst_among_inv_mixed_lookup"].max(skipna=True),
            "max_fst_chi_mixed_lookup":       g["fst_chi_mixed_lookup"].max(skipna=True),
        })

    gene_agg = (joined.groupby("gene_id", sort=False)
                     .apply(agg_gene, include_groups=False)
                     .reset_index())

    # Left-join onto foreground so the output row set is exactly the 2,648 genes
    # in foreground_genes.tsv, in its order. Every foreground gene should match by
    # construction (it was defined as a gene with >=1 two-method SNP).
    out_df = fg.merge(gene_agg, on="gene_id", how="left")
    n_missing = int(out_df["n_2plus_method_snps"].isna().sum())
    if n_missing:
        print(f"[warn] {n_missing} foreground genes have no matching 2+ method SNP "
              f"in snp_to_gene_mapping (unexpected)", file=sys.stderr)
    # Zero-fill integer count columns for any unmatched genes
    count_cols = ["n_2plus_method_snps", "n_snps_among_inv", "n_snps_chi", "n_snps_mixed"]
    for c in count_cols:
        out_df[c] = out_df[c].fillna(0).astype(int)

    # --- Derive gene_state and selection_block (single source of truth) ---
    sb = out_df.apply(derive_state_and_block, axis=1, result_type="expand")
    sb.columns = ["gene_state", "selection_block"]
    out_df = pd.concat([out_df, sb], axis=1)

    # Reorder to spec column order (aux mixed-lookup columns placed right
    # after the primary max_fst columns so they're easy to diff against;
    # gene_state and selection_block placed near the top so downstream
    # consumers can see the classification at a glance).
    out_df = out_df[[
        "gene_id", "chrom", "gene_start", "gene_end",
        "gene_state", "selection_block",
        "consensus_pattern",
        "n_outlier_snps", "n_2plus_method_snps", "max_n_methods",
        "n_snps_among_inv", "n_snps_chi", "n_snps_mixed",
        "max_fst_among_invasive", "max_fst_chi",
        "max_fst_among_inv_mixed_lookup", "max_fst_chi_mixed_lookup",
        "detecting_methods", "location_type",
    ]]

    # --- Step 4: write outputs ---
    out_tsv = args.out_dir / "candidate_gene_fst_by_axis.tsv"
    out_df.to_csv(out_tsv, sep="\t", index=False, na_rep="NA")
    print(f"[write] {out_tsv} ({len(out_df)} rows)", file=sys.stderr)

    # Summary statistics
    n_genes = len(out_df)
    n_with_among = int(out_df["max_fst_among_invasive"].notna().sum())
    n_with_chi = int(out_df["max_fst_chi"].notna().sum())
    n_with_both = int((out_df["max_fst_among_invasive"].notna()
                       & out_df["max_fst_chi"].notna()).sum())
    n_with_any_snp = int((out_df["n_2plus_method_snps"] > 0).sum())

    n_with_mixed_lookup = int(
        out_df["max_fst_among_inv_mixed_lookup"].notna().sum()
        | out_df["max_fst_chi_mixed_lookup"].notna().sum()
    )

    summary = [
        "# Task 4.4: FST axis separation summary",
        "# (Task 4.7 Step 7e: pure-SNP-only convention for primary max FSTs;",
        "#  mixed-axis SNPs' true per-axis FSTs reported in aux columns.)",
        f"input master rows (all): {n_master_total}",
        f"input master rows (n_methods >= 2): {n_2plus}",
        f"input foreground genes: {n_genes}",
        "",
        "axis classification of 2+ method SNPs:",
        f"  among_invasive: {axis_counts.get('among_invasive', 0)}",
        f"  chi_vs_invasive: {axis_counts.get('chi_vs_invasive', 0)}",
        f"  mixed:           {axis_counts.get('mixed', 0)}",
        f"  no_percentile:   {axis_counts.get('no_percentile', 0)}",
        "",
        "per-gene coverage (primary, pure-SNP-only):",
        f"  genes with any 2+ method SNP in flank: {n_with_any_snp}/{n_genes}",
        f"  genes with max_fst_among_invasive:     {n_with_among}/{n_genes}",
        f"  genes with max_fst_chi:                {n_with_chi}/{n_genes}",
        f"  genes with both axes populated:        {n_with_both}/{n_genes}",
        f"  genes with any mixed-lookup value:     {n_with_mixed_lookup}/{n_genes}",
        "",
        "# Acceptance check: POST_INVASION genic genes with n_outlier_snps >= 5",
    ]

    post_check = out_df[(out_df["consensus_pattern"] == "POST_INVASION")
                       & (out_df["location_type"] == "genic")
                       & (out_df["n_outlier_snps"] >= 5)]
    n_post = len(post_check)
    n_post_with_among = int(post_check["max_fst_among_invasive"].notna().sum())
    summary.append(f"  POST_INVASION genic n>=5 total: {n_post}")
    summary.append(f"  of those with max_fst_among_invasive populated: "
                   f"{n_post_with_among}/{n_post}")

    summary_path = args.out_dir / "fst_axis_separation_summary.txt"
    summary_path.write_text("\n".join(summary) + "\n")
    print(f"[write] {summary_path}", file=sys.stderr)
    print("\n".join(summary), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
