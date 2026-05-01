#!/usr/bin/env Rscript

# title: build_candidate_tables_and_figures.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-04-20
# last modified: 2026-04-29 (Threshold cleanup pass: dropped both
#                            horizontal top-1% F_ST threshold lines from
#                            Figure 1. The native-axis threshold (~0.93)
#                            is a sample-size artifact of CHI n=5
#                            saturating F_ST near 1, so rendering it as a
#                            visual cutoff was misleading. The bottom-
#                            panel threshold (~0.19) is biologically
#                            meaningful but dropped for visual symmetry.
#                            The figure-specific layered rule still uses
#                            both thresholds as its mechanical criterion;
#                            the rule is now described in prose in the
#                            caption rather than visually anchored.)
#                2026-04-29 (Figure 1 rebuilt as a mirrored two-panel
#                            Manhattan: top panel = per-SNP max F_ST across
#                            the 3 CHI-vs-X comparisons (CHI-filtered),
#                            bottom panel = per-SNP max F_ST across the 3
#                            invasive-vs-invasive comparisons, y-axis
#                            reversed so peaks grow downward and chromosome
#                            labels sit in the gap between panels. Per-SNP
#                            inputs are produced by
#                            scripts/manhattan_mirrored_dataprep.R, which
#                            must be run before --step 5c. Figure-specific
#                            "effective LAYERED" rule applied for visual
#                            logic only: a featured gene gets the dashed
#                            vertical connector + dual-panel label iff its
#                            representative SNPs cross the panel-specific
#                            top-1% F_ST threshold on BOTH panels. Featured
#                            points colored by functional_theme (5-level
#                            collapse of functional_category, added to
#                            featured_candidate_table.tsv by the dataprep
#                            script). Caption written to
#                            manhattan_candidates_caption.txt as a draft.
#                            Plotdata TSV is now long-format (one row per
#                            featured gene per panel).)
#                2026-04-28 (Manhattan refresh: removed NC_058196.1 region
#                            rect, added gold halos behind each featured
#                            candidate so outliers are distinguishable from
#                            the background scatter, added dotted vertical
#                            chromosome separators, and dropped the inline
#                            geom_text_repel labels in favor of relying on
#                            the functional-category fill legend.)
#                2026-04-24 (Task 4.7 Step 7f rewrite: 4-block axis-based
#                            categorization replacing the old 4-block
#                            per-SNP consensus_pattern scheme; new blocks
#                            are AMONG_INVASIVE, NATIVE_VS_INVASIVE,
#                            LAYERED_SELECTION, POPULATION_SPECIFIC;
#                            NC_058196.1 cluster rule dropped; LOC123319966
#                            demoted; LOC123318563 FST corrected to 0.314
#                            via Step 7e pure-SNP-only convention; axis-
#                            category color palette added.)
#
# purpose:
#   Task 4.5 of the ladybug popgen pipeline. Builds the full supplementary
#   candidate gene table, the manuscript featured candidate table, and the
#   two accompanying figures (annotated Manhattan plot, functional-category
#   dot plot). All reported effect sizes draw on the axis-separated FST
#   from Task 4.4 (Step 7e pure-SNP-only convention) so that AMONG_INVASIVE
#   candidates are benchmarked against the among-invasive genome-wide
#   background (~0.010) rather than the inflated CHI-vs-invasive axis
#   (~0.06), and mixed-axis SNPs do not inflate primary effect sizes.
#
#   Candidate selection blocks (Task 4.7 Step 7f; see handoff Section 3.1
#   and resolved question 2026-04-24):
#     AMONG_INVASIVE      State A: gene has pure among-invasive SNPs only
#                         (max_fst_among_invasive populated, max_fst_chi NA)
#     NATIVE_VS_INVASIVE  State B: gene has pure CHI SNPs only
#                         (max_fst_chi populated, max_fst_among_invasive NA)
#     LAYERED_SELECTION   State C: pure SNPs on both axes (two independent
#                         selection signals at different sites in the gene)
#     POPULATION_SPECIFIC State D: only mixed-axis SNPs (both primary NA,
#                         mixed_lookup populated; conceptually PBS-style
#                         single-branch selection -- one population odd out)
#
#   Filter (decision 2026-04-24): all four blocks use a single criterion --
#   location_type == "genic". The foreground gene set itself (n_methods >= 2
#   from cross-method concordance) is the quality control; n_outlier_snps
#   and max_n_methods thresholds were dropped from block filters because
#   they were excluding curated featured candidates with strong 2-method
#   support but no 3-method-concordant SNP. Magnitude-based ranking via the
#   per-gene max_FST column replaces threshold-based filtering, letting
#   readers and reviewers choose their own stringency.
#
# inputs:
#   - results/enrichment/foreground_genes.tsv
#       Columns: gene_id, chrom, gene_start, gene_end, strand, gene_name,
#       ncbi_gene_id, n_outlier_snps, consensus_pattern, detecting_methods,
#       location_type. One row per foreground gene (n = 2,648).
#   - results/enrichment/candidate_gene_fst_by_axis.tsv
#       Task 4.4 output with axis-separated max FST per gene
#       (max_fst_among_invasive, max_fst_chi, n_snps_among_inv, etc.).
#   - results/enrichment/gene2go_foreground.tsv
#       Columns: gene_id, go_id, go_term, go_ontology, evidence_source,
#       ortholog_id, blast_identity. Long format (one row per GO term).
#
# outputs:
#   - results/enrichment/supplementary_candidate_table.tsv
#   - results/enrichment/featured_candidate_table.tsv
#   - results/enrichment/plots/manhattan_candidates.pdf
#   - results/enrichment/plots/manhattan_candidates.png
#   - results/enrichment/plots/manhattan_candidates_plotdata.tsv
#       Row-for-row plot data (gene-level): the 2,039 background points plus
#       a `featured` boolean and the functional descriptor/category where
#       applicable, so downstream readers (and chat-style LLMs that can't
#       inspect the raster) can interrogate the plot numerically.
#   - results/enrichment/plots/candidate_summary_dotplot.pdf
#   - results/enrichment/plots/candidate_summary_dotplot.png
#   - results/enrichment/plots/candidate_summary_dotplot_plotdata.tsv
#       The 10 featured rows with the exact (x, y, size, color) mapping used
#       by the dot plot, plus the plotted FST source axis.
#   - results/enrichment/plots/block_distribution_barplot.pdf
#   - results/enrichment/plots/block_distribution_barplot.png
#   - results/enrichment/plots/block_distribution_barplot_plotdata.tsv
#       4-block axis-separation distribution of the genic foreground gene set
#       (AMONG_INVASIVE / NATIVE_VS_INVASIVE / LAYERED_SELECTION /
#       POPULATION_SPECIFIC). Makes the axis-separation framework legible at a
#       glance.
#   - results/enrichment/logs/task4_5_summary.txt
#
# required packages:
#   - tidyverse (dplyr, readr, stringr, tibble, tidyr)
#   - ggplot2
#   - ggrepel
#
# usage examples:
#   # Run all three sub-steps end to end:
#   Rscript scripts/build_candidate_tables_and_figures.R --step all
#
#   # Run just Step 5a (supplementary table) for QC inspection:
#   Rscript scripts/build_candidate_tables_and_figures.R --step 5a
#
#   # After validating 5a, build the featured table:
#   Rscript scripts/build_candidate_tables_and_figures.R --step 5b
#
#   # After validating 5b, build the figures:
#   Rscript scripts/build_candidate_tables_and_figures.R --step 5c

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(ggrastr)
  library(hexbin)
  library(ggnewscale)
})

# --- CLI argument parsing ---------------------------------------------------

parse_args <- function(argv) {
  step <- "all"
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("--step", "-s")) {
      step <- argv[i + 1]
      i <- i + 2
    } else if (a %in% c("--help", "-h")) {
      cat(
        "Usage: Rscript build_candidate_tables_and_figures.R --step {5a|5b|5c|all}\n"
      )
      quit(status = 0)
    } else {
      stop(sprintf("unknown argument: %s", a))
    }
  }
  if (!step %in% c("5a", "5b", "5c", "all")) {
    stop(sprintf("--step must be one of 5a, 5b, 5c, all (got '%s')", step))
  }
  list(step = step)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

# --- Paths ------------------------------------------------------------------

project_root <- if (file.exists("results/enrichment/foreground_genes.tsv")) {
  "."
} else if (file.exists("../results/enrichment/foreground_genes.tsv")) {
  ".."
} else {
  stop("Cannot locate results/enrichment/. Run from project root or scripts/.")
}

in_foreground  <- file.path(project_root, "results/enrichment/foreground_genes.tsv")
in_fst_by_axis <- file.path(project_root, "results/enrichment/candidate_gene_fst_by_axis.tsv")
in_gene2go     <- file.path(project_root, "results/enrichment/gene2go_foreground.tsv")

out_supp       <- file.path(project_root, "results/enrichment/supplementary_candidate_table.tsv")
out_featured   <- file.path(project_root, "results/enrichment/featured_candidate_table.tsv")
out_plot_dir   <- file.path(project_root, "results/enrichment/plots")
out_log_dir    <- file.path(project_root, "results/enrichment/logs")
out_summary    <- file.path(out_log_dir, "task4_5_summary.txt")

dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir,  recursive = TRUE, showWarnings = FALSE)

# --- Logging helper ---------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

# --- Step 5a: Supplementary table -------------------------------------------

build_supplementary_table <- function() {
  log_msg("\n=== Step 5a: Supplementary candidate table ===")

  foreground <- read_tsv(in_foreground, show_col_types = FALSE)
  fst_axis   <- read_tsv(in_fst_by_axis, show_col_types = FALSE)
  gene2go    <- read_tsv(in_gene2go, show_col_types = FALSE)

  log_msg(sprintf("  foreground genes: %d", nrow(foreground)))
  log_msg(sprintf("  axis-FST rows:    %d", nrow(fst_axis)))
  log_msg(sprintf("  gene2go rows:     %d  (unique genes: %d)",
                  nrow(gene2go), dplyr::n_distinct(gene2go$gene_id)))

  # Collapse gene2go to one row per gene: separate BP and MF term strings,
  # max blast_identity, and a single evidence source label.
  go_collapsed <- gene2go %>%
    group_by(gene_id) %>%
    summarise(
      GO_BP_terms = paste(sort(unique(go_term[go_ontology == "BP"])), collapse = "; "),
      GO_MF_terms = paste(sort(unique(go_term[go_ontology == "MF"])), collapse = "; "),
      blast_identity = suppressWarnings(max(blast_identity, na.rm = TRUE)),
      evidence_source = {
        sources <- unique(evidence_source)
        if ("direct" %in% sources) "direct"
        else if (any(grepl("orthology", sources))) "orthology"
        else paste(sources, collapse = ";")
      },
      .groups = "drop"
    ) %>%
    mutate(
      GO_BP_terms    = if_else(GO_BP_terms == "", NA_character_, GO_BP_terms),
      GO_MF_terms    = if_else(GO_MF_terms == "", NA_character_, GO_MF_terms),
      blast_identity = if_else(is.finite(blast_identity), blast_identity, NA_real_)
    )

  # The axis-FST table (Task 4.4 / Step 7e, 19 columns) is the upstream
  # source of truth for gene_state and selection_block. classification and
  # block assignment now happen in `scripts/fst_axis_separation.py`, so this
  # script just reads them. Block filter (handoff §3.1 Step 7f, 2026-04-24
  # decision): location_type == "genic" only -- threshold-based filtering
  # replaced with magnitude-based ranking via max_FST.
  joined <- fst_axis %>%
    select(
      gene_id, chrom, gene_start, gene_end,
      gene_state, selection_block,
      consensus_pattern, n_outlier_snps,
      n_2plus_method_snps, max_n_methods,
      n_snps_among_inv, n_snps_chi, n_snps_mixed,
      max_fst_among_invasive, max_fst_chi,
      max_fst_among_inv_mixed_lookup, max_fst_chi_mixed_lookup,
      location_type, detecting_methods
    ) %>%
    left_join(go_collapsed, by = "gene_id")

  # Per-gene max_FST: the maximum FST across all four axis columns (primary
  # and mixed-lookup). This is the magnitude-based ranking key for the
  # supplementary table -- readers can sort by max_FST to surface the
  # strongest signals regardless of axis. max_FST_source records which of
  # the four columns supplied the max, so the value is interpretable.
  joined <- joined %>%
    rowwise() %>%
    mutate(
      max_FST = suppressWarnings(max(c(
        max_fst_among_invasive, max_fst_chi,
        max_fst_among_inv_mixed_lookup, max_fst_chi_mixed_lookup
      ), na.rm = TRUE)),
      max_FST_source = {
        vals <- c(
          among_invasive          = max_fst_among_invasive,
          chi                     = max_fst_chi,
          among_inv_mixed_lookup  = max_fst_among_inv_mixed_lookup,
          chi_mixed_lookup        = max_fst_chi_mixed_lookup
        )
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) NA_character_ else names(vals)[which.max(vals)]
      }
    ) %>%
    ungroup() %>%
    mutate(max_FST = if_else(is.finite(max_FST), max_FST, NA_real_))

  # Log full-foreground state census (for QC / methods section reporting).
  state_census <- joined %>%
    count(gene_state, name = "n_genes") %>%
    arrange(gene_state)
  log_msg("  foreground gene-state census (pre-filter):")
  for (i in seq_len(nrow(state_census))) {
    st <- state_census$gene_state[i]
    st_key <- if (is.na(st)) "NA" else st
    label <- switch(st_key,
                    A = "A: among-invasive only",
                    B = "B: native-vs-invasive only",
                    C = "C: both pure axes (layered)",
                    D = "D: mixed-axis only",
                    "(unclassified)")
    log_msg(sprintf("    %-30s %d", label, state_census$n_genes[i]))
  }

  supplementary <- joined %>%
    filter(!is.na(selection_block)) %>%
    mutate(selection_block = factor(
      selection_block,
      levels = c("AMONG_INVASIVE", "NATIVE_VS_INVASIVE",
                 "LAYERED_SELECTION", "POPULATION_SPECIFIC")
    )) %>%
    # Sort: block first (preserves grouping for table layout), then by
    # max_FST descending (magnitude-based ranking within each block), then
    # by n_outlier_snps as a secondary tiebreaker.
    arrange(selection_block, desc(max_FST), desc(n_outlier_snps)) %>%
    select(
      gene_id, chrom, gene_start, gene_end,
      selection_block, gene_state,
      consensus_pattern,           # cross-reference to Task 4.3 enrichment stratification
      max_FST, max_FST_source,     # magnitude-based ranking key
      n_outlier_snps, max_n_methods,
      n_snps_among_inv, n_snps_chi, n_snps_mixed,
      max_fst_among_invasive, max_fst_chi,
      max_fst_among_inv_mixed_lookup, max_fst_chi_mixed_lookup,
      location_type, detecting_methods,
      GO_BP_terms, GO_MF_terms, evidence_source, blast_identity
    )

  write_tsv(supplementary, out_supp, na = "NA")

  log_msg(sprintf("  wrote %s (%d rows)", out_supp, nrow(supplementary)))
  block_counts <- supplementary %>%
    count(selection_block, name = "n") %>%
    arrange(selection_block)
  log_msg("  block counts:")
  for (i in seq_len(nrow(block_counts))) {
    log_msg(sprintf("    %-18s %d",
                    as.character(block_counts$selection_block[i]),
                    block_counts$n[i]))
  }

  go_cov <- supplementary %>%
    summarise(
      any_go = sum(!is.na(GO_BP_terms) | !is.na(GO_MF_terms)),
      bp     = sum(!is.na(GO_BP_terms)),
      mf     = sum(!is.na(GO_MF_terms))
    )
  log_msg(sprintf("  GO coverage (any | BP | MF): %d | %d | %d out of %d",
                  go_cov$any_go, go_cov$bp, go_cov$mf, nrow(supplementary)))

  invisible(supplementary)
}

# --- Step 5b: Featured table ------------------------------------------------

# Featured candidate set for the manuscript. Each row is one candidate with
# a short functional descriptor and a functional category used for figure
# grouping.
#
# Task 4.7 Step 7f notes (2026-04-24):
#   - LOC123319966 was NEVER in this tribble (verified against archived
#     featured_candidate_table.tsv); handoff resolved question's "drop from
#     featured table" was imprecise wording for "do not promote" -- the gene
#     was discussed as a featured candidate during Task 4.5 spec'ing but
#     never made the published 10-gene set. Tribble stays at 10 rows.
#   - LOC123318563 label updated: under Step 7e pure-SNP-only convention
#     this gene is AMONG_INVASIVE (State A) driven by 16 pure EEU_vs_USA
#     SNPs at max 0.314; the 0.930 CHI signal lives only on 2 mixed-axis
#     SNPs and is now reported in the mixed-lookup columns, not the primary
#     axis. The "DISCORDANT cluster" framing in the original label is
#     replaced with "EEU-USA differentiation" to match the corrected biology.
#   - Replacement candidates for any future expansion or substitution are
#     listed in handoff §3.1 Task 4.7 Step 7f (rebuilt 2026-04-24
#     supplementary table). Top picks: LOC123307937 (formaldehyde
#     catabolism / detoxification) and LOC123310273 (oxidative stress
#     response).
featured_spec <- tibble::tribble(
  ~gene_id,       ~functional_annotation,                        ~functional_category,
  "LOC123318554", "phospholipase C",                             "Lipid metabolism / membrane",
  "LOC123316249", "G protein-coupled receptor",                  "Sensory / chemoreception",
  "LOC123311185", "EGF receptor signaling",                      "Growth / proliferation",
  "LOC123314071", "inositol polyphosphate kinase",               "Inositol / phospholipid signaling",
  "LOC123315519", "cation channel",                              "Ion homeostasis",
  "LOC123318563", "TWiK K+ leak channel (EEU-USA differentiation)", "Ion channel / K+ transport",
  "LOC123318056", "Toll pathway component",                      "Immune / Toll signaling",
  "LOC123316739", "dynein, germ cell development",               "Cytoskeletal / developmental",
  "LOC123317320", "monocarboxylic acid transporter",             "Metabolic transport",
  "LOC123312215", "fatty acid-CoA ligase",                       "Fatty acid metabolism"
)

build_featured_table <- function(supplementary = NULL) {
  log_msg("\n=== Step 5b: Featured candidate table ===")

  if (is.null(supplementary)) {
    if (!file.exists(out_supp)) {
      stop("Supplementary table not found. Run Step 5a first.")
    }
    supplementary <- read_tsv(out_supp, show_col_types = FALSE)
  }

  # Verify every proposed featured gene is in the supplementary table. Any
  # gene that fell out (e.g., didn't meet block inclusion criteria) is a
  # flag to the user, not a silent drop.
  missing <- setdiff(featured_spec$gene_id, supplementary$gene_id)
  if (length(missing) > 0) {
    log_msg(sprintf("  WARNING: %d featured gene(s) not in supplementary table: %s",
                    length(missing), paste(missing, collapse = ", ")))
  }

  featured <- supplementary %>%
    filter(gene_id %in% featured_spec$gene_id) %>%
    left_join(featured_spec, by = "gene_id") %>%
    mutate(
      # Display FST priority: primary among-invasive -> primary CHI ->
      # mixed-axis lookup (if primaries are both NA, e.g. a POPULATION_SPECIFIC
      # candidate). Report the axis source explicitly so the manuscript can
      # cite the right value.
      featured_fst = coalesce(
        max_fst_among_invasive, max_fst_chi,
        max_fst_among_inv_mixed_lookup, max_fst_chi_mixed_lookup
      ),
      featured_fst_source = case_when(
        !is.na(max_fst_among_invasive)         ~ "among_invasive",
        !is.na(max_fst_chi)                    ~ "chi",
        !is.na(max_fst_among_inv_mixed_lookup) ~ "mixed_lookup_among_inv",
        !is.na(max_fst_chi_mixed_lookup)       ~ "mixed_lookup_chi",
        TRUE                                   ~ NA_character_
      )
    ) %>%
    select(
      gene_id,
      chromosome            = chrom,
      gene_start,
      gene_end,
      selection_block,
      gene_state,
      consensus_pattern,              # cross-reference to Task 4.3
      max_FST,
      max_FST_source,
      n_outlier_SNPs        = n_outlier_snps,
      max_methods           = max_n_methods,
      among_invasive_max_FST            = max_fst_among_invasive,
      chi_max_FST                       = max_fst_chi,
      among_inv_max_FST_mixed_lookup    = max_fst_among_inv_mixed_lookup,
      chi_max_FST_mixed_lookup          = max_fst_chi_mixed_lookup,
      featured_FST_display              = featured_fst,
      featured_FST_source               = featured_fst_source,
      functional_annotation,
      functional_category
    ) %>%
    arrange(selection_block, desc(max_FST), gene_id)

  write_tsv(featured, out_featured, na = "NA")

  log_msg(sprintf("  wrote %s (%d rows)", out_featured, nrow(featured)))
  log_msg(sprintf("  functional categories represented: %d (%s)",
                  dplyr::n_distinct(featured$functional_category),
                  paste(sort(unique(featured$functional_category)), collapse = "; ")))
  log_msg("  per-gene summary:")
  for (i in seq_len(nrow(featured))) {
    ai_fst <- featured$among_invasive_max_FST[i]
    chi_fst <- featured$chi_max_FST[i]
    log_msg(sprintf(
      "    %-14s block=%-20s state=%s  n_SNPs=%d  among_inv=%s  chi=%s  category=%s",
      featured$gene_id[i],
      as.character(featured$selection_block[i]),
      as.character(featured$gene_state[i]),
      featured$n_outlier_SNPs[i],
      if (is.na(ai_fst))  "NA   " else format(round(ai_fst, 3), nsmall = 3),
      if (is.na(chi_fst)) "NA   " else format(round(chi_fst, 3), nsmall = 3),
      featured$functional_category[i]
    ))
  }

  invisible(featured)
}

# --- Step 5c: Figures -------------------------------------------------------

# Genome-wide among-invasive FST from Section 5.2 expected ranges. Falls in
# 0.006-0.013; using a conservative midpoint of 0.010 for the reference line.
GW_AMONG_INV_FST <- 0.010

# --- Color palettes (see handoff Section 5.3) -------------------------------

# Population palette: IBM colorblind-safe, mandatory for any population-level
# distinction. Not used in this script's candidate figures directly (no pop
# stratification), but declared here for consistency if the script is ever
# extended with pop-level overlays.
pop_colors <- c(
  CHI = "#648FFF",
  EEU = "#785EF0",
  WEU = "#DC267F",
  USA = "#FFB000"
)

# Axis/block palette: 4 colors distinct from the 4 pop colors, used for
# selection_block stratification in the Manhattan and dot-plot legends.
# Choices:
#   AMONG_INVASIVE      teal    (#009E73) -- Okabe-Ito green, good contrast to
#                                            the orange USA color
#   NATIVE_VS_INVASIVE  magenta (#CC79A7) -- Okabe-Ito pink, distinct from WEU
#   LAYERED_SELECTION   amber   (#E69F00) -- Okabe-Ito orange, darker than USA
#   POPULATION_SPECIFIC black   (#000000) -- maximum contrast for the rarest
#                                            category so it stands out
block_colors <- c(
  AMONG_INVASIVE      = "#009E73",
  NATIVE_VS_INVASIVE  = "#CC79A7",
  LAYERED_SELECTION   = "#E69F00",
  POPULATION_SPECIFIC = "#000000"
)

# Functional-theme palette for the mirrored Manhattan (Figure 1). Okabe-Ito
# subset chosen to be colorblind-safe and visually distinct from `pop_colors`.
# Note that "Cell signaling" (#009E73) shares a hue with block_colors
# AMONG_INVASIVE (#009E73), and "Growth & development" (#E69F00) shares a hue
# with block_colors LAYERED_SELECTION (#E69F00). Acceptable because the two
# palettes are never co-displayed (block_colors only appears in the dot plot
# and block barplot, theme_colors only in the Manhattan).
theme_colors <- c(
  "Metabolism"               = "#D55E00",
  "Membrane / ion transport" = "#56B4E9",
  "Cell signaling"           = "#009E73",
  "Growth & development"     = "#E69F00",
  "Sensory"                  = "#000000"
)

build_figures <- function() {
  log_msg("\n=== Step 5c: Figures ===")

  if (!file.exists(out_supp) || !file.exists(out_featured)) {
    stop("Supplementary and featured tables must exist. Run 5a and 5b first.")
  }

  supplementary <- read_tsv(out_supp, show_col_types = FALSE)
  featured      <- read_tsv(out_featured, show_col_types = FALSE)
  fst_axis      <- read_tsv(in_fst_by_axis, show_col_types = FALSE)

  # --- Manhattan: mirrored two-panel design (Figure 1) --------------------
  # Top panel: per-SNP max FST across the 3 CHI-vs-X comparisons (CHI filter
  # applied). Bottom panel: per-SNP max FST across the 3 invasive-vs-invasive
  # comparisons (no CHI filter). Bottom panel y-axis is reversed so peaks
  # grow downward and the chromosome labels sit in the gap between panels.
  # Y-axes are scaled independently because native-vs-invasive FST reaches
  # systematically higher values than among-invasive FST -- caller-defined
  # bounds: top panel 0..1.0, bottom panel 0..0.60.
  #
  # Per-SNP inputs and per-axis top-1% thresholds come from
  # scripts/manhattan_mirrored_dataprep.R. That script must be run first.
  #
  # Figure-specific LAYERED rule (replaces the gene-level Step 7e block call
  # for visual logic only): a featured gene gets the dashed vertical
  # connector + dual-panel label iff its representative SNPs cross the
  # panel-specific top-1% threshold on BOTH panels. This keeps the visual
  # logic tied to what the reader sees on the figure rather than to the
  # gene-level pure-SNP-only block call (which can route mixed-axis SNPs
  # away from the primary axis columns and thereby hide a strong cross-axis
  # SNP-level signal). The gene-level Step 7e block call is preserved
  # everywhere else in the manuscript (supplementary table, dot plot,
  # block-distribution barplot, GO enrichment stratification).
  dataprep_dir <- file.path(out_plot_dir, "manhattan_axis_dataprep")
  per_snp_in   <- file.path(dataprep_dir, "per_snp_axis_fst.tsv")
  rep_in       <- file.path(dataprep_dir, "featured_rep_snps.tsv")
  thresh_in    <- file.path(dataprep_dir, "axis_thresholds.tsv")
  offsets_in   <- file.path(dataprep_dir, "chrom_offsets.tsv")

  required_in <- c(per_snp_in, rep_in, thresh_in, offsets_in)
  if (!all(file.exists(required_in))) {
    missing <- required_in[!file.exists(required_in)]
    stop(sprintf(
      "Manhattan dataprep outputs missing (%s). Run scripts/manhattan_mirrored_dataprep.R first.",
      paste(basename(missing), collapse = ", ")
    ))
  }

  per_snp    <- read_tsv(per_snp_in,  show_col_types = FALSE)
  rep_df     <- read_tsv(rep_in,      show_col_types = FALSE)
  thresh_df  <- read_tsv(thresh_in,   show_col_types = FALSE)
  offsets_df <- read_tsv(offsets_in,  show_col_types = FALSE)

  native_top1   <- thresh_df$top1pct_threshold[thresh_df$axis == "native_vs_invasive"]
  invasive_top1 <- thresh_df$top1pct_threshold[thresh_df$axis == "among_invasive"]

  log_msg(sprintf("  per-SNP background: %d rows", nrow(per_snp)))
  log_msg(sprintf("  native-axis top-1%% threshold:  %.4f", native_top1))
  log_msg(sprintf("  invasive-axis top-1%% threshold: %.4f", invasive_top1))

  # Apply figure-specific LAYERED rule.
  rep_df <- rep_df %>%
    mutate(
      crosses_top = !is.na(top_rep_native_fst)   & top_rep_native_fst   >= native_top1,
      crosses_bot = !is.na(bot_rep_invasive_fst) & bot_rep_invasive_fst >= invasive_top1,
      figure_treatment = case_when(
        crosses_top & crosses_bot ~ "layered",
        crosses_bot               ~ "bottom_only",
        crosses_top               ~ "top_only",
        TRUE                      ~ "neither"
      ),
      label_top = figure_treatment %in% c("layered", "top_only"),
      label_bot = figure_treatment %in% c("layered", "bottom_only"),
      is_layered_guide_target = figure_treatment == "layered"
    )

  log_msg(sprintf("  figure_treatment: layered=%d  bottom_only=%d  top_only=%d  neither=%d",
                  sum(rep_df$figure_treatment == "layered"),
                  sum(rep_df$figure_treatment == "bottom_only"),
                  sum(rep_df$figure_treatment == "top_only"),
                  sum(rep_df$figure_treatment == "neither")))

  # Alternating chromosome bands: subtle gray for odd chromosomes (chr1, 3,
  # 5, 7, 9). Lays under both panels at the same x-extent.
  band_df <- offsets_df %>%
    filter((chrom_idx %% 2) == 1) %>%
    select(CHROM, min_cum, max_cum)

  # Y-axis bounds chosen for the manuscript figure: native panel runs 0..1
  # (FST is bounded above by 1; the small CHI sample drives saturation up
  # there and the threshold lives near the ceiling), invasive panel runs
  # 0..0.60 (max observed ~0.59, threshold ~0.19, plenty of room to read).
  y_top_max <- 1.00
  y_bot_max <- 0.60

  # Color/size constants. The cloud layer is a hex-bin density heatmap
  # (geom_hex) with TWO alternating fill gradients keyed to chromosome
  # parity (odd vs even chrom_idx), so that adjacent chromosomes are
  # visually distinct -- a colored take on the classic alternating-chrom
  # Manhattan convention while still encoding per-bin SNP count via fill
  # lightness. Both gradients share a light base near the band color so
  # low-density bins read as "background" rather than as additional
  # signal. Outlier highlight is drawn as solid black points on top of
  # the hex bins for a clean accent.
  HEX_LOW_ODD    <- "#DCE5EE"   # light blue tint, low count, odd chroms
  HEX_HIGH_ODD   <- "#1F4F7A"   # deep navy, high count, odd chroms
  HEX_LOW_EVEN   <- "#EFE3D2"   # light warm tint, low count, even chroms
  HEX_HIGH_EVEN  <- "#7A4818"   # deep amber/brown, high count, even chroms
  OUTLIER_COLOR  <- "#000000"
  OUTLIER_SIZE   <- 0.9
  OUTLIER_ALPHA  <- 0.85
  FEATURED_SIZE  <- 2.6
  FEATURED_STROKE <- 0.5
  GUIDE_COLOR    <- "#3A3A3A"   # darker than before so the dashed cross-
                                # panel connectors stay readable on top of
                                # the hex density heatmap
  GUIDE_LINEWIDTH <- 0.5
  GUIDE_ALPHA     <- 0.95
  LABEL_FONT_SIZE <- 2.4        # ~7pt
  # Hex bin sizes in data units. binwidth_x ≈ 4 Mb gives ~150 bins across
  # the genome (cum_pos spans ~600 Mb). y binwidths sized so visual hex
  # heights match between the 0..1 top panel and the 0..0.6 bottom panel.
  HEX_BINWIDTH_X     <- 4e6
  HEX_BINWIDTH_Y_TOP <- 0.020
  HEX_BINWIDTH_Y_BOT <- 0.012

  # Tag every per-SNP row with its chromosome parity (odd vs even index in
  # the autosome ordering) so the two-gradient hex layers below have a
  # clean filter key. offsets_df has chrom_idx 1..9.
  chrom_parity <- offsets_df %>%
    transmute(CHROM, chrom_idx, alt = if_else(chrom_idx %% 2 == 1, "odd", "even"))
  per_snp <- per_snp %>% left_join(chrom_parity, by = "CHROM")

  base_axis_breaks_top <- seq(0, y_top_max, by = 0.20)
  base_axis_breaks_bot <- seq(0, y_bot_max, by = 0.10)

  # Build a small data frame for the alternating bands so geom_rect inherits
  # nothing extra.
  band_layer_top <- geom_rect(
    data = band_df,
    aes(xmin = min_cum, xmax = max_cum, ymin = -Inf, ymax = Inf),
    fill = "#F1EFE8", color = NA, inherit.aes = FALSE
  )
  band_layer_bot <- geom_rect(
    data = band_df,
    aes(xmin = min_cum, xmax = max_cum, ymin = -Inf, ymax = Inf),
    fill = "#F1EFE8", color = NA, inherit.aes = FALSE
  )

  # --- Top panel: native-vs-invasive ----------------------------------------
  per_snp_top <- per_snp %>% filter(!is.na(native_max_fst))
  per_snp_top_odd  <- per_snp_top %>% filter(alt == "odd")
  per_snp_top_even <- per_snp_top %>% filter(alt == "even")
  per_snp_top_outliers <- per_snp_top %>% filter(is_outlier_2plus)
  rep_top_pts <- rep_df %>%
    filter(!is.na(top_rep_native_fst)) %>%
    mutate(y_plot = pmin(top_rep_native_fst, y_top_max))
  guide_top <- rep_df %>% filter(is_layered_guide_target,
                                 !is.na(top_rep_cum_pos),
                                 !is.na(top_rep_native_fst))
  # Label every featured gene that has a top-panel rep SNP, regardless of
  # threshold-crossing status. The dashed connector still encodes the
  # figure-specific layered classification; the label set is no longer
  # gated on it.
  label_top_df <- rep_df %>% filter(!is.na(top_rep_cum_pos),
                                    !is.na(top_rep_native_fst))

  p_top <- ggplot() +
    band_layer_top +
    # Odd-chromosome hex layer (chr1, 3, 5, 7, 9): blue gradient.
    geom_hex(
      data = per_snp_top_odd,
      aes(x = cum_pos, y = native_max_fst, fill = after_stat(count)),
      binwidth = c(HEX_BINWIDTH_X, HEX_BINWIDTH_Y_TOP),
      color = NA
    ) +
    scale_fill_gradient(
      low = HEX_LOW_ODD, high = HEX_HIGH_ODD,
      trans = "log10", name = NULL, guide = "none",
      na.value = "transparent"
    ) +
    ggnewscale::new_scale_fill() +
    # Even-chromosome hex layer (chr2, 4, 6, 8): amber gradient.
    geom_hex(
      data = per_snp_top_even,
      aes(x = cum_pos, y = native_max_fst, fill = after_stat(count)),
      binwidth = c(HEX_BINWIDTH_X, HEX_BINWIDTH_Y_TOP),
      color = NA
    ) +
    scale_fill_gradient(
      low = HEX_LOW_EVEN, high = HEX_HIGH_EVEN,
      trans = "log10", name = NULL, guide = "none",
      na.value = "transparent"
    ) +
    geom_point(
      data = per_snp_top_outliers,
      aes(x = cum_pos, y = native_max_fst),
      size = OUTLIER_SIZE, alpha = OUTLIER_ALPHA, color = OUTLIER_COLOR,
      shape = 16, inherit.aes = FALSE
    ) +
    # No horizontal threshold line drawn. The native-axis top-1% threshold
    # (~0.93) is a sample-size artifact of CHI n=5 saturating F_ST at 1, so
    # rendering it would mislead readers into reading the line as a
    # biological cutoff. The figure-specific layered rule (which determines
    # the LAYERED dashed vertical guides below) still uses both panel
    # thresholds as its mechanical criterion; the rule is described in
    # the caption rather than visually anchored on the panel.
    # LAYERED guide stubs (drawn AFTER the cloud and outliers so they sit
    # on top): vertical dashed segment from rep point down to the panel
    # x-axis (y = 0). Pairs with a matching stub in the bottom panel; the
    # eye fills in the gap across the chromosome-label area to read this
    # as a single cross-panel connector for genes with strong signal on
    # both selection axes.
    geom_segment(
      data = guide_top,
      aes(x = top_rep_cum_pos, xend = top_rep_cum_pos,
          y = top_rep_native_fst, yend = 0),
      linetype = "dashed", color = GUIDE_COLOR,
      linewidth = GUIDE_LINEWIDTH, alpha = GUIDE_ALPHA,
      inherit.aes = FALSE
    ) +
    # Reset the fill scale so featured points get their own theme palette
    # rather than inheriting the hex density gradient.
    ggnewscale::new_scale_fill() +
    geom_point(
      data = rep_top_pts,
      aes(x = top_rep_cum_pos, y = y_plot, fill = functional_theme),
      shape = 21, color = "black", size = FEATURED_SIZE, stroke = FEATURED_STROKE,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = label_top_df,
      aes(x = top_rep_cum_pos, y = top_rep_native_fst,
          label = functional_annotation),
      size = LABEL_FONT_SIZE, fontface = "plain", color = "black",
      fill = scales::alpha("white", 0.85),
      label.size = 0.15, label.r = unit(0.08, "lines"),
      label.padding = unit(0.12, "lines"),
      min.segment.length = 0, segment.size = 0.25, segment.color = "grey30",
      box.padding = 0.45, point.padding = 0.25, max.overlaps = Inf,
      seed = 42, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = theme_colors, drop = FALSE,
                      name = "Functional theme") +
    scale_x_continuous(
      breaks = offsets_df$center, labels = offsets_df$label,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    scale_y_continuous(
      limits = c(0, y_top_max), breaks = base_axis_breaks_top,
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = expression(F[ST] ~ "(CHI vs invasive, max)")
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x  = element_text(size = 7),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      legend.position = "bottom",
      legend.key.size = unit(3, "mm"),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8)
    )

  # --- Bottom panel: among-invasive (mirrored) ------------------------------
  per_snp_bot <- per_snp %>% filter(!is.na(invasive_max_fst))
  per_snp_bot_odd  <- per_snp_bot %>% filter(alt == "odd")
  per_snp_bot_even <- per_snp_bot %>% filter(alt == "even")
  per_snp_bot_outliers <- per_snp_bot %>% filter(is_outlier_2plus)
  rep_bot_pts <- rep_df %>%
    filter(!is.na(bot_rep_invasive_fst)) %>%
    mutate(y_plot = pmin(bot_rep_invasive_fst, y_bot_max))
  guide_bot <- rep_df %>% filter(is_layered_guide_target,
                                 !is.na(bot_rep_cum_pos),
                                 !is.na(bot_rep_invasive_fst))
  label_bot_df <- rep_df %>% filter(!is.na(bot_rep_cum_pos),
                                    !is.na(bot_rep_invasive_fst))

  p_bot <- ggplot() +
    band_layer_bot +
    geom_hex(
      data = per_snp_bot_odd,
      aes(x = cum_pos, y = invasive_max_fst, fill = after_stat(count)),
      binwidth = c(HEX_BINWIDTH_X, HEX_BINWIDTH_Y_BOT),
      color = NA
    ) +
    scale_fill_gradient(
      low = HEX_LOW_ODD, high = HEX_HIGH_ODD,
      trans = "log10", name = NULL, guide = "none",
      na.value = "transparent"
    ) +
    ggnewscale::new_scale_fill() +
    geom_hex(
      data = per_snp_bot_even,
      aes(x = cum_pos, y = invasive_max_fst, fill = after_stat(count)),
      binwidth = c(HEX_BINWIDTH_X, HEX_BINWIDTH_Y_BOT),
      color = NA
    ) +
    scale_fill_gradient(
      low = HEX_LOW_EVEN, high = HEX_HIGH_EVEN,
      trans = "log10", name = NULL, guide = "none",
      na.value = "transparent"
    ) +
    geom_point(
      data = per_snp_bot_outliers,
      aes(x = cum_pos, y = invasive_max_fst),
      size = OUTLIER_SIZE, alpha = OUTLIER_ALPHA, color = OUTLIER_COLOR,
      shape = 16, inherit.aes = FALSE
    ) +
    # No horizontal threshold line drawn (see top panel comment). The
    # bottom-panel threshold (~0.190) is biologically meaningful (well-
    # sampled among-invasive comparisons), but is omitted here for visual
    # symmetry with the top panel. The figure-specific layered rule still
    # uses it mechanically; the rule is described in the caption.
    geom_segment(
      data = guide_bot,
      aes(x = bot_rep_cum_pos, xend = bot_rep_cum_pos,
          y = 0, yend = bot_rep_invasive_fst),
      linetype = "dashed", color = GUIDE_COLOR,
      linewidth = GUIDE_LINEWIDTH, alpha = GUIDE_ALPHA,
      inherit.aes = FALSE
    ) +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = rep_bot_pts,
      aes(x = bot_rep_cum_pos, y = y_plot, fill = functional_theme),
      shape = 21, color = "black", size = FEATURED_SIZE, stroke = FEATURED_STROKE,
      inherit.aes = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = label_bot_df,
      aes(x = bot_rep_cum_pos, y = bot_rep_invasive_fst,
          label = functional_annotation),
      size = LABEL_FONT_SIZE, fontface = "plain", color = "black",
      fill = scales::alpha("white", 0.85),
      label.size = 0.15, label.r = unit(0.08, "lines"),
      label.padding = unit(0.12, "lines"),
      min.segment.length = 0, segment.size = 0.25, segment.color = "grey30",
      box.padding = 0.45, point.padding = 0.25, max.overlaps = Inf,
      seed = 42, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = theme_colors, drop = FALSE,
                      name = "Functional theme") +
    scale_x_continuous(
      breaks = offsets_df$center, labels = offsets_df$label,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    # scale_y_reverse: data y is in natural FST units, but the visual
    # mapping is flipped so y=0 sits at the TOP of the panel (against the
    # gap-with-labels) and the panel maximum sits at the bottom. Peaks
    # grow downward. Limits are given in data order; ggplot reverses.
    scale_y_reverse(
      limits = c(y_bot_max, 0), breaks = base_axis_breaks_bot,
      expand = expansion(mult = c(0.02, 0))
    ) +
    labs(
      x = NULL,
      y = expression(F[ST] ~ "(among invasive, max)")
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x  = element_blank(),     # chromosome labels live on top panel
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y  = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      legend.position = "bottom",
      legend.key.size = unit(3, "mm"),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 8)
    )

  # --- Combine with patchwork ----------------------------------------------
  # Equal panel heights; `guides = "collect"` merges the theme-color fill
  # legends from each panel into a single shared one anchored below both
  # panels. Hex-density gradients are not collected (guide = "none").
  # Embedded "how to read" caption. Each element of `caption_lines` becomes
  # one rendered line; ggplot's plot.caption does NOT auto-wrap, so explicit
  # line breaks are required. Keep each line under ~130 characters at 6.5pt
  # to fit a 183mm-wide figure without horizontal clipping.
  caption_lines <- c(
    "How to read",
    "Top panel: per-SNP F_ST across the 3 CHI-vs-invasive comparisons (max of CHI-EEU, CHI-WEU, CHI-USA).",
    "        CHI genotype-count filter applied (>=3 of 5 CHI individuals genotyped per site).",
    "Bottom panel: per-SNP F_ST across the 3 among-invasive comparisons (max of EEU-WEU, EEU-USA, WEU-USA), y-axis flipped.",
    "Y-axes scaled independently (top 0..1, bottom 0..0.6); peak heights are NOT directly comparable across panels.",
    "Hex shading: per-bin SNP density (log scale); blue on odd chromosomes (1, 3, 5, 7, 9), amber on even (2, 4, 6, 8).",
    "Black dots: SNPs flagged by 2+ of 3 selection-scan methods (per-SNP F_ST percentile, OutFLANK, pcadapt; n = 4,183).",
    "Colored points: 10 featured candidate genes at their per-axis max-F_ST SNP, colored by functional theme.",
    "Dashed vertical lines: connect featured genes whose representative SNPs sit in the panel-specific top 1% of",
    "        per-SNP F_ST on BOTH panels (figure-specific layered classification).",
    "Sex chromosome (NC_058198.1) excluded throughout."
  )
  on_figure_caption <- paste(caption_lines, collapse = "\n")

  p_manhattan <- (p_top / p_bot) +
    plot_layout(heights = c(1, 1), guides = "collect") +
    plot_annotation(
      caption = on_figure_caption,
      theme = theme(plot.caption = element_text(
        size = 6.5, hjust = 0, lineheight = 1.15,
        margin = margin(t = 4, b = 0)
      ))
    ) &
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.box.just = "center") &
    guides(fill = guide_legend(title = "Functional theme",
                               nrow = 1, byrow = TRUE,
                               override.aes = list(size = 3.5)))

  pdf_path <- file.path(out_plot_dir, "manhattan_candidates.pdf")
  png_path <- file.path(out_plot_dir, "manhattan_candidates.png")
  ggsave(pdf_path, p_manhattan, width = 183, height = 175, units = "mm")
  ggsave(png_path, p_manhattan, width = 183, height = 175, units = "mm",
         dpi = 300)
  pdf_kb <- round(file.info(pdf_path)$size / 1024, 1)
  log_msg(sprintf("  wrote %s (%.1f KB)", pdf_path, pdf_kb))
  log_msg(sprintf("  wrote %s",            png_path))

  # --- Caption file --------------------------------------------------------
  caption_path <- file.path(out_plot_dir, "manhattan_candidates_caption.txt")
  caption <- paste(
    "Figure 1. Mirrored Manhattan plot showing per-SNP F_ST across the",
    "Coccinella septempunctata genome on two selection axes. Top panel:",
    "maximum per-SNP F_ST across the three CHI-vs-invasive comparisons",
    "(CHI-EEU, CHI-WEU, CHI-USA), with the CHI genotype-count filter",
    "applied (>= 3 of 5 CHI individuals genotyped per site). Bottom panel:",
    "maximum per-SNP F_ST across the three among-invasive comparisons",
    "(EEU-WEU, EEU-USA, WEU-USA), with the y-axis reversed so peaks grow",
    "downward. Y-axes are scaled independently because native-vs-invasive",
    "F_ST reaches systematically higher values than among-invasive F_ST",
    "(driven in part by the small CHI sample, n=5); tall peaks in the two",
    "panels are NOT directly comparable in absolute magnitude. Background",
    "shading: hex-bin density of all SNPs passing per-axis filtering",
    "(monochromatic gradient, log-scaled count; lighter = sparse, darker",
    "= dense). Black points: SNPs flagged as outliers by 2 or more of 3",
    "selection-scan methods (per-SNP F_ST percentile, OutFLANK, pcadapt) in",
    "the cross-method concordance analysis (n = 4,183). Note that black",
    "dots span the full F_ST range on each panel because method-flagging",
    "is orthogonal to per-axis max F_ST: pcadapt's multivariate signal can",
    "flag SNPs that are not in the per-axis top-1% F_ST tail, and a SNP",
    "flagged by methods on one axis may sit at low max F_ST on the other.",
    "Colored points: the 10 featured candidate genes plotted at their per-",
    "axis maximum-F_ST SNP, colored by functional theme. Vertical dashed",
    "lines connect featured genes whose representative SNPs sit in the",
    sprintf("panel-specific top 1%% of per-SNP F_ST on BOTH panels (top: %.3f;",
            native_top1),
    sprintf("bottom: %.3f; figure-specific layered classification, applied for",
            invasive_top1),
    "visual logic only). The gene-level Step 7e block classification used",
    "in the supplementary table and downstream enrichment stratification",
    "is preserved as a separate, gene-level scheme that may differ from",
    "this SNP-level visual rule by individual genes.",
    "Sex chromosome (NC_058198.1) excluded throughout. Chromosomes 1-9",
    "correspond to NC_058189.1 through NC_058197.1. See Supplementary",
    "Table SX for the full candidate gene list and Figure 2 for the",
    "per-axis gene-level dot plot.",
    sep = " "
  )
  writeLines(caption, caption_path)
  log_msg(sprintf("  wrote %s", caption_path))

  # --- Plot-data TSV -------------------------------------------------------
  # One row per (gene x panel) for the 10 featured genes. Captures every
  # featured point on the figure with the value plotted, the panel it sits
  # on, the figure-specific layered classification, and whether the gene
  # gets a label / a dashed connector. Background-cloud and outlier-layer
  # SNP-level data are already in
  # results/enrichment/plots/manhattan_axis_dataprep/per_snp_axis_fst.tsv.
  plotdata_top <- rep_df %>%
    transmute(
      gene_id, chromosome, gene_start, gene_end,
      selection_block_gene_level = selection_block,
      figure_treatment, is_layered_guide_target,
      panel = "top",
      rep_pos       = top_rep_pos,
      cum_pos       = top_rep_cum_pos,
      axis_max_fst  = top_rep_native_fst,
      labeled_in_panel = label_top,
      crosses_panel_threshold = crosses_top,
      functional_theme, functional_annotation, functional_category
    )
  plotdata_bot <- rep_df %>%
    transmute(
      gene_id, chromosome, gene_start, gene_end,
      selection_block_gene_level = selection_block,
      figure_treatment, is_layered_guide_target,
      panel = "bottom",
      rep_pos       = bot_rep_pos,
      cum_pos       = bot_rep_cum_pos,
      axis_max_fst  = bot_rep_invasive_fst,
      labeled_in_panel = label_bot,
      crosses_panel_threshold = crosses_bot,
      functional_theme, functional_annotation, functional_category
    )
  manhattan_plotdata <- bind_rows(plotdata_top, plotdata_bot) %>%
    arrange(gene_id, panel)
  write_tsv(manhattan_plotdata,
            file.path(out_plot_dir, "manhattan_candidates_plotdata.tsv"),
            na = "NA")
  log_msg(sprintf("  wrote %s (%d rows; %d featured genes x 2 panels)",
                  file.path(out_plot_dir, "manhattan_candidates_plotdata.tsv"),
                  nrow(manhattan_plotdata),
                  nrow(rep_df)))

  # --- Dot plot: featured candidates grouped by functional category -------
  dot_df <- featured %>%
    mutate(
      y_fst = featured_FST_display,  # prioritized: pure among-inv -> pure CHI
                                     # -> mixed_lookup; computed in Step 5b
      gene_label = paste0(gene_id, "\n", functional_annotation)
    ) %>%
    arrange(functional_category, desc(y_fst)) %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(unique(gene_label))),
      functional_category = factor(functional_category,
                                   levels = unique(functional_category)),
      selection_block = factor(selection_block,
                               levels = names(block_colors))
    )

  # Human-readable legend labels for selection_block. Keys must match
  # names(block_colors) so the labels swap in via scale_color_manual without
  # changing the underlying factor levels (which other code filters on).
  block_legend_labels <- c(
    AMONG_INVASIVE      = "A: among-invasive only",
    NATIVE_VS_INVASIVE  = "B: native-vs-invasive only",
    LAYERED_SELECTION   = "C: layered (both axes)",
    POPULATION_SPECIFIC = "D: population-specific"
  )

  p_dot <- ggplot(dot_df,
                  aes(x = y_fst, y = gene_label,
                      size = n_outlier_SNPs, color = selection_block)) +
    geom_vline(xintercept = GW_AMONG_INV_FST, linetype = "dashed",
               color = "steelblue", linewidth = 0.4) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = block_colors,
                       labels = block_legend_labels,
                       drop = FALSE) +
    scale_size_continuous(range = c(2, 6)) +
    # Force every color-legend key to render at a fixed visible size,
    # otherwise factor levels with no points in the data inherit a missing
    # size aesthetic and their swatches disappear from the legend. Wrap to
    # 2 rows so the longest block label fits within the 183mm figure width.
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1),
                                nrow = 2, byrow = TRUE, order = 1),
           size  = guide_legend(order = 2)) +
    facet_grid(functional_category ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    labs(
      x = expression(max~F[ST]~"(priority: among-invasive > CHI > mixed-lookup)"),
      y = NULL,
      color = "Selection block",
      size = "n outlier SNPs",
      title = "Featured candidate genes by functional category",
      subtitle = sprintf(
        paste0("%d candidates across %d categories\n",
               "dashed line: among-invasive axis genome-wide ",
               "F_ST (~%.3f, null)"),
        nrow(dot_df), dplyr::n_distinct(dot_df$functional_category),
        GW_AMONG_INV_FST
      )
    ) +
    theme_bw(base_size = 9) +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.major.y = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.size = unit(3, "mm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8)
    )

  ggsave(file.path(out_plot_dir, "candidate_summary_dotplot.pdf"),
         p_dot, width = 183, height = 160, units = "mm")
  ggsave(file.path(out_plot_dir, "candidate_summary_dotplot.png"),
         p_dot, width = 183, height = 160, units = "mm", dpi = 300)
  log_msg(sprintf("  wrote %s (+ .png)",
                  file.path(out_plot_dir, "candidate_summary_dotplot.pdf")))

  # Plot-data TSV for the dot plot: the 10 featured rows with the exact
  # aesthetic mapping (x = plotted FST, y = gene label, size = n_SNPs,
  # color = consensus_pattern) used by the figure, plus which axis the x
  # value came from.
  dotplot_plotdata <- dot_df %>%
    transmute(
      gene_id,
      functional_category,
      functional_annotation,
      selection_block,
      gene_state,
      consensus_pattern,
      n_outlier_SNPs,
      plotted_fst = y_fst,
      plotted_fst_source = featured_FST_source,
      among_invasive_max_FST,
      chi_max_FST,
      among_inv_max_FST_mixed_lookup,
      chi_max_FST_mixed_lookup
    )
  write_tsv(dotplot_plotdata,
            file.path(out_plot_dir, "candidate_summary_dotplot_plotdata.tsv"),
            na = "NA")
  log_msg(sprintf("  wrote %s (%d rows)",
                  file.path(out_plot_dir, "candidate_summary_dotplot_plotdata.tsv"),
                  nrow(dotplot_plotdata)))

  # --- 4-block distribution bar plot --------------------------------------
  # Visualizes the axis-separation block census across the genic foreground
  # gene set. The supplementary table is already filtered to
  # location_type == "genic" (Step 5a), so its row count *is* the genic
  # foreground n. Same `block_colors` palette as the dot plot, kept in the
  # canonical AMONG_INVASIVE > NATIVE_VS_INVASIVE > LAYERED_SELECTION >
  # POPULATION_SPECIFIC order.
  block_df <- supplementary %>%
    count(selection_block, name = "n_genes") %>%
    mutate(
      selection_block = factor(selection_block, levels = names(block_colors)),
      # Wrap long block labels onto two lines so the x-axis doesn't collide
      # (LAYERED_SELECTION / POPULATION_SPECIFIC overlap otherwise).
      block_label_wrapped = factor(
        gsub("_", "\n", as.character(selection_block)),
        levels = gsub("_", "\n", names(block_colors))
      ),
      pct   = 100 * n_genes / sum(n_genes),
      label = sprintf("%d\n(%.1f%%)", n_genes, pct)
    ) %>%
    arrange(selection_block)

  total_genic <- sum(block_df$n_genes)

  p_blocks <- ggplot(block_df,
                     aes(x = block_label_wrapped, y = n_genes,
                         fill = selection_block)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    geom_text(aes(label = label), vjust = -0.25, size = 3, lineheight = 0.9) +
    scale_fill_manual(values = block_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(
      x = NULL,
      y = "Foreground genes (genic)",
      title = "Axis-separation framework: genic foreground gene distribution",
      subtitle = sprintf(
        "n = %s genic foreground genes; pure-SNP-only convention (Step 7e)",
        format(total_genic, big.mark = ",")
      )
    ) +
    theme_bw(base_size = 9) +
    theme(
      legend.position    = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 7.5, lineheight = 0.85),
      plot.title         = element_text(size = 10, face = "bold"),
      plot.subtitle      = element_text(size = 8)
    )

  ggsave(file.path(out_plot_dir, "block_distribution_barplot.pdf"),
         p_blocks, width = 130, height = 90, units = "mm")
  ggsave(file.path(out_plot_dir, "block_distribution_barplot.png"),
         p_blocks, width = 130, height = 90, units = "mm", dpi = 300)
  log_msg(sprintf("  wrote %s (+ .png)",
                  file.path(out_plot_dir, "block_distribution_barplot.pdf")))

  write_tsv(block_df,
            file.path(out_plot_dir, "block_distribution_barplot_plotdata.tsv"),
            na = "NA")
  log_msg(sprintf("  block census: %s",
                  paste(sprintf("%s=%d", block_df$selection_block, block_df$n_genes),
                        collapse = "; ")))

  invisible(list(manhattan = p_manhattan, dotplot = p_dot, blocks = p_blocks))
}

# --- Run --------------------------------------------------------------------

log_msg(sprintf("Task 4.5: step = %s", args$step))

if (args$step %in% c("5a", "all")) {
  supp <- build_supplementary_table()
} else {
  supp <- NULL
}

if (args$step %in% c("5b", "all")) {
  build_featured_table(supp)
}

if (args$step %in% c("5c", "all")) {
  build_figures()
}

writeLines(log_lines, out_summary)
cat(sprintf("\nSummary log written to %s\n", out_summary))
