#!/usr/bin/env Rscript

# title: diagnostic_fst_histograms.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-01
# last modified: 2026-05-01
#
# purpose:
#   Manhattan-scatter diagnostic, Task 1. For each of the 6 pairwise
#   comparisons, produces (a) the per-SNP Weir & Cockerham F_ST distribution
#   as a histogram on linear and log10 y-axes, and (b) a summary table of
#   first four moments (mean, median, SD, skewness), the 99th percentile,
#   and the fraction of negative F_ST values (vcftools can return negative
#   F_ST when within-pop variance exceeds between-pop variance). The shape
#   diagnostic distinguishes a healthy right-skewed distribution with a
#   long upper tail (consistent with selection acting in a structured but
#   well-behaved background) from a broad/flat distribution centered well
#   above zero with no clear tail (consistent with structure leakage
#   and/or small-sample noise inflating the bulk of the distribution).
#
#   Skewness is computed manually (third standardized central moment)
#   because the `moments` package is not installed in the project's R
#   library. The bias-corrected adjusted Fisher-Pearson form is reported
#   so it matches what most statistics texts and the moments package
#   return.
#
# inputs:
#   - results/fst/{COMP}/fst_genomewide/{CHROM}.filtered.weir.fst
#       Per-SNP Weir & Cockerham F_ST. Columns: CHROM, POS,
#       WEIR_AND_COCKERHAM_FST. One file per (comparison, chromosome).
#       Comparisons: CHI_vs_EEU, CHI_vs_USA, CHI_vs_WEU, EEU_vs_USA,
#                    EEU_vs_WEU, USA_vs_WEU.
#       Sex chromosome (NC_058198.1) excluded; 9 autosomes used.
#
# outputs:
#   - figures/diagnostic/fst_histograms/
#       {COMP}_linear.png       : Per-comparison histogram, linear y-axis.
#       {COMP}_logy.png         : Per-comparison histogram, log10 y-axis.
#       all_comparisons_linear.png  : 2x3 facet, CHI-axis row vs invasive-axis row.
#       all_comparisons_logy.png    : Same layout, log10 y-axis.
#   - results/diagnostic/fst_distribution_summary.tsv
#       One row per comparison with: comparison, comp_type, n_snps,
#       n_negative, frac_negative, mean, median, sd, skewness, q99.
#
# usage example:
#   Rscript scripts/diagnostic_fst_histograms.R
#   Rscript scripts/diagnostic_fst_histograms.R --keep-negatives FALSE

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

# --- CLI argument parsing ---------------------------------------------------

parse_args <- function(argv) {
  keep_negatives <- TRUE
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("--keep-negatives")) {
      keep_negatives <- as.logical(argv[i + 1])
      if (is.na(keep_negatives)) {
        stop("--keep-negatives must be TRUE or FALSE")
      }
      i <- i + 2
    } else if (a %in% c("--help", "-h")) {
      cat("Usage: Rscript diagnostic_fst_histograms.R [--keep-negatives TRUE|FALSE]\n")
      cat("  --keep-negatives  Include negative F_ST in histograms and stats (default TRUE).\n")
      cat("                    Negative F_ST is a real noise-floor signal for vcftools\n")
      cat("                    Weir-Cockerham; dropping them masks how broad the bulk is.\n")
      quit(status = 0)
    } else {
      stop(sprintf("unknown argument: %s", a))
    }
    i <- i + 1 - 1   # no-op, kept for symmetry with future flags
  }
  list(keep_negatives = keep_negatives)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

# --- Paths ------------------------------------------------------------------

project_root <- if (file.exists("results/fst")) {
  "."
} else if (file.exists("../results/fst")) {
  ".."
} else {
  stop("Cannot locate results/fst/. Run from project root or scripts/.")
}

fst_base   <- file.path(project_root, "results/fst")
fig_dir    <- file.path(project_root, "figures/diagnostic/fst_histograms")
out_tsv    <- file.path(project_root, "results/diagnostic/fst_distribution_summary.tsv")

dir.create(fig_dir,        recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

# --- Configuration ----------------------------------------------------------

chi_comparisons      <- c("CHI_vs_EEU", "CHI_vs_USA", "CHI_vs_WEU")
invasive_comparisons <- c("EEU_vs_USA", "EEU_vs_WEU", "USA_vs_WEU")
all_comparisons      <- c(chi_comparisons, invasive_comparisons)

# 9 autosomes; sex chromosome NC_058198.1 excluded throughout, consistent
# with cross_method_concordance.R and manhattan_mirrored_dataprep.R.
autosomes <- sprintf("NC_0581%02d.1", 89:97)

# Comparison-type lookup for the summary table and facet layout.
comp_type <- function(comp) {
  if (comp %in% chi_comparisons) "CHI_vs_invasive" else "Among_invasive"
}

# --- Logging ----------------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

log_msg("=== Task 1: Per-comparison F_ST histograms ===")
log_msg(sprintf("project root: %s", normalizePath(project_root)))
log_msg(sprintf("keep_negatives: %s", args$keep_negatives))
log_msg(sprintf("autosomes (n=%d): %s", length(autosomes),
                paste(autosomes, collapse = ", ")))

# --- Step 1: Load per-SNP F_ST for all 6 comparisons ------------------------

log_msg("\n[1/3] Loading per-SNP F_ST for 6 comparisons x 9 autosomes")

read_fst_one <- function(comp, chrom) {
  f <- file.path(fst_base, comp, "fst_genomewide",
                 sprintf("%s.filtered.weir.fst", chrom))
  if (!file.exists(f)) {
    stop(sprintf("missing per-SNP F_ST file: %s", f))
  }
  dt <- fread(f, sep = "\t", header = TRUE,
              colClasses = c(CHROM = "character",
                             POS = "integer",
                             WEIR_AND_COCKERHAM_FST = "double"))
  # vcftools writes "-nan"/"nan" for monomorphic sites; fread converts those
  # to NA. Drop NA F_ST so they don't bias the bulk-shape diagnostic.
  dt <- dt[!is.na(WEIR_AND_COCKERHAM_FST)]
  dt[, comparison := comp]
  dt
}

fst_list <- list()
for (comp in all_comparisons) {
  comp_chunks <- lapply(autosomes, read_fst_one, comp = comp)
  fst_list[[comp]] <- rbindlist(comp_chunks)
  log_msg(sprintf("  %-12s %d SNPs (n_negative=%d, %.2f%%)",
                  comp,
                  nrow(fst_list[[comp]]),
                  sum(fst_list[[comp]]$WEIR_AND_COCKERHAM_FST < 0),
                  100 * mean(fst_list[[comp]]$WEIR_AND_COCKERHAM_FST < 0)))
}

# --- Step 2: Summary statistics ---------------------------------------------

log_msg("\n[2/3] Computing summary statistics")

# Manual skewness: bias-corrected adjusted Fisher-Pearson g1, matching the
# default of the `moments::skewness` function. Formula:
#   g1 = (n / ((n-1)*(n-2))) * sum( ((x - mean) / sd)^3 )
# with sd computed using the (n-1) denominator (i.e. base R sd()).
skewness_fp <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  mu <- mean(x)
  s  <- sd(x)
  if (s == 0) return(NA_real_)
  (n / ((n - 1) * (n - 2))) * sum(((x - mu) / s)^3)
}

summarize_one <- function(comp) {
  vals_all <- fst_list[[comp]]$WEIR_AND_COCKERHAM_FST
  vals <- if (args$keep_negatives) vals_all else vals_all[vals_all >= 0]
  tibble(
    comparison    = comp,
    comp_type     = comp_type(comp),
    n_snps        = length(vals),
    n_negative    = sum(vals_all < 0),
    frac_negative = mean(vals_all < 0),
    mean          = mean(vals),
    median        = median(vals),
    sd            = sd(vals),
    skewness      = skewness_fp(vals),
    q99           = as.numeric(quantile(vals, 0.99, na.rm = TRUE))
  )
}

summary_df <- bind_rows(lapply(all_comparisons, summarize_one)) %>%
  mutate(comp_type = factor(comp_type,
                            levels = c("CHI_vs_invasive", "Among_invasive"))) %>%
  arrange(comp_type, comparison)

write_tsv(summary_df, out_tsv, na = "NA")
log_msg(sprintf("  wrote summary table: %s", out_tsv))

# Pretty-print the summary so it lands in the run log too.
log_msg("\n  per-comparison summary:")
fmt_row <- function(r) {
  sprintf("    %-12s [%s]  n=%d  mean=%6.4f  median=%6.4f  sd=%6.4f  skew=%6.3f  q99=%6.4f  frac_neg=%5.3f",
          r$comparison, r$comp_type, r$n_snps,
          r$mean, r$median, r$sd, r$skewness, r$q99, r$frac_negative)
}
for (i in seq_len(nrow(summary_df))) {
  log_msg(fmt_row(summary_df[i, ]))
}

# --- Step 3: Per-comparison histograms --------------------------------------

log_msg("\n[3/3] Plotting histograms")

# Common x-axis breaks. F_ST is mathematically bounded above by 1; vcftools
# can return small negatives, so we extend the lower end to -0.2 to capture
# the noise floor. Bin width 0.01 gives ~120 bins, fine enough to see tail
# shape without aliasing.
X_MIN  <- -0.2
X_MAX  <- 1.0
BINW   <- 0.01

# Color the two comparison types so the diagnostic distinction is visible at
# a glance. CHI-vs-invasive in coral, among-invasive in steelblue, matching
# the convention in scripts/unified_selection_scan.R.
type_colors <- c(CHI_vs_invasive = "#E07B7B",
                 Among_invasive  = "#5B8AB7")

make_hist <- function(comp, log_y) {
  vals_all <- fst_list[[comp]]$WEIR_AND_COCKERHAM_FST
  vals <- if (args$keep_negatives) vals_all else vals_all[vals_all >= 0]
  df <- tibble(fst = vals)

  this_summary <- summary_df %>% filter(comparison == comp)
  fill_col <- type_colors[[as.character(this_summary$comp_type)]]

  subtitle <- sprintf("n=%s  mean=%.4f  median=%.4f  sd=%.4f  skew=%.2f  q99=%.4f  frac_neg=%.3f",
                      format(this_summary$n_snps, big.mark = ","),
                      this_summary$mean,
                      this_summary$median,
                      this_summary$sd,
                      this_summary$skewness,
                      this_summary$q99,
                      this_summary$frac_negative)

  p <- ggplot(df, aes(x = fst)) +
    geom_histogram(binwidth = BINW, boundary = 0,
                   fill = fill_col, color = "grey20", linewidth = 0.1) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    geom_vline(xintercept = this_summary$median,
               linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_vline(xintercept = this_summary$q99,
               linetype = "dashed", color = "darkred", linewidth = 0.4) +
    labs(
      title    = sprintf("%s -- per-SNP F_ST distribution", comp),
      subtitle = subtitle,
      x = expression(F[ST] ~ "(Weir & Cockerham, per SNP)"),
      y = if (log_y) "Count (log10 scale)" else "Count",
      caption = "Dotted = 0; dashed black = median; dashed red = 99th percentile."
    ) +
    coord_cartesian(xlim = c(X_MIN, X_MAX)) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, color = "grey20"),
          plot.caption  = element_text(size = 7, color = "grey40", hjust = 0))

  if (log_y) {
    p <- p + scale_y_continuous(trans = "log10",
                                labels = scales::label_number(big.mark = ","))
  } else {
    p <- p + scale_y_continuous(labels = scales::label_number(big.mark = ","))
  }
  p
}

# Save individual figures.
for (comp in all_comparisons) {
  for (log_y in c(FALSE, TRUE)) {
    p <- make_hist(comp, log_y)
    suffix <- if (log_y) "logy" else "linear"
    fpath  <- file.path(fig_dir, sprintf("%s_%s.png", comp, suffix))
    ggsave(fpath, p, width = 160, height = 100, units = "mm", dpi = 200)
  }
  log_msg(sprintf("  wrote: %s_{linear,logy}.png", comp))
}

# Combined 2x3 panel: CHI-axis comparisons on top row, among-invasive on
# bottom row, so the two distribution shapes can be compared at a glance.
build_combined <- function(log_y) {
  panels <- lapply(all_comparisons, function(comp) {
    p <- make_hist(comp, log_y)
    # Strip subtitle/caption from individual panels in the combined layout
    # to keep the per-panel area readable; the per-comparison figures keep
    # the full annotation for inspection.
    p + labs(title = comp, subtitle = NULL, caption = NULL) +
      theme(plot.title = element_text(size = 9, face = "bold"))
  })
  combined <- (panels[[1]] | panels[[2]] | panels[[3]]) /
              (panels[[4]] | panels[[5]] | panels[[6]]) +
    plot_annotation(
      title    = sprintf("Per-SNP F_ST distributions across 6 comparisons (%s)",
                         if (log_y) "log10 y" else "linear y"),
      subtitle = "Top row: CHI-vs-invasive. Bottom row: among-invasive.",
      caption  = sprintf("Bin width = %.3f. x-axis clipped to [%.1f, %.1f]. Dotted = 0; dashed black = median; dashed red = q99.",
                         BINW, X_MIN, X_MAX),
      theme = theme(plot.title = element_text(face = "bold", size = 11),
                    plot.caption = element_text(size = 7, color = "grey40",
                                                hjust = 0))
    )
  combined
}

for (log_y in c(FALSE, TRUE)) {
  p <- build_combined(log_y)
  suffix <- if (log_y) "logy" else "linear"
  fpath  <- file.path(fig_dir, sprintf("all_comparisons_%s.png", suffix))
  ggsave(fpath, p, width = 320, height = 180, units = "mm", dpi = 200)
  log_msg(sprintf("  wrote: %s", fpath))
}

# --- Run log ----------------------------------------------------------------

log_path <- file.path(project_root, "results/diagnostic/fst_histograms_run.log")
writeLines(log_lines, log_path)
log_msg(sprintf("\nRun log: %s", log_path))
