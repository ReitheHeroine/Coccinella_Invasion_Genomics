#!/usr/bin/env Rscript

# title: diagnostic_genomic_inflation.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-01
# last modified: 2026-05-01
#
# purpose:
#   Manhattan-scatter diagnostic, Task 2. Computes the genomic inflation
#   factor lambda for each of the 6 pairwise per-SNP F_ST comparisons
#   following the Devlin & Roeder convention adapted to F_ST. Per the
#   project spec: filter to F_ST >= 0 (vcftools W&C can return negatives
#   when within-pop variance > between-pop variance), standardize the
#   surviving F_ST values to a z distribution, convert |z| to a one-df
#   chi-squared statistic via the two-sided p-value identity
#   (qchisq(1 - 2*pnorm(-|z|), df=1) is mathematically z^2 under df=1),
#   and compute lambda = median(chi^2) / qchisq(0.5, 1).
#
#   Reading the output: lambda close to 1 means the *standardized* F_ST
#   distribution's median sits where a chi-squared(1) would predict.
#   Lambda well above 1 means the bulk of the distribution is heavier-
#   tailed than chi-squared after standardization. Note that this is the
#   shape of the standardized distribution, not raw inflation: a F_ST
#   distribution centered higher (e.g. CHI vs invasive at mean ~0.07)
#   but with the same standardized shape as a comparison centered near
#   zero will give a similar lambda. The standardization removes the
#   raw mean shift. So lambda here detects *tail-shape deviation from
#   chi-squared* rather than mean-level inflation per se. See the
#   discussion at the end of the run log for caveats and an alternative
#   formulation that does pick up raw inflation.
#
# inputs:
#   - results/fst/{COMP}/fst_genomewide/{CHROM}.filtered.weir.fst
#       Per-SNP Weir & Cockerham F_ST. Same source as Task 1.
#       Sex chromosome (NC_058198.1) excluded; 9 autosomes used.
#
# outputs:
#   - results/diagnostic/genomic_inflation_summary.tsv
#       Per-comparison row with: comparison, comp_type, n_total,
#       n_negative_dropped, n_used, mean_pos, sd_pos, median_chi2,
#       lambda, interpretation.
#   - results/diagnostic/genomic_inflation_summary.txt
#       Human-readable version with comparison columns aligned for the
#       handoff doc, plus an interpretation block at the end.
#
# usage example:
#   Rscript scripts/diagnostic_genomic_inflation.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# --- Paths ------------------------------------------------------------------

project_root <- if (file.exists("results/fst")) {
  "."
} else if (file.exists("../results/fst")) {
  ".."
} else {
  stop("Cannot locate results/fst/. Run from project root or scripts/.")
}

fst_base   <- file.path(project_root, "results/fst")
out_dir    <- file.path(project_root, "results/diagnostic")
out_tsv    <- file.path(out_dir, "genomic_inflation_summary.tsv")
out_txt    <- file.path(out_dir, "genomic_inflation_summary.txt")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Configuration ----------------------------------------------------------

chi_comparisons      <- c("CHI_vs_EEU", "CHI_vs_USA", "CHI_vs_WEU")
invasive_comparisons <- c("EEU_vs_USA", "EEU_vs_WEU", "USA_vs_WEU")
all_comparisons      <- c(chi_comparisons, invasive_comparisons)
autosomes <- sprintf("NC_0581%02d.1", 89:97)

CHISQ_DF1_MEDIAN <- qchisq(0.5, df = 1)   # ~0.4549

# --- Logging ----------------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

log_msg("=== Task 2: Genomic inflation factor (lambda) per comparison ===")
log_msg(sprintf("project root: %s", normalizePath(project_root)))
log_msg(sprintf("chi-squared(df=1) median: %.4f", CHISQ_DF1_MEDIAN))
log_msg("")

# --- Step 1: Load per-SNP F_ST ---------------------------------------------

read_fst_one <- function(comp, chrom) {
  f <- file.path(fst_base, comp, "fst_genomewide",
                 sprintf("%s.filtered.weir.fst", chrom))
  if (!file.exists(f)) stop(sprintf("missing per-SNP F_ST file: %s", f))
  dt <- fread(f, sep = "\t", header = TRUE,
              colClasses = c(CHROM = "character",
                             POS = "integer",
                             WEIR_AND_COCKERHAM_FST = "double"))
  dt <- dt[!is.na(WEIR_AND_COCKERHAM_FST)]
  dt
}

comp_type_of <- function(comp) {
  if (comp %in% chi_comparisons) "CHI_vs_invasive" else "Among_invasive"
}

# --- Step 2: Per-comparison lambda -----------------------------------------

# Implements the spec's recipe:
#   1. Drop F_ST < 0 (vcftools noise floor).
#   2. z = scale(F_ST_pos)   -- standardize.
#   3. chi^2 = qchisq(1 - 2*pnorm(-|z|), df=1)
#      [mathematically equivalent to z^2 for df=1; written this way to
#      stay one-to-one with the test-statistic-from-p-value convention.]
#   4. lambda = median(chi^2) / qchisq(0.5, df=1)
compute_lambda_one <- function(comp) {
  dt <- rbindlist(lapply(autosomes, read_fst_one, comp = comp))
  fst_all <- dt$WEIR_AND_COCKERHAM_FST
  n_total <- length(fst_all)

  fst_pos <- fst_all[fst_all >= 0]
  n_neg   <- n_total - length(fst_pos)

  z       <- as.numeric(scale(fst_pos))      # mean 0, sd 1
  pvals_2 <- 2 * pnorm(-abs(z))              # 2-sided p
  # Floor extreme p-values to avoid p == 0 producing chi^2 = Inf at the
  # very tails. machine eps ~2.2e-16; clamp at that scale.
  pvals_2 <- pmax(pvals_2, .Machine$double.eps)
  chi_sq  <- qchisq(1 - pvals_2, df = 1)

  lambda  <- median(chi_sq, na.rm = TRUE) / CHISQ_DF1_MEDIAN

  interp <- if (is.na(lambda)) "NA" else
            if (lambda < 1.05) "well_calibrated" else
            if (lambda < 1.5)  "mild_inflation"  else
                                "severe_inflation"

  tibble(
    comparison         = comp,
    comp_type          = comp_type_of(comp),
    n_total            = n_total,
    n_negative_dropped = n_neg,
    n_used             = length(fst_pos),
    mean_pos_fst       = mean(fst_pos),
    sd_pos_fst         = sd(fst_pos),
    median_chi2        = median(chi_sq, na.rm = TRUE),
    lambda             = lambda,
    interpretation     = interp
  )
}

log_msg("[1/2] Computing lambda per comparison")
summary_rows <- vector("list", length(all_comparisons))
for (i in seq_along(all_comparisons)) {
  comp <- all_comparisons[i]
  log_msg(sprintf("  loading %s...", comp))
  summary_rows[[i]] <- compute_lambda_one(comp)
}
summary_df <- bind_rows(summary_rows) %>%
  mutate(comp_type = factor(comp_type,
                            levels = c("CHI_vs_invasive", "Among_invasive"))) %>%
  arrange(comp_type, comparison)

write_tsv(summary_df, out_tsv, na = "NA")
log_msg(sprintf("\n  wrote summary table: %s", out_tsv))

# --- Step 3: Pretty text output --------------------------------------------

log_msg("\n[2/2] Writing human-readable summary")

txt_lines <- character()
add <- function(s) txt_lines <<- c(txt_lines, s)

add("Genomic inflation factor (lambda) per pairwise comparison")
add("==========================================================")
add("")
add(sprintf("Method: Devlin & Roeder, applied to per-SNP Weir & Cockerham F_ST."))
add(sprintf("  1. Drop F_ST < 0 (vcftools noise-floor convention)."))
add(sprintf("  2. Standardize remaining F_ST (z = scale(F_ST))."))
add(sprintf("  3. chi^2 = qchisq(1 - 2*pnorm(-|z|), df=1)   [equivalent to z^2 for df=1]."))
add(sprintf("  4. lambda = median(chi^2) / qchisq(0.5, df=1) = median(chi^2) / %.4f.",
            CHISQ_DF1_MEDIAN))
add("")
add("Interpretation thresholds (per spec):")
add("  lambda < 1.05  : well_calibrated     (no inflation)")
add("  1.05 - 1.50    : mild_inflation")
add("  > 1.50         : severe_inflation")
add("")
add(sprintf("%-12s  %-16s  %10s  %10s  %10s  %8s  %8s  %10s  %8s  %s",
            "comparison", "comp_type", "n_total", "n_neg_drop", "n_used",
            "mean_pos", "sd_pos", "median_chi2", "lambda", "interp"))
add(strrep("-", 132))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  add(sprintf("%-12s  %-16s  %10d  %10d  %10d  %8.4f  %8.4f  %10.4f  %8.4f  %s",
              r$comparison, as.character(r$comp_type),
              r$n_total, r$n_negative_dropped, r$n_used,
              r$mean_pos_fst, r$sd_pos_fst,
              r$median_chi2, r$lambda, r$interpretation))
}
add("")

# Group means.
chi_lambda_mean <- mean(summary_df$lambda[summary_df$comp_type == "CHI_vs_invasive"])
inv_lambda_mean <- mean(summary_df$lambda[summary_df$comp_type == "Among_invasive"])
add(sprintf("Mean lambda, CHI-vs-invasive comparisons: %.4f", chi_lambda_mean))
add(sprintf("Mean lambda, among-invasive comparisons:  %.4f", inv_lambda_mean))
add("")

add("Methodological caveat (read before interpreting):")
add("  This is lambda on STANDARDIZED F_ST. Standardization (subtracting the")
add("  mean and dividing by SD) removes the raw mean shift. So a comparison")
add("  with a much higher genome-wide F_ST (e.g. CHI-vs-invasive, mean ~0.07)")
add("  is not penalized by lambda for its higher mean -- it is penalized only")
add("  if its STANDARDIZED tail is heavier than chi-squared(1). The user's")
add("  framing 'lambda > 1.5 means structure is inflating everything' is")
add("  therefore not a clean test of mean-level structure inflation under this")
add("  formulation; it is a test of tail-shape deviation from a chi-squared(1)")
add("  null after z-scoring.")
add("")
add("  An alternative formulation that does pick up raw-FST inflation is to")
add("  convert F_ST to chi-squared via the relationship for a 2-pop allele-")
add("  count contrast: chi^2 ~= (n1+n2) * 2 * F_ST  (Workman & Niswander 1970,")
add("  approximately). Under this scaling, lambda would diagnose whether the")
add("  median F_ST exceeds what neutral drift among samples of size n1, n2")
add("  would predict, which is a more direct structure-leakage test. If the")
add("  current lambda values do not match the spec's interpretive framework")
add("  cleanly, the alternative formulation is a one-line modification.")

writeLines(txt_lines, out_txt)
log_msg(sprintf("  wrote: %s", out_txt))

cat("\n", paste(txt_lines, collapse = "\n"), "\n", sep = "")
