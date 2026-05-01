# title: outflank_ldpruned.R
# project: BIOL624 Final Project -- C. septempunctata population genomics
# author: Reina Hastings <reinahastings13@gmail.com>
# date created: 2026-04-23
# last modified: 2026-04-23
#
# purpose:
#   Re-run OutFLANK on an LD-pruned SNP set for all 6 pairwise population
#   comparisons (CHI/EEU/WEU/USA), per the method's independence requirement
#   (Whitlock & Lotterhos 2015). Writes per-comparison diagnostics, outlier
#   tables, and full-result RDS files; plus a combined summary across all
#   comparisons. This is Task 4.7 Step 7b in the project handoff. Outputs
#   feed the cross-method concordance rebuild in Step 7c.
#
#   Design note: differs from scripts/outflank_diagnostic.R only at the input
#   layer (single merged VCF of LD-pruned autosomal SNPs, produced by
#   scripts/ld_prune.sh, instead of per-chromosome filtered VCFs). The
#   OutFLANK call, He filter, and diagnostic outputs match the original so
#   results are directly comparable. Plots are intentionally skipped here;
#   outflank_full_results.rds preserves enough state to regenerate them
#   later via the helpers in outflank_diagnostic.R.
#
# inputs:
#   - data/plink/ladybug_snps_ldpruned.vcf.gz (bgzipped, tabix-indexed;
#       454,840 autosomal SNPs x 88 samples; NC_058198.1 excluded at pruning)
#   - metadata/pop_{CHI,EEU,WEU,USA}.txt (one sample ID per line)
#
# outputs:
#   results/outflank_ldpruned/<comparison>/
#     - diagnostic_stats.tsv       (dfInferred, FSTbar, outlier counts)
#     - outliers_q05.tsv           (top-100 SNPs at qvalue < 0.05)
#     - outliers_q01.tsv           (all SNPs at qvalue < 0.01)
#     - outflank_full_results.rds  (outflank object + FST data + diagnostics)
#   results/outflank_ldpruned/diagnostic_summary_all.tsv
#
# usage example:
#   conda activate pop_gen
#   cd scripts
#   Rscript outflank_ldpruned.R
#
# copy/paste: Rscript outflank_ldpruned.R

################################################################################
# LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(vcfR)
  library(OutFLANK)
  library(dplyr)
  library(tidyr)
  library(readr)
})

################################################################################
# CONFIGURATION
################################################################################

# Paths (relative to scripts/ directory)
VCF_PATH   <- "../data/plink/ladybug_snps_ldpruned.vcf.gz"
POP_DIR    <- "../metadata"
OUTPUT_DIR <- "../results/outflank_ldpruned"

# OutFLANK parameters (match outflank_diagnostic.R for apples-to-apples
# comparison with the Jan 2026 run; only the input SNP set has changed)
LEFT_TRIM_FRACTION  <- 0.05
RIGHT_TRIM_FRACTION <- 0.05
H_MIN               <- 0.10
Q_THRESHOLD         <- 0.05
MISSING_THRESHOLD   <- 0.20   # per-SNP max missing rate
# Pairwise comparisons (same 6 as the original diagnostic run)
COMPARISONS <- list(
  "CHI_vs_EEU" = c("CHI", "EEU"),
  "CHI_vs_USA" = c("CHI", "USA"),
  "CHI_vs_WEU" = c("CHI", "WEU"),
  "EEU_vs_USA" = c("EEU", "USA"),
  "EEU_vs_WEU" = c("EEU", "WEU"),
  "USA_vs_WEU" = c("USA", "WEU")
)

################################################################################
# SETUP
################################################################################

cat("\n===============================================================================\n")
cat("OUTFLANK RE-RUN ON LD-PRUNED DATA (Task 4.7, Step 7b)\n")
cat("===============================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("VCF:       ", VCF_PATH, "\n")
cat("Out dir:   ", OUTPUT_DIR, "\n")
cat("===============================================================================\n\n")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# UTILITIES
################################################################################

# --- Read a population file: one sample ID per line, skip blank lines ---
read_population_file <- function(pop_name) {
  f <- file.path(POP_DIR, paste0("pop_", pop_name, ".txt"))
  if (!file.exists(f)) stop("Population file not found: ", f)
  s <- trimws(readLines(f))
  s[nchar(s) > 0]
}

# --- Vectorized GT -> 0/1/2/NA matrix. Handles phased | and unphased /.
# Non-biallelic genotypes (e.g., "1/2") stay NA. ~100x faster than the
# nested loop in outflank_diagnostic.R; this matters on 454K x 88 cells.
convert_to_biallelic <- function(gt_matrix) {
  result <- matrix(NA_integer_,
                   nrow = nrow(gt_matrix), ncol = ncol(gt_matrix),
                   dimnames = dimnames(gt_matrix))
  result[gt_matrix %in% c("0/0", "0|0")] <- 0L
  result[gt_matrix %in% c("0/1", "0|1", "1/0", "1|0")] <- 1L
  result[gt_matrix %in% c("1/1", "1|1")] <- 2L
  result
}

# --- Mean-imputation of missing calls, rounded to {0, 1, 2}. Matches the
# behaviour of outflank_diagnostic.R so the two runs are comparable.
impute_rowwise <- function(geno_mat) {
  missing_rows <- which(rowSums(is.na(geno_mat)) > 0)
  for (i in missing_rows) {
    m <- mean(geno_mat[i, ], na.rm = TRUE)
    if (is.finite(m)) {
      geno_mat[i, is.na(geno_mat[i, ])] <- round(pmin(2, pmax(0, m)))
    }
  }
  geno_mat
}

################################################################################
# LOAD VCF AND BUILD GENOTYPE MATRIX
################################################################################

cat("Loading VCF...\n")
t0 <- Sys.time()
vcf <- read.vcfR(VCF_PATH, verbose = FALSE)
cat("  Loaded in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

sample_names <- colnames(vcf@gt)[-1]
cat("  Samples:", length(sample_names), "\n")

cat("  Extracting genotypes...\n")
gt <- extract.gt(vcf)
geno <- convert_to_biallelic(gt)
variants <- data.frame(CHROM = getCHROM(vcf), POS = getPOS(vcf),
                       stringsAsFactors = FALSE)

cat("  Raw variant count:", nrow(geno), "\n")

# --- Filter SNPs with too much missingness, then mean-impute the rest ---
miss_rate <- rowMeans(is.na(geno))
keep <- miss_rate <= MISSING_THRESHOLD
cat("  SNPs dropped for >", 100 * MISSING_THRESHOLD, "% missing: ",
    sum(!keep), " (", round(100 * mean(!keep), 2), "%)\n", sep = "")

geno     <- geno[keep, , drop = FALSE]
variants <- variants[keep, , drop = FALSE]

cat("  Imputing remaining missing calls (per-SNP mean, rounded)...\n")
geno <- impute_rowwise(geno)
cat("  Ready: ", nrow(geno), " SNPs x ", ncol(geno), " samples\n\n", sep = "")

# Free the raw VCF now that we have what we need -- saves ~500 MB on large sets
rm(vcf, gt); invisible(gc(verbose = FALSE))

################################################################################
# PER-COMPARISON OUTFLANK
################################################################################

run_comparison <- function(comp_name, pop1, pop2) {

  cat("-------------------------------------------------------------------------------\n")
  cat("COMPARISON:", comp_name, "\n")
  cat("-------------------------------------------------------------------------------\n")

  comp_dir <- file.path(OUTPUT_DIR, comp_name)
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Resolve sample indices for the two pops --------------------------------
  s1 <- read_population_file(pop1)
  s2 <- read_population_file(pop2)
  # Defense: fail loud if pop-file samples are missing from the VCF. Silent
  # truncation here would shrink n (esp. for CHI, n=5) without any warning
  # in the diagnostic_stats.tsv output.
  missing1 <- setdiff(s1, sample_names)
  missing2 <- setdiff(s2, sample_names)
  if (length(missing1) > 0 || length(missing2) > 0) {
    stop(sprintf(
      "Pop-file samples absent from VCF. %s missing: %s. %s missing: %s.",
      pop1, ifelse(length(missing1), paste(missing1, collapse = ","), "(none)"),
      pop2, ifelse(length(missing2), paste(missing2, collapse = ","), "(none)")
    ))
  }
  i1 <- which(sample_names %in% s1)
  i2 <- which(sample_names %in% s2)
  if (length(i1) == 0 || length(i2) == 0) {
    cat("  ERROR: pop sample match failed (pop1=", length(i1),
        ", pop2=", length(i2), ")\n", sep = "")
    return(NULL)
  }
  cat("  ", pop1, ": ", length(i1), "   ", pop2, ": ", length(i2), "\n", sep = "")

  geno_sub <- geno[, c(i1, i2), drop = FALSE]
  pop_vec  <- factor(c(rep(pop1, length(i1)), rep(pop2, length(i2))))
  loci     <- paste0(variants$CHROM, ":", variants$POS)

  # --- MakeDiploidFSTMat expects samples-as-rows ------------------------------
  cat("  MakeDiploidFSTMat...\n")
  fst_data <- tryCatch(
    MakeDiploidFSTMat(SNPmat = t(geno_sub), locusNames = loci,
                      popNames = pop_vec),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fst_data)) return(NULL)
  cat("    FST matrix rows:", nrow(fst_data), "\n")

  # --- He filter (same threshold as the original diagnostic run) --------------
  fst_filtered <- fst_data[!is.na(fst_data$He) & is.finite(fst_data$He) &
                             fst_data$He >= H_MIN, , drop = FALSE]
  cat("    After He >= ", H_MIN, ": ", nrow(fst_filtered), "\n", sep = "")
  if (nrow(fst_filtered) < 1000) {
    cat("  SKIPPED: <1000 SNPs after filter\n")
    return(NULL)
  }

  # --- OutFLANK ----------------------------------------------------------------
  cat("  Running OutFLANK...\n")
  of_obj <- tryCatch(
    OutFLANK(FstDataFrame = fst_filtered,
             LeftTrimFraction  = LEFT_TRIM_FRACTION,
             RightTrimFraction = RIGHT_TRIM_FRACTION,
             Hmin              = H_MIN,
             NumberOfSamples   = 2,
             qthreshold        = Q_THRESHOLD),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(of_obj)) return(NULL)

  results_df <- of_obj$results
  df_inf     <- if (is.null(of_obj$dfInferred))   NA else of_obj$dfInferred
  fst_bar    <- if (is.null(of_obj$FSTbar))       NA else of_obj$FSTbar
  fst_no_c   <- if (is.null(of_obj$FSTNoCorrbar)) NA else of_obj$FSTNoCorrbar
  n_low      <- if (is.null(of_obj$numberLowFSTloci))  NA else of_obj$numberLowFSTloci
  n_high     <- if (is.null(of_obj$numberHighFSTloci)) NA else of_obj$numberHighFSTloci

  fit_quality <- if (is.na(df_inf) || !is.finite(df_inf)) "FAILED"
                 else if (df_inf <= 10) "GOOD"
                 else if (df_inf <= 20) "ACCEPTABLE"
                 else if (df_inf <= 50) "POOR"
                 else "FAILED"

  n_q05 <- sum(results_df$qvalues < 0.05, na.rm = TRUE)
  n_q01 <- sum(results_df$qvalues < 0.01, na.rm = TRUE)

  cat(sprintf("    dfInferred=%.2f  FSTbar=%.4f  outliers(q<.05)=%d  (%s)\n",
              df_inf, fst_bar, n_q05, fit_quality))

  # --- Write diagnostics -------------------------------------------------------
  diag_stats <- data.frame(
    Comparison       = comp_name,
    Pop1             = pop1,   Pop2 = pop2,
    N_Pop1           = length(i1),  N_Pop2 = length(i2),
    Total_SNPs       = nrow(results_df),
    dfInferred       = df_inf,
    FSTbar           = fst_bar,
    FSTNoCorrbar     = fst_no_c,
    numberLowFSTloci = n_low,
    numberHighFSTloci = n_high,
    Outliers_q05     = n_q05,
    Outliers_q01     = n_q01,
    Outlier_pct_q05  = round(100 * n_q05 / nrow(results_df), 4),
    Fit_Quality      = fit_quality
  )
  write_tsv(diag_stats, file.path(comp_dir, "diagnostic_stats.tsv"))

  if (n_q05 > 0) {
    outl <- results_df %>% filter(qvalues < 0.05) %>% arrange(qvalues) %>% head(100)
    write_tsv(outl, file.path(comp_dir, "outliers_q05.tsv"))
  }
  if (n_q01 > 0) {
    outl <- results_df %>% filter(qvalues < 0.01) %>% arrange(qvalues)
    write_tsv(outl, file.path(comp_dir, "outliers_q01.tsv"))
  }

  saveRDS(
    list(outflank_obj = of_obj, results = results_df, diagnostics = diag_stats),
    file.path(comp_dir, "outflank_full_results.rds")
  )

  diag_stats
}

################################################################################
# MAIN LOOP
################################################################################

all_stats <- list()
for (comp_name in names(COMPARISONS)) {
  pops <- COMPARISONS[[comp_name]]
  res  <- run_comparison(comp_name, pops[1], pops[2])
  if (!is.null(res)) all_stats[[comp_name]] <- res
}

################################################################################
# COMBINED SUMMARY
################################################################################

cat("\n===============================================================================\n")
cat("COMBINED SUMMARY\n")
cat("===============================================================================\n")

summary_df <- bind_rows(all_stats)
summary_df$Recommendation <- case_when(
  summary_df$Fit_Quality == "GOOD"       ~ "Use OutFLANK",
  summary_df$Fit_Quality == "ACCEPTABLE" ~ "Use OutFLANK with caution",
  TRUE                                    ~ "Use percentile method"
)

summary_file <- file.path(OUTPUT_DIR, "diagnostic_summary_all.tsv")
write_tsv(summary_df, summary_file)

print(summary_df[, c("Comparison", "dfInferred", "FSTbar",
                     "Outliers_q05", "Fit_Quality")])

cat("\nSummary written to:", summary_file, "\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("===============================================================================\n")
