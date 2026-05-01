# title: outflank_diagnostic.R
# project: BIOL624 Final Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-01-09
# last modified: 2026-01-12
#
# purpose:
#   Run comprehensive OutFLANK diagnostic analysis on ALL population comparisons
#   to evaluate whether OutFLANK can be applied uniformly across high-FST and
#   low-FST comparisons, or if the dual-method approach is necessary.
#
#   Key diagnostic outputs:
#   1. dfInferred - degrees of freedom for chi-squared null distribution fit
#   2. FSTbar - mean neutral FST estimated by OutFLANK
#   3. FSTNoCorrbar - uncorrected mean FST
#   4. numberLowFSTloci - loci used to infer null distribution
#   5. numberHighFSTloci - potential outlier loci in right tail
#
#   Interpretation guide for dfInferred:
#   - df 2-10: Good fit, results reliable
#   - df 10-20: Acceptable, use with caution
#   - df >20: Poor fit, OutFLANK unreliable
#   - df >50 or Inf: Failed fit, use alternative method
#
#   UPDATE 2026-01-12: Now preserves CHROM:POS format for SNP IDs to enable
#   concordance validation with percentile outliers. Previous version used
#   generic SNP_# IDs which prevented cross-method validation.
#
# inputs:
#   - Per-chromosome VCF files in data/VARIANTS_BY_CHR/
#   - Population files in metadata/pop_*.txt
#
# outputs:
#   results/outflank_diagnostic/<comparison>/
#     - diagnostic_stats.tsv
#     - plot_fst_distribution.png/pdf
#     - plot_fst_vs_het.png/pdf
#     - plot_qvalue_distribution.png/pdf
#     - outliers_q05.tsv
#     - outliers_q01.tsv
#     - outflank_full_results.rds
#   results/outflank_diagnostic/diagnostic_summary_all.tsv
#
# usage example:
#   conda activate pop_gen
#   cd scripts
#   Rscript outflank_diagnostic.R
#
# copy/paste: Rscript outflank_diagnostic.R

################################################################################
# LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(vcfR)
  library(OutFLANK)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

################################################################################
# CONFIGURATION
################################################################################

# Paths (relative to scripts/ directory)
VCF_DIR <- '../data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS'
POP_DIR <- '../metadata'
OUTPUT_DIR <- '../results/outflank_diagnostic'

# OutFLANK parameters
LEFT_TRIM_FRACTION <- 0.05
RIGHT_TRIM_FRACTION <- 0.05
H_MIN <- 0.1
Q_THRESHOLD <- 0.05
MISSING_THRESHOLD <- 0.20
MIN_MAF <- 0.05

# Define all comparisons
COMPARISONS <- list(
  'CHI_vs_EEU' = c('CHI', 'EEU'),
  'CHI_vs_USA' = c('CHI', 'USA'),
  'CHI_vs_WEU' = c('CHI', 'WEU'),
  'EEU_vs_USA' = c('EEU', 'USA'),
  'EEU_vs_WEU' = c('EEU', 'WEU'),
  'USA_vs_WEU' = c('USA', 'WEU')
)

################################################################################
# SETUP
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('OUTFLANK COMPREHENSIVE DIAGNOSTIC ANALYSIS\n')
cat('===============================================================================\n')
cat('Start time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')
cat('===============================================================================\n\n')

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# UTILITY FUNCTIONS
################################################################################

# Function to read population file
read_population_file <- function(pop_name) {
  file_path <- file.path(POP_DIR, paste0('pop_', pop_name, '.txt'))
  if (!file.exists(file_path)) {
    stop(paste('Population file not found:', file_path))
  }
  samples <- readLines(file_path)
  samples <- samples[nchar(trimws(samples)) > 0]
  return(trimws(samples))
}

# Function to convert genotypes to biallelic format
convert_to_biallelic <- function(gt_matrix) {
  geno_mat <- matrix(NA, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))

  for (i in seq_len(nrow(gt_matrix))) {
    for (j in seq_len(ncol(gt_matrix))) {
      g <- gt_matrix[i, j]
      if (is.na(g) || g == './.') {
        geno_mat[i, j] <- NA
      } else {
        alleles <- as.numeric(unlist(strsplit(g, '[/|]')))
        if (length(alleles) == 2) {
          ref_count <- sum(alleles == 0)
          if (ref_count == 2) geno_mat[i, j] <- 0
          else if (ref_count == 1) geno_mat[i, j] <- 1
          else geno_mat[i, j] <- 2
        } else {
          geno_mat[i, j] <- NA
        }
      }
    }
  }

  return(geno_mat)
}

# Function to create diagnostic plots
create_diagnostic_plots <- function(outflank_obj, results_df, comparison_name, comp_dir) {
  cat('  Creating diagnostic plots...\n')

  # Get expected chi-squared distribution parameters
  df_inferred <- ifelse(is.null(outflank_obj$dfInferred), NA, outflank_obj$dfInferred)
  fst_bar <- ifelse(is.null(outflank_obj$FSTbar), NA, outflank_obj$FSTbar)

  # -------------------------------------------------------------------------
  # Plot 1: Observed vs Expected FST Distribution
  # -------------------------------------------------------------------------

  # Create histogram of observed FST
  fst_values <- results_df$FST[!is.na(results_df$FST) & is.finite(results_df$FST)]
  fst_positive <- fst_values[fst_values >= 0]

  if (length(fst_positive) > 0 && !is.na(df_inferred) && is.finite(df_inferred) && df_inferred > 0) {
    # Compute expected chi-squared distribution scaled to FST
    # OutFLANK models FST ~ (df * FSTbar / chi-sq(df))
    x_seq <- seq(0, max(fst_positive) * 1.1, length.out = 200)

    # The transformation from chi-squared to FST
    # FST = FSTbar * (chi-sq value) / df
    # So chi-sq = FST * df / FSTbar
    if (!is.na(fst_bar) && fst_bar > 0) {
      chi_sq_values <- x_seq * df_inferred / fst_bar
      expected_density <- dchisq(chi_sq_values, df = df_inferred) * df_inferred / fst_bar
      expected_df <- data.frame(x = x_seq, density = expected_density)
    } else {
      expected_df <- data.frame(x = numeric(0), density = numeric(0))
    }

    # Create plot
    p1 <- ggplot() +
      geom_histogram(data = data.frame(FST = fst_positive),
                    aes(x = FST, y = after_stat(density)),
                    bins = 80, fill = 'steelblue', alpha = 0.7, color = 'white') +
      labs(
        title = paste('FST Distribution:', comparison_name),
        subtitle = paste0('dfInferred = ', round(df_inferred, 2),
                         ' | FSTbar = ', round(fst_bar, 4)),
        x = 'FST',
        y = 'Density'
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = 'bold'),
        plot.subtitle = element_text(size = 11)
      )

    # Add expected distribution if available
    if (nrow(expected_df) > 0 && any(expected_df$density > 0)) {
      p1 <- p1 +
        geom_line(data = expected_df, aes(x = x, y = density),
                 color = 'red', linewidth = 1.2, linetype = 'solid') +
        annotate('text', x = max(fst_positive) * 0.7, y = max(expected_df$density) * 0.8,
                label = 'Expected (chi-sq)', color = 'red', size = 4)
    }

    # Add interpretation
    interpretation <- if (is.na(df_inferred) || !is.finite(df_inferred)) {
      'FAILED: Could not fit distribution'
    } else if (df_inferred <= 10) {
      'GOOD FIT: df 2-10 indicates reliable results'
    } else if (df_inferred <= 20) {
      'ACCEPTABLE: df 10-20, use with caution'
    } else if (df_inferred <= 50) {
      'POOR FIT: df 20-50, results unreliable'
    } else {
      'FAILED: df > 50, use alternative method'
    }

    p1 <- p1 + labs(caption = interpretation)

  } else {
    # Fallback plot without expected distribution
    p1 <- ggplot(data.frame(FST = fst_values), aes(x = FST)) +
      geom_histogram(bins = 80, fill = 'steelblue', alpha = 0.7) +
      labs(
        title = paste('FST Distribution:', comparison_name),
        subtitle = 'Warning: Could not compute expected distribution',
        x = 'FST',
        y = 'Count'
      ) +
      theme_minimal()
  }

  # Save Plot 1
  ggsave(file.path(comp_dir, 'plot_fst_distribution.png'), p1,
         width = 10, height = 7, dpi = 300)
  ggsave(file.path(comp_dir, 'plot_fst_distribution.pdf'), p1,
         width = 10, height = 7)

  # -------------------------------------------------------------------------
  # Plot 2: FST vs Heterozygosity
  # -------------------------------------------------------------------------

  if ('He' %in% colnames(results_df) && 'FST' %in% colnames(results_df)) {
    plot_data <- results_df %>%
      filter(!is.na(FST) & !is.na(He) & is.finite(FST) & is.finite(He))

    if (nrow(plot_data) > 0) {
      # Determine outlier status
      if ('OutlierFlag' %in% colnames(plot_data)) {
        plot_data$Outlier <- ifelse(plot_data$OutlierFlag == TRUE, 'Outlier', 'Non-outlier')
      } else {
        plot_data$Outlier <- 'Non-outlier'
      }

      p2 <- ggplot(plot_data, aes(x = He, y = FST, color = Outlier)) +
        geom_point(alpha = 0.5, size = 1) +
        scale_color_manual(values = c('Non-outlier' = 'gray60', 'Outlier' = 'red')) +
        labs(
          title = paste('FST vs Heterozygosity:', comparison_name),
          subtitle = 'Outliers should be distributed across heterozygosity range',
          x = 'Heterozygosity (He)',
          y = 'FST',
          color = ''
        ) +
        theme_minimal() +
        theme(legend.position = 'bottom')

      ggsave(file.path(comp_dir, 'plot_fst_vs_het.png'), p2,
             width = 10, height = 7, dpi = 300)
      ggsave(file.path(comp_dir, 'plot_fst_vs_het.pdf'), p2,
             width = 10, height = 7)
    }
  }

  # -------------------------------------------------------------------------
  # Plot 3: Q-value Distribution
  # -------------------------------------------------------------------------

  if ('qvalues' %in% colnames(results_df)) {
    qvals <- results_df$qvalues[!is.na(results_df$qvalues) & is.finite(results_df$qvalues)]

    if (length(qvals) > 0) {
      p3 <- ggplot(data.frame(qvalue = qvals), aes(x = qvalue)) +
        geom_histogram(bins = 50, fill = 'steelblue', alpha = 0.7, color = 'white') +
        geom_vline(xintercept = 0.05, color = 'red', linetype = 'dashed', linewidth = 1) +
        geom_vline(xintercept = 0.01, color = 'darkred', linetype = 'dashed', linewidth = 1) +
        annotate('text', x = 0.07, y = Inf, label = 'q=0.05', vjust = 2, color = 'red') +
        annotate('text', x = 0.03, y = Inf, label = 'q=0.01', vjust = 2, color = 'darkred') +
        labs(
          title = paste('Q-value Distribution:', comparison_name),
          subtitle = 'Peak near 1 = neutral; tail near 0 = outliers',
          x = 'Q-value (FDR-adjusted)',
          y = 'Count'
        ) +
        theme_minimal()

      ggsave(file.path(comp_dir, 'plot_qvalue_distribution.png'), p3,
             width = 10, height = 7, dpi = 300)
      ggsave(file.path(comp_dir, 'plot_qvalue_distribution.pdf'), p3,
             width = 10, height = 7)
    }
  }

  cat('    Plots saved to:', comp_dir, '\n')
}

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

run_diagnostic_for_comparison <- function(comparison_name, pop1_name, pop2_name,
                                          combined_geno, combined_variants, all_samples) {

  cat('\n===============================================================================\n')
  cat('COMPARISON:', comparison_name, '\n')
  cat('===============================================================================\n')

  # Create output directory
  comp_dir <- file.path(OUTPUT_DIR, comparison_name)
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)

  # Get population samples
  pop1_samples <- read_population_file(pop1_name)
  pop2_samples <- read_population_file(pop2_name)

  cat('Population 1 (', pop1_name, '):', length(pop1_samples), 'samples\n')
  cat('Population 2 (', pop2_name, '):', length(pop2_samples), 'samples\n')

  # Find indices
  pop1_idx <- which(all_samples %in% pop1_samples)
  pop2_idx <- which(all_samples %in% pop2_samples)

  if (length(pop1_idx) == 0 || length(pop2_idx) == 0) {
    cat('ERROR: Could not match samples to VCF\n')
    return(NULL)
  }

  cat('Matched:', length(pop1_idx), '+', length(pop2_idx), '=',
      length(pop1_idx) + length(pop2_idx), 'samples\n')

  # Subset genotypes
  pop_idx <- c(pop1_idx, pop2_idx)
  geno_subset <- combined_geno[, pop_idx]

  # Create population vector
  pop_vector <- factor(c(rep(pop1_name, length(pop1_idx)),
                         rep(pop2_name, length(pop2_idx))))

  cat('\nPreparing OutFLANK input...\n')
  cat('  SNPs:', nrow(geno_subset), '\n')
  cat('  Samples:', ncol(geno_subset), '\n')

  # Transpose for OutFLANK (samples as rows, SNPs as columns)
  geno_transposed <- t(geno_subset)

  # Create locus names using CHROM:POS format for concordance validation
  # This enables mapping OutFLANK outliers to genomic windows for comparison
  # with percentile-based outliers
  locus_names <- paste0(combined_variants$CHROM, ':', combined_variants$POS)
  colnames(geno_transposed) <- locus_names

  # Run OutFLANK
  cat('\nRunning OutFLANK...\n')

  outflank_result <- tryCatch({
    # Step 1: Create FST data frame
    cat('  Step 1: Computing FST matrix...\n')
    fst_data <- MakeDiploidFSTMat(
      SNPmat = geno_transposed,
      locusNames = locus_names,
      popNames = pop_vector
    )

    cat('  FST matrix created:', nrow(fst_data), 'SNPs\n')

    # Step 2: Filter by heterozygosity
    cat('  Step 2: Filtering by heterozygosity (He >=', H_MIN, ')...\n')
    valid_he <- !is.na(fst_data$He) & is.finite(fst_data$He) & fst_data$He >= H_MIN
    fst_filtered <- fst_data[valid_he, ]
    cat('  SNPs after He filter:', nrow(fst_filtered), '\n')

    if (nrow(fst_filtered) < 1000) {
      stop('Too few SNPs after filtering')
    }

    # Step 3: Run OutFLANK
    cat('  Step 3: Running OutFLANK with parameters:\n')
    cat('    LeftTrimFraction:', LEFT_TRIM_FRACTION, '\n')
    cat('    RightTrimFraction:', RIGHT_TRIM_FRACTION, '\n')
    cat('    Hmin:', H_MIN, '\n')
    cat('    qthreshold:', Q_THRESHOLD, '\n')

    outflank_obj <- OutFLANK(
      FstDataFrame = fst_filtered,
      LeftTrimFraction = LEFT_TRIM_FRACTION,
      RightTrimFraction = RIGHT_TRIM_FRACTION,
      Hmin = H_MIN,
      NumberOfSamples = 2,
      qthreshold = Q_THRESHOLD
    )

    cat('  OutFLANK completed successfully!\n')

    list(
      success = TRUE,
      outflank_obj = outflank_obj,
      fst_data = fst_filtered,
      results = outflank_obj$results
    )

  }, error = function(e) {
    cat('  ERROR:', e$message, '\n')
    list(
      success = FALSE,
      error = e$message
    )
  })

  # Process results
  if (!outflank_result$success) {
    # Create error diagnostic file
    error_df <- data.frame(
      Comparison = comparison_name,
      Status = 'FAILED',
      Error = outflank_result$error
    )
    write.table(error_df, file.path(comp_dir, 'diagnostic_stats.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)

    return(list(
      comparison = comparison_name,
      success = FALSE,
      error = outflank_result$error
    ))
  }

  # Extract diagnostic statistics
  outflank_obj <- outflank_result$outflank_obj
  results_df <- outflank_result$results

  # Key diagnostics
  df_inferred <- ifelse(is.null(outflank_obj$dfInferred), NA, outflank_obj$dfInferred)
  fst_bar <- ifelse(is.null(outflank_obj$FSTbar), NA, outflank_obj$FSTbar)
  fst_no_corr <- ifelse(is.null(outflank_obj$FSTNoCorrbar), NA, outflank_obj$FSTNoCorrbar)
  n_low_fst <- ifelse(is.null(outflank_obj$numberLowFSTloci), NA, outflank_obj$numberLowFSTloci)
  n_high_fst <- ifelse(is.null(outflank_obj$numberHighFSTloci), NA, outflank_obj$numberHighFSTloci)

  cat('\n--- DIAGNOSTIC STATISTICS ---\n')
  cat('dfInferred:', df_inferred, '\n')
  cat('FSTbar:', fst_bar, '\n')
  cat('FSTNoCorrbar:', fst_no_corr, '\n')
  cat('numberLowFSTloci:', n_low_fst, '\n')
  cat('numberHighFSTloci:', n_high_fst, '\n')

  # Determine fit quality
  fit_quality <- if (is.na(df_inferred) || !is.finite(df_inferred)) {
    'FAILED'
  } else if (df_inferred <= 10) {
    'GOOD'
  } else if (df_inferred <= 20) {
    'ACCEPTABLE'
  } else if (df_inferred <= 50) {
    'POOR'
  } else {
    'FAILED'
  }

  cat('Fit quality:', fit_quality, '\n')

  # Count outliers
  n_outliers_q05 <- sum(results_df$qvalues < 0.05, na.rm = TRUE)
  n_outliers_q01 <- sum(results_df$qvalues < 0.01, na.rm = TRUE)

  cat('\nOutliers at q < 0.05:', n_outliers_q05, '\n')
  cat('Outliers at q < 0.01:', n_outliers_q01, '\n')

  # Save diagnostic statistics
  diag_stats <- data.frame(
    Comparison = comparison_name,
    Pop1 = pop1_name,
    Pop2 = pop2_name,
    N_Pop1 = length(pop1_idx),
    N_Pop2 = length(pop2_idx),
    Total_SNPs = nrow(results_df),
    dfInferred = df_inferred,
    FSTbar = fst_bar,
    FSTNoCorrbar = fst_no_corr,
    numberLowFSTloci = n_low_fst,
    numberHighFSTloci = n_high_fst,
    Outliers_q05 = n_outliers_q05,
    Outliers_q01 = n_outliers_q01,
    Outlier_pct_q05 = round(100 * n_outliers_q05 / nrow(results_df), 4),
    Fit_Quality = fit_quality
  )

  write.table(diag_stats, file.path(comp_dir, 'diagnostic_stats.tsv'),
              sep = '\t', row.names = FALSE, quote = FALSE)

  # Save outlier lists
  if (n_outliers_q05 > 0) {
    outliers_q05 <- results_df %>%
      filter(qvalues < 0.05) %>%
      arrange(qvalues) %>%
      head(100)  # Top 100
    write.table(outliers_q05, file.path(comp_dir, 'outliers_q05.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
  }

  if (n_outliers_q01 > 0) {
    outliers_q01 <- results_df %>%
      filter(qvalues < 0.01) %>%
      arrange(qvalues)
    write.table(outliers_q01, file.path(comp_dir, 'outliers_q01.tsv'),
                sep = '\t', row.names = FALSE, quote = FALSE)
  }

  # Save full results as RDS
  saveRDS(list(
    outflank_obj = outflank_obj,
    results = results_df,
    diagnostics = diag_stats
  ), file.path(comp_dir, 'outflank_full_results.rds'))

  # Create diagnostic plots
  create_diagnostic_plots(outflank_obj, results_df, comparison_name, comp_dir)

  cat('\nResults saved to:', comp_dir, '\n')

  return(diag_stats)
}

################################################################################
# LOAD AND PROCESS VCF DATA
################################################################################

cat('Loading VCF data...\n')

# Find VCF files
vcf_files <- list.files(VCF_DIR, pattern = '\\.vcf\\.gz$', full.names = TRUE)
cat('Found', length(vcf_files), 'VCF files\n')

if (length(vcf_files) == 0) {
  stop('No VCF files found in:', VCF_DIR)
}

# Process each chromosome
all_geno_list <- list()
all_var_list <- list()
sample_names <- NULL

for (vcf_file in vcf_files) {
  chrom <- gsub('\\.vcf\\.gz$', '', basename(vcf_file))
  cat('\n  Processing', chrom, '...')

  # Read VCF
  vcf <- read.vcfR(vcf_file, verbose = FALSE)

  # Get samples
  if (is.null(sample_names)) {
    sample_names <- colnames(vcf@gt)[-1]
  }

  # Extract genotypes
  gt <- extract.gt(vcf)

  # Convert to biallelic
  geno_mat <- convert_to_biallelic(gt)

  # Filter by missing data
  missing_rate <- rowMeans(is.na(geno_mat))
  keep <- missing_rate <= MISSING_THRESHOLD
  geno_mat <- geno_mat[keep, ]

  # Impute remaining missing
  for (i in seq_len(nrow(geno_mat))) {
    if (any(is.na(geno_mat[i, ]))) {
      snp_mean <- mean(geno_mat[i, ], na.rm = TRUE)
      if (!is.na(snp_mean)) {
        geno_mat[i, is.na(geno_mat[i, ])] <- round(max(0, min(2, snp_mean)))
      }
    }
  }

  # Store
  all_geno_list[[chrom]] <- geno_mat
  all_var_list[[chrom]] <- data.frame(
    CHROM = getCHROM(vcf)[keep],
    POS = getPOS(vcf)[keep]
  )

  cat(' retained', nrow(geno_mat), 'SNPs')
}

# Combine all chromosomes
cat('\n\nCombining all chromosomes...\n')
combined_geno <- do.call(rbind, all_geno_list)
combined_variants <- do.call(rbind, all_var_list)

cat('Total SNPs:', nrow(combined_geno), '\n')
cat('Total samples:', ncol(combined_geno), '\n')

################################################################################
# RUN DIAGNOSTICS FOR ALL COMPARISONS
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('RUNNING OUTFLANK DIAGNOSTICS FOR ALL COMPARISONS\n')
cat('===============================================================================\n')

all_diagnostics <- list()

for (comp_name in names(COMPARISONS)) {
  pops <- COMPARISONS[[comp_name]]

  result <- run_diagnostic_for_comparison(
    comparison_name = comp_name,
    pop1_name = pops[1],
    pop2_name = pops[2],
    combined_geno = combined_geno,
    combined_variants = combined_variants,
    all_samples = sample_names
  )

  if (!is.null(result)) {
    all_diagnostics[[comp_name]] <- result
  }
}

################################################################################
# CREATE COMBINED SUMMARY
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('GENERATING COMBINED SUMMARY\n')
cat('===============================================================================\n')

# Combine all diagnostic stats
summary_df <- do.call(rbind, all_diagnostics)
rownames(summary_df) <- NULL

# Add recommendation column
summary_df$Recommendation <- ifelse(
  summary_df$Fit_Quality == 'GOOD', 'Use OutFLANK',
  ifelse(summary_df$Fit_Quality == 'ACCEPTABLE', 'Use OutFLANK with caution',
         'Use percentile method')
)

# Save combined summary
summary_file <- file.path(OUTPUT_DIR, 'diagnostic_summary_all.tsv')
write.table(summary_df, summary_file, sep = '\t', row.names = FALSE, quote = FALSE)
cat('Combined summary saved to:', summary_file, '\n')

# Print summary table
cat('\n--- DIAGNOSTIC SUMMARY ---\n\n')
print(summary_df[, c('Comparison', 'dfInferred', 'FSTbar', 'Outliers_q05', 'Fit_Quality', 'Recommendation')])

################################################################################
# FINAL REPORT
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('ANALYSIS COMPLETE\n')
cat('===============================================================================\n')
cat('End time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

# Summary of fit quality
good_fits <- sum(summary_df$Fit_Quality == 'GOOD')
acceptable_fits <- sum(summary_df$Fit_Quality == 'ACCEPTABLE')
poor_fits <- sum(summary_df$Fit_Quality == 'POOR')
failed_fits <- sum(summary_df$Fit_Quality == 'FAILED')

cat('FIT QUALITY SUMMARY:\n')
cat('  GOOD (df 2-10):', good_fits, 'comparisons\n')
cat('  ACCEPTABLE (df 10-20):', acceptable_fits, 'comparisons\n')
cat('  POOR (df 20-50):', poor_fits, 'comparisons\n')
cat('  FAILED (df > 50 or error):', failed_fits, 'comparisons\n')

cat('\nOUTPUT DIRECTORIES:\n')
for (comp_name in names(COMPARISONS)) {
  cat(' ', file.path(OUTPUT_DIR, comp_name), '\n')
}
cat('\nSUMMARY FILE:', summary_file, '\n')

# Overall recommendation
cat('\n--- OVERALL RECOMMENDATION ---\n')
if (good_fits == length(COMPARISONS)) {
  cat('OutFLANK can be applied uniformly to all comparisons.\n')
} else if (good_fits + acceptable_fits >= 4) {
  cat('OutFLANK works well for most comparisons.\n')
  cat('Consider dual-method approach for comparisons with POOR/FAILED fits.\n')
} else {
  cat('Dual-method approach recommended.\n')
  cat('Use OutFLANK only for comparisons with GOOD/ACCEPTABLE fits.\n')
}

cat('\n===============================================================================\n')
