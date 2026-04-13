# title: calculate_outflank_concordance.R
# project: BIOL624 Final Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-01-12
# last modified: 2026-01-12
#
# purpose:
#   Calculate concordance between OutFLANK outliers and percentile-based outliers.
#   Maps OutFLANK SNP-level outliers (CHROM:POS format) to 5kb genomic windows
#   and calculates overlap with percentile outlier windows.
#
#   This script validates the unified percentile approach by showing agreement
#   between the two methods where OutFLANK performs well (among-invasive comparisons).
#
# inputs:
#   - results/outflank_diagnostic/<comparison>/outflank_full_results.rds
#   - results/unified_selection_scan/percentile_outliers_*.tsv
#
# outputs:
#   results/outflank_diagnostic/concordance/
#     - concordance_<comparison>.tsv
#     - concordance_summary.txt
#     - concordance_plot.png
#
# usage example:
#   conda activate pop_gen
#   cd scripts
#   Rscript calculate_outflank_concordance.R
#
# copy/paste: Rscript calculate_outflank_concordance.R

################################################################################
# LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

################################################################################
# CONFIGURATION
################################################################################

# Paths (relative to scripts/ directory)
OUTFLANK_DIR <- '../results/outflank_diagnostic'
PERCENTILE_DIR <- '../results/unified_selection_scan'
OUTPUT_DIR <- '../results/outflank_diagnostic/concordance'

# Window size for mapping SNPs to genomic windows (must match percentile analysis)
WINDOW_SIZE <- 5000

# Comparisons with OutFLANK outliers (from diagnostic analysis)
VALID_COMPARISONS <- c('EEU_vs_WEU', 'EEU_vs_USA')

# Q-value threshold for OutFLANK outliers
Q_THRESHOLD <- 0.05

# Percentile thresholds to test
PERCENTILE_THRESHOLDS <- c('0.5pct', '1pct', '2pct')

################################################################################
# SETUP
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('OUTFLANK-PERCENTILE CONCORDANCE CALCULATION\n')
cat('===============================================================================\n')
cat('Start time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')
cat('===============================================================================\n\n')

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# UTILITY FUNCTIONS
################################################################################

# Function to map SNP position to 5kb window
map_snp_to_window <- function(chrom, pos, window_size = 5000) {
  bin_start <- floor((pos - 1) / window_size) * window_size + 1
  bin_end <- bin_start + window_size - 1
  window_id <- paste0(chrom, '_', bin_start, '_', bin_end)
  return(data.frame(
    CHROM = chrom,
    BIN_START = bin_start,
    BIN_END = bin_end,
    Window_ID = window_id
  ))
}

# Function to parse CHROM:POS locus name
parse_locus_name <- function(locus_name) {
  parts <- strsplit(locus_name, ':')[[1]]
  if (length(parts) == 2) {
    return(data.frame(
      CHROM = parts[1],
      POS = as.numeric(parts[2])
    ))
  } else {
    return(data.frame(
      CHROM = NA_character_,
      POS = NA_real_
    ))
  }
}

################################################################################
# LOAD PERCENTILE OUTLIERS
################################################################################

cat('Loading percentile outliers...\n')

percentile_outliers <- list()

for (threshold in PERCENTILE_THRESHOLDS) {
  file_path <- file.path(PERCENTILE_DIR, paste0('percentile_outliers_', threshold, '.tsv'))

  if (file.exists(file_path)) {
    df <- read.table(file_path, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    percentile_outliers[[threshold]] <- df
    cat('  Loaded', threshold, ':', nrow(df), 'outlier windows\n')
  } else {
    cat('  WARNING: File not found:', file_path, '\n')
  }
}

################################################################################
# PROCESS EACH COMPARISON
################################################################################

cat('\nProcessing OutFLANK results...\n')

all_concordance_results <- list()

for (comparison in VALID_COMPARISONS) {
  cat('\n--- Processing', comparison, '---\n')

  # Load OutFLANK results
  rds_path <- file.path(OUTFLANK_DIR, comparison, 'outflank_full_results.rds')

  if (!file.exists(rds_path)) {
    cat('  WARNING: RDS file not found:', rds_path, '\n')
    next
  }

  outflank_data <- readRDS(rds_path)
  results_df <- outflank_data$results

  cat('  Total SNPs in OutFLANK results:', nrow(results_df), '\n')

  # Get outliers at q < 0.05
  outliers <- results_df %>%
    filter(!is.na(qvalues) & qvalues < Q_THRESHOLD)

  n_outliers <- nrow(outliers)
  cat('  OutFLANK outliers (q <', Q_THRESHOLD, '):', n_outliers, '\n')

  if (n_outliers == 0) {
    cat('  No outliers to process\n')
    next
  }

  # Check if locus names are in CHROM:POS format
  sample_locus <- outliers$LocusName[1]
  if (!grepl(':', sample_locus)) {
    cat('  ERROR: Locus names are not in CHROM:POS format\n')
    cat('  Sample locus name:', sample_locus, '\n')
    cat('  Please re-run outflank_diagnostic.R with the updated version\n')
    next
  }

  # Parse locus names and map to windows
  cat('  Mapping SNPs to windows...\n')

  outlier_windows <- outliers %>%
    rowwise() %>%
    mutate(
      parsed = list(parse_locus_name(LocusName))
    ) %>%
    unnest(parsed) %>%
    filter(!is.na(CHROM) & !is.na(POS)) %>%
    rowwise() %>%
    mutate(
      window_info = list(map_snp_to_window(CHROM, POS, WINDOW_SIZE))
    ) %>%
    unnest(window_info, names_sep = '_') %>%
    ungroup()

  # Get unique windows containing OutFLANK outliers
  outflank_window_ids <- unique(outlier_windows$window_info_Window_ID)
  cat('  Unique windows with OutFLANK outliers:', length(outflank_window_ids), '\n')

  # Calculate concordance for each percentile threshold
  for (threshold in PERCENTILE_THRESHOLDS) {
    if (!(threshold %in% names(percentile_outliers))) next

    # Get percentile outliers for this comparison
    pct_df <- percentile_outliers[[threshold]] %>%
      filter(Comparison == comparison)

    if (nrow(pct_df) == 0) {
      cat('    No percentile outliers for', comparison, 'at', threshold, '\n')
      next
    }

    percentile_window_ids <- unique(pct_df$Window_ID)

    # Calculate overlap
    both <- intersect(outflank_window_ids, percentile_window_ids)
    outflank_only <- setdiff(outflank_window_ids, percentile_window_ids)
    percentile_only <- setdiff(percentile_window_ids, outflank_window_ids)

    # Concordance metrics
    n_both <- length(both)
    n_outflank_only <- length(outflank_only)
    n_percentile_only <- length(percentile_only)
    n_outflank_total <- length(outflank_window_ids)
    n_percentile_total <- length(percentile_window_ids)

    # Concordance percentage (what fraction of OutFLANK windows are also percentile outliers)
    concordance_pct <- ifelse(n_outflank_total > 0,
                               round(100 * n_both / n_outflank_total, 2), NA)

    # Jaccard index (intersection over union)
    union_size <- length(union(outflank_window_ids, percentile_window_ids))
    jaccard <- ifelse(union_size > 0, round(n_both / union_size, 4), NA)

    cat('    ', threshold, ': Concordance =', concordance_pct, '% |',
        'Jaccard =', jaccard, '|',
        'Both =', n_both, '|',
        'OutFLANK only =', n_outflank_only, '|',
        'Percentile only =', n_percentile_only, '\n')

    # Store results
    all_concordance_results[[paste(comparison, threshold, sep = '_')]] <- data.frame(
      Comparison = comparison,
      Threshold = threshold,
      OutFLANK_SNPs = n_outliers,
      OutFLANK_Windows = n_outflank_total,
      Percentile_Windows = n_percentile_total,
      Both_Methods = n_both,
      OutFLANK_Only = n_outflank_only,
      Percentile_Only = n_percentile_only,
      Concordance_Pct = concordance_pct,
      Jaccard_Index = jaccard
    )

    # Save detailed concordance for this comparison and threshold
    if (n_both > 0) {
      both_windows_df <- data.frame(
        Window_ID = both,
        Comparison = comparison,
        Threshold = threshold,
        Method = 'BOTH'
      )

      outfile <- file.path(OUTPUT_DIR,
                           paste0('concordant_windows_', comparison, '_', threshold, '.tsv'))
      write.table(both_windows_df, outfile, sep = '\t', row.names = FALSE, quote = FALSE)
    }
  }
}

################################################################################
# GENERATE SUMMARY
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('CONCORDANCE SUMMARY\n')
cat('===============================================================================\n\n')

if (length(all_concordance_results) > 0) {
  summary_df <- do.call(rbind, all_concordance_results)
  rownames(summary_df) <- NULL

  # Save summary table
  summary_file <- file.path(OUTPUT_DIR, 'concordance_summary.tsv')
  write.table(summary_df, summary_file, sep = '\t', row.names = FALSE, quote = FALSE)
  cat('Summary table saved to:', summary_file, '\n\n')

  # Print summary
  print(summary_df)

  # Create summary plot
  if (nrow(summary_df) > 0) {
    p <- ggplot(summary_df, aes(x = Threshold, y = Concordance_Pct, fill = Comparison)) +
      geom_bar(stat = 'identity', position = 'dodge', alpha = 0.8) +
      geom_text(aes(label = paste0(Concordance_Pct, '%')),
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 3) +
      labs(
        title = 'OutFLANK-Percentile Concordance',
        subtitle = 'Percentage of OutFLANK outlier windows also identified by percentile method',
        x = 'Percentile Threshold',
        y = 'Concordance (%)',
        fill = 'Comparison'
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = 'bold'),
        legend.position = 'bottom'
      ) +
      scale_fill_brewer(palette = 'Set2') +
      ylim(0, 100)

    ggsave(file.path(OUTPUT_DIR, 'concordance_plot.png'), p,
           width = 10, height = 7, dpi = 300)
    ggsave(file.path(OUTPUT_DIR, 'concordance_plot.pdf'), p,
           width = 10, height = 7)
    cat('\nConcordance plot saved\n')
  }

  # Generate text summary report
  report_file <- file.path(OUTPUT_DIR, 'concordance_summary.txt')
  sink(report_file)

  cat('================================================================================\n')
  cat('OUTFLANK-PERCENTILE CONCORDANCE VALIDATION REPORT\n')
  cat('================================================================================\n')
  cat('Generated:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

  cat('METHODOLOGY\n')
  cat('-----------\n')
  cat('1. OutFLANK identifies SNP-level outliers using q-value < 0.05\n')
  cat('2. SNP positions are mapped to 5kb genomic windows\n')
  cat('3. Windows are compared with percentile-based outlier windows\n')
  cat('4. Concordance = (Windows identified by both methods) / (OutFLANK windows)\n\n')

  cat('RESULTS BY COMPARISON\n')
  cat('---------------------\n\n')

  for (comparison in VALID_COMPARISONS) {
    comp_data <- summary_df %>% filter(Comparison == comparison)
    if (nrow(comp_data) > 0) {
      cat('### ', comparison, ' ###\n')
      cat('OutFLANK outlier SNPs:', unique(comp_data$OutFLANK_SNPs), '\n')
      cat('OutFLANK outlier windows:', unique(comp_data$OutFLANK_Windows), '\n\n')

      for (i in 1:nrow(comp_data)) {
        cat('  ', comp_data$Threshold[i], ':\n')
        cat('    Percentile windows:', comp_data$Percentile_Windows[i], '\n')
        cat('    Both methods:', comp_data$Both_Methods[i], '\n')
        cat('    Concordance:', comp_data$Concordance_Pct[i], '%\n')
        cat('    Jaccard index:', comp_data$Jaccard_Index[i], '\n\n')
      }
    }
  }

  cat('INTERPRETATION\n')
  cat('--------------\n')

  # Calculate overall concordance at 1% threshold
  pct1_data <- summary_df %>% filter(Threshold == '1pct')
  if (nrow(pct1_data) > 0) {
    mean_concordance <- mean(pct1_data$Concordance_Pct, na.rm = TRUE)
    cat('Mean concordance at 1% threshold:', round(mean_concordance, 1), '%\n\n')

    if (mean_concordance >= 70) {
      cat('EXCELLENT AGREEMENT: >70% of OutFLANK outlier windows are also\n')
      cat('identified by the percentile method, validating the unified approach.\n')
    } else if (mean_concordance >= 50) {
      cat('GOOD AGREEMENT: 50-70% concordance suggests both methods identify\n')
      cat('similar regions, with expected differences due to methodology.\n')
    } else if (mean_concordance >= 30) {
      cat('MODERATE AGREEMENT: 30-50% concordance indicates partial overlap.\n')
      cat('Methods may be identifying different aspects of differentiation.\n')
    } else {
      cat('LOW AGREEMENT: <30% concordance suggests substantial methodological\n')
      cat('differences. Consider investigating discordant regions.\n')
    }
  }

  cat('\n\nNOTE: Low concordance does not invalidate either method. OutFLANK uses\n')
  cat('formal statistical testing while percentile is empirical. Discordant\n')
  cat('windows may represent true biological differences in detection sensitivity.\n')

  sink()
  cat('\nReport saved to:', report_file, '\n')

} else {
  cat('No concordance results generated.\n')
  cat('This may indicate that OutFLANK results still use generic SNP_# format.\n')
  cat('Please re-run outflank_diagnostic.R with the updated version.\n')
}

################################################################################
# FINAL OUTPUT
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('ANALYSIS COMPLETE\n')
cat('===============================================================================\n')
cat('End time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

cat('Output files:\n')
cat('  ', file.path(OUTPUT_DIR, 'concordance_summary.tsv'), '\n')
cat('  ', file.path(OUTPUT_DIR, 'concordance_summary.txt'), '\n')
cat('  ', file.path(OUTPUT_DIR, 'concordance_plot.png'), '\n')
cat('  ', file.path(OUTPUT_DIR, 'concordance_plot.pdf'), '\n')

cat('\n===============================================================================\n')
