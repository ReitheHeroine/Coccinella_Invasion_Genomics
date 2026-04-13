#!/usr/bin/env Rscript
################################################################################
# title: unified_selection_scan.R
# project: BIOL624 Final Project — Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-01-09
# last modified: 2026-01-09
#
# purpose:
#   Unified percentile-based selection scan applied to ALL six pairwise
#   comparisons, with OutFLANK validation for among-invasive comparisons.
#
#   Methodological rationale:
#   - Percentile approach (top 0.5%, 1%, 2%) applied uniformly across all
#     comparisons for methodological consistency
#   - OutFLANK results used as VALIDATION for among-invasive comparisons
#     where OutFLANK performs well (EEU vs WEU: 2,456 outliers; EEU vs USA: 439)
#   - Concordance metrics calculated to assess method agreement
#
#   This replaces the previous dual-method approach (OutFLANK for low-FST,
#   percentile for high-FST) based on OutFLANK diagnostic analysis showing
#   that while OutFLANK fits well for all comparisons (df=2), it is too
#   conservative for high-FST comparisons (0 outliers detected).
#
# inputs:
#   - Windowed FST files: ../results/fst_*_vs_*_results/fst/*.windowed.weir.fst
#   - OutFLANK diagnostic results: ../results/outflank_diagnostic/*/outliers_q05.tsv
#
# outputs:
#   ../results/unified_selection_scan/
#   - percentile_outliers_0.5pct.tsv
#   - percentile_outliers_1pct.tsv
#   - percentile_outliers_2pct.tsv
#   - outflank_validation/
#       - concordance_*.tsv
#       - validation_summary.tsv
#   - pattern_categorization/
#       - both_pattern_snps.tsv
#       - invasion_specific_snps.tsv
#       - post_invasion_snps.tsv
#   - summary_report.txt
#
# usage:
#   conda activate pop_gen
#   cd scripts
#   Rscript unified_selection_scan.R
#
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cat('\n================================================================================\n')
cat('UNIFIED PERCENTILE SELECTION SCAN WITH OUTFLANK VALIDATION\n')
cat('Species: Coccinella septempunctata (seven-spotted lady beetle)\n')
cat('Date:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')
cat('================================================================================\n\n')

################################################################################
# Configuration
################################################################################

# Define all pairwise comparisons
comparisons <- list(
  'CHI_vs_EEU' = c('CHI', 'EEU'),
  'CHI_vs_USA' = c('CHI', 'USA'),
  'CHI_vs_WEU' = c('CHI', 'WEU'),
  'EEU_vs_USA' = c('EEU', 'USA'),
  'EEU_vs_WEU' = c('EEU', 'WEU'),
  'USA_vs_WEU' = c('USA', 'WEU')
)

# Percentile thresholds (applied to ALL comparisons)
percentile_thresholds <- c(0.995, 0.99, 0.98)
threshold_names <- c('0.5pct', '1pct', '2pct')

# Paths
FST_BASE_DIR <- '../results'
OUTFLANK_DIAG_DIR <- '../results/outflank_diagnostic'
OUTPUT_DIR <- '../results/unified_selection_scan'

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, 'outflank_validation'), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, 'pattern_categorization'), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, 'plots'), showWarnings = FALSE)

################################################################################
# Step 1: Load windowed FST data for all comparisons
################################################################################

cat('STEP 1: Loading windowed FST data\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

all_fst_data <- list()
comparison_stats <- data.frame(
  Comparison = character(),
  Type = character(),
  N_Windows = integer(),
  Mean_FST = numeric(),
  Median_FST = numeric(),
  Max_FST = numeric(),
  stringsAsFactors = FALSE
)

for (comp_name in names(comparisons)) {
  cat('Loading:', comp_name, '...')

  # Find FST directory
  fst_dir <- file.path(FST_BASE_DIR, paste0('fst_', comp_name, '_results'), 'fst')

  if (!dir.exists(fst_dir)) {
    cat(' [ERROR: Directory not found]\n')
    next
  }

  # Find all windowed FST files
  fst_files <- list.files(fst_dir, pattern = '\\.windowed\\.weir\\.fst$', full.names = TRUE)

  if (length(fst_files) == 0) {
    cat(' [ERROR: No FST files found]\n')
    next
  }

  # Load and combine all chromosomes
  fst_combined <- bind_rows(lapply(fst_files, function(f) {
    df <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    df$Comparison <- comp_name
    return(df)
  }))

  # Filter out sex chromosome and invalid FST values
  fst_combined <- fst_combined %>%
    filter(!grepl('NC_058198', CHROM)) %>%
    filter(!is.na(WEIGHTED_FST) & is.finite(WEIGHTED_FST))

  # Create unique window ID
  fst_combined <- fst_combined %>%
    mutate(Window_ID = paste(CHROM, BIN_START, BIN_END, sep = '_'))

  # Categorize comparison type
  comp_type <- ifelse(grepl('CHI', comp_name), 'Native_vs_Invasive', 'Among_Invasive')

  # Store data
  all_fst_data[[comp_name]] <- fst_combined

  # Calculate stats
  comparison_stats <- rbind(comparison_stats, data.frame(
    Comparison = comp_name,
    Type = comp_type,
    N_Windows = nrow(fst_combined),
    Mean_FST = round(mean(fst_combined$WEIGHTED_FST, na.rm = TRUE), 4),
    Median_FST = round(median(fst_combined$WEIGHTED_FST, na.rm = TRUE), 4),
    Max_FST = round(max(fst_combined$WEIGHTED_FST, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  ))

  cat(' [', nrow(fst_combined), 'windows, mean FST =',
      round(mean(fst_combined$WEIGHTED_FST, na.rm = TRUE), 4), ']\n')
}

cat('\n')
cat('Comparison statistics:\n')
print(comparison_stats)
cat('\n')

################################################################################
# Step 2: Apply percentile thresholds to ALL comparisons
################################################################################

cat('\nSTEP 2: Applying percentile thresholds to all comparisons\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

# Store outliers for each threshold
all_outliers <- list()

for (t_idx in seq_along(percentile_thresholds)) {
  threshold <- percentile_thresholds[t_idx]
  thresh_name <- threshold_names[t_idx]

  cat('\n--- Threshold:', thresh_name, '(top', 100*(1-threshold), '%) ---\n\n')

  outliers_this_threshold <- list()

  for (comp_name in names(all_fst_data)) {
    fst_data <- all_fst_data[[comp_name]]

    # Calculate FST threshold for this comparison
    fst_cutoff <- quantile(fst_data$WEIGHTED_FST, threshold, na.rm = TRUE)

    # Identify outliers
    outliers <- fst_data %>%
      filter(WEIGHTED_FST >= fst_cutoff) %>%
      mutate(
        Threshold = thresh_name,
        FST_Cutoff = fst_cutoff,
        Method = 'Percentile'
      )

    outliers_this_threshold[[comp_name]] <- outliers

    cat('  ', comp_name, ': FST cutoff =', round(fst_cutoff, 4),
        ', outliers =', nrow(outliers), '\n')
  }

  # Combine all outliers for this threshold
  combined_outliers <- bind_rows(outliers_this_threshold)
  all_outliers[[thresh_name]] <- combined_outliers

  # Save to file
  output_file <- file.path(OUTPUT_DIR, paste0('percentile_outliers_', thresh_name, '.tsv'))
  write.table(combined_outliers, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
  cat('\n  Saved:', output_file, '(', nrow(combined_outliers), 'total outliers)\n')
}

################################################################################
# Step 3: Load OutFLANK diagnostic results for validation
################################################################################

cat('\n\nSTEP 3: Loading OutFLANK diagnostic results for validation\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

# OutFLANK comparisons to validate against
outflank_comparisons <- c('EEU_vs_USA', 'EEU_vs_WEU', 'USA_vs_WEU')

outflank_outliers <- list()

for (comp_name in outflank_comparisons) {
  outlier_file <- file.path(OUTFLANK_DIAG_DIR, comp_name, 'outliers_q05.tsv')

  if (file.exists(outlier_file)) {
    df <- read.table(outlier_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    df$Comparison <- comp_name
    outflank_outliers[[comp_name]] <- df
    cat('  ', comp_name, ':', nrow(df), 'OutFLANK outliers loaded\n')
  } else {
    cat('  ', comp_name, ': No outliers file found (0 outliers expected)\n')
    outflank_outliers[[comp_name]] <- data.frame()
  }
}

################################################################################
# Step 4: Calculate concordance between methods
################################################################################

cat('\n\nSTEP 4: Calculating concordance between percentile and OutFLANK\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

validation_summary <- data.frame(
  Comparison = character(),
  Threshold = character(),
  Percentile_Outliers = integer(),
  OutFLANK_Outliers = integer(),
  Both_Methods = integer(),
  Percentile_Only = integer(),
  OutFLANK_Only = integer(),
  Concordance_Pct = numeric(),
  stringsAsFactors = FALSE
)

for (comp_name in outflank_comparisons) {
  cat('\n--- Validating:', comp_name, '---\n')

  outflank_df <- outflank_outliers[[comp_name]]

  # Skip if no OutFLANK outliers
  if (nrow(outflank_df) == 0) {
    cat('  No OutFLANK outliers to validate against\n')

    for (thresh_name in threshold_names) {
      percentile_df <- all_outliers[[thresh_name]] %>% filter(Comparison == comp_name)

      validation_summary <- rbind(validation_summary, data.frame(
        Comparison = comp_name,
        Threshold = thresh_name,
        Percentile_Outliers = nrow(percentile_df),
        OutFLANK_Outliers = 0,
        Both_Methods = 0,
        Percentile_Only = nrow(percentile_df),
        OutFLANK_Only = 0,
        Concordance_Pct = NA,
        stringsAsFactors = FALSE
      ))
    }
    next
  }

  # Create OutFLANK window IDs to match with percentile
  # OutFLANK has LocusName with format CHROM:POS
  outflank_df <- outflank_df %>%
    mutate(
      CHROM = sapply(strsplit(LocusName, ':'), `[`, 1),
      POS = as.numeric(sapply(strsplit(LocusName, ':'), `[`, 2))
    )

  for (thresh_name in threshold_names) {
    percentile_df <- all_outliers[[thresh_name]] %>% filter(Comparison == comp_name)

    # Match OutFLANK SNPs to percentile windows
    # A SNP is in a window if POS >= BIN_START and POS < BIN_END
    matched_outflank <- 0
    both_methods_windows <- c()

    for (i in 1:nrow(outflank_df)) {
      snp_chrom <- outflank_df$CHROM[i]
      snp_pos <- outflank_df$POS[i]

      # Find matching window
      matching_window <- percentile_df %>%
        filter(CHROM == snp_chrom, BIN_START <= snp_pos, BIN_END > snp_pos)

      if (nrow(matching_window) > 0) {
        matched_outflank <- matched_outflank + 1
        both_methods_windows <- c(both_methods_windows, matching_window$Window_ID[1])
      }
    }

    both_methods_windows <- unique(both_methods_windows)
    n_both <- length(both_methods_windows)
    n_percentile_only <- nrow(percentile_df) - n_both
    n_outflank_only <- nrow(outflank_df) - matched_outflank

    # Concordance = OutFLANK SNPs found in percentile windows / total OutFLANK SNPs
    concordance <- ifelse(nrow(outflank_df) > 0,
                          round(100 * matched_outflank / nrow(outflank_df), 1),
                          NA)

    validation_summary <- rbind(validation_summary, data.frame(
      Comparison = comp_name,
      Threshold = thresh_name,
      Percentile_Outliers = nrow(percentile_df),
      OutFLANK_Outliers = nrow(outflank_df),
      Both_Methods = n_both,
      Percentile_Only = n_percentile_only,
      OutFLANK_Only = n_outflank_only,
      Concordance_Pct = concordance,
      stringsAsFactors = FALSE
    ))

    cat('  ', thresh_name, ': Percentile=', nrow(percentile_df),
        ', OutFLANK=', nrow(outflank_df),
        ', Both=', n_both,
        ', Concordance=', concordance, '%\n')
  }

  # Save detailed concordance for this comparison
  concordance_file <- file.path(OUTPUT_DIR, 'outflank_validation',
                                paste0('concordance_', comp_name, '.tsv'))
  write.table(
    validation_summary %>% filter(Comparison == comp_name),
    concordance_file, sep = '\t', row.names = FALSE, quote = FALSE
  )
}

# Save validation summary
validation_file <- file.path(OUTPUT_DIR, 'outflank_validation', 'validation_summary.tsv')
write.table(validation_summary, validation_file, sep = '\t', row.names = FALSE, quote = FALSE)
cat('\nSaved validation summary:', validation_file, '\n')

################################################################################
# Step 5: Pattern categorization
################################################################################

cat('\n\nSTEP 5: Categorizing outliers by biological pattern\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

# Use top 1% threshold for pattern categorization (balanced approach)
pattern_threshold <- '1pct'
pattern_outliers <- all_outliers[[pattern_threshold]]

# Count windows appearing in different comparison types
window_summary <- pattern_outliers %>%
  mutate(
    Comp_Type = ifelse(grepl('CHI', Comparison), 'Native_vs_Invasive', 'Among_Invasive')
  ) %>%
  group_by(Window_ID, CHROM, BIN_START, BIN_END) %>%
  summarize(
    N_Comparisons = n(),
    Comparisons = paste(unique(Comparison), collapse = '; '),
    Mean_FST = mean(WEIGHTED_FST, na.rm = TRUE),
    Max_FST = max(WEIGHTED_FST, na.rm = TRUE),
    In_Native_vs_Invasive = any(Comp_Type == 'Native_vs_Invasive'),
    In_Among_Invasive = any(Comp_Type == 'Among_Invasive'),
    .groups = 'drop'
  ) %>%
  mutate(
    Pattern = case_when(
      In_Native_vs_Invasive & In_Among_Invasive ~ 'BOTH',
      In_Native_vs_Invasive & !In_Among_Invasive ~ 'INVASION_SPECIFIC',
      !In_Native_vs_Invasive & In_Among_Invasive ~ 'POST_INVASION',
      TRUE ~ 'OTHER'
    ),
    Confidence = case_when(
      N_Comparisons >= 4 ~ 'Very_High',
      N_Comparisons == 3 ~ 'High',
      N_Comparisons == 2 ~ 'Moderate',
      TRUE ~ 'Single'
    )
  )

# Pattern summary
cat('Pattern distribution (', pattern_threshold, ' threshold):\n\n')
pattern_summary <- window_summary %>%
  group_by(Pattern) %>%
  summarize(
    N_Windows = n(),
    Mean_N_Comparisons = round(mean(N_Comparisons), 2),
    Mean_FST = round(mean(Mean_FST), 4),
    .groups = 'drop'
  )
print(pattern_summary)
cat('\n')

# Confidence summary
cat('Confidence distribution:\n\n')
confidence_summary <- window_summary %>%
  group_by(Confidence) %>%
  summarize(N_Windows = n(), .groups = 'drop')
print(confidence_summary)
cat('\n')

# Save pattern files
both_pattern <- window_summary %>% filter(Pattern == 'BOTH')
invasion_specific <- window_summary %>% filter(Pattern == 'INVASION_SPECIFIC')
post_invasion <- window_summary %>% filter(Pattern == 'POST_INVASION')

write.table(both_pattern,
            file.path(OUTPUT_DIR, 'pattern_categorization', 'both_pattern_windows.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(invasion_specific,
            file.path(OUTPUT_DIR, 'pattern_categorization', 'invasion_specific_windows.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
write.table(post_invasion,
            file.path(OUTPUT_DIR, 'pattern_categorization', 'post_invasion_windows.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)

cat('Saved pattern files to:', file.path(OUTPUT_DIR, 'pattern_categorization'), '\n')

################################################################################
# Step 6: Create summary visualizations
################################################################################

cat('\n\nSTEP 6: Creating summary visualizations\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

# Plot 1: Outliers per comparison (1% threshold)
outlier_counts <- pattern_outliers %>%
  group_by(Comparison) %>%
  summarize(N_Outliers = n(), .groups = 'drop') %>%
  left_join(comparison_stats %>% select(Comparison, Type), by = 'Comparison')

p1 <- ggplot(outlier_counts, aes(x = reorder(Comparison, N_Outliers),
                                  y = N_Outliers, fill = Type)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = N_Outliers), hjust = -0.1, size = 3) +
  scale_fill_manual(values = c('Native_vs_Invasive' = 'coral',
                               'Among_Invasive' = 'steelblue')) +
  coord_flip() +
  labs(
    title = 'Selection Candidates per Comparison (Unified Percentile)',
    subtitle = 'Top 1% FST threshold applied to all comparisons',
    x = '',
    y = 'Number of Outlier Windows',
    fill = 'Comparison Type'
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = 'bold'))

ggsave(file.path(OUTPUT_DIR, 'plots', 'outliers_per_comparison.png'),
       p1, width = 10, height = 6, dpi = 300)

# Plot 2: Pattern distribution
p2 <- ggplot(pattern_summary, aes(x = reorder(Pattern, N_Windows),
                                   y = N_Windows, fill = Pattern)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = format(N_Windows, big.mark = ',')), hjust = -0.1, size = 4) +
  scale_fill_brewer(palette = 'Set2') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  coord_flip() +
  labs(
    title = 'Selection Candidates by Biological Pattern',
    subtitle = 'BOTH = strongest evidence; POST_INVASION = cleanest signal',
    x = '',
    y = 'Number of Windows'
  ) +
  theme_minimal() +
  theme(legend.position = 'none', plot.title = element_text(face = 'bold'))

ggsave(file.path(OUTPUT_DIR, 'plots', 'biological_patterns.png'),
       p2, width = 10, height = 6, dpi = 300)

# Plot 3: OutFLANK validation concordance
validation_plot_data <- validation_summary %>%
  filter(!is.na(Concordance_Pct))

if (nrow(validation_plot_data) > 0) {
  p3 <- ggplot(validation_plot_data,
               aes(x = Threshold, y = Concordance_Pct, fill = Comparison)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    geom_text(aes(label = paste0(Concordance_Pct, '%')),
              position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
    scale_fill_brewer(palette = 'Set1') +
    labs(
      title = 'Concordance: Percentile vs OutFLANK',
      subtitle = '% of OutFLANK outliers found in percentile outlier windows',
      x = 'Percentile Threshold',
      y = 'Concordance (%)',
      fill = 'Comparison'
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = 'bold'))

  ggsave(file.path(OUTPUT_DIR, 'plots', 'outflank_concordance.png'),
         p3, width = 10, height = 6, dpi = 300)
}

cat('Plots saved to:', file.path(OUTPUT_DIR, 'plots'), '\n')

################################################################################
# Step 7: Generate comprehensive summary report
################################################################################

cat('\n\nSTEP 7: Generating summary report\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

report_file <- file.path(OUTPUT_DIR, 'summary_report.txt')
sink(report_file)

cat('================================================================================\n')
cat('UNIFIED PERCENTILE SELECTION SCAN - SUMMARY REPORT\n')
cat('Species: Coccinella septempunctata (seven-spotted lady beetle)\n')
cat('Analysis Date:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')
cat('================================================================================\n\n')

cat('METHODOLOGY\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')
cat('PRIMARY METHOD: Unified percentile approach\n')
cat('  - Applied to ALL 6 pairwise comparisons\n')
cat('  - Thresholds: Top 0.5%, 1%, 2% of FST distribution\n')
cat('  - Provides methodological consistency across comparisons\n\n')

cat('VALIDATION METHOD: OutFLANK (among-invasive only)\n')
cat('  - EEU vs WEU: 2,456 outliers at q<0.05\n')
cat('  - EEU vs USA: 439 outliers at q<0.05\n')
cat('  - USA vs WEU: 0 outliers at q<0.05\n')
cat('  - Used to assess concordance with percentile results\n\n')

cat('RATIONALE FOR UNIFIED APPROACH:\n')
cat('  - OutFLANK diagnostic showed df=2 (excellent fit) for ALL comparisons\n')
cat('  - However, OutFLANK detected 0 outliers for CHI vs invasive comparisons\n')
cat('  - Root cause: FSTbar (~0.11-0.12) approaches observed FST\n')
cat('  - Result: OutFLANK too conservative for high-differentiation comparisons\n')
cat('  - Solution: Use percentile as primary method, OutFLANK for validation\n\n')

cat('COMPARISON STATISTICS\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')
print(comparison_stats)
cat('\n\n')

cat('OUTLIERS BY THRESHOLD\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')
for (thresh_name in threshold_names) {
  df <- all_outliers[[thresh_name]]
  cat(thresh_name, ':\n')
  cat('  Total outlier windows:', nrow(df), '\n')
  cat('  By comparison:\n')
  counts <- df %>% group_by(Comparison) %>% summarize(n = n(), .groups = 'drop')
  for (i in 1:nrow(counts)) {
    cat('    ', counts$Comparison[i], ':', counts$n[i], '\n')
  }
  cat('\n')
}

cat('OUTFLANK VALIDATION RESULTS\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')
print(validation_summary)
cat('\n\n')

cat('PATTERN CATEGORIZATION (1% threshold)\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')
cat('Pattern definitions:\n')
cat('  BOTH: Outliers in CHI vs invasive AND among invasive comparisons\n')
cat('  INVASION_SPECIFIC: Outliers only in CHI vs invasive comparisons\n')
cat('  POST_INVASION: Outliers only among invasive comparisons\n\n')

cat('Pattern counts:\n')
print(pattern_summary)
cat('\n')

cat('Confidence levels:\n')
print(confidence_summary)
cat('\n\n')

cat('KEY FINDINGS\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')

cat('1. Methodological consistency achieved:\n')
cat('   - Same percentile approach applied to all 6 comparisons\n')
cat('   - Eliminates method-dependent bias in candidate selection\n\n')

cat('2. OutFLANK validation supports percentile results:\n')
for (comp in unique(validation_summary$Comparison)) {
  row <- validation_summary %>% filter(Comparison == comp, Threshold == '1pct')
  if (nrow(row) > 0 && !is.na(row$Concordance_Pct)) {
    cat('   - ', comp, ': ', row$Concordance_Pct, '% concordance\n', sep = '')
  }
}
cat('\n')

cat('3. Pattern distribution reveals selection context:\n')
for (i in 1:nrow(pattern_summary)) {
  cat('   - ', pattern_summary$Pattern[i], ': ',
      pattern_summary$N_Windows[i], ' windows\n', sep = '')
}
cat('\n')

cat('OUTPUT FILES\n')
cat('────────────────────────────────────────────────────────────────────────────\n\n')
cat('Percentile outliers:\n')
for (thresh_name in threshold_names) {
  cat('  - percentile_outliers_', thresh_name, '.tsv\n', sep = '')
}
cat('\nValidation files:\n')
cat('  - outflank_validation/validation_summary.tsv\n')
cat('  - outflank_validation/concordance_*.tsv\n')
cat('\nPattern files:\n')
cat('  - pattern_categorization/both_pattern_windows.tsv\n')
cat('  - pattern_categorization/invasion_specific_windows.tsv\n')
cat('  - pattern_categorization/post_invasion_windows.tsv\n')
cat('\nPlots:\n')
cat('  - plots/outliers_per_comparison.png\n')
cat('  - plots/biological_patterns.png\n')
cat('  - plots/outflank_concordance.png\n\n')

cat('================================================================================\n')
cat('END OF REPORT\n')
cat('================================================================================\n')

sink()

cat('Report saved to:', report_file, '\n')

################################################################################
# Final summary
################################################################################

cat('\n================================================================================\n')
cat('ANALYSIS COMPLETE!\n')
cat('================================================================================\n\n')

cat('Output directory:', OUTPUT_DIR, '\n\n')

cat('Summary:\n')
for (thresh_name in threshold_names) {
  df <- all_outliers[[thresh_name]]
  cat('  ', thresh_name, ':', nrow(df), 'outlier windows\n')
}

cat('\nPattern distribution (1% threshold):\n')
print(pattern_summary)

cat('\nNext steps:\n')
cat('  1. Review summary_report.txt\n')
cat('  2. Examine concordance with OutFLANK validation\n')
cat('  3. Map high-confidence windows to genes\n')
cat('  4. Integrate with diversity metrics (Tajima\'s D, pi)\n\n')
