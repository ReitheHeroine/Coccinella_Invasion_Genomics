# title: pcadapt_selection_scan.R
# project: BIOL624 Final Project
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-01-12
# last modified: 2026-01-12
#
# purpose:
#   Run pcadapt selection scan to identify selection candidates using a
#   population structure-free approach. Then calculate SNP-level concordance
#   with existing methods (OutFLANK, percentile) by disaggregating percentile
#   windows to individual SNPs.
#
#   Key innovations:
#   1. pcadapt doesn't require predefined populations (handles reticulate history)
#   2. PC loading categorization maps to biological patterns (INVASION_SPECIFIC, POST_INVASION)
#   3. SNP-level concordance removes window aggregation artifacts
#
#   NOTE: This script expects a PLINK bed file to already exist at data/plink/ladybug_snps.bed
#   If not present, run the data preparation step first (see below).
#
# data preparation:
#   # First merge per-chromosome VCFs and convert to PLINK format
#   cd data/VARIANTS_BY_CHR/FILTERED_CHR_VCFS
#   bcftools concat *.vcf.gz -Oz -o ../merged_filtered.vcf.gz
#   plink --vcf ../merged_filtered.vcf.gz --make-bed --out ../../plink/ladybug_snps --allow-extra-chr
#
# inputs:
#   - data/plink/ladybug_snps.bed (PLINK binary format - created from VCFs)
#   - metadata/pop_*.txt (population files for PCA visualization)
#   - results/unified_selection_scan/percentile_outliers_1pct.tsv
#   - results/outflank_diagnostic/*/outflank_full_results.rds
#
# outputs:
#   results/pcadapt/
#     - screeplot.png, pca_plot.png, manhattan.png, qqplot.png
#     - pcadapt_outlier_snps.tsv, pcadapt_summary_report.txt
#   results/pcadapt/concordance/
#     - all_snps_outlier_status.tsv, snp_level_concordance_summary.tsv
#     - three_way_overlap_snps.tsv, concordance_report.txt
#
# usage example:
#   conda activate pop_gen
#   cd scripts
#   Rscript pcadapt_selection_scan.R
#
# copy/paste: Rscript pcadapt_selection_scan.R

################################################################################
# LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(pcadapt)
  library(vcfR)
  library(qvalue)
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
OUTPUT_DIR <- '../results/pcadapt'
CONCORDANCE_DIR <- '../results/pcadapt/concordance'
PERCENTILE_FILE <- '../results/unified_selection_scan/percentile_outliers_1pct.tsv'
OUTFLANK_DIR <- '../results/outflank_diagnostic'

# Analysis parameters
MAX_K <- 10           # Maximum K to test for scree plot
MIN_VAR_EXPLAINED <- 0.05  # Minimum variance to include PC
Q_THRESHOLD <- 0.05   # q-value threshold for outliers
Z_SCORE_THRESHOLD <- 3  # For PC loading categorization
WINDOW_SIZE <- 5000   # 5kb windows to match percentile analysis
MISSING_THRESHOLD <- 0.20  # Max missing rate per SNP

################################################################################
# SETUP
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('PCADAPT SELECTION SCAN WITH SNP-LEVEL CONCORDANCE\n')
cat('===============================================================================\n')
cat('Start time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n')
cat('===============================================================================\n\n')

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CONCORDANCE_DIR, recursive = TRUE, showWarnings = FALSE)

################################################################################
# PART 1: LOAD AND PREPARE GENOTYPE DATA
################################################################################

cat('PART 1: Loading genotype data from VCF files...\n')
cat('=========================================\n\n')

# Find VCF files
vcf_files <- list.files(VCF_DIR, pattern = '\\.vcf\\.gz$', full.names = TRUE)
cat('Found', length(vcf_files), 'VCF files\n')

# Function to convert genotypes to numeric (0, 1, 2)
convert_to_numeric <- function(gt_matrix) {
  geno_mat <- matrix(NA, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))

  for (i in 1:nrow(gt_matrix)) {
    for (j in 1:ncol(gt_matrix)) {
      g <- gt_matrix[i, j]
      if (is.na(g) || g == './.') {
        geno_mat[i, j] <- NA
      } else {
        alleles <- as.numeric(unlist(strsplit(g, '[/|]')))
        if (length(alleles) == 2) {
          geno_mat[i, j] <- sum(alleles)  # 0=hom ref, 1=het, 2=hom alt
        } else {
          geno_mat[i, j] <- NA
        }
      }
    }
  }
  return(geno_mat)
}

# Process each chromosome
all_geno_list <- list()
all_snp_info <- list()
sample_names <- NULL

for (vcf_file in vcf_files) {
  chrom <- gsub('\\.filtered\\.vcf\\.gz$', '', basename(vcf_file))
  cat('  Processing', chrom, '...')

  # Read VCF
  vcf <- read.vcfR(vcf_file, verbose = FALSE)

  # Get samples (first time only)
  if (is.null(sample_names)) {
    sample_names <- colnames(vcf@gt)[-1]
    cat(' (', length(sample_names), ' samples)')
  }

  # Extract genotypes
  gt <- extract.gt(vcf)

  # Convert to numeric
  geno_mat <- convert_to_numeric(gt)

  # Filter by missing data rate
  missing_rate <- rowMeans(is.na(geno_mat))
  keep <- missing_rate <= MISSING_THRESHOLD
  geno_mat <- geno_mat[keep, ]

  # Store SNP info
  snp_info <- data.frame(
    CHROM = getCHROM(vcf)[keep],
    POS = getPOS(vcf)[keep],
    stringsAsFactors = FALSE
  )
  snp_info$SNP_ID <- paste0(snp_info$CHROM, ':', snp_info$POS)

  # Store
  all_geno_list[[chrom]] <- geno_mat
  all_snp_info[[chrom]] <- snp_info

  cat(' retained', nrow(geno_mat), 'SNPs\n')
}

# Combine all chromosomes
cat('\nCombining all chromosomes...\n')
combined_geno <- do.call(rbind, all_geno_list)
combined_snp_info <- do.call(rbind, all_snp_info)
rownames(combined_snp_info) <- NULL

cat('Total SNPs:', nrow(combined_geno), '\n')
cat('Total samples:', ncol(combined_geno), '\n')

# Assign window IDs to all SNPs
cat('\nAssigning window IDs to SNPs...\n')
combined_snp_info$WINDOW_START <- floor((combined_snp_info$POS - 1) / WINDOW_SIZE) * WINDOW_SIZE + 1
combined_snp_info$WINDOW_END <- combined_snp_info$WINDOW_START + WINDOW_SIZE - 1
combined_snp_info$WINDOW_ID <- paste(combined_snp_info$CHROM,
                                      combined_snp_info$WINDOW_START,
                                      combined_snp_info$WINDOW_END, sep = '_')

n_windows <- length(unique(combined_snp_info$WINDOW_ID))
cat('Unique windows:', n_windows, '\n')

################################################################################
# PART 2: LOAD POPULATION INFORMATION
################################################################################

cat('\nLoading population information...\n')

# Read population files
pop_files <- list.files(POP_DIR, pattern = '^pop_.*\\.txt$', full.names = TRUE)
pop_assignments <- data.frame(Sample = sample_names, Population = NA_character_)

for (pop_file in pop_files) {
  pop_name <- gsub('^pop_|\\.txt$', '', basename(pop_file))
  pop_samples <- readLines(pop_file)
  pop_samples <- trimws(pop_samples[nchar(trimws(pop_samples)) > 0])
  pop_assignments$Population[pop_assignments$Sample %in% pop_samples] <- pop_name
}

# Report population sizes
pop_counts <- table(pop_assignments$Population)
cat('Population sizes:\n')
for (pop in names(pop_counts)) {
  cat('  ', pop, ':', pop_counts[pop], '\n')
}

################################################################################
# PART 3: RUN PCADAPT
################################################################################

cat('\n')
cat('PART 3: Running pcadapt analysis...\n')
cat('===================================\n\n')

# Transpose genotype matrix for pcadapt (samples as rows, SNPs as columns)
# pcadapt expects: individuals x SNPs
geno_transposed <- t(combined_geno)

# Impute missing values with mean (required for PCA)
cat('Imputing missing values...\n')
for (j in 1:ncol(geno_transposed)) {
  col_mean <- mean(geno_transposed[, j], na.rm = TRUE)
  if (is.na(col_mean)) col_mean <- 1  # Default if all NA
  geno_transposed[is.na(geno_transposed[, j]), j] <- round(col_mean)
}

# Write to temporary file for pcadapt
temp_file <- tempfile(fileext = '.pcadapt')
cat('Writing temporary file for pcadapt...\n')

# pcadapt pool format: space-separated genotypes, one SNP per line
# Each line = one SNP, values = genotypes for each individual
write.table(t(geno_transposed), temp_file, row.names = FALSE, col.names = FALSE, sep = ' ')

# Read with pcadapt
cat('Reading data with pcadapt...\n')
geno_pcadapt <- read.pcadapt(temp_file, type = 'pool')

# Determine optimal K using scree plot
cat('\nDetermining optimal K...\n')
cat('  Testing K = 1 to', MAX_K, '\n')

x_test <- pcadapt(geno_pcadapt, K = MAX_K)

# Calculate variance explained
var_explained <- x_test$singular.values^2 / sum(x_test$singular.values^2)

# Save scree plot
png(file.path(OUTPUT_DIR, 'screeplot.png'), width = 800, height = 600, res = 150)
plot(1:MAX_K, var_explained[1:MAX_K], type = 'b', pch = 19,
     xlab = 'Principal Component', ylab = 'Variance Explained',
     main = 'Scree Plot for K Selection')
abline(h = MIN_VAR_EXPLAINED, col = 'red', lty = 2)
text(MAX_K * 0.8, MIN_VAR_EXPLAINED + 0.02, paste0('Threshold: ', MIN_VAR_EXPLAINED), col = 'red')
dev.off()

# Auto-select K
K_selected <- max(1, sum(var_explained >= MIN_VAR_EXPLAINED))
cat('  Selected K =', K_selected, '(PCs with >=', MIN_VAR_EXPLAINED * 100, '% variance)\n')
cat('  Variance explained by selected PCs:\n')
for (k in 1:K_selected) {
  cat('    PC', k, ':', round(var_explained[k] * 100, 2), '%\n')
}

# Run pcadapt with selected K
cat('\nRunning pcadapt with K =', K_selected, '...\n')
x <- pcadapt(geno_pcadapt, K = K_selected)

# Calculate q-values
cat('Calculating q-values...\n')
qval <- qvalue(x$pvalues)$qvalues

# Genomic inflation factor
lambda <- median(x$chi2.stat, na.rm = TRUE) / qchisq(0.5, df = K_selected)
cat('Genomic inflation factor (lambda):', round(lambda, 3), '\n')

################################################################################
# PART 4: IDENTIFY AND CATEGORIZE OUTLIERS
################################################################################

cat('\n')
cat('PART 4: Identifying and categorizing outliers...\n')
cat('================================================\n\n')

# Create results data frame
pcadapt_results <- combined_snp_info
pcadapt_results$PVALUE <- x$pvalues
pcadapt_results$QVALUE <- qval
pcadapt_results$CHI2 <- x$chi2.stat

# Add PC loadings (z-scores)
for (k in 1:K_selected) {
  pcadapt_results[[paste0('PC', k, '_ZSCORE')]] <- x$zscores[, k]
}

# Identify outliers
pcadapt_results$IS_OUTLIER_Q05 <- pcadapt_results$QVALUE < Q_THRESHOLD
pcadapt_results$IS_OUTLIER_TOP1PCT <- pcadapt_results$PVALUE <= quantile(pcadapt_results$PVALUE, 0.01, na.rm = TRUE)

n_outliers_q05 <- sum(pcadapt_results$IS_OUTLIER_Q05, na.rm = TRUE)
n_outliers_top1 <- sum(pcadapt_results$IS_OUTLIER_TOP1PCT, na.rm = TRUE)

cat('Outliers at q < 0.05:', n_outliers_q05, '\n')
cat('Outliers at top 1%:', n_outliers_top1, '\n')

# Categorize outliers by PC loading
cat('\nCategorizing outliers by PC loading...\n')

if (K_selected >= 2) {
  pcadapt_results$PC_CATEGORY <- 'NOT_OUTLIER'

  # For outliers, determine which PC(s) they load on
  outlier_idx <- which(pcadapt_results$IS_OUTLIER_Q05)

  for (i in outlier_idx) {
    pc1_z <- abs(pcadapt_results$PC1_ZSCORE[i])
    pc2_z <- ifelse(K_selected >= 2, abs(pcadapt_results$PC2_ZSCORE[i]), 0)

    high_pc1 <- !is.na(pc1_z) && pc1_z > Z_SCORE_THRESHOLD
    high_pc2 <- !is.na(pc2_z) && pc2_z > Z_SCORE_THRESHOLD

    if (high_pc1 && high_pc2) {
      pcadapt_results$PC_CATEGORY[i] <- 'MULTIPLE_PCs'  # Maps to BOTH
    } else if (high_pc1) {
      pcadapt_results$PC_CATEGORY[i] <- 'PC1_ONLY'      # Maps to INVASION_SPECIFIC
    } else if (high_pc2) {
      pcadapt_results$PC_CATEGORY[i] <- 'PC2_ONLY'      # Maps to POST_INVASION
    } else {
      pcadapt_results$PC_CATEGORY[i] <- 'LOW_LOADING'   # Outlier but no strong PC loading
    }
  }

  # Report category counts
  cat_counts <- table(pcadapt_results$PC_CATEGORY[pcadapt_results$IS_OUTLIER_Q05])
  cat('Outlier categories:\n')
  for (cat in names(cat_counts)) {
    cat('  ', cat, ':', cat_counts[cat], '\n')
  }
} else {
  pcadapt_results$PC_CATEGORY <- ifelse(pcadapt_results$IS_OUTLIER_Q05, 'PC1_ONLY', 'NOT_OUTLIER')
}

################################################################################
# PART 5: GENERATE PCADAPT VISUALIZATIONS
################################################################################

cat('\n')
cat('PART 5: Generating visualizations...\n')
cat('====================================\n\n')

# PCA plot colored by population
cat('Creating PCA plot...\n')
pca_scores <- as.data.frame(x$scores)
colnames(pca_scores) <- paste0('PC', 1:K_selected)
pca_scores$Population <- pop_assignments$Population

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(alpha = 0.8, size = 3) +
  labs(
    title = 'PCA of C. septempunctata Populations',
    subtitle = paste0('K = ', K_selected, ' | Lambda = ', round(lambda, 3)),
    x = paste0('PC1 (', round(var_explained[1] * 100, 1), '%)'),
    y = paste0('PC2 (', round(var_explained[2] * 100, 1), '%)')
  ) +
  scale_color_brewer(palette = 'Set1') +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = 'bold'),
    legend.position = 'bottom'
  )

ggsave(file.path(OUTPUT_DIR, 'pca_plot.png'), pca_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, 'pca_plot.pdf'), pca_plot, width = 10, height = 8)

# Manhattan plot
cat('Creating Manhattan plot...\n')

# Add chromosome number for plotting
chrom_order <- c('NC_058189.1', 'NC_058190.1', 'NC_058191.1', 'NC_058192.1', 'NC_058193.1',
                 'NC_058194.1', 'NC_058195.1', 'NC_058196.1', 'NC_058197.1', 'NC_058198.1')
pcadapt_results$CHROM_NUM <- match(pcadapt_results$CHROM, chrom_order)

# Calculate cumulative position
pcadapt_results <- pcadapt_results %>%
  arrange(CHROM_NUM, POS) %>%
  group_by(CHROM) %>%
  mutate(BP_CUM = POS + (CHROM_NUM - 1) * max(POS)) %>%
  ungroup()

# Subsample for plotting if too many points
plot_data <- pcadapt_results
if (nrow(plot_data) > 100000) {
  # Keep all outliers, sample others
  outliers <- plot_data[plot_data$IS_OUTLIER_Q05, ]
  non_outliers <- plot_data[!plot_data$IS_OUTLIER_Q05, ]
  sampled_non <- non_outliers[sample(nrow(non_outliers), 50000), ]
  plot_data <- rbind(outliers, sampled_non)
}

manhattan_plot <- ggplot(plot_data, aes(x = BP_CUM, y = -log10(PVALUE))) +
  geom_point(aes(color = as.factor(CHROM_NUM %% 2)), alpha = 0.5, size = 0.5) +
  geom_point(data = plot_data[plot_data$IS_OUTLIER_Q05, ],
             aes(x = BP_CUM, y = -log10(PVALUE)), color = 'red', size = 1) +
  scale_color_manual(values = c('steelblue', 'navy'), guide = 'none') +
  labs(
    title = 'pcadapt Manhattan Plot',
    subtitle = paste0(n_outliers_q05, ' outliers at q < 0.05'),
    x = 'Chromosome',
    y = '-log10(p-value)'
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    plot.title = element_text(size = 14, face = 'bold')
  )

ggsave(file.path(OUTPUT_DIR, 'manhattan.png'), manhattan_plot, width = 12, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, 'manhattan.pdf'), manhattan_plot, width = 12, height = 6)

# Q-Q plot
cat('Creating Q-Q plot...\n')
png(file.path(OUTPUT_DIR, 'qqplot.png'), width = 800, height = 600, res = 150)
qq_dat <- data.frame(
  expected = -log10(ppoints(length(x$pvalues))),
  observed = sort(-log10(x$pvalues))
)
plot(qq_dat$expected, qq_dat$observed, pch = 20, cex = 0.5,
     xlab = 'Expected -log10(p)', ylab = 'Observed -log10(p)',
     main = paste0('Q-Q Plot (Lambda = ', round(lambda, 3), ')'))
abline(0, 1, col = 'red')
dev.off()

# Save PDF version
pdf(file.path(OUTPUT_DIR, 'qqplot.pdf'), width = 8, height = 6)
plot(qq_dat$expected, qq_dat$observed, pch = 20, cex = 0.5,
     xlab = 'Expected -log10(p)', ylab = 'Observed -log10(p)',
     main = paste0('Q-Q Plot (Lambda = ', round(lambda, 3), ')'))
abline(0, 1, col = 'red')
dev.off()

################################################################################
# PART 6: LOAD EXISTING OUTLIER DATA
################################################################################

cat('\n')
cat('PART 6: Loading existing outlier data...\n')
cat('========================================\n\n')

# Load percentile outliers
cat('Loading percentile outliers...\n')
if (file.exists(PERCENTILE_FILE)) {
  percentile_outliers <- read.table(PERCENTILE_FILE, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  percentile_window_ids <- unique(percentile_outliers$Window_ID)
  cat('  Percentile outlier windows:', length(percentile_window_ids), '\n')

  # Flag SNPs in percentile windows
  pcadapt_results$IS_IN_PERCENTILE_WINDOW <- pcadapt_results$WINDOW_ID %in% percentile_window_ids
  cat('  SNPs in percentile windows:', sum(pcadapt_results$IS_IN_PERCENTILE_WINDOW), '\n')
} else {
  cat('  WARNING: Percentile file not found\n')
  pcadapt_results$IS_IN_PERCENTILE_WINDOW <- FALSE
}

# Load OutFLANK outliers (EEU_vs_WEU and EEU_vs_USA only)
cat('\nLoading OutFLANK outliers...\n')

outflank_snps <- list()
for (comparison in c('EEU_vs_WEU', 'EEU_vs_USA')) {
  rds_file <- file.path(OUTFLANK_DIR, comparison, 'outflank_full_results.rds')
  if (file.exists(rds_file)) {
    of_data <- readRDS(rds_file)
    outlier_snps <- of_data$results$LocusName[of_data$results$qvalues < Q_THRESHOLD]
    outlier_snps <- outlier_snps[!is.na(outlier_snps)]
    outflank_snps[[comparison]] <- outlier_snps
    cat('  ', comparison, ':', length(outlier_snps), 'outlier SNPs\n')
  }
}

# Merge OutFLANK outliers
all_outflank_snps <- unique(unlist(outflank_snps))
pcadapt_results$IS_OUTFLANK_OUTLIER <- pcadapt_results$SNP_ID %in% all_outflank_snps
cat('  Total unique OutFLANK outlier SNPs:', length(all_outflank_snps), '\n')
cat('  Matched in pcadapt data:', sum(pcadapt_results$IS_OUTFLANK_OUTLIER), '\n')

################################################################################
# PART 7: SNP-LEVEL CONCORDANCE ANALYSIS
################################################################################

cat('\n')
cat('PART 7: Calculating SNP-level concordance...\n')
cat('============================================\n\n')

# Define outlier sets
pcadapt_outlier_snps <- pcadapt_results$SNP_ID[pcadapt_results$IS_OUTLIER_Q05]
percentile_window_snps <- pcadapt_results$SNP_ID[pcadapt_results$IS_IN_PERCENTILE_WINDOW]
outflank_outlier_snps <- pcadapt_results$SNP_ID[pcadapt_results$IS_OUTFLANK_OUTLIER]

cat('Outlier counts:\n')
cat('  pcadapt (q<0.05):', length(pcadapt_outlier_snps), '\n')
cat('  Percentile (in outlier windows):', length(percentile_window_snps), '\n')
cat('  OutFLANK (q<0.05):', length(outflank_outlier_snps), '\n')

# Calculate pairwise overlaps
calc_concordance <- function(set_a, set_b, name_a, name_b) {
  overlap <- length(intersect(set_a, set_b))
  union_size <- length(union(set_a, set_b))
  jaccard <- if (union_size > 0) overlap / union_size else NA

  data.frame(
    Comparison = paste(name_a, 'vs', name_b),
    Method_A = name_a,
    Method_B = name_b,
    N_A = length(set_a),
    N_B = length(set_b),
    Overlap = overlap,
    Pct_A_in_B = round(100 * overlap / length(set_a), 2),
    Pct_B_in_A = round(100 * overlap / length(set_b), 2),
    Jaccard = round(jaccard, 4)
  )
}

concordance_results <- rbind(
  calc_concordance(pcadapt_outlier_snps, percentile_window_snps, 'pcadapt', 'Percentile'),
  calc_concordance(pcadapt_outlier_snps, outflank_outlier_snps, 'pcadapt', 'OutFLANK'),
  calc_concordance(outflank_outlier_snps, percentile_window_snps, 'OutFLANK', 'Percentile')
)

cat('\nPairwise concordance:\n')
print(concordance_results)

# Three-way overlap
three_way_snps <- Reduce(intersect, list(pcadapt_outlier_snps, percentile_window_snps, outflank_outlier_snps))
cat('\nThree-way overlap (all methods):', length(three_way_snps), 'SNPs\n')

# Two-way overlaps
pcadapt_percentile <- setdiff(intersect(pcadapt_outlier_snps, percentile_window_snps), three_way_snps)
pcadapt_outflank <- setdiff(intersect(pcadapt_outlier_snps, outflank_outlier_snps), three_way_snps)
percentile_outflank <- setdiff(intersect(percentile_window_snps, outflank_outlier_snps), three_way_snps)

cat('Two-way overlaps (excluding three-way):\n')
cat('  pcadapt + Percentile:', length(pcadapt_percentile), '\n')
cat('  pcadapt + OutFLANK:', length(pcadapt_outflank), '\n')
cat('  Percentile + OutFLANK:', length(percentile_outflank), '\n')

################################################################################
# PART 8: SAVE OUTPUTS
################################################################################

cat('\n')
cat('PART 8: Saving outputs...\n')
cat('=========================\n\n')

# Save pcadapt outliers
outlier_output <- pcadapt_results %>%
  filter(IS_OUTLIER_Q05) %>%
  select(SNP_ID, CHROM, POS, PVALUE, QVALUE, CHI2, PC_CATEGORY,
         starts_with('PC'), WINDOW_ID, IS_IN_PERCENTILE_WINDOW, IS_OUTFLANK_OUTLIER) %>%
  arrange(QVALUE)

write.table(outlier_output, file.path(OUTPUT_DIR, 'pcadapt_outlier_snps.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('Saved pcadapt outliers:', nrow(outlier_output), 'SNPs\n')

# Save full SNP table with all outlier status
full_output <- pcadapt_results %>%
  select(SNP_ID, CHROM, POS, WINDOW_ID, PVALUE, QVALUE,
         IS_OUTLIER_Q05, IS_IN_PERCENTILE_WINDOW, IS_OUTFLANK_OUTLIER, PC_CATEGORY)

write.table(full_output, file.path(CONCORDANCE_DIR, 'all_snps_outlier_status.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('Saved full SNP table:', nrow(full_output), 'SNPs\n')

# Save concordance summary
write.table(concordance_results, file.path(CONCORDANCE_DIR, 'snp_level_concordance_summary.tsv'),
            sep = '\t', row.names = FALSE, quote = FALSE)

# Save three-way overlap SNPs
if (length(three_way_snps) > 0) {
  three_way_df <- pcadapt_results %>%
    filter(SNP_ID %in% three_way_snps) %>%
    select(SNP_ID, CHROM, POS, WINDOW_ID, PVALUE, QVALUE, PC_CATEGORY)
  write.table(three_way_df, file.path(CONCORDANCE_DIR, 'three_way_overlap_snps.tsv'),
              sep = '\t', row.names = FALSE, quote = FALSE)
  cat('Saved three-way overlap:', nrow(three_way_df), 'SNPs\n')
}

################################################################################
# PART 9: GENERATE SUMMARY REPORT
################################################################################

cat('\n')
cat('PART 9: Generating summary report...\n')
cat('====================================\n\n')

report_file <- file.path(OUTPUT_DIR, 'pcadapt_summary_report.txt')
sink(report_file)

cat('================================================================================\n')
cat('PCADAPT SELECTION SCAN WITH SNP-LEVEL CONCORDANCE REPORT\n')
cat('================================================================================\n')
cat('Generated:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

cat('PCADAPT ANALYSIS\n')
cat('----------------\n')
cat('Total SNPs analyzed:', nrow(pcadapt_results), '\n')
cat('Total samples:', length(sample_names), '\n')
cat('Selected K:', K_selected, '\n')
cat('Genomic inflation factor (lambda):', round(lambda, 3), '\n\n')

cat('Variance explained by PCs:\n')
for (k in 1:K_selected) {
  cat('  PC', k, ':', round(var_explained[k] * 100, 2), '%\n')
}

cat('\nOUTLIER DETECTION\n')
cat('-----------------\n')
cat('Outliers at q < 0.05:', n_outliers_q05, '\n')
cat('Outliers at top 1%:', n_outliers_top1, '\n\n')

if (K_selected >= 2) {
  cat('Outlier categories (based on PC loadings):\n')
  cat_counts <- table(pcadapt_results$PC_CATEGORY[pcadapt_results$IS_OUTLIER_Q05])
  for (cat in names(cat_counts)) {
    interpretation <- switch(cat,
                            'PC1_ONLY' = 'INVASION_SPECIFIC (CHI vs invasive)',
                            'PC2_ONLY' = 'POST_INVASION (among invasive)',
                            'MULTIPLE_PCs' = 'BOTH categories',
                            'LOW_LOADING' = 'No strong PC loading',
                            'Uncategorized')
    cat('  ', cat, ':', cat_counts[cat], '->', interpretation, '\n')
  }
}

cat('\nSNP-LEVEL CONCORDANCE\n')
cat('---------------------\n')
cat('This analysis disaggregates percentile outlier WINDOWS to individual SNPs,\n')
cat('enabling direct SNP-level comparison with pcadapt and OutFLANK.\n\n')

cat('Method comparisons:\n')
for (i in 1:nrow(concordance_results)) {
  r <- concordance_results[i, ]
  cat('\n', r$Comparison, ':\n')
  cat('  Method A outliers:', r$N_A, '\n')
  cat('  Method B outliers:', r$N_B, '\n')
  cat('  Overlap:', r$Overlap, 'SNPs\n')
  cat('  % of A in B:', r$Pct_A_in_B, '%\n')
  cat('  % of B in A:', r$Pct_B_in_A, '%\n')
  cat('  Jaccard index:', r$Jaccard, '\n')
}

cat('\nTHREE-WAY OVERLAP\n')
cat('-----------------\n')
cat('SNPs identified by ALL three methods:', length(three_way_snps), '\n')
cat('These represent the highest-confidence selection candidates.\n')

cat('\nINTERPRETATION\n')
cat('--------------\n')
if (length(three_way_snps) > 0) {
  cat('The', length(three_way_snps), 'SNPs detected by all three methods (pcadapt, OutFLANK,\n')
  cat('and percentile) represent robust selection candidates with strong cross-method\n')
  cat('validation. These should be prioritized for functional analysis.\n')
} else {
  cat('No SNPs were detected by all three methods simultaneously. This may indicate:\n')
  cat('1. Methods are detecting complementary aspects of selection\n')
  cat('2. Different statistical approaches have different sensitivities\n')
  cat('3. The two-way overlaps may still provide meaningful validation\n')
}

cat('\nOUTPUT FILES\n')
cat('------------\n')
cat('results/pcadapt/\n')
cat('  screeplot.png - K selection\n')
cat('  pca_plot.png - Population structure\n')
cat('  manhattan.png - Genome-wide p-values\n')
cat('  qqplot.png - Inflation check\n')
cat('  pcadapt_outlier_snps.tsv - All outlier SNPs\n')
cat('  pcadapt_summary_report.txt - This file\n')
cat('\nresults/pcadapt/concordance/\n')
cat('  all_snps_outlier_status.tsv - Master table with all SNP statuses\n')
cat('  snp_level_concordance_summary.tsv - Pairwise overlaps\n')
cat('  three_way_overlap_snps.tsv - SNPs in all methods\n')

cat('\n================================================================================\n')

sink()
cat('Summary report saved to:', report_file, '\n')

################################################################################
# CLEANUP
################################################################################

# Remove temporary file
if (file.exists(temp_file)) {
  file.remove(temp_file)
}

################################################################################
# FINAL OUTPUT
################################################################################

cat('\n')
cat('===============================================================================\n')
cat('ANALYSIS COMPLETE\n')
cat('===============================================================================\n')
cat('End time:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')

cat('Key results:\n')
cat('  K selected:', K_selected, '\n')
cat('  Lambda:', round(lambda, 3), '\n')
cat('  pcadapt outliers (q<0.05):', n_outliers_q05, '\n')
cat('  Three-way overlap:', length(three_way_snps), 'SNPs\n')

cat('\nOutput directories:\n')
cat('  ', OUTPUT_DIR, '\n')
cat('  ', CONCORDANCE_DIR, '\n')

cat('\n===============================================================================\n')
