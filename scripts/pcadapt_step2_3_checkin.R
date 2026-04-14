# title: pcadapt_step2_3_checkin.R
# project: BIOL624 Final Project
# author: Reina Hastings (reinahastings13@gmail.com)
# date created: 2026-04-13
# last modified: 2026-04-13
#
# purpose:
#   Steps 2-3 of pcadapt selection scan: generate scree plot for K selection
#   and inspect PCA for missingness artifacts per Yi (2022) protocol.
#   This is a check-in script -- review outputs before proceeding to full scan.
#
# inputs:
#   - data/plink/ladybug_snps_pcadapt.bed (81 individuals, post-exclusion)
#   - data/plink/ladybug_snps_pcadapt.imiss (per-sample missingness)
#   - metadata/pop_CHI.txt, pop_EEU.txt, pop_WEU.txt, pop_USA.txt
#
# outputs:
#   results/pcadapt/plots/screeplot.png
#   results/pcadapt/plots/pca_scores_by_pop.png
#   results/pcadapt/plots/pca_scores_by_missingness.png
#   results/pcadapt/step2_3_checkin.txt (console summary)
#
# usage example:
#   conda activate pop_gen
#   Rscript scripts/pcadapt_step2_3_checkin.R

# --- Libraries ---

suppressPackageStartupMessages({
  library(pcadapt)
  library(ggplot2)
  library(dplyr)
})

# --- Configuration ---

PLINK_PREFIX <- 'data/plink/ladybug_snps_pcadapt'
IMISS_FILE   <- 'data/plink/ladybug_snps_pcadapt.imiss'
POP_DIR      <- 'metadata'
PLOT_DIR     <- 'results/pcadapt/plots'
OUTPUT_DIR   <- 'results/pcadapt'
MAX_K        <- 10

# Samples to exclude (already removed from PLINK file)
EXCLUDED_SAMPLES <- c('SUS5', 'MUS113', 'CZE10', 'MUS114', 'DEU8', 'CHN11', 'SUS66')

# Borderline samples to inspect (20-50% missing)
BORDERLINE_SAMPLES <- c('WUS7', 'NLD3', 'DEU12', 'SVK1')

# Population color palette
POP_COLORS <- c('CHI' = '#E41A1C', 'EEU' = '#377EB8', 'WEU' = '#4DAF4A', 'USA' = '#984EA3')

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load Population Assignments ---

cat('\n========================================================================\n')
cat('PCADAPT STEPS 2-3: SCREE PLOT AND PCA MISSINGNESS INSPECTION\n')
cat('========================================================================\n\n')

# Read .fam to get sample order in PLINK file
fam <- read.table(paste0(PLINK_PREFIX, '.fam'), stringsAsFactors = FALSE)
sample_ids <- fam$V2
cat('Samples in PLINK file:', length(sample_ids), '\n')

# Load population files and assign
pop_chi <- readLines(file.path(POP_DIR, 'pop_CHI.txt'))
pop_eeu <- readLines(file.path(POP_DIR, 'pop_EEU.txt'))
pop_weu <- readLines(file.path(POP_DIR, 'pop_WEU.txt'))
pop_usa <- readLines(file.path(POP_DIR, 'pop_USA.txt'))

pop_chi <- trimws(pop_chi[nchar(trimws(pop_chi)) > 0])
pop_eeu <- trimws(pop_eeu[nchar(trimws(pop_eeu)) > 0])
pop_weu <- trimws(pop_weu[nchar(trimws(pop_weu)) > 0])
pop_usa <- trimws(pop_usa[nchar(trimws(pop_usa)) > 0])

# Remove excluded samples from pop lists
pop_chi <- setdiff(pop_chi, EXCLUDED_SAMPLES)
pop_eeu <- setdiff(pop_eeu, EXCLUDED_SAMPLES)
pop_weu <- setdiff(pop_weu, EXCLUDED_SAMPLES)
pop_usa <- setdiff(pop_usa, EXCLUDED_SAMPLES)

# Assign populations to sample_ids in .fam order
pop_labels <- rep(NA_character_, length(sample_ids))
pop_labels[sample_ids %in% pop_chi] <- 'CHI'
pop_labels[sample_ids %in% pop_eeu] <- 'EEU'
pop_labels[sample_ids %in% pop_weu] <- 'WEU'
pop_labels[sample_ids %in% pop_usa] <- 'USA'

cat('Population sizes (post-exclusion):\n')
cat('  CHI:', sum(pop_labels == 'CHI', na.rm = TRUE), '\n')
cat('  EEU:', sum(pop_labels == 'EEU', na.rm = TRUE), '\n')
cat('  WEU:', sum(pop_labels == 'WEU', na.rm = TRUE), '\n')
cat('  USA:', sum(pop_labels == 'USA', na.rm = TRUE), '\n')
cat('  Unassigned:', sum(is.na(pop_labels)), '\n')

if (any(is.na(pop_labels))) {
  cat('  WARNING: Unassigned samples:', sample_ids[is.na(pop_labels)], '\n')
}

# --- Load Missingness Data ---

cat('\nLoading per-sample missingness...\n')
imiss <- read.table(IMISS_FILE, header = TRUE, stringsAsFactors = FALSE)
# Match to .fam order
miss_rate <- imiss$F_MISS[match(sample_ids, imiss$IID)]
cat('  Mean missingness:', round(mean(miss_rate) * 100, 2), '%\n')
cat('  Range:', round(min(miss_rate) * 100, 2), '-', round(max(miss_rate) * 100, 2), '%\n')

# Report borderline samples
cat('\nBorderline samples (20-50% missing, inspect for PCA drift):\n')
for (s in BORDERLINE_SAMPLES) {
  idx <- which(sample_ids == s)
  if (length(idx) > 0) {
    pop <- pop_labels[idx]
    mr <- round(miss_rate[idx] * 100, 2)
    cat('  ', s, '(', pop, '):', mr, '% missing\n')
  } else {
    cat('  ', s, ': NOT FOUND in PLINK file (already excluded?)\n')
  }
}

# --- Step 2: Read PLINK and Generate Scree Plot ---

cat('\n--- Step 2: Scree Plot for K Selection ---\n\n')

bed_file <- paste0(PLINK_PREFIX, '.bed')
cat('Reading PLINK bed file:', bed_file, '\n')
geno <- read.pcadapt(bed_file, type = 'bed')

cat('Running pcadapt with K =', MAX_K, '(overspecified for scree)...\n')
x_scree <- pcadapt(geno, K = MAX_K)

# Variance explained
singular_vals <- x_scree$singular.values
var_explained <- singular_vals^2 / sum(singular_vals^2)

cat('\nVariance explained by each PC:\n')
for (k in 1:MAX_K) {
  cat('  PC', k, ':', round(var_explained[k] * 100, 2), '%\n')
}

# Cattell's rule: K = number of PCs before the plateau (elbow)
# Compute differences in variance explained
diffs <- diff(var_explained[1:MAX_K])
# The elbow is where the rate of decrease slows dramatically
# Look for the first point where the drop is less than 1/3 of the previous drop
K_cattell <- 1
for (k in 2:length(diffs)) {
  if (abs(diffs[k]) < abs(diffs[k-1]) / 3) {
    K_cattell <- k
    break
  }
}
# Ensure K >= 2 (we need at least 2 for pattern categorization)
K_cattell <- max(2, K_cattell)
cat('\nCattell rule suggests K =', K_cattell, '\n')

# Scree plot
png(file.path(PLOT_DIR, 'screeplot.png'), width = 800, height = 600, res = 150)
par(mar = c(5, 5, 4, 2))
plot(1:MAX_K, var_explained[1:MAX_K] * 100, type = 'b', pch = 19, cex = 1.2,
     xlab = 'Principal Component', ylab = 'Variance Explained (%)',
     main = 'Scree Plot for K Selection (pcadapt)',
     las = 1, xaxt = 'n')
axis(1, at = 1:MAX_K)
abline(v = K_cattell + 0.5, col = 'red', lty = 2, lwd = 1.5)
text(K_cattell + 0.7, max(var_explained[1:MAX_K] * 100) * 0.9,
     paste0('Suggested K = ', K_cattell), col = 'red', adj = 0, cex = 0.8)
dev.off()
cat('Scree plot saved to:', file.path(PLOT_DIR, 'screeplot.png'), '\n')

# --- Step 3: PCA Score Plots with Population and Missingness Coloring ---

cat('\n--- Step 3: PCA Missingness Inspection (Yi 2022 Protocol) ---\n\n')

# Extract PCA scores
scores <- as.data.frame(x_scree$scores[, 1:min(4, MAX_K)])
colnames(scores) <- paste0('PC', 1:ncol(scores))
scores$Sample <- sample_ids
scores$Population <- factor(pop_labels, levels = c('CHI', 'EEU', 'WEU', 'USA'))
scores$Missingness <- miss_rate * 100
scores$Is_Borderline <- sample_ids %in% BORDERLINE_SAMPLES

# Plot 1: PCA colored by population
p_pop <- ggplot(scores, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_point(data = filter(scores, Is_Borderline),
             aes(x = PC1, y = PC2), shape = 1, size = 5, color = 'black', stroke = 1.2) +
  geom_text(data = filter(scores, Is_Borderline),
            aes(x = PC1, y = PC2, label = Sample),
            vjust = -1.2, hjust = 0.5, size = 2.5, color = 'black') +
  scale_color_manual(values = POP_COLORS) +
  labs(
    title = 'PCA Score Plot - Colored by Population',
    subtitle = paste0('81 individuals | K = ', MAX_K, ' (overspecified) | Circled = borderline samples'),
    x = paste0('PC1 (', round(var_explained[1] * 100, 1), '%)'),
    y = paste0('PC2 (', round(var_explained[2] * 100, 1), '%)')
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'bottom')

ggsave(file.path(PLOT_DIR, 'pca_scores_by_pop.png'), p_pop, width = 9, height = 7, dpi = 300)
cat('PCA population plot saved.\n')

# Plot 2: PCA colored by missingness (Yi 2022 protocol)
p_miss <- ggplot(scores, aes(x = PC1, y = PC2, color = Missingness)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_point(data = filter(scores, Is_Borderline),
             aes(x = PC1, y = PC2), shape = 1, size = 5, color = 'black', stroke = 1.2) +
  geom_text(data = filter(scores, Is_Borderline),
            aes(x = PC1, y = PC2, label = paste0(Sample, '\n(', round(Missingness, 1), '%)')),
            vjust = -1.5, hjust = 0.5, size = 2.2, color = 'black') +
  scale_color_gradient(low = 'navy', high = 'red', name = 'Missing (%)') +
  labs(
    title = 'PCA Score Plot - Colored by Missingness (Yi 2022)',
    subtitle = 'Check: do borderline samples drift toward origin or cluster with populations?',
    x = paste0('PC1 (', round(var_explained[1] * 100, 1), '%)'),
    y = paste0('PC2 (', round(var_explained[2] * 100, 1), '%)')
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'right')

ggsave(file.path(PLOT_DIR, 'pca_scores_by_missingness.png'), p_miss, width = 10, height = 7, dpi = 300)
cat('PCA missingness plot saved.\n')

# --- Report on Borderline Samples ---

cat('\nBorderline sample PCA coordinates:\n')
borderline_scores <- scores[scores$Is_Borderline, ]
if (nrow(borderline_scores) > 0) {
  for (i in 1:nrow(borderline_scores)) {
    row <- borderline_scores[i, ]
    # Find population centroid
    pop_scores <- scores[scores$Population == row$Population & !scores$Is_Borderline, ]
    centroid_pc1 <- mean(pop_scores$PC1, na.rm = TRUE)
    centroid_pc2 <- mean(pop_scores$PC2, na.rm = TRUE)
    dist_to_centroid <- sqrt((row$PC1 - centroid_pc1)^2 + (row$PC2 - centroid_pc2)^2)
    pop_sd <- mean(c(sd(pop_scores$PC1, na.rm = TRUE), sd(pop_scores$PC2, na.rm = TRUE)))
    dist_in_sd <- dist_to_centroid / pop_sd

    cat(sprintf('  %s (%s, %.1f%% missing): PC1=%.3f, PC2=%.3f | %.1f SD from %s centroid\n',
                row$Sample, as.character(row$Population), row$Missingness,
                row$PC1, row$PC2, dist_in_sd, as.character(row$Population)))
  }
}

# Overall origin distance check
origin_dist <- sqrt(scores$PC1^2 + scores$PC2^2)
scores$Origin_Dist <- origin_dist
cat('\nDistance from PCA origin (mean by population):\n')
for (pop in c('CHI', 'EEU', 'WEU', 'USA')) {
  pop_dists <- origin_dist[pop_labels == pop]
  cat(sprintf('  %s: mean=%.3f, sd=%.3f\n', pop, mean(pop_dists, na.rm = TRUE), sd(pop_dists, na.rm = TRUE)))
}
cat('Borderline samples:\n')
for (s in BORDERLINE_SAMPLES) {
  idx <- which(sample_ids == s)
  if (length(idx) > 0) {
    cat(sprintf('  %s: origin_dist=%.3f\n', s, origin_dist[idx]))
  }
}

# --- Summary ---

cat('\n========================================================================\n')
cat('STEP 2-3 CHECK-IN SUMMARY\n')
cat('========================================================================\n')
cat('PLINK file: 803,821 variants x 81 individuals\n')
cat('Post-exclusion genotyping rate: 93.7%\n')
cat('Cattell rule K:', K_cattell, '\n')
cat('Variance: PC1=', round(var_explained[1]*100, 1), '%, PC2=', round(var_explained[2]*100, 1), '%\n')
cat('\nPlots to review:\n')
cat('  1. ', file.path(PLOT_DIR, 'screeplot.png'), '\n')
cat('  2. ', file.path(PLOT_DIR, 'pca_scores_by_pop.png'), '\n')
cat('  3. ', file.path(PLOT_DIR, 'pca_scores_by_missingness.png'), '\n')
cat('\nDECISIONS NEEDED:\n')
cat('  1. Confirm K (Cattell suggests K =', K_cattell, ')\n')
cat('  2. Keep or exclude borderline samples based on PCA plots\n')
cat('========================================================================\n')

# Save checkin report
sink(file.path(OUTPUT_DIR, 'step2_3_checkin.txt'))
cat('Step 2-3 Check-in Report\n')
cat('Generated:', format(Sys.time(), '%Y-%m-%d %H:%M:%S'), '\n\n')
cat('Variance explained:\n')
for (k in 1:MAX_K) {
  cat('  PC', k, ':', round(var_explained[k] * 100, 2), '%\n')
}
cat('\nCattell rule K:', K_cattell, '\n')
cat('\nBorderline sample distances from population centroids (in SD):\n')
if (nrow(borderline_scores) > 0) {
  for (i in 1:nrow(borderline_scores)) {
    row <- borderline_scores[i, ]
    pop_scores <- scores[scores$Population == row$Population & !scores$Is_Borderline, ]
    centroid_pc1 <- mean(pop_scores$PC1, na.rm = TRUE)
    centroid_pc2 <- mean(pop_scores$PC2, na.rm = TRUE)
    dist_to_centroid <- sqrt((row$PC1 - centroid_pc1)^2 + (row$PC2 - centroid_pc2)^2)
    pop_sd <- mean(c(sd(pop_scores$PC1, na.rm = TRUE), sd(pop_scores$PC2, na.rm = TRUE)))
    dist_in_sd <- dist_to_centroid / pop_sd
    cat(sprintf('  %s (%s, %.1f%% missing): %.1f SD from centroid\n',
                row$Sample, as.character(row$Population), row$Missingness, dist_in_sd))
  }
}
sink()
cat('\nCheck-in report saved to:', file.path(OUTPUT_DIR, 'step2_3_checkin.txt'), '\n')
