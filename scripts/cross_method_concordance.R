# title: cross_method_concordance.R
# project: BIOL624 Final Project
# author: Reina Hastings (reinahastings13@gmail.com)
# date created: 2026-04-13
# last modified: 2026-04-13
#
# purpose:
#   Integrate outlier SNPs from three detection methods (per-SNP percentile,
#   OutFLANK, pcadapt) into a unified concordance analysis. Computes pairwise
#   and three-way overlap metrics (Jaccard, directional containment), assigns
#   confidence tiers based on cross-method support, reconciles pattern categories,
#   and produces a master candidate table for downstream annotation.
#
# inputs:
#   - results/per_snp_selection_scan/snp_summary_1pct_chifiltered.tsv
#   - results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv
#   - results/outflank_diagnostic/EEU_vs_WEU/outflank_full_results.rds
#   - results/outflank_diagnostic/EEU_vs_USA/outflank_full_results.rds
#   - results/pcadapt/pcadapt_outliers_q05.tsv
#   - results/pcadapt/pcadapt_outliers_q01.tsv
#   - data/plink/ladybug_snps.bim (88-sample SNP universe)
#   - data/plink/ladybug_snps_pcadapt.bim (81-sample SNP universe)
#
# outputs:
#   - results/cross_method_concordance/master_candidate_table.tsv
#   - results/cross_method_concordance/concordance_summary.txt
#   - results/cross_method_concordance/plots/venn_post_invasion.png
#   - results/cross_method_concordance/plots/venn_invasion_specific.png
#   - results/cross_method_concordance/plots/manhattan_tier1.png
#
# usage example:
#   cd scripts && Rscript cross_method_concordance.R

################################################################################
# LIBRARIES
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(VennDiagram)
  library(futile.logger)  # suppress VennDiagram logging
})

################################################################################
# CONFIGURATION
################################################################################

# --- Paths (relative to scripts/ directory) ---
BASE_DIR <- ".."

# Per-SNP percentile inputs (CHI-filtered)
PERCENTILE_SUMMARY <- file.path(BASE_DIR, "results/per_snp_selection_scan/snp_summary_1pct_chifiltered.tsv")
PERCENTILE_DETAIL  <- file.path(BASE_DIR, "results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv")

# OutFLANK RDS files (full results, not truncated TSVs)
OUTFLANK_RDS <- list(
  EEU_vs_WEU = file.path(BASE_DIR, "results/outflank_diagnostic/EEU_vs_WEU/outflank_full_results.rds"),
  EEU_vs_USA = file.path(BASE_DIR, "results/outflank_diagnostic/EEU_vs_USA/outflank_full_results.rds")
)

# pcadapt inputs
PCADAPT_Q05 <- file.path(BASE_DIR, "results/pcadapt/pcadapt_outliers_q05.tsv")
PCADAPT_Q01 <- file.path(BASE_DIR, "results/pcadapt/pcadapt_outliers_q01.tsv")

# PLINK BIM files for SNP universe
BIM_FULL    <- file.path(BASE_DIR, "data/plink/ladybug_snps.bim")
BIM_PCADAPT <- file.path(BASE_DIR, "data/plink/ladybug_snps_pcadapt.bim")

# Output directory
OUT_DIR  <- file.path(BASE_DIR, "results/cross_method_concordance")
PLOT_DIR <- file.path(OUT_DIR, "plots")

# Sex chromosome to exclude from concordance
SEX_CHROM <- "NC_058189.1"

# Thresholds for Tier 2 "strong support" criteria
PCADAPT_STRONG_Q <- 0.001
PERCENTILE_STRONG_PCT <- 0.001  # top 0.1%

################################################################################
# SETUP
################################################################################

cat("\n")
cat("========================================================================\n")
cat("CROSS-METHOD CONCORDANCE ANALYSIS\n")
cat("========================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Suppress VennDiagram log file creation
flog.threshold(ERROR, name = "VennDiagramLogger")

################################################################################
# STEP 1: LOAD AND VALIDATE INPUTS
################################################################################

cat("--- Step 1: Loading inputs ---\n\n")

# --- 1a. Per-SNP percentile (CHI-filtered, 1% threshold) ---
cat("Loading per-SNP percentile outliers (CHI-filtered, 1%)...\n")
pct_summary <- read_tsv(PERCENTILE_SUMMARY, show_col_types = FALSE)
pct_detail  <- read_tsv(PERCENTILE_DETAIL, show_col_types = FALSE)

# Normalize column names to lowercase
names(pct_summary) <- tolower(names(pct_summary))

cat("  Unique outlier SNPs:", nrow(pct_summary), "\n")
cat("  Detail rows (SNP x comparison):", nrow(pct_detail), "\n")
cat("  Patterns:", paste(table(pct_summary$pattern), collapse = ", "), "\n")
cat("  Pattern names:", paste(names(table(pct_summary$pattern)), collapse = ", "), "\n\n")

# --- 1b. OutFLANK (extract from RDS, q < 0.05) ---
cat("Loading OutFLANK results from RDS files...\n")

outflank_outliers <- list()
for (comp_name in names(OUTFLANK_RDS)) {
  rds_path <- OUTFLANK_RDS[[comp_name]]
  if (!file.exists(rds_path)) {
    cat("  WARNING: RDS not found:", rds_path, "\n")
    next
  }
  rds_data <- readRDS(rds_path)
  results_df <- rds_data$results

  # Extract outliers at q < 0.05
  outliers <- results_df %>%
    filter(!is.na(qvalues) & qvalues < 0.05) %>%
    mutate(comparison = comp_name)

  outflank_outliers[[comp_name]] <- outliers
  cat("  ", comp_name, ":", nrow(outliers), "outliers at q<0.05\n")
}

outflank_all <- bind_rows(outflank_outliers)
cat("  Total OutFLANK outliers:", nrow(outflank_all), "\n")

# Parse OutFLANK SNP IDs (LocusName = CHROM:POS)
outflank_all <- outflank_all %>%
  mutate(
    chrom = sub(":.*", "", LocusName),
    pos   = as.integer(sub(".*:", "", LocusName)),
    snp_id = LocusName
  )

# Deduplicate: a SNP can be outlier in both comparisons
outflank_snps <- outflank_all %>%
  group_by(snp_id, chrom, pos) %>%
  summarise(
    outflank_comparisons = paste(comparison, collapse = ";"),
    outflank_q_min = min(qvalues),
    outflank_n_comparisons = n(),
    .groups = "drop"
  )
cat("  Unique OutFLANK outlier SNPs:", nrow(outflank_snps), "\n\n")

# --- 1c. pcadapt (q < 0.05) ---
cat("Loading pcadapt outliers (q<0.05)...\n")
pcadapt_q05 <- read_tsv(PCADAPT_Q05, show_col_types = FALSE)
pcadapt_q01 <- read_tsv(PCADAPT_Q01, show_col_types = FALSE)
cat("  Outliers at q<0.05:", nrow(pcadapt_q05), "\n")
cat("  Outliers at q<0.01:", nrow(pcadapt_q01), "\n")
cat("  PC categories:", paste(table(pcadapt_q05$pc_category), collapse = ", "), "\n")
cat("  Category names:", paste(names(table(pcadapt_q05$pc_category)), collapse = ", "), "\n\n")

# --- 1d. SNP universes ---
cat("Loading SNP universes from BIM files...\n")
bim_cols <- c("chrom", "snp_id", "cm", "pos", "a1", "a2")

bim_full <- read_tsv(BIM_FULL, col_names = bim_cols, show_col_types = FALSE) %>%
  mutate(snp_id = paste0(chrom, ":", pos))
bim_pcadapt <- read_tsv(BIM_PCADAPT, col_names = bim_cols, show_col_types = FALSE) %>%
  mutate(snp_id = paste0(chrom, ":", pos))

cat("  Full BIM (88 samples):", nrow(bim_full), "SNPs\n")
cat("  pcadapt BIM (81 samples):", nrow(bim_pcadapt), "SNPs\n")

# SNP universe intersection (should be identical since same VCFs, just different samples)
universe_intersect <- intersect(bim_full$snp_id, bim_pcadapt$snp_id)
cat("  SNP universe intersection:", length(universe_intersect), "\n\n")

################################################################################
# STEP 2: EXCLUDE SEX CHROMOSOME
################################################################################

cat("--- Step 2: Excluding sex chromosome (", SEX_CHROM, ") ---\n\n")

# Count sex chromosome outliers before exclusion
pct_sex <- sum(pct_summary$chrom == SEX_CHROM)
outflank_sex <- sum(outflank_snps$chrom == SEX_CHROM)
pcadapt_sex <- sum(pcadapt_q05$chrom == SEX_CHROM)

cat("  Sex chromosome outliers:\n")
cat("    Percentile:", pct_sex, "\n")
cat("    OutFLANK:", outflank_sex, "\n")
cat("    pcadapt:", pcadapt_sex, "\n")

# Save pcadapt sex-chromosome outliers separately (as noted in spec)
pcadapt_sex_outliers <- pcadapt_q05 %>% filter(chrom == SEX_CHROM)
if (nrow(pcadapt_sex_outliers) > 0) {
  write_tsv(pcadapt_sex_outliers, file.path(OUT_DIR, "pcadapt_sex_chrom_outliers.tsv"))
  cat("    pcadapt sex-chrom outliers saved separately\n")
}

# Exclude sex chromosome
pct_summary   <- pct_summary %>% filter(chrom != SEX_CHROM)
pct_detail    <- pct_detail %>% filter(chrom != SEX_CHROM)
outflank_snps <- outflank_snps %>% filter(chrom != SEX_CHROM)
pcadapt_q05   <- pcadapt_q05 %>% filter(chrom != SEX_CHROM)
pcadapt_q01   <- pcadapt_q01 %>% filter(chrom != SEX_CHROM)

# Also exclude from universe
universe_autosomal <- universe_intersect[!grepl(paste0("^", SEX_CHROM, ":"), universe_intersect)]

cat("\n  After exclusion (autosomes only):\n")
cat("    Percentile:", nrow(pct_summary), "SNPs\n")
cat("    OutFLANK:", nrow(outflank_snps), "SNPs\n")
cat("    pcadapt:", nrow(pcadapt_q05), "SNPs\n")
cat("    Autosomal universe:", length(universe_autosomal), "SNPs\n\n")

################################################################################
# STEP 3: BUILD MASTER CANDIDATE TABLE
################################################################################

cat("--- Step 3: Building master candidate table ---\n\n")

# --- 3a. Prepare per-method data frames with consistent columns ---

# Percentile: one row per unique SNP
pct_for_join <- pct_summary %>%
  select(snp_id, chrom, pos,
         percentile_comparisons = comparisons,
         percentile_fst_max = max_fst,
         percentile_pattern = pattern,
         percentile_n_comparisons = n_comparisons)

# Determine strong percentile support (top 0.1%)
# Need detail file for percentile_rank to check top 0.1%
pct_strong <- pct_detail %>%
  group_by(snp_id) %>%
  summarise(percentile_rank_min = min(percentile_rank), .groups = "drop") %>%
  mutate(percentile_strong = percentile_rank_min >= (1 - PERCENTILE_STRONG_PCT))

pct_for_join <- pct_for_join %>%
  left_join(pct_strong, by = "snp_id") %>%
  mutate(percentile_strong = coalesce(percentile_strong, FALSE))

# OutFLANK: one row per unique SNP
outflank_for_join <- outflank_snps %>%
  select(snp_id, chrom, pos,
         outflank_comparisons,
         outflank_q_min)

# pcadapt: one row per SNP
pcadapt_for_join <- pcadapt_q05 %>%
  select(snp_id, chrom, pos,
         pcadapt_q = qvalue,
         pcadapt_z_PC1 = z_PC1,
         pcadapt_z_PC2 = z_PC2,
         pcadapt_pc_category = pc_category)

# Flag strong pcadapt support
pcadapt_for_join <- pcadapt_for_join %>%
  mutate(pcadapt_strong = pcadapt_q < PCADAPT_STRONG_Q)

# --- 3b. Full outer join on snp_id ---
master <- pct_for_join %>%
  select(snp_id, chrom, pos) %>%
  full_join(outflank_for_join %>% select(snp_id, chrom, pos), by = c("snp_id", "chrom", "pos")) %>%
  full_join(pcadapt_for_join %>% select(snp_id, chrom, pos), by = c("snp_id", "chrom", "pos")) %>%
  distinct()

# Join method-specific columns
master <- master %>%
  left_join(pct_for_join, by = c("snp_id", "chrom", "pos")) %>%
  left_join(outflank_for_join, by = c("snp_id", "chrom", "pos")) %>%
  left_join(pcadapt_for_join, by = c("snp_id", "chrom", "pos"))

# --- 3c. Add outlier flags ---
master <- master %>%
  mutate(
    percentile_outlier = !is.na(percentile_comparisons),
    outflank_outlier   = !is.na(outflank_comparisons),
    pcadapt_outlier    = !is.na(pcadapt_q)
  )

# --- 3d. Count methods and classify method families ---
master <- master %>%
  mutate(
    n_methods = as.integer(percentile_outlier) + as.integer(outflank_outlier) + as.integer(pcadapt_outlier),
    method_families = case_when(
      percentile_outlier & pcadapt_outlier & outflank_outlier ~ "both",
      (percentile_outlier | outflank_outlier) & pcadapt_outlier ~ "both",
      percentile_outlier | outflank_outlier ~ "FST",
      pcadapt_outlier ~ "PCA",
      TRUE ~ NA_character_
    )
  )

# --- 3e. Determine methods_available ---
# POST_INVASION SNPs (among-invasive) have 3 methods available
# INVASION_SPECIFIC SNPs (CHI-vs-invasive) have 2 methods available
# Need to determine this from available pattern info
# A SNP's pattern context determines availability:
# - If any method classifies it as POST_INVASION or BOTH, OutFLANK could detect it -> 3
# - If only INVASION_SPECIFIC, OutFLANK cannot detect it -> 2
# Use the consensus pattern (computed below) to assign, but we need a preliminary pass

# Preliminary pattern from whatever methods flagged the SNP
master <- master %>%
  mutate(
    # Collect all pattern assignments
    pct_pat = ifelse(percentile_outlier, percentile_pattern, NA_character_),
    pca_pat = ifelse(pcadapt_outlier, pcadapt_pc_category, NA_character_),
    # OutFLANK is implicitly POST_INVASION
    of_pat  = ifelse(outflank_outlier, "POST_INVASION", NA_character_)
  )

# --- 3f. Pattern reconciliation ---
# Consensus pattern when methods agree; DISCORDANT when they disagree
master <- master %>%
  rowwise() %>%
  mutate(
    patterns_present = list(na.omit(c(pct_pat, pca_pat, of_pat))),
    unique_patterns  = list(unique(na.omit(c(pct_pat, pca_pat, of_pat))))
  ) %>%
  ungroup() %>%
  mutate(
    consensus_pattern = case_when(
      lengths(unique_patterns) == 0 ~ NA_character_,
      lengths(unique_patterns) == 1 ~ sapply(unique_patterns, `[`, 1),
      # BOTH is compatible with either specific pattern
      lengths(unique_patterns) == 2 & sapply(unique_patterns, function(x) "BOTH" %in% x) ~ {
        sapply(unique_patterns, function(x) {
          other <- setdiff(x, "BOTH")
          if (length(other) == 1) other else "BOTH"
        })
      },
      TRUE ~ "DISCORDANT"
    )
  )

# --- 3g. Assign methods_available ---
master <- master %>%
  mutate(
    methods_available = case_when(
      consensus_pattern == "INVASION_SPECIFIC" ~ 2L,
      consensus_pattern %in% c("POST_INVASION", "BOTH") ~ 3L,
      consensus_pattern == "DISCORDANT" ~ 3L,  # conservative: if discordant, some method saw among-invasive
      # For SNPs with no clear pattern (shouldn't happen), use 2 as conservative default
      TRUE ~ 2L
    )
  )

# --- 3h. Confidence tier assignment ---
master <- master %>%
  mutate(
    confidence_tier = case_when(
      # Tier 1: Cross-family concordance (PCA + FST-based)
      pcadapt_outlier & (percentile_outlier | outflank_outlier) ~ 1L,

      # Tier 2a: Two FST-based methods agree (percentile + OutFLANK)
      percentile_outlier & outflank_outlier ~ 2L,

      # Tier 2b: Single method with very strong support
      pcadapt_outlier & pcadapt_strong ~ 2L,
      percentile_outlier & percentile_strong ~ 2L,

      # Tier 3: Single method at standard threshold
      TRUE ~ 3L
    )
  )

# --- 3i. Clean up and select final columns ---
master_out <- master %>%
  select(
    chrom, pos, snp_id,
    percentile_outlier, percentile_comparisons, percentile_fst_max, percentile_pattern,
    outflank_outlier, outflank_comparisons, outflank_q_min,
    pcadapt_outlier, pcadapt_q, pcadapt_z_PC1, pcadapt_z_PC2, pcadapt_pc_category,
    n_methods, method_families, methods_available,
    confidence_tier, consensus_pattern
  ) %>%
  arrange(chrom, pos)

cat("Master candidate table built:\n")
cat("  Total unique outlier SNPs:", nrow(master_out), "\n")
cat("  By n_methods:\n")
print(table(master_out$n_methods))
cat("  By confidence tier:\n")
print(table(master_out$confidence_tier))
cat("  By consensus pattern:\n")
print(table(master_out$consensus_pattern, useNA = "ifany"))
cat("  By method families:\n")
print(table(master_out$method_families))
cat("\n")

# Write master candidate table
write_tsv(master_out, file.path(OUT_DIR, "master_candidate_table.tsv"))
cat("  Written to:", file.path(OUT_DIR, "master_candidate_table.tsv"), "\n\n")

################################################################################
# STEP 4: COMPUTE CONCORDANCE METRICS
################################################################################

cat("--- Step 4: Computing concordance metrics ---\n\n")

# Define outlier SNP sets (as character vectors of snp_id)
set_pct     <- pct_summary$snp_id
set_outflank <- outflank_snps$snp_id
set_pcadapt <- pcadapt_q05$snp_id

# Helper function for pairwise concordance
compute_concordance <- function(set_a, set_b, name_a, name_b) {
  overlap   <- length(intersect(set_a, set_b))
  union_ab  <- length(union(set_a, set_b))
  jaccard   <- if (union_ab > 0) overlap / union_ab else 0
  contain_a_in_b <- if (length(set_a) > 0) overlap / length(set_a) else 0
  contain_b_in_a <- if (length(set_b) > 0) overlap / length(set_b) else 0

  data.frame(
    method_A = name_a,
    method_B = name_b,
    n_A = length(set_a),
    n_B = length(set_b),
    overlap = overlap,
    union = union_ab,
    jaccard = round(jaccard, 4),
    pct_A_in_B = round(100 * contain_a_in_b, 2),
    pct_B_in_A = round(100 * contain_b_in_a, 2),
    stringsAsFactors = FALSE
  )
}

# Pairwise concordance
conc_pct_of  <- compute_concordance(set_pct, set_outflank, "Percentile", "OutFLANK")
conc_pct_pca <- compute_concordance(set_pct, set_pcadapt, "Percentile", "pcadapt")
conc_of_pca  <- compute_concordance(set_outflank, set_pcadapt, "OutFLANK", "pcadapt")

concordance_table <- bind_rows(conc_pct_of, conc_pct_pca, conc_of_pca)

cat("Pairwise concordance:\n")
print(concordance_table, row.names = FALSE)
cat("\n")

# Three-way overlap (POST_INVASION SNPs where all 3 methods can contribute)
three_way <- Reduce(intersect, list(set_pct, set_outflank, set_pcadapt))
cat("Three-way overlap (all methods):", length(three_way), "SNPs\n")

# Also compute: SNPs in exactly 2 of 3 methods, and exactly 1
two_of_three <- setdiff(
  union(intersect(set_pct, set_outflank),
        union(intersect(set_pct, set_pcadapt),
              intersect(set_outflank, set_pcadapt))),
  three_way
)
cat("Two-of-three overlap:", length(two_of_three), "SNPs\n")
cat("Single-method only:", nrow(master_out) - length(three_way) - length(two_of_three), "SNPs\n\n")

################################################################################
# STEP 5: WRITE CONCORDANCE SUMMARY
################################################################################

cat("--- Step 5: Writing concordance summary ---\n\n")

sink(file.path(OUT_DIR, "concordance_summary.txt"))

cat("========================================================================\n")
cat("CROSS-METHOD CONCORDANCE SUMMARY\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================================================\n\n")

cat("--- SNP Universe ---\n")
cat("Full BIM (88 samples):", nrow(bim_full), "SNPs\n")
cat("pcadapt BIM (81 samples):", nrow(bim_pcadapt), "SNPs\n")
cat("Universe intersection:", length(universe_intersect), "SNPs\n")
cat("Autosomal universe (excl.", SEX_CHROM, "):", length(universe_autosomal), "SNPs\n\n")

cat("--- Sex Chromosome Exclusion ---\n")
cat("Sex chromosome:", SEX_CHROM, "\n")
cat("Percentile outliers on sex chrom:", pct_sex, "\n")
cat("OutFLANK outliers on sex chrom:", outflank_sex, "\n")
cat("pcadapt outliers on sex chrom:", pcadapt_sex, "\n")
cat("(pcadapt sex-chrom outliers saved to pcadapt_sex_chrom_outliers.tsv)\n\n")

cat("--- Method Summary (autosomes only) ---\n")
cat("Per-SNP percentile (CHI-filtered, 1%):", nrow(pct_summary), "outlier SNPs\n")
cat("  Note: CHI genotype count filter applied (>=3 of 4 CHI individuals genotyped)\n")
cat("OutFLANK (q<0.05):", nrow(outflank_snps), "outlier SNPs\n")
cat("  EEU_vs_WEU:", sum(grepl("EEU_vs_WEU", outflank_snps$outflank_comparisons)), "\n")
cat("  EEU_vs_USA:", sum(grepl("EEU_vs_USA", outflank_snps$outflank_comparisons)), "\n")
cat("  (Only among-invasive comparisons produce OutFLANK outliers)\n")
cat("pcadapt (q<0.05):", nrow(pcadapt_q05), "outlier SNPs\n\n")

cat("--- Pairwise Concordance ---\n\n")
for (i in 1:nrow(concordance_table)) {
  row <- concordance_table[i, ]
  cat(sprintf("%s vs %s:\n", row$method_A, row$method_B))
  cat(sprintf("  %s: %d outliers | %s: %d outliers\n", row$method_A, row$n_A, row$method_B, row$n_B))
  cat(sprintf("  Overlap: %d SNPs\n", row$overlap))
  cat(sprintf("  Jaccard index: %.4f\n", row$jaccard))
  cat(sprintf("  Directional containment:\n"))
  cat(sprintf("    %.1f%% of %s outliers also in %s\n", row$pct_A_in_B, row$method_A, row$method_B))
  cat(sprintf("    %.1f%% of %s outliers also in %s\n", row$pct_B_in_A, row$method_B, row$method_A))
  cat("\n")
}

cat("--- Three-Way Overlap ---\n")
cat("SNPs detected by all 3 methods:", length(three_way), "\n")
cat("SNPs detected by exactly 2 methods:", length(two_of_three), "\n")
cat("SNPs detected by exactly 1 method:", nrow(master_out) - length(three_way) - length(two_of_three), "\n\n")

cat("--- Confidence Tier Distribution ---\n\n")
tier_table <- master_out %>%
  count(confidence_tier, consensus_pattern) %>%
  arrange(confidence_tier, consensus_pattern)
cat("Tier x Pattern:\n")
print(as.data.frame(tier_table), row.names = FALSE)
cat("\n")

tier_summary <- master_out %>% count(confidence_tier)
cat("Tier totals:\n")
for (i in 1:nrow(tier_summary)) {
  cat(sprintf("  Tier %d: %d SNPs (%.1f%%)\n",
              tier_summary$confidence_tier[i],
              tier_summary$n[i],
              100 * tier_summary$n[i] / nrow(master_out)))
}
cat("\n")

cat("--- Pattern Distribution ---\n\n")
pattern_summary <- master_out %>% count(consensus_pattern)
for (i in 1:nrow(pattern_summary)) {
  cat(sprintf("  %s: %d SNPs (%.1f%%)\n",
              pattern_summary$consensus_pattern[i],
              pattern_summary$n[i],
              100 * pattern_summary$n[i] / nrow(master_out)))
}
cat("\n")

# Discordant SNPs detail
discordant <- master_out %>% filter(consensus_pattern == "DISCORDANT")
if (nrow(discordant) > 0) {
  cat("--- Discordant Pattern SNPs ---\n")
  cat("Count:", nrow(discordant), "\n")
  cat("These SNPs were assigned different patterns by different methods.\n")
  cat("First 20:\n")
  print(as.data.frame(head(discordant %>%
    select(snp_id, percentile_pattern, pcadapt_pc_category, n_methods, confidence_tier), 20)),
    row.names = FALSE)
  cat("\n")
}

cat("--- Methods Available Distribution ---\n")
cat("  methods_available=2 (INVASION_SPECIFIC, no OutFLANK coverage):",
    sum(master_out$methods_available == 2), "\n")
cat("  methods_available=3 (POST_INVASION/BOTH, all methods can contribute):",
    sum(master_out$methods_available == 3), "\n\n")

cat("========================================================================\n")
cat("END OF CONCORDANCE SUMMARY\n")
cat("========================================================================\n")

sink()

cat("  Concordance summary written to:", file.path(OUT_DIR, "concordance_summary.txt"), "\n\n")

################################################################################
# STEP 6: VENN DIAGRAMS
################################################################################

cat("--- Step 6: Generating Venn diagrams ---\n\n")

# --- 6a. Three-way Venn for POST_INVASION SNPs ---
# (where all 3 methods can contribute)
# Filter to SNPs with POST_INVASION or BOTH pattern from any method
post_inv_pct <- pct_summary %>%
  filter(pattern %in% c("POST_INVASION", "BOTH")) %>%
  pull(snp_id)
post_inv_of <- outflank_snps$snp_id  # all OutFLANK outliers are implicitly POST_INVASION
post_inv_pca <- pcadapt_q05 %>%
  filter(pc_category %in% c("POST_INVASION", "BOTH")) %>%
  pull(snp_id)

cat("POST_INVASION Venn (3-way):\n")
cat("  Percentile (POST_INVASION/BOTH):", length(post_inv_pct), "\n")
cat("  OutFLANK:", length(post_inv_of), "\n")
cat("  pcadapt (POST_INVASION/BOTH):", length(post_inv_pca), "\n")

# Generate 3-way Venn (use VennDiagram::venn.diagram for file-based output)
venn_3way_file <- file.path(PLOT_DIR, "venn_post_invasion.png")

# Suppress the default log file
invisible(venn.diagram(
  x = list(
    Percentile = post_inv_pct,
    OutFLANK = post_inv_of,
    pcadapt = post_inv_pca
  ),
  filename = venn_3way_file,
  imagetype = "png",
  height = 2400, width = 2400, resolution = 300,
  fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
  alpha = 0.4,
  cex = 1.2,
  cat.cex = 1.1,
  cat.fontface = "bold",
  main = "POST_INVASION Outlier Concordance (3-way)",
  main.cex = 1.3
))
cat("  3-way Venn saved to:", venn_3way_file, "\n")

# --- 6b. Two-way Venn for INVASION_SPECIFIC SNPs ---
inv_spec_pct <- pct_summary %>%
  filter(pattern %in% c("INVASION_SPECIFIC", "BOTH")) %>%
  pull(snp_id)
inv_spec_pca <- pcadapt_q05 %>%
  filter(pc_category %in% c("INVASION_SPECIFIC", "BOTH")) %>%
  pull(snp_id)

cat("\nINVASION_SPECIFIC Venn (2-way):\n")
cat("  Percentile (INVASION_SPECIFIC/BOTH):", length(inv_spec_pct), "\n")
cat("  pcadapt (INVASION_SPECIFIC/BOTH):", length(inv_spec_pca), "\n")

venn_2way_file <- file.path(PLOT_DIR, "venn_invasion_specific.png")

invisible(venn.diagram(
  x = list(
    Percentile = inv_spec_pct,
    pcadapt = inv_spec_pca
  ),
  filename = venn_2way_file,
  imagetype = "png",
  height = 1800, width = 2400, resolution = 300,
  fill = c("#E41A1C", "#4DAF4A"),
  alpha = 0.4,
  cex = 1.3,
  cat.cex = 1.2,
  cat.fontface = "bold",
  main = "INVASION_SPECIFIC Outlier Concordance (2-way)",
  main.cex = 1.3
))
cat("  2-way Venn saved to:", venn_2way_file, "\n\n")

# --- 6c. Compute and save Venn partition counts ---
# These are the exact numbers shown in each Venn region, saved as text
# so they are readable without viewing the PNG images.

cat("Computing Venn partition counts...\n")

# POST_INVASION 3-way partitions
pi_pct_only  <- setdiff(setdiff(post_inv_pct, post_inv_of), post_inv_pca)
pi_of_only   <- setdiff(setdiff(post_inv_of, post_inv_pct), post_inv_pca)
pi_pca_only  <- setdiff(setdiff(post_inv_pca, post_inv_pct), post_inv_of)
pi_pct_of    <- setdiff(intersect(post_inv_pct, post_inv_of), post_inv_pca)
pi_pct_pca   <- setdiff(intersect(post_inv_pct, post_inv_pca), post_inv_of)
pi_of_pca    <- setdiff(intersect(post_inv_of, post_inv_pca), post_inv_pct)
pi_all_three <- Reduce(intersect, list(post_inv_pct, post_inv_of, post_inv_pca))

# INVASION_SPECIFIC 2-way partitions
is_pct_only <- setdiff(inv_spec_pct, inv_spec_pca)
is_pca_only <- setdiff(inv_spec_pca, inv_spec_pct)
is_both     <- intersect(inv_spec_pct, inv_spec_pca)

# ALL METHODS 3-way partitions (full outlier sets, not pattern-filtered)
all_pct_only  <- setdiff(setdiff(set_pct, set_outflank), set_pcadapt)
all_of_only   <- setdiff(setdiff(set_outflank, set_pct), set_pcadapt)
all_pca_only  <- setdiff(setdiff(set_pcadapt, set_pct), set_outflank)
all_pct_of    <- setdiff(intersect(set_pct, set_outflank), set_pcadapt)
all_pct_pca   <- setdiff(intersect(set_pct, set_pcadapt), set_outflank)
all_of_pca    <- setdiff(intersect(set_outflank, set_pcadapt), set_pct)
all_all_three <- Reduce(intersect, list(set_pct, set_outflank, set_pcadapt))

# Append to concordance summary
sink(file.path(OUT_DIR, "concordance_summary.txt"), append = TRUE)

cat("\n\n--- Venn Diagram Partition Counts ---\n\n")

cat("POST_INVASION 3-way Venn (pattern-filtered: POST_INVASION/BOTH from each method):\n")
cat("  Percentile only:", length(pi_pct_only), "\n")
cat("  OutFLANK only:", length(pi_of_only), "\n")
cat("  pcadapt only:", length(pi_pca_only), "\n")
cat("  Percentile + OutFLANK only:", length(pi_pct_of), "\n")
cat("  Percentile + pcadapt only:", length(pi_pct_pca), "\n")
cat("  OutFLANK + pcadapt only:", length(pi_of_pca), "\n")
cat("  All three methods:", length(pi_all_three), "\n")
cat("  Total:", length(pi_pct_only) + length(pi_of_only) + length(pi_pca_only) +
    length(pi_pct_of) + length(pi_pct_pca) + length(pi_of_pca) + length(pi_all_three), "\n\n")

cat("INVASION_SPECIFIC 2-way Venn (pattern-filtered: INVASION_SPECIFIC/BOTH from each method):\n")
cat("  Percentile only:", length(is_pct_only), "\n")
cat("  pcadapt only:", length(is_pca_only), "\n")
cat("  Both methods:", length(is_both), "\n")
cat("  Total:", length(is_pct_only) + length(is_pca_only) + length(is_both), "\n\n")

cat("ALL METHODS 3-way partitions (full outlier sets, no pattern filter):\n")
cat("  Percentile only:", length(all_pct_only), "\n")
cat("  OutFLANK only:", length(all_of_only), "\n")
cat("  pcadapt only:", length(all_pca_only), "\n")
cat("  Percentile + OutFLANK only:", length(all_pct_of), "\n")
cat("  Percentile + pcadapt only:", length(all_pct_pca), "\n")
cat("  OutFLANK + pcadapt only:", length(all_of_pca), "\n")
cat("  All three methods:", length(all_all_three), "\n")
cat("  Total:", length(all_pct_only) + length(all_of_only) + length(all_pca_only) +
    length(all_pct_of) + length(all_pct_pca) + length(all_of_pca) + length(all_all_three), "\n")

sink()

cat("  Venn partitions appended to concordance_summary.txt\n\n")

################################################################################
# STEP 7: MANHATTAN PLOT HIGHLIGHTING TIER 1 SNPS
################################################################################

cat("--- Step 7: Manhattan plot (Tier 1 highlights) ---\n\n")

# Use pcadapt -log10(q) as the y-axis for all SNPs, highlight Tier 1
# Need chromosome ordering for Manhattan layout
chrom_order <- sort(unique(master_out$chrom))
chrom_numeric <- setNames(seq_along(chrom_order), chrom_order)

# Build Manhattan data from pcadapt (largest SNP set for background)
manhattan_data <- pcadapt_q05 %>%
  select(chrom, pos, snp_id, qvalue) %>%
  mutate(
    chrom_num = chrom_numeric[chrom],
    neg_log_q = -log10(pmax(qvalue, .Machine$double.xmin))
  )

# Assign cumulative positions
chrom_lengths <- manhattan_data %>%
  group_by(chrom_num) %>%
  summarise(max_pos = max(pos), .groups = "drop") %>%
  arrange(chrom_num) %>%
  mutate(offset = cumsum(lag(max_pos, default = 0)))

manhattan_data <- manhattan_data %>%
  left_join(chrom_lengths %>% select(chrom_num, offset), by = "chrom_num") %>%
  mutate(cumulative_pos = pos + offset)

# Identify Tier 1 SNPs
tier1_snps <- master_out %>% filter(confidence_tier == 1) %>% pull(snp_id)
manhattan_data <- manhattan_data %>%
  mutate(tier = ifelse(snp_id %in% tier1_snps, "Tier 1", "Other outlier"))

# Axis labels
axis_labels <- manhattan_data %>%
  group_by(chrom_num, chrom) %>%
  summarise(center = median(cumulative_pos), .groups = "drop") %>%
  mutate(label = gsub("NC_0581", "", chrom))  # shorten chromosome names

p_manhattan <- ggplot(manhattan_data, aes(x = cumulative_pos, y = neg_log_q)) +
  geom_point(data = manhattan_data %>% filter(tier == "Other outlier"),
             aes(color = factor(chrom_num %% 2)),
             size = 0.5, alpha = 0.3) +
  geom_point(data = manhattan_data %>% filter(tier == "Tier 1"),
             color = "red", size = 1.2, alpha = 0.8) +
  scale_color_manual(values = c("0" = "grey50", "1" = "grey70"), guide = "none") +
  scale_x_continuous(breaks = axis_labels$center, labels = axis_labels$label) +
  labs(
    title = "pcadapt Outliers with Tier 1 Concordant SNPs Highlighted",
    subtitle = paste0("Tier 1 (cross-family concordance): ", length(tier1_snps),
                      " SNPs (red) out of ", nrow(pcadapt_q05), " pcadapt outliers"),
    x = "Chromosome",
    y = expression(-log[10](q))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

manhattan_file <- file.path(PLOT_DIR, "manhattan_tier1.png")
ggsave(manhattan_file, p_manhattan, width = 14, height = 5, dpi = 300)
cat("  Manhattan plot saved to:", manhattan_file, "\n\n")

################################################################################
# FINAL SUMMARY
################################################################################

cat("========================================================================\n")
cat("CONCORDANCE ANALYSIS COMPLETE\n")
cat("========================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("Output files:\n")
cat("  ", file.path(OUT_DIR, "master_candidate_table.tsv"), "\n")
cat("  ", file.path(OUT_DIR, "concordance_summary.txt"), "\n")
cat("  ", file.path(OUT_DIR, "pcadapt_sex_chrom_outliers.tsv"), "\n")
cat("  ", venn_3way_file, "\n")
cat("  ", venn_2way_file, "\n")
cat("  ", manhattan_file, "\n\n")

cat("Key results:\n")
cat("  Total unique outlier SNPs:", nrow(master_out), "\n")
cat("  Tier 1 (cross-family):", sum(master_out$confidence_tier == 1), "\n")
cat("  Tier 2 (same-family or strong single):", sum(master_out$confidence_tier == 2), "\n")
cat("  Tier 3 (single method, standard):", sum(master_out$confidence_tier == 3), "\n")
n_discordant <- sum(master_out$consensus_pattern == "DISCORDANT", na.rm = TRUE)
if (n_discordant > 0) {
  cat("  Discordant pattern SNPs:", n_discordant, "(flagged for review)\n")
}
cat("\n========================================================================\n")
