#!/usr/bin/env Rscript
################################################################################
# title: per_snp_selection_scan.R
# project: BIOL624 Final Project - Selection Detection in Coccinella septempunctata
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-04-13
# last modified: 2026-04-13
#
# purpose:
#   Per-SNP percentile-based selection scan across all six pairwise population
#   comparisons. Replaces the windowed unified_selection_scan.R approach to put
#   FST-based detection on the same SNP-level scale as OutFLANK and pcadapt,
#   enabling meaningful cross-method concordance.
#
#   For each comparison, applies genome-wide percentile thresholds (top 0.5%,
#   1%, 2%) to per-SNP Weir & Cockerham FST values. Categorizes outlier SNPs
#   by biological pattern (BOTH, INVASION_SPECIFIC, POST_INVASION) and assigns
#   confidence tiers. Clusters nearby outlier SNPs into genomic regions for
#   downstream gene annotation.
#
# inputs:
#   - Per-SNP FST files: ../results/fst/{POP1}_vs_{POP2}/fst_genomewide/*.weir.fst
#     Schema: CHROM, POS, WEIR_AND_COCKERHAM_FST
#   - OutFLANK outlier files (optional, for cross-referencing):
#     ../results/outflank_diagnostic/{comparison}/outliers_q05.tsv
#
# outputs:
#   ../results/per_snp_selection_scan/
#     - outliers_0.5pct.tsv, outliers_1pct.tsv, outliers_2pct.tsv
#     - pattern_categorization/both_pattern_snps.tsv
#     - pattern_categorization/invasion_specific_snps.tsv
#     - pattern_categorization/post_invasion_snps.tsv
#     - clustered_regions/outlier_regions_1pct.tsv
#     - plots/manhattan_{comparison}.png
#     - plots/manhattan_all_comparisons.png
#     - summary_report.txt
#
# usage:
#   cd scripts
#   Rscript per_snp_selection_scan.R
#   Rscript per_snp_selection_scan.R --cluster-dist 5000   # custom cluster distance
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
})

# --- Argument Parsing --------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# Default: 10kb clustering distance for post-hoc region merging
CLUSTER_DIST <- 10000

if (length(args) > 0) {
  if ("--cluster-dist" %in% args) {
    idx <- which(args == "--cluster-dist")
    if (idx < length(args)) {
      CLUSTER_DIST <- as.integer(args[idx + 1])
    }
  }
  if ("--help" %in% args || "-h" %in% args) {
    cat(
      "Usage: Rscript per_snp_selection_scan.R [--cluster-dist N]\n\n",
      "Options:\n",
      "  --cluster-dist N   Distance (bp) for merging nearby outlier SNPs into\n",
      "                     regions (default: 10000)\n",
      "  --help, -h         Show this help\n",
      sep = ""
    )
    quit(save = "no", status = 0)
  }
}

cat("\n================================================================================\n")
cat("PER-SNP PERCENTILE SELECTION SCAN\n")
cat("Species: Coccinella septempunctata (seven-spotted lady beetle)\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Cluster distance:", CLUSTER_DIST, "bp\n")
cat("================================================================================\n\n")

# --- Configuration -----------------------------------------------------------

comparisons <- c(
  "CHI_vs_EEU", "CHI_vs_USA", "CHI_vs_WEU",
  "EEU_vs_USA", "EEU_vs_WEU", "USA_vs_WEU"
)

chi_comparisons <- c("CHI_vs_EEU", "CHI_vs_USA", "CHI_vs_WEU")
invasive_comparisons <- c("EEU_vs_USA", "EEU_vs_WEU", "USA_vs_WEU")

percentile_thresholds <- c(0.995, 0.99, 0.98)
threshold_labels <- c("0.5pct", "1pct", "2pct")

FST_BASE <- "../results/fst"
OUTFLANK_DIR <- "../results/outflank_diagnostic"
OUTPUT_DIR <- "../results/per_snp_selection_scan"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "pattern_categorization"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "clustered_regions"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE)

################################################################################
# Step 1: Load per-SNP FST data for all comparisons
################################################################################

cat("STEP 1: Loading per-SNP FST data\n")
cat("------------------------------------------------------------------------\n\n")

all_fst <- list()
comparison_stats <- tibble(
  comparison = character(),
  n_snps = integer(),
  mean_fst = double(),
  median_fst = double(),
  max_fst = double(),
  n_negative = integer(),
  pct_negative = double()
)

for (comp in comparisons) {
  gw_dir <- file.path(FST_BASE, comp, "fst_genomewide")

  if (!dir.exists(gw_dir)) {
    cat("  WARNING: Directory not found:", gw_dir, "\n")
    next
  }

  fst_files <- list.files(gw_dir, pattern = "\\.weir\\.fst$", full.names = TRUE)

  if (length(fst_files) == 0) {
    cat("  WARNING: No .weir.fst files in:", gw_dir, "\n")
    next
  }

  # Read and combine all chromosomes
  fst_combined <- map_dfr(fst_files, function(f) {
    read_tsv(f, col_types = cols(
      CHROM = col_character(),
      POS = col_integer(),
      WEIR_AND_COCKERHAM_FST = col_double()
    ))
  })

  # Remove rows where FST is NA (vcftools outputs nan for monomorphic sites)
  n_before <- nrow(fst_combined)
  fst_combined <- fst_combined %>% filter(!is.na(WEIR_AND_COCKERHAM_FST))
  n_removed <- n_before - nrow(fst_combined)

  # Add SNP ID and comparison label
  fst_combined <- fst_combined %>%
    mutate(
      snp_id = paste0(CHROM, ":", POS),
      comparison = comp
    )

  all_fst[[comp]] <- fst_combined

  # Collect stats
  n_neg <- sum(fst_combined$WEIR_AND_COCKERHAM_FST < 0)
  comparison_stats <- bind_rows(comparison_stats, tibble(
    comparison = comp,
    n_snps = nrow(fst_combined),
    mean_fst = mean(fst_combined$WEIR_AND_COCKERHAM_FST),
    median_fst = median(fst_combined$WEIR_AND_COCKERHAM_FST),
    max_fst = max(fst_combined$WEIR_AND_COCKERHAM_FST),
    n_negative = n_neg,
    pct_negative = 100 * n_neg / nrow(fst_combined)
  ))

  cat(sprintf("  %-12s  %7d SNPs  (mean=%.4f, median=%.4f, max=%.4f, %d NA removed)\n",
              comp, nrow(fst_combined),
              mean(fst_combined$WEIR_AND_COCKERHAM_FST),
              median(fst_combined$WEIR_AND_COCKERHAM_FST),
              max(fst_combined$WEIR_AND_COCKERHAM_FST),
              n_removed))
}

cat("\nComparison summary:\n")
print(as.data.frame(comparison_stats %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))
cat("\n")

if (length(all_fst) == 0) {
  stop("No FST data loaded. Check that results/fst/*/fst_genomewide/ directories exist.")
}

################################################################################
# Step 2: Apply genome-wide percentile thresholds
################################################################################

cat("\nSTEP 2: Applying genome-wide percentile thresholds\n")
cat("------------------------------------------------------------------------\n\n")

all_outliers <- list()  # keyed by threshold label

for (t_idx in seq_along(percentile_thresholds)) {
  threshold <- percentile_thresholds[t_idx]
  label <- threshold_labels[t_idx]

  cat(sprintf("--- Threshold: %s (top %.1f%%) ---\n", label, 100 * (1 - threshold)))

  outliers_list <- list()

  for (comp in names(all_fst)) {
    fst_data <- all_fst[[comp]]

    cutoff <- quantile(fst_data$WEIR_AND_COCKERHAM_FST, threshold)

    outliers <- fst_data %>%
      filter(WEIR_AND_COCKERHAM_FST >= cutoff) %>%
      mutate(
        threshold = label,
        fst_cutoff = cutoff,
        percentile_rank = percent_rank(WEIR_AND_COCKERHAM_FST)
      )

    # Recalculate percentile_rank from the full distribution, not just outliers
    full_ecdf <- ecdf(fst_data$WEIR_AND_COCKERHAM_FST)
    outliers$percentile_rank <- full_ecdf(outliers$WEIR_AND_COCKERHAM_FST)

    outliers_list[[comp]] <- outliers

    cat(sprintf("  %-12s  cutoff=%.4f  outliers=%d\n",
                comp, cutoff, nrow(outliers)))
  }

  combined <- bind_rows(outliers_list)
  all_outliers[[label]] <- combined

  # Save outlier table
  out_file <- file.path(OUTPUT_DIR, paste0("outliers_", label, ".tsv"))
  write_tsv(combined %>% select(
    chrom = CHROM, pos = POS, snp_id, comparison, fst = WEIR_AND_COCKERHAM_FST,
    percentile_rank, threshold, fst_cutoff
  ), out_file)
  cat(sprintf("  Saved: %s (%d total outlier SNP-comparison pairs)\n\n",
              out_file, nrow(combined)))
}

################################################################################
# Step 3: Pattern categorization at the SNP level
################################################################################

cat("\nSTEP 3: Categorizing outlier SNPs by biological pattern\n")
cat("------------------------------------------------------------------------\n\n")

# Use top 1% as the primary threshold for pattern categorization
pattern_label <- "1pct"
pattern_data <- all_outliers[[pattern_label]]

cat("Using", pattern_label, "threshold for pattern categorization\n\n")

# Summarize each SNP across comparisons
snp_summary <- pattern_data %>%
  mutate(
    comp_type = if_else(comparison %in% chi_comparisons,
                        "native_vs_invasive", "among_invasive")
  ) %>%
  group_by(snp_id, CHROM, POS) %>%
  summarize(
    n_comparisons = n(),
    comparisons = paste(sort(unique(comparison)), collapse = ";"),
    mean_fst = mean(WEIR_AND_COCKERHAM_FST),
    max_fst = max(WEIR_AND_COCKERHAM_FST),
    in_native_vs_invasive = any(comp_type == "native_vs_invasive"),
    in_among_invasive = any(comp_type == "among_invasive"),
    n_native = sum(comp_type == "native_vs_invasive"),
    n_among = sum(comp_type == "among_invasive"),
    .groups = "drop"
  ) %>%
  mutate(
    pattern = case_when(
      in_native_vs_invasive & in_among_invasive ~ "BOTH",
      in_native_vs_invasive & !in_among_invasive ~ "INVASION_SPECIFIC",
      !in_native_vs_invasive & in_among_invasive ~ "POST_INVASION"
    ),
    confidence = case_when(
      n_comparisons >= 4 ~ "Very_High",
      n_comparisons == 3 ~ "High",
      n_comparisons == 2 ~ "Moderate",
      TRUE ~ "Single"
    )
  )

# Pattern summary
pattern_counts <- snp_summary %>%
  count(pattern, name = "n_snps") %>%
  arrange(desc(n_snps))

confidence_counts <- snp_summary %>%
  count(confidence, name = "n_snps") %>%
  arrange(match(confidence, c("Very_High", "High", "Moderate", "Single")))

cat("Pattern distribution:\n")
print(as.data.frame(pattern_counts))
cat("\nConfidence distribution:\n")
print(as.data.frame(confidence_counts))
cat("\n")

# Cross-tabulation: pattern x confidence
cat("Pattern x Confidence:\n")
print(as.data.frame(
  snp_summary %>% count(pattern, confidence) %>% pivot_wider(names_from = confidence, values_from = n, values_fill = 0)
))
cat("\n")

# Save pattern files
for (pat in c("BOTH", "INVASION_SPECIFIC", "POST_INVASION")) {
  pat_data <- snp_summary %>% filter(pattern == pat)
  out_file <- file.path(OUTPUT_DIR, "pattern_categorization",
                        paste0(tolower(pat), "_snps.tsv"))
  write_tsv(pat_data, out_file)
  cat("  Saved:", out_file, "(", nrow(pat_data), "SNPs)\n")
}

# Save the full SNP summary (all patterns combined)
write_tsv(snp_summary, file.path(OUTPUT_DIR, "snp_summary_1pct.tsv"))
cat("  Saved: snp_summary_1pct.tsv (", nrow(snp_summary), "unique SNPs)\n\n")

################################################################################
# Step 4: Post-hoc clustering of outlier SNPs into regions
################################################################################

cat("\nSTEP 4: Clustering nearby outlier SNPs into regions\n")
cat("------------------------------------------------------------------------\n")
cat("  Merge distance:", CLUSTER_DIST, "bp\n\n")

# Cluster SNPs within CLUSTER_DIST bp on the same chromosome
cluster_snps <- function(df, max_dist) {
  if (nrow(df) == 0) return(tibble())

  df <- df %>% arrange(CHROM, POS)

  regions <- list()
  region_id <- 1
  current_chrom <- df$CHROM[1]
  current_start <- df$POS[1]
  current_end <- df$POS[1]
  current_snps <- 1

  for (i in seq_len(nrow(df))) {
    if (i == 1) next

    same_chrom <- df$CHROM[i] == current_chrom
    within_dist <- (df$POS[i] - current_end) <= max_dist

    if (same_chrom && within_dist) {
      # Extend current region
      current_end <- df$POS[i]
      current_snps <- current_snps + 1
    } else {
      # Close current region, start new one
      regions[[region_id]] <- tibble(
        region_id = region_id,
        chrom = current_chrom,
        start = current_start,
        end = current_end,
        n_snps = current_snps,
        span_bp = current_end - current_start
      )
      region_id <- region_id + 1
      current_chrom <- df$CHROM[i]
      current_start <- df$POS[i]
      current_end <- df$POS[i]
      current_snps <- 1
    }
  }

  # Close final region
  regions[[region_id]] <- tibble(
    region_id = region_id,
    chrom = current_chrom,
    start = current_start,
    end = current_end,
    n_snps = current_snps,
    span_bp = current_end - current_start
  )

  bind_rows(regions)
}

clustered <- cluster_snps(snp_summary, CLUSTER_DIST)

# Annotate regions with pattern and confidence info by joining back to SNPs
if (nrow(clustered) > 0) {
  # For each region, find which SNPs fall in it and summarize
  region_annotations <- map_dfr(seq_len(nrow(clustered)), function(r) {
    reg <- clustered[r, ]
    snps_in_region <- snp_summary %>%
      filter(CHROM == reg$chrom, POS >= reg$start, POS <= reg$end)

    tibble(
      region_id = reg$region_id,
      patterns = paste(sort(unique(snps_in_region$pattern)), collapse = ";"),
      max_confidence = snps_in_region$confidence[which.max(match(
        snps_in_region$confidence,
        c("Very_High", "High", "Moderate", "Single")
      ))],
      mean_fst = mean(snps_in_region$mean_fst),
      max_fst = max(snps_in_region$max_fst),
      all_comparisons = paste(sort(unique(unlist(
        strsplit(snps_in_region$comparisons, ";")
      ))), collapse = ";")
    )
  })

  clustered <- left_join(clustered, region_annotations, by = "region_id")

  cat("  Total regions:", nrow(clustered), "\n")
  cat("  Singleton SNPs (1-SNP regions):", sum(clustered$n_snps == 1), "\n")
  cat("  Multi-SNP regions:", sum(clustered$n_snps > 1), "\n")
  cat("  Largest region:", max(clustered$n_snps), "SNPs,",
      max(clustered$span_bp), "bp span\n\n")

  cat("  Region size distribution:\n")
  print(as.data.frame(
    clustered %>%
      mutate(size_class = case_when(
        n_snps == 1 ~ "1 SNP",
        n_snps <= 3 ~ "2-3 SNPs",
        n_snps <= 10 ~ "4-10 SNPs",
        TRUE ~ ">10 SNPs"
      )) %>%
      count(size_class, name = "n_regions") %>%
      arrange(match(size_class, c("1 SNP", "2-3 SNPs", "4-10 SNPs", ">10 SNPs")))
  ))
  cat("\n")

  out_file <- file.path(OUTPUT_DIR, "clustered_regions",
                        paste0("outlier_regions_", pattern_label, ".tsv"))
  write_tsv(clustered, out_file)
  cat("  Saved:", out_file, "\n\n")
} else {
  cat("  No regions to cluster.\n\n")
}

################################################################################
# Step 5: Manhattan plots
################################################################################

cat("\nSTEP 5: Generating Manhattan plots\n")
cat("------------------------------------------------------------------------\n\n")

# Define chromosome order for Manhattan x-axis
chrom_order <- paste0("NC_0581", 89:98, ".1")

# Assign alternating colors and cumulative positions for Manhattan layout
# We need a genome-wide position offset per chromosome
chrom_lengths <- map_dfr(names(all_fst), function(comp) {
  all_fst[[comp]] %>%
    group_by(CHROM) %>%
    summarize(max_pos = max(POS), .groups = "drop")
}) %>%
  group_by(CHROM) %>%
  summarize(max_pos = max(max_pos), .groups = "drop") %>%
  arrange(match(CHROM, chrom_order))

chrom_lengths <- chrom_lengths %>%
  mutate(
    offset = cumsum(lag(max_pos, default = 0)),
    mid = offset + max_pos / 2
  )

# --- Per-comparison Manhattan plots ---

for (comp in names(all_fst)) {
  cat("  Plotting:", comp, "...")

  fst_data <- all_fst[[comp]] %>%
    left_join(chrom_lengths %>% select(CHROM, offset), by = "CHROM") %>%
    mutate(genome_pos = POS + offset)

  # Mark outliers at 1% threshold
  outlier_snps <- all_outliers[["1pct"]] %>%
    filter(comparison == comp) %>%
    pull(snp_id)

  fst_data <- fst_data %>%
    mutate(is_outlier = snp_id %in% outlier_snps)

  # Chromosome alternating color index
  fst_data <- fst_data %>%
    mutate(chrom_idx = match(CHROM, chrom_order) %% 2)

  # Subsample non-outliers for plotting speed (keep all outliers)
  n_total <- nrow(fst_data)
  max_plot_points <- 200000
  if (n_total > max_plot_points) {
    non_outlier <- fst_data %>% filter(!is_outlier)
    outlier <- fst_data %>% filter(is_outlier)
    frac <- max_plot_points / nrow(non_outlier)
    non_outlier <- non_outlier %>% slice_sample(prop = min(frac, 1))
    plot_data <- bind_rows(non_outlier, outlier) %>% arrange(genome_pos)
  } else {
    plot_data <- fst_data
  }

  # Determine comparison type for color
  is_chi_comp <- comp %in% chi_comparisons
  outlier_color <- if (is_chi_comp) "firebrick" else "dodgerblue3"

  p <- ggplot() +
    # Non-outliers: alternating gray shades
    geom_point(
      data = plot_data %>% filter(!is_outlier),
      aes(x = genome_pos, y = WEIR_AND_COCKERHAM_FST,
          color = factor(chrom_idx)),
      size = 0.3, alpha = 0.3
    ) +
    scale_color_manual(values = c("0" = "gray60", "1" = "gray30"), guide = "none") +
    # Outliers on top
    geom_point(
      data = plot_data %>% filter(is_outlier),
      aes(x = genome_pos, y = WEIR_AND_COCKERHAM_FST),
      color = outlier_color, size = 0.6, alpha = 0.7
    ) +
    # Chromosome labels
    scale_x_continuous(
      breaks = chrom_lengths$mid,
      labels = gsub("NC_0581", "", chrom_lengths$CHROM) %>% gsub("\\.1$", "", .),
      expand = c(0.01, 0)
    ) +
    labs(
      title = paste("Per-SNP FST -", comp),
      subtitle = paste0("Top 1% outliers highlighted (n=",
                        sum(fst_data$is_outlier), ")"),
      x = "Chromosome",
      y = expression(F[ST]~"(Weir & Cockerham)")
    ) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 8)
    )

  ggsave(file.path(OUTPUT_DIR, "plots", paste0("manhattan_", comp, ".png")),
         p, width = 10, height = 4, dpi = 300)
  cat(" done\n")
}

# --- Combined faceted Manhattan (all 6 comparisons) ---

cat("\n  Plotting combined Manhattan (all comparisons)...")

# Build a combined dataset with subsampling
combined_plot_data <- map_dfr(names(all_fst), function(comp) {
  fst_data <- all_fst[[comp]] %>%
    left_join(chrom_lengths %>% select(CHROM, offset), by = "CHROM") %>%
    mutate(genome_pos = POS + offset)

  outlier_snps <- all_outliers[["1pct"]] %>%
    filter(comparison == comp) %>%
    pull(snp_id)

  fst_data <- fst_data %>%
    mutate(
      is_outlier = snp_id %in% outlier_snps,
      chrom_idx = match(CHROM, chrom_order) %% 2
    )

  # Subsample per comparison for combined plot
  non_outlier <- fst_data %>% filter(!is_outlier)
  outlier <- fst_data %>% filter(is_outlier)
  max_per_comp <- 80000
  if (nrow(non_outlier) > max_per_comp) {
    non_outlier <- non_outlier %>% slice_sample(n = max_per_comp)
  }
  bind_rows(non_outlier, outlier)
})

# Color outliers by comparison type
combined_plot_data <- combined_plot_data %>%
  mutate(
    comp_type = if_else(comparison %in% chi_comparisons,
                        "Native vs Invasive", "Among Invasive"),
    point_color = case_when(
      is_outlier & comp_type == "Native vs Invasive" ~ "chi_outlier",
      is_outlier & comp_type == "Among Invasive" ~ "inv_outlier",
      TRUE ~ as.character(chrom_idx)
    )
  )

p_combined <- ggplot() +
  geom_point(
    data = combined_plot_data %>% filter(!is_outlier),
    aes(x = genome_pos, y = WEIR_AND_COCKERHAM_FST,
        color = factor(chrom_idx)),
    size = 0.15, alpha = 0.2
  ) +
  scale_color_manual(
    values = c("0" = "gray60", "1" = "gray30"),
    guide = "none"
  ) +
  ggnewscale::new_scale_color() +
  geom_point(
    data = combined_plot_data %>% filter(is_outlier),
    aes(x = genome_pos, y = WEIR_AND_COCKERHAM_FST, color = comp_type),
    size = 0.4, alpha = 0.6
  ) +
  scale_color_manual(
    values = c("Native vs Invasive" = "firebrick",
               "Among Invasive" = "dodgerblue3"),
    name = "Outlier type"
  ) +
  scale_x_continuous(
    breaks = chrom_lengths$mid,
    labels = gsub("NC_0581", "", chrom_lengths$CHROM) %>% gsub("\\.1$", "", .),
    expand = c(0.01, 0)
  ) +
  facet_wrap(~ comparison, ncol = 2, scales = "free_y") +
  labs(
    title = "Per-SNP FST across all pairwise comparisons",
    subtitle = "Top 1% outliers highlighted",
    x = "Chromosome",
    y = expression(F[ST])
  ) +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 6),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "plots", "manhattan_all_comparisons.png"),
       p_combined, width = 12, height = 10, dpi = 300)
cat(" done\n\n")

################################################################################
# Step 6: Cross-reference with OutFLANK (among-invasive only)
################################################################################

cat("STEP 6: Cross-referencing with OutFLANK outliers\n")
cat("------------------------------------------------------------------------\n\n")

outflank_available <- FALSE

for (comp in invasive_comparisons) {
  of_file <- file.path(OUTFLANK_DIR, comp, "outliers_q05.tsv")

  if (!file.exists(of_file)) {
    cat("  ", comp, ": OutFLANK file not found, skipping\n")
    next
  }

  of_data <- read_tsv(of_file, col_types = cols(.default = col_character())) %>%
    mutate(
      of_snp_id = LocusName,
      of_fst = as.double(FST),
      of_qvalue = as.double(qvalues)
    ) %>%
    select(of_snp_id, of_fst, of_qvalue)

  if (nrow(of_data) == 0) {
    cat("  ", comp, ": 0 OutFLANK outliers\n")
    next
  }

  outflank_available <- TRUE

  # Compare at each threshold
  for (label in threshold_labels) {
    pct_outliers <- all_outliers[[label]] %>%
      filter(comparison == comp) %>%
      pull(snp_id)

    overlap <- intersect(pct_outliers, of_data$of_snp_id)
    pct_only <- setdiff(pct_outliers, of_data$of_snp_id)
    of_only <- setdiff(of_data$of_snp_id, pct_outliers)
    union_size <- length(union(pct_outliers, of_data$of_snp_id))

    jaccard <- if (union_size > 0) length(overlap) / union_size else NA

    cat(sprintf("  %-12s %s: percentile=%d, OutFLANK=%d, overlap=%d (Jaccard=%.3f)\n",
                comp, label, length(pct_outliers), nrow(of_data),
                length(overlap), jaccard))
  }
  cat("\n")
}

if (!outflank_available) {
  cat("  No OutFLANK results available for cross-referencing.\n\n")
}

################################################################################
# Step 7: Summary report
################################################################################

cat("\nSTEP 7: Writing summary report\n")
cat("------------------------------------------------------------------------\n\n")

report_file <- file.path(OUTPUT_DIR, "summary_report.txt")
sink(report_file)

cat("================================================================================\n")
cat("PER-SNP PERCENTILE SELECTION SCAN - SUMMARY REPORT\n")
cat("Species: Coccinella septempunctata (seven-spotted lady beetle)\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

cat("METHODOLOGY\n")
cat("------------------------------------------------------------------------\n\n")
cat("Per-SNP Weir & Cockerham FST values from vcftools (no windowing).\n")
cat("Genome-wide percentile thresholds applied per comparison.\n")
cat("Negative FST values retained in the distribution (indicate low differentiation).\n")
cat("Post-hoc clustering distance:", CLUSTER_DIST, "bp\n\n")

cat("This REPLACES the windowed unified percentile approach (unified_selection_scan.R).\n")
cat("Rationale: puts FST-based detection on the same SNP-level scale as pcadapt and\n")
cat("OutFLANK, enabling direct cross-method concordance at the SNP level.\n\n")

cat("INPUT DATA\n")
cat("------------------------------------------------------------------------\n\n")
print(as.data.frame(comparison_stats %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))
cat("\n\n")

cat("OUTLIERS BY THRESHOLD\n")
cat("------------------------------------------------------------------------\n\n")
for (label in threshold_labels) {
  df <- all_outliers[[label]]
  cat(label, ":\n")
  cat("  Total outlier SNP-comparison pairs:", nrow(df), "\n")
  counts <- df %>% count(comparison)
  for (i in seq_len(nrow(counts))) {
    cat("    ", counts$comparison[i], ":", counts$n[i], "\n")
  }
  cat("\n")
}

cat("UNIQUE OUTLIER SNPS (1% threshold)\n")
cat("------------------------------------------------------------------------\n\n")
cat("Total unique SNPs:", nrow(snp_summary), "\n\n")

cat("Pattern distribution:\n")
print(as.data.frame(pattern_counts))
cat("\n")

cat("Confidence distribution:\n")
print(as.data.frame(confidence_counts))
cat("\n\n")

cat("CLUSTERED REGIONS (1% threshold, ", CLUSTER_DIST, " bp merge distance)\n")
cat("------------------------------------------------------------------------\n\n")
if (nrow(clustered) > 0) {
  cat("Total regions:", nrow(clustered), "\n")
  cat("Singleton SNPs (1-SNP regions):", sum(clustered$n_snps == 1), "\n")
  cat("Multi-SNP regions:", sum(clustered$n_snps > 1), "\n")
  cat("Largest region:", max(clustered$n_snps), "SNPs,",
      max(clustered$span_bp), "bp\n\n")
}

cat("OUTPUT FILES\n")
cat("------------------------------------------------------------------------\n\n")
cat("Outlier tables:\n")
for (label in threshold_labels) {
  cat("  - outliers_", label, ".tsv\n", sep = "")
}
cat("\nPattern categorization:\n")
cat("  - pattern_categorization/both_snps.tsv\n")
cat("  - pattern_categorization/invasion_specific_snps.tsv\n")
cat("  - pattern_categorization/post_invasion_snps.tsv\n")
cat("  - snp_summary_1pct.tsv (all patterns, full annotation)\n")
cat("\nClustered regions:\n")
cat("  - clustered_regions/outlier_regions_1pct.tsv\n")
cat("\nPlots:\n")
cat("  - plots/manhattan_{comparison}.png (6 per-comparison Manhattans)\n")
cat("  - plots/manhattan_all_comparisons.png (combined faceted)\n\n")

cat("NOTES\n")
cat("------------------------------------------------------------------------\n\n")
cat("1. CHI genotype count filter NOT yet applied. The per-SNP FST files from\n")
cat("   vcftools do not include per-site sample counts. A separate vcftools run\n")
cat("   with --counts or extraction from the VCF is needed to implement the\n")
cat("   CHI >=3 genotype filter for CHI-vs-invasive comparisons.\n\n")
cat("2. Sample exclusion decision still pending. These results use the full\n")
cat("   88-sample dataset. Re-run after exclusions are decided.\n\n")
cat("3. Cross-method concordance (percentile + OutFLANK + pcadapt) will be\n")
cat("   calculated in a separate script after pcadapt is complete.\n\n")

cat("================================================================================\n")
cat("END OF REPORT\n")
cat("================================================================================\n")

sink()
cat("Report saved:", report_file, "\n")

################################################################################
# Done
################################################################################

cat("\n================================================================================\n")
cat("PER-SNP SELECTION SCAN COMPLETE\n")
cat("================================================================================\n\n")

cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("Quick summary (1% threshold):\n")
cat("  Unique outlier SNPs:", nrow(snp_summary), "\n")
for (i in seq_len(nrow(pattern_counts))) {
  cat("    ", pattern_counts$pattern[i], ":", pattern_counts$n_snps[i], "\n")
}
cat("  Clustered regions:", nrow(clustered), "\n\n")

cat("Next steps:\n")
cat("  1. Review Manhattan plots and summary report\n")
cat("  2. Implement CHI genotype count filter (requires per-site sample counts)\n")
cat("  3. Run pcadapt selection scan\n")
cat("  4. Calculate cross-method concordance at SNP level\n\n")
