#!/usr/bin/env Rscript

# title: diagnostic_chrom_distribution.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-01
# last modified: 2026-05-01
#
# purpose:
#   Manhattan-scatter diagnostic, Task 3. For each of the 6 pairwise
#   comparisons, computes the per-chromosome density of top-1% per-SNP
#   percentile-method outliers (outlier count / SNP count per chromosome)
#   and produces (a) one bar chart per comparison plus a combined 2x3
#   panel, and (b) a chi-squared goodness-of-fit test of the null that
#   per-chromosome outlier counts are proportional to per-chromosome SNP
#   counts. Roughly uniform density across chromosomes is consistent
#   with noise floor / genome-wide structure / polygenic adaptation
#   (explanations 1, 2, 4); concentration on a few chromosomes is
#   consistent with locus-specific selection blunted at SNP resolution.
#   Sex chromosome NC_058198.1 is excluded throughout (autosomes only),
#   matching the convention used in cross_method_concordance.R; the
#   number of sex-chrom outliers in the upstream file is reported as a
#   sanity-check side note in the run log.
#
#   Outlier source: results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv,
#   which is the canonical input that feeds cross_method_concordance.R
#   (CHI genotype-count filter applied to CHI-vs-X comparisons; no filter
#   on among-invasive comparisons).
#
# inputs:
#   - results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv
#       Columns: chrom, pos, snp_id, comparison, fst, percentile_rank,
#       threshold, fst_cutoff. One row per outlier SNP per comparison.
#   - results/fst/{COMP}/fst_genomewide/{CHROM}.filtered.weir.fst
#       Used as the denominator (per-chromosome SNP count, post any
#       CHI-filter that was already applied upstream of the percentile
#       call). NA F_ST rows are dropped to match the per-SNP scan's
#       denominator.
#
# outputs:
#   - figures/diagnostic/chrom_distribution/
#       {COMP}.png                      : Per-comparison bar chart.
#       all_comparisons.png             : 2x3 panel (CHI row / invasive row).
#   - results/diagnostic/chrom_uniformity_test.tsv
#       One row per comparison: comparison, comp_type, n_chroms,
#       n_outliers, chi_sq, df, p_value, max_density_chrom,
#       max_over_min_density_ratio, sex_chrom_outliers (sanity check).
#   - results/diagnostic/chrom_distribution_run.log
#
# usage example:
#   Rscript scripts/diagnostic_chrom_distribution.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

# --- Paths ------------------------------------------------------------------

project_root <- if (file.exists("results/fst")) {
  "."
} else if (file.exists("../results/fst")) {
  ".."
} else {
  stop("Cannot locate results/fst/. Run from project root or scripts/.")
}

fst_base    <- file.path(project_root, "results/fst")
outliers_in <- file.path(project_root,
                         "results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv")

fig_dir <- file.path(project_root, "figures/diagnostic/chrom_distribution")
out_tsv <- file.path(project_root, "results/diagnostic/chrom_uniformity_test.tsv")
out_log <- file.path(project_root, "results/diagnostic/chrom_distribution_run.log")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

# --- Configuration ----------------------------------------------------------

chi_comparisons      <- c("CHI_vs_EEU", "CHI_vs_USA", "CHI_vs_WEU")
invasive_comparisons <- c("EEU_vs_USA", "EEU_vs_WEU", "USA_vs_WEU")
all_comparisons      <- c(chi_comparisons, invasive_comparisons)

autosomes <- sprintf("NC_0581%02d.1", 89:97)   # 9 autosomes
sex_chrom <- "NC_058198.1"

# Map NCBI accession to short chromosome label (1..9) for x-axis ticks.
chrom_labels <- setNames(seq_along(autosomes), autosomes)

# --- Logging ----------------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

log_msg("=== Task 3: Per-chromosome outlier-density diagnostic ===")
log_msg(sprintf("project root: %s", normalizePath(project_root)))

# --- Step 1: Load outlier file ---------------------------------------------

log_msg("\n[1/4] Loading outliers")
out_dt <- fread(outliers_in, sep = "\t", header = TRUE)
log_msg(sprintf("  loaded %d outlier rows from %s",
                nrow(out_dt), basename(outliers_in)))

# Sanity check: how many outliers landed on the sex chromosome? The per-SNP
# percentile scan does not exclude it; downstream concordance does. Report
# the count per comparison and then drop them.
sex_counts <- out_dt[chrom == sex_chrom, .N, by = comparison][order(comparison)]
if (nrow(sex_counts) > 0) {
  log_msg("  sex-chromosome outliers in upstream file (will be dropped here):")
  for (i in seq_len(nrow(sex_counts))) {
    log_msg(sprintf("    %-12s n=%d", sex_counts$comparison[i], sex_counts$N[i]))
  }
} else {
  log_msg("  no sex-chromosome outliers in upstream file (good).")
}

out_auto <- out_dt[chrom %in% autosomes]
log_msg(sprintf("  retained %d outlier rows on the 9 autosomes", nrow(out_auto)))

# --- Step 2: Per-comparison per-chromosome SNP counts ---------------------

log_msg("\n[2/4] Computing per-chromosome SNP counts (denominators)")

count_snps_one <- function(comp, chrom) {
  f <- file.path(fst_base, comp, "fst_genomewide",
                 sprintf("%s.filtered.weir.fst", chrom))
  if (!file.exists(f)) stop(sprintf("missing per-SNP F_ST file: %s", f))
  dt <- fread(f, sep = "\t", header = TRUE,
              colClasses = c(CHROM = "character",
                             POS = "integer",
                             WEIR_AND_COCKERHAM_FST = "double"))
  sum(!is.na(dt$WEIR_AND_COCKERHAM_FST))
}

# Loop in nested order: comparison x chromosome. About 50 iters total.
denom_rows <- vector("list", length(all_comparisons) * length(autosomes))
k <- 1
for (comp in all_comparisons) {
  for (chrom in autosomes) {
    n <- count_snps_one(comp, chrom)
    denom_rows[[k]] <- tibble(comparison = comp, chrom = chrom, n_snps = n)
    k <- k + 1
  }
  log_msg(sprintf("  %-12s denominators loaded", comp))
}
denom_df <- bind_rows(denom_rows)

# --- Step 3: Density and chi-squared test ----------------------------------

log_msg("\n[3/4] Computing per-chromosome outlier density and uniformity test")

# Per-chrom outlier counts.
out_counts <- out_auto[, .(n_outliers = .N), by = .(comparison, chrom)]
out_counts <- as_tibble(out_counts)

# Full grid (every comparison x autosome combination), to make sure
# zero-outlier chromosomes are not silently dropped.
grid <- crossing(comparison = all_comparisons, chrom = autosomes)
density_df <- grid %>%
  left_join(denom_df,    by = c("comparison", "chrom")) %>%
  left_join(out_counts,  by = c("comparison", "chrom")) %>%
  mutate(
    n_outliers = replace_na(n_outliers, 0L),
    density    = n_outliers / n_snps,
    chrom_idx  = chrom_labels[chrom],
    comp_type  = if_else(comparison %in% chi_comparisons,
                         "CHI_vs_invasive", "Among_invasive")
  ) %>%
  arrange(comparison, chrom_idx)

# Chi-squared test per comparison: null is outliers proportional to SNP count.
test_rows <- vector("list", length(all_comparisons))
for (i in seq_along(all_comparisons)) {
  comp <- all_comparisons[i]
  d <- density_df %>% filter(comparison == comp)
  observed <- d$n_outliers
  total    <- sum(observed)
  expected <- total * d$n_snps / sum(d$n_snps)
  chi_sq   <- sum((observed - expected)^2 / expected)
  df       <- length(observed) - 1L     # 9 - 1 = 8
  p_val    <- pchisq(chi_sq, df = df, lower.tail = FALSE)

  max_dens   <- max(d$density)
  min_dens   <- min(d$density[d$density > 0])
  max_chrom  <- d$chrom[which.max(d$density)]
  ratio      <- if (min_dens == 0) Inf else max_dens / min_dens
  sex_n      <- if (nrow(sex_counts) == 0) 0L else
                sex_counts$N[sex_counts$comparison == comp]
  if (length(sex_n) == 0) sex_n <- 0L

  test_rows[[i]] <- tibble(
    comparison      = comp,
    comp_type       = if (comp %in% chi_comparisons)
                        "CHI_vs_invasive" else "Among_invasive",
    n_chroms        = length(observed),
    n_outliers      = total,
    chi_sq          = chi_sq,
    df              = df,
    p_value         = p_val,
    max_density_chrom = max_chrom,
    max_density       = max_dens,
    max_over_min_density_ratio = ratio,
    sex_chrom_outliers_in_upstream_file = sex_n
  )
}
test_df <- bind_rows(test_rows) %>%
  mutate(comp_type = factor(comp_type,
                            levels = c("CHI_vs_invasive", "Among_invasive"))) %>%
  arrange(comp_type, comparison)

write_tsv(test_df, out_tsv, na = "NA")
log_msg(sprintf("  wrote: %s", out_tsv))
log_msg("\n  uniformity test results:")
for (i in seq_len(nrow(test_df))) {
  r <- test_df[i, ]
  log_msg(sprintf(
    "    %-12s [%s]  n=%d  chi^2=%.1f (df=%d)  p=%.3g  max_dens=%.4f on chr%d  ratio=%s",
    r$comparison, as.character(r$comp_type), r$n_outliers,
    r$chi_sq, r$df, r$p_value, r$max_density,
    chrom_labels[r$max_density_chrom],
    if (is.infinite(r$max_over_min_density_ratio)) "Inf"
      else sprintf("%.1fx", r$max_over_min_density_ratio)
  ))
}

# --- Step 4: Bar charts -----------------------------------------------------

log_msg("\n[4/4] Plotting bar charts")

type_colors <- c(CHI_vs_invasive = "#E07B7B",
                 Among_invasive  = "#5B8AB7")

make_bar <- function(comp) {
  d <- density_df %>% filter(comparison == comp)
  ct <- unique(d$comp_type)
  fill_col <- type_colors[[ct]]
  this_test <- test_df %>% filter(comparison == comp)

  # Reference line: genome-wide density under perfect uniformity =
  # total_outliers / total_snps for this comparison. If observed densities
  # match the line, distribution is uniform.
  total_out  <- sum(d$n_outliers)
  total_snps <- sum(d$n_snps)
  ref_dens   <- total_out / total_snps

  subtitle <- sprintf(
    "n_outliers=%s  uniformity chi^2=%.1f (df=%d)  p=%.3g",
    format(total_out, big.mark = ","),
    this_test$chi_sq, this_test$df, this_test$p_value
  )

  ggplot(d, aes(x = factor(chrom_idx), y = density)) +
    geom_col(fill = fill_col, color = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = ref_dens, linetype = "dashed",
               color = "grey20", linewidth = 0.4) +
    geom_text(aes(label = format(n_outliers, big.mark = ",")),
              vjust = -0.4, size = 2.5, color = "grey20") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                       labels = scales::label_percent(accuracy = 0.1)) +
    labs(
      title    = sprintf("%s -- per-chromosome outlier density", comp),
      subtitle = subtitle,
      x = "Chromosome (1-9; autosomes only)",
      y = "Outlier density (outliers / SNPs on chromosome)",
      caption = sprintf(
        "Dashed line: genome-wide outlier rate under uniformity (=%.3f%%). Bar labels: outlier counts.",
        100 * ref_dens
      )
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, color = "grey20"),
          plot.caption  = element_text(size = 7, color = "grey40", hjust = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())
}

# Individual figures.
for (comp in all_comparisons) {
  p <- make_bar(comp)
  fpath <- file.path(fig_dir, sprintf("%s.png", comp))
  ggsave(fpath, p, width = 160, height = 100, units = "mm", dpi = 200)
  log_msg(sprintf("  wrote: %s.png", comp))
}

# Combined 2x3 panel: CHI on top row, invasive on bottom row.
panels <- lapply(all_comparisons, function(comp) {
  p <- make_bar(comp)
  p + labs(title = comp, subtitle = NULL, caption = NULL) +
    theme(plot.title = element_text(size = 9, face = "bold"))
})
combined <- (panels[[1]] | panels[[2]] | panels[[3]]) /
            (panels[[4]] | panels[[5]] | panels[[6]]) +
  plot_annotation(
    title    = "Per-chromosome outlier density across 6 comparisons",
    subtitle = "Top row: CHI-vs-invasive. Bottom row: among-invasive. Dashed line per panel = genome-wide rate under uniformity.",
    caption  = "Bar labels: outlier counts. Bar height: outlier count / per-chromosome SNP count.",
    theme = theme(plot.title    = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 8.5, color = "grey20"),
                  plot.caption  = element_text(size = 7, color = "grey40", hjust = 0))
  )
fpath <- file.path(fig_dir, "all_comparisons.png")
ggsave(fpath, combined, width = 320, height = 180, units = "mm", dpi = 200)
log_msg(sprintf("  wrote: %s", fpath))

# --- Run log ----------------------------------------------------------------

writeLines(log_lines, out_log)
log_msg(sprintf("\nRun log: %s", out_log))
