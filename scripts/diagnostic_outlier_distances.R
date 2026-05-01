#!/usr/bin/env Rscript

# title: diagnostic_outlier_distances.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-01
# last modified: 2026-05-01
#
# purpose:
#   Manhattan-scatter diagnostic, Task 4 Part A (stratified). For each of
#   the 6 pairwise comparisons, computes the distribution of distances
#   (in bp) between consecutive outlier SNPs along the genome, stratified
#   into three groups:
#
#     1. Chr-1 28 Mb peak    : outliers on chr 1 in 26 <= pos < 30 Mb.
#                              Centred on the cross-comparison-shared
#                              peak identified in the chr-1 follow-up.
#     2. Chr-1 control       : a 4 Mb window elsewhere on chr 1 (not
#                              overlapping 24-32 Mb buffer) selected per
#                              comparison to have outlier count closest
#                              to the peak window. This is "matched-count
#                              elsewhere on chr 1" -- if peak distances
#                              are tighter than control distances, the
#                              peak is unusually concentrated relative to
#                              count-matched chr-1 background.
#     3. Other autosomes     : outliers on chr 2-9. Distances computed
#                              within each chromosome separately and
#                              pooled across chromosomes (no cross-chrom
#                              distances).
#
#   For each stratum and each comparison, the script computes adjacent-
#   outlier distance vectors, plots log-x histograms (one panel per
#   comparison with the three strata overlaid), summarizes per-stratum
#   distance distributions (n, median, IQR), and runs pairwise
#   Kolmogorov-Smirnov tests across strata. Strata are confounded with
#   different upper bounds on possible distance (peak/control are
#   capped at 4 Mb; other autosomes can run to ~70 Mb), so the
#   informative comparison is the LOWER tail of the log-distance
#   distribution -- KS and median capture that.
#
#   Predictions (per the user's spec):
#     If peak distances are tighter than control distances => sustained
#     signal in a low-recombination context, consistent with a structural
#     variant or strong sweep at the peak.
#     If peak distances match control distances => locus is real but not
#     structurally distinguished; multiple independent SNPs around 28 Mb
#     are tagging the same selection signal at LD-decay scale.
#
# inputs:
#   - results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv
#       Outlier set used as the SNP universe for the distance analysis.
#       Columns: chrom, pos, snp_id, comparison, fst, percentile_rank,
#                threshold, fst_cutoff.
#
# outputs:
#   - figures/diagnostic/outlier_distances/
#       {COMP}_distance_hist.png        : Per-comparison distance histogram,
#                                         3 strata overlaid (log-x).
#       all_comparisons.png             : 2x3 facet panel.
#       control_window_choices.png      : Bar chart of which 4 Mb window
#                                         on chr 1 was picked as control
#                                         for each comparison.
#   - results/diagnostic/outlier_distances_summary.tsv
#       One row per (comparison, stratum) with n_outliers, n_distances,
#       median_distance, q25, q75, mean_distance, control_window_start_mb.
#   - results/diagnostic/outlier_distances_ks_tests.tsv
#       Pairwise KS-test p-values per comparison (peak vs control,
#       peak vs other, control vs other).
#   - results/diagnostic/outlier_distances_run.log
#
# usage example:
#   Rscript scripts/diagnostic_outlier_distances.R

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

outliers_in <- file.path(project_root,
                         "results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv")
fig_dir     <- file.path(project_root, "figures/diagnostic/outlier_distances")
out_dir     <- file.path(project_root, "results/diagnostic")
out_summary <- file.path(out_dir, "outlier_distances_summary.tsv")
out_ks      <- file.path(out_dir, "outlier_distances_ks_tests.tsv")
out_log     <- file.path(out_dir, "outlier_distances_run.log")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Configuration ----------------------------------------------------------

chi_comparisons      <- c("CHI_vs_EEU", "CHI_vs_USA", "CHI_vs_WEU")
invasive_comparisons <- c("EEU_vs_USA", "EEU_vs_WEU", "USA_vs_WEU")
all_comparisons      <- c(chi_comparisons, invasive_comparisons)
autosomes <- sprintf("NC_0581%02d.1", 89:97)
chr1      <- "NC_058189.1"
chrom_labels <- setNames(seq_along(autosomes), autosomes)

# Chr-1 peak window (Mb) and buffer for excluding overlap when picking
# the count-matched control.
peak_start_mb <- 26
peak_end_mb   <- 30
buffer_start_mb <- 24
buffer_end_mb   <- 32
window_size_mb  <- 4

# Strata levels in canonical plotting order.
strat_levels <- c("chr1_peak_26_30Mb",
                  "chr1_control",
                  "other_autosomes_2to9")
strat_colors <- c(chr1_peak_26_30Mb     = "#D55E00",   # orange
                  chr1_control          = "#009E73",   # teal
                  other_autosomes_2to9  = "#5B8AB7")   # blue

# --- Logging ----------------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

log_msg("=== Task 4 Part A: outlier-outlier distance distributions (stratified) ===")
log_msg(sprintf("peak window: chr 1 [%d, %d) Mb", peak_start_mb, peak_end_mb))
log_msg(sprintf("buffer:      chr 1 [%d, %d) Mb (excluded from control choice)",
                buffer_start_mb, buffer_end_mb))
log_msg(sprintf("window size: %d Mb", window_size_mb))

# --- Step 1: Load outliers ---------------------------------------------------

log_msg("\n[1] Loading outliers")
out_dt <- fread(outliers_in, sep = "\t", header = TRUE,
                colClasses = c(chrom = "character", pos = "integer"))
out_dt <- out_dt[chrom %in% autosomes]
log_msg(sprintf("  retained %d autosomal outlier rows across %d comparisons",
                nrow(out_dt), length(unique(out_dt$comparison))))

# --- Step 2: Pick chr-1 control window per comparison -----------------------

log_msg("\n[2] Selecting count-matched chr-1 control window per comparison")

# Chr-1 length (use max observed pos as a defensive upper bound; matches
# the chr-1 follow-up). 71.18 Mb known from VCF header.
chr1_max_pos_mb <- ceiling(max(out_dt$pos[out_dt$chrom == chr1]) / 1e6)
chr1_len_mb     <- 72  # rounded up from 71.18 Mb so windows up to 68 Mb start are valid

# Sliding-window over chr 1 at 1 Mb step, 4 Mb wide. For each comparison,
# find the window NOT overlapping the buffer with outlier count closest to
# the peak window count.
window_starts_mb <- 0:(chr1_len_mb - window_size_mb)
window_choices <- vector("list", length(all_comparisons))
for (i in seq_along(all_comparisons)) {
  comp <- all_comparisons[i]
  chr1_pos_mb <- out_dt$pos[out_dt$chrom == chr1 & out_dt$comparison == comp] / 1e6

  # Peak count.
  n_peak <- sum(chr1_pos_mb >= peak_start_mb & chr1_pos_mb < peak_end_mb)

  # Window counts.
  win_counts <- sapply(window_starts_mb, function(s) {
    sum(chr1_pos_mb >= s & chr1_pos_mb < s + window_size_mb)
  })

  # Mark windows overlapping the buffer.
  overlap_buffer <- (window_starts_mb < buffer_end_mb) &
                    (window_starts_mb + window_size_mb > buffer_start_mb)
  eligible <- !overlap_buffer

  if (!any(eligible)) {
    stop(sprintf("[%s] no eligible chr-1 window outside buffer; check config.",
                 comp))
  }

  diffs <- abs(win_counts[eligible] - n_peak)
  best_idx_in_eligible <- which.min(diffs)
  best_start <- window_starts_mb[eligible][best_idx_in_eligible]
  best_n     <- win_counts[eligible][best_idx_in_eligible]

  window_choices[[i]] <- tibble(
    comparison         = comp,
    n_peak             = n_peak,
    control_start_mb   = best_start,
    control_end_mb     = best_start + window_size_mb,
    n_control          = best_n,
    abs_diff_in_count  = abs(best_n - n_peak)
  )
  log_msg(sprintf("  %-12s peak n=%d   chosen control [%d, %d) Mb n=%d   |diff|=%d",
                  comp, n_peak, best_start, best_start + window_size_mb,
                  best_n, abs(best_n - n_peak)))
}
window_choices_df <- bind_rows(window_choices)

# --- Step 3: Stratify outliers and compute distances ------------------------

log_msg("\n[3] Stratifying outliers and computing adjacent-outlier distances")

# Helper: compute distances between consecutive outliers within each chrom,
# pooled across chromosomes within the stratum.
distances_within_chroms <- function(df) {
  # df: chrom, pos. Sort and diff within each chromosome.
  if (nrow(df) < 2) return(integer(0))
  df <- df %>% arrange(chrom, pos)
  res <- df %>%
    group_by(chrom) %>%
    arrange(pos, .by_group = TRUE) %>%
    mutate(d = c(NA_integer_, diff(pos))) %>%
    ungroup() %>%
    pull(d)
  res[!is.na(res)]
}

stratify_one <- function(comp) {
  o <- out_dt[comparison == comp]
  ctrl <- window_choices_df %>% filter(comparison == comp)

  is_peak <- (o$chrom == chr1) &
             (o$pos >= peak_start_mb * 1e6) &
             (o$pos <  peak_end_mb   * 1e6)
  is_ctrl <- (o$chrom == chr1) &
             (o$pos >= ctrl$control_start_mb * 1e6) &
             (o$pos <  ctrl$control_end_mb   * 1e6)
  is_oth  <- o$chrom != chr1

  list(
    peak    = as_tibble(o[is_peak, .(chrom, pos)]),
    control = as_tibble(o[is_ctrl, .(chrom, pos)]),
    other   = as_tibble(o[is_oth,  .(chrom, pos)])
  )
}

dist_long_rows <- vector("list", length(all_comparisons) * 3)
summary_rows   <- vector("list", length(all_comparisons) * 3)
ki <- 1
for (i in seq_along(all_comparisons)) {
  comp   <- all_comparisons[i]
  parts  <- stratify_one(comp)
  ctrl   <- window_choices_df %>% filter(comparison == comp)

  for (s in c("peak", "control", "other")) {
    d <- distances_within_chroms(parts[[s]])
    stratum_name <- switch(s,
                           peak    = "chr1_peak_26_30Mb",
                           control = "chr1_control",
                           other   = "other_autosomes_2to9")
    dist_long_rows[[ki]] <- tibble(
      comparison = comp,
      stratum    = stratum_name,
      distance   = d
    )
    summary_rows[[ki]] <- tibble(
      comparison              = comp,
      stratum                 = stratum_name,
      n_outliers              = nrow(parts[[s]]),
      n_distances             = length(d),
      median_distance_bp      = if (length(d) == 0) NA_real_ else median(d),
      q25_distance_bp         = if (length(d) == 0) NA_real_ else quantile(d, 0.25),
      q75_distance_bp         = if (length(d) == 0) NA_real_ else quantile(d, 0.75),
      mean_distance_bp        = if (length(d) == 0) NA_real_ else mean(d),
      control_window_start_mb = if (s == "control") ctrl$control_start_mb else NA_real_
    )
    ki <- ki + 1
  }
}
dist_long  <- bind_rows(dist_long_rows)
summary_df <- bind_rows(summary_rows) %>%
  mutate(stratum = factor(stratum, levels = strat_levels)) %>%
  arrange(comparison, stratum)
write_tsv(summary_df, out_summary, na = "NA")
log_msg(sprintf("  wrote: %s", out_summary))

# Console summary.
log_msg("\n  per-comparison median distance (bp) by stratum:")
for (comp in all_comparisons) {
  rows <- summary_df %>% filter(comparison == comp)
  log_msg(sprintf(
    "    %-12s peak: med=%6.0f (n_d=%d)   control: med=%6.0f (n_d=%d)   other: med=%6.0f (n_d=%d)",
    comp,
    rows$median_distance_bp[rows$stratum == "chr1_peak_26_30Mb"],
    rows$n_distances[rows$stratum == "chr1_peak_26_30Mb"],
    rows$median_distance_bp[rows$stratum == "chr1_control"],
    rows$n_distances[rows$stratum == "chr1_control"],
    rows$median_distance_bp[rows$stratum == "other_autosomes_2to9"],
    rows$n_distances[rows$stratum == "other_autosomes_2to9"]
  ))
}

# --- Step 4: Pairwise KS tests ---------------------------------------------

log_msg("\n[4] Pairwise Kolmogorov-Smirnov tests per comparison")

ks_rows <- vector("list", length(all_comparisons) * 3)
ki <- 1
ks_pairs <- list(
  c("chr1_peak_26_30Mb", "chr1_control"),
  c("chr1_peak_26_30Mb", "other_autosomes_2to9"),
  c("chr1_control",      "other_autosomes_2to9")
)
for (comp in all_comparisons) {
  for (pp in ks_pairs) {
    a <- dist_long %>% filter(comparison == comp, stratum == pp[1]) %>% pull(distance)
    b <- dist_long %>% filter(comparison == comp, stratum == pp[2]) %>% pull(distance)
    if (length(a) >= 2 && length(b) >= 2) {
      kt <- suppressWarnings(ks.test(a, b))
      ks_stat <- as.numeric(kt$statistic)
      p_val   <- as.numeric(kt$p.value)
    } else {
      ks_stat <- NA_real_; p_val <- NA_real_
    }
    ks_rows[[ki]] <- tibble(
      comparison    = comp,
      stratum_a     = pp[1],
      stratum_b     = pp[2],
      n_a           = length(a),
      n_b           = length(b),
      median_a_bp   = if (length(a) == 0) NA_real_ else median(a),
      median_b_bp   = if (length(b) == 0) NA_real_ else median(b),
      ks_statistic  = ks_stat,
      ks_p_value    = p_val
    )
    ki <- ki + 1
  }
}
ks_df <- bind_rows(ks_rows)
write_tsv(ks_df, out_ks, na = "NA")
log_msg(sprintf("  wrote: %s", out_ks))

log_msg("\n  KS test results (peak-vs-control most directly tests the user's prediction):")
for (i in seq_len(nrow(ks_df))) {
  r <- ks_df[i, ]
  log_msg(sprintf("    %-12s %-22s vs %-22s  D=%.3f  p=%.3g",
                  r$comparison, r$stratum_a, r$stratum_b,
                  r$ks_statistic, r$ks_p_value))
}

# --- Step 5: Plotting -------------------------------------------------------

log_msg("\n[5] Plotting distance histograms")

# Common log-x bins: 1 bp to 100 Mb (10^0 to 10^8) in steps of 0.1 in log10.
log_breaks <- 10^seq(0, 8, by = 0.1)

stratum_label <- function(x) {
  recode(x,
         chr1_peak_26_30Mb   = "chr 1 peak (26-30 Mb)",
         chr1_control        = "chr 1 control (count-matched)",
         other_autosomes_2to9 = "other autosomes (chr 2-9)")
}

make_hist <- function(comp) {
  d <- dist_long %>% filter(comparison == comp, !is.na(distance))
  ctrl <- window_choices_df %>% filter(comparison == comp)
  ctrl_lbl <- sprintf("[%d, %d) Mb", ctrl$control_start_mb, ctrl$control_end_mb)
  meds <- summary_df %>% filter(comparison == comp)

  subtitle <- sprintf("peak n=%d  control %s n=%d  other n=%d",
                      meds$n_distances[meds$stratum == "chr1_peak_26_30Mb"],
                      ctrl_lbl,
                      meds$n_distances[meds$stratum == "chr1_control"],
                      meds$n_distances[meds$stratum == "other_autosomes_2to9"])

  ggplot(d, aes(x = distance, fill = stratum, color = stratum)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 50, alpha = 0.5,
                   position = "identity", linewidth = 0.15) +
    geom_vline(data = meds %>% filter(!is.na(median_distance_bp)),
               aes(xintercept = median_distance_bp, color = stratum),
               linetype = "dashed", linewidth = 0.4) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10",
                                                scales::math_format(10^.x))) +
    scale_fill_manual(values = strat_colors, labels = stratum_label,
                      drop = FALSE) +
    scale_color_manual(values = strat_colors, labels = stratum_label,
                       drop = FALSE) +
    labs(
      title    = sprintf("%s -- adjacent-outlier distance distribution", comp),
      subtitle = subtitle,
      x        = "Distance to next outlier within stratum (bp; log10)",
      y        = "Density",
      fill     = "Stratum",
      color    = "Stratum",
      caption  = "Dashed verticals = per-stratum medians."
    ) +
    theme_bw(base_size = 9) +
    theme(plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8),
          plot.caption  = element_text(size = 7, color = "grey40", hjust = 0),
          legend.position = "bottom",
          legend.text   = element_text(size = 7),
          legend.title  = element_text(size = 8))
}

# Per-comparison figures.
for (comp in all_comparisons) {
  p <- make_hist(comp)
  ggsave(file.path(fig_dir, sprintf("%s_distance_hist.png", comp)),
         p, width = 200, height = 120, units = "mm", dpi = 200)
  log_msg(sprintf("  wrote: %s_distance_hist.png", comp))
}

# Combined panel.
panels <- lapply(all_comparisons, function(comp) {
  make_hist(comp) +
    labs(title = comp, subtitle = NULL, caption = NULL) +
    theme(plot.title = element_text(size = 9, face = "bold"))
})
combined <- (panels[[1]] | panels[[2]] | panels[[3]]) /
            (panels[[4]] | panels[[5]] | panels[[6]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Adjacent-outlier distance distributions across 6 comparisons",
    subtitle = "Top row: CHI-vs-invasive. Bottom row: among-invasive. Three strata overlaid per panel.",
    caption = "Dashed verticals = per-stratum medians. Distance is bp to next outlier within the same stratum, computed within chromosomes.",
    theme = theme(plot.title    = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(size = 9),
                  plot.caption  = element_text(size = 7, color = "grey40", hjust = 0),
                  legend.position = "bottom")
  )
ggsave(file.path(fig_dir, "all_comparisons.png"),
       combined, width = 320, height = 200, units = "mm", dpi = 200)
log_msg(sprintf("  wrote: all_comparisons.png"))

# Control window choices figure.
p_ctrl <- ggplot(window_choices_df,
                 aes(x = comparison, y = control_start_mb, fill = comparison)) +
  geom_col(color = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = c(peak_start_mb, peak_end_mb),
             linetype = "dashed", color = "grey20") +
  geom_text(aes(label = sprintf("[%d, %d) Mb\nn=%d (peak n=%d)",
                                control_start_mb, control_end_mb,
                                n_control, n_peak)),
            vjust = -0.4, size = 2.5) +
  annotate("text", x = 0.6, y = (peak_start_mb + peak_end_mb) / 2,
           label = "peak window\n26-30 Mb", hjust = 0, size = 2.8,
           color = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    title    = "Count-matched chr-1 control window selected per comparison",
    subtitle = "Each bar shows the start of the chosen 4 Mb control window. Dashed lines bracket the peak.",
    x = NULL, y = "Control window start position (Mb on chr 1)",
    fill = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(fig_dir, "control_window_choices.png"),
       p_ctrl, width = 200, height = 110, units = "mm", dpi = 200)
log_msg("  wrote: control_window_choices.png")

writeLines(log_lines, out_log)
log_msg(sprintf("\nRun log: %s", out_log))
