#!/usr/bin/env Rscript

# title: diagnostic_ld_decay.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-05-01
# last modified: 2026-05-01
#
# purpose:
#   Manhattan-scatter diagnostic, Task 5 (final task in the chain).
#   Per-chromosome LD-decay analysis to test whether chr 1 has a
#   distinct recombination regime versus the other 8 autosomes, or
#   whether the chr-1 outlier-clustering pattern from Tasks 3-4 is
#   purely a concentration of selection signal.
#
#   Two interpretations of the chr-1 pattern carry different manuscript
#   stories:
#     (a) Chr 1 has shorter LD half-life than other autosomes
#         => chr 1 has a different recombination regime; outlier
#         clustering at <500 bp scale is partly explained by faster LD
#         decay (more independent SNPs per Mb under selection).
#     (b) Chr 1 has similar LD half-life to other autosomes
#         => the outlier clustering is from denser true selection
#         targets, not faster recombination. This is the "concentrated
#         selection signal" story.
#
#   Approach: compute pairwise r^2 within USA samples (n=40, the
#   largest single population) per autosome, restricted to SNPs with
#   MAF >= 0.05 and within 100 kb. Bin pairs by physical distance on
#   a log scale, compute mean r^2 per bin per chromosome, and report
#   the empirical distance at which mean r^2 drops below 0.20 (a
#   conventional LD threshold) as the half-decay metric.
#
#   USA-only is chosen over genome-wide LD because mixing populations
#   inflates baseline r^2 via population structure, which masks
#   per-chromosome differences in true recombination. USA is the most
#   recently established invasive population (largest n among the four
#   pops; also the most uniformly sampled North American set), so its
#   LD reflects modern recombination rather than historical
#   bottleneck-driven LD.
#
# inputs:
#   - data/plink/chr_{CHROM}.{bed,bim,fam}    : per-chromosome PLINK files
#                                                (one set per autosome)
#   - metadata/pop_USA.txt                     : 40 USA sample IDs
#
# outputs:
#   - /tmp/ld_decay/{CHROM}.ld                : per-chrom PLINK r^2 output
#                                                (cached; regenerated only
#                                                if absent)
#   - figures/diagnostic/ld_decay/
#       per_chrom_decay_curves.png              : decay curves per chrom
#                                                  (linear and log-x panels)
#       chr1_vs_others_decay.png                : chr 1 highlighted vs the
#                                                  other 8 autosomes pooled
#       half_decay_distances.png                : bar chart of per-chrom
#                                                  half-decay distance
#   - results/diagnostic/ld_decay_summary.tsv
#       Per-chrom row: chrom, n_pairs, mean_r2_at_<100bp, mean_r2_at_5kb,
#       mean_r2_at_25kb, mean_r2_at_100kb, half_decay_distance_bp.
#   - results/diagnostic/ld_decay_run.log
#
# usage example:
#   Rscript scripts/diagnostic_ld_decay.R [--population USA]

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

# --- CLI arg parsing --------------------------------------------------------

parse_args <- function(argv) {
  pop <- "USA"
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("--population", "-p")) {
      pop <- argv[i + 1]
      i <- i + 2
    } else if (a %in% c("--help", "-h")) {
      cat("Usage: Rscript diagnostic_ld_decay.R [--population CHI|EEU|WEU|USA]\n")
      cat("  Default: USA (n=40, largest single population).\n")
      quit(status = 0)
    } else {
      stop(sprintf("unknown argument: %s", a))
    }
  }
  list(population = pop)
}
args <- parse_args(commandArgs(trailingOnly = TRUE))

# --- Paths ------------------------------------------------------------------

project_root <- if (file.exists("results/fst")) {
  "."
} else if (file.exists("../results/fst")) {
  ".."
} else {
  stop("Cannot locate results/fst/. Run from project root or scripts/.")
}

plink_bin   <- "/Users/reina/miniconda3/envs/pop_gen/bin/plink"
plink_dir   <- file.path(project_root, "data/plink")
pop_file_in <- file.path(project_root, sprintf("metadata/pop_%s.txt", args$population))

ld_cache_dir <- file.path("/tmp", sprintf("ld_decay_%s", args$population))
fig_dir      <- file.path(project_root, "figures/diagnostic/ld_decay")
out_dir      <- file.path(project_root, "results/diagnostic")
out_summary  <- file.path(out_dir, sprintf("ld_decay_summary_%s.tsv", args$population))
out_log      <- file.path(out_dir, sprintf("ld_decay_run_%s.log", args$population))
keep_file    <- file.path(ld_cache_dir, "keep_ids.txt")

dir.create(ld_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Configuration ----------------------------------------------------------

autosomes <- sprintf("NC_0581%02d.1", 89:97)
chr1      <- "NC_058189.1"
chrom_labels <- setNames(seq_along(autosomes), autosomes)

# PLINK r^2 parameters.
LD_WINDOW_KB    <- 100   # max physical distance for pairs (kb)
LD_WINDOW_SNPS  <- 99999 # effectively unbounded SNP count
LD_R2_MIN       <- 0     # keep all r^2 values (no threshold)
MAF_MIN         <- 0.05

# Half-decay metric: distance at which mean r^2 first drops below this value.
HALF_DECAY_THRESHOLD <- 0.20

# --- Logging ----------------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

log_msg(sprintf("=== Task 5: LD decay per chromosome (population = %s) ===",
                args$population))
log_msg(sprintf("max physical distance: %d kb", LD_WINDOW_KB))
log_msg(sprintf("MAF cutoff: %.2f", MAF_MIN))
log_msg(sprintf("half-decay threshold: r^2 < %.2f", HALF_DECAY_THRESHOLD))

# --- Step 1: Build PLINK --keep file ----------------------------------------

log_msg("\n[1] Building PLINK --keep file")
ids <- readLines(pop_file_in)
ids <- ids[nzchar(ids)]
log_msg(sprintf("  population %s: n=%d samples", args$population, length(ids)))
# PLINK expects FID IID (two columns). Per fam inspection, FID == IID for this
# dataset, so we duplicate the column.
keep_lines <- paste(ids, ids, sep = "\t")
writeLines(keep_lines, keep_file)
log_msg(sprintf("  wrote: %s", keep_file))

# --- Step 2: PLINK --r2 per chromosome --------------------------------------

log_msg("\n[2] Running PLINK --r2 per chromosome (cached in /tmp/)")

run_plink_r2 <- function(chrom) {
  bfile <- file.path(plink_dir, sprintf("chr_%s", chrom))
  out_prefix <- file.path(ld_cache_dir, chrom)
  ld_file <- paste0(out_prefix, ".ld")
  if (file.exists(ld_file) && file.size(ld_file) > 0) {
    log_msg(sprintf("  %-12s already cached (%s, %s MB)", chrom,
                    basename(ld_file),
                    format(round(file.size(ld_file) / 1e6, 1))))
    return(ld_file)
  }
  cmd <- sprintf(
    "%s --bfile %s --keep %s --maf %g --r2 --ld-window %d --ld-window-kb %d --ld-window-r2 %g --out %s --allow-extra-chr 2>&1",
    shQuote(plink_bin), shQuote(bfile), shQuote(keep_file),
    MAF_MIN, LD_WINDOW_SNPS, LD_WINDOW_KB, LD_R2_MIN,
    shQuote(out_prefix)
  )
  log_msg(sprintf("  running plink for %s ...", chrom))
  out <- system(cmd, intern = TRUE)
  # Capture plink's "X variants and Y people pass filters" line for the log.
  pass_line <- grep("variants and .* people pass filters", out, value = TRUE)
  if (length(pass_line) > 0) {
    log_msg(sprintf("    %s", trimws(pass_line[1])))
  }
  if (!file.exists(ld_file)) {
    stop(sprintf("PLINK did not produce %s", ld_file))
  }
  log_msg(sprintf("    wrote %s (%s MB)", basename(ld_file),
                  format(round(file.size(ld_file) / 1e6, 1))))
  ld_file
}

ld_files <- setNames(vapply(autosomes, run_plink_r2, character(1)), autosomes)

# --- Step 3: Read .ld files and bin distances -------------------------------

log_msg("\n[3] Reading .ld files and binning by distance")

# Logarithmic distance bins (bp). 11 bins from 1e1 to 1e5 bp, log-spaced.
bin_edges_bp <- 10 ^ seq(1, log10(LD_WINDOW_KB * 1000), length.out = 12)
bin_centers_bp <- sqrt(head(bin_edges_bp, -1) * tail(bin_edges_bp, -1))
n_bins <- length(bin_edges_bp) - 1

# Process one chromosome at a time to keep memory bounded.
process_one <- function(chrom) {
  ld_path <- ld_files[[chrom]]
  ld <- fread(ld_path, header = TRUE)
  setnames(ld, gsub("\\.", "_", names(ld)))   # safety for any dots
  if (!"R2" %in% names(ld)) {
    # PLINK 1.9 output: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
    if ("R^2" %in% names(ld)) setnames(ld, "R^2", "R2")
  }
  ld[, distance := abs(BP_B - BP_A)]
  ld <- ld[distance > 0 & distance <= LD_WINDOW_KB * 1000]
  ld[, bin_idx := findInterval(distance, bin_edges_bp,
                                rightmost.closed = TRUE)]
  ld <- ld[bin_idx >= 1 & bin_idx <= n_bins]
  binned <- ld[, .(n_pairs = .N,
                    mean_r2  = mean(R2, na.rm = TRUE),
                    median_r2 = median(R2, na.rm = TRUE)),
                by = bin_idx]
  binned[, chrom := chrom]
  binned[, bin_center_bp := bin_centers_bp[bin_idx]]
  binned[]
}

binned_rows <- vector("list", length(autosomes))
for (i in seq_along(autosomes)) {
  chrom <- autosomes[i]
  binned_rows[[i]] <- process_one(chrom)
  log_msg(sprintf("  %-12s %s pairs across %d bins",
                  chrom,
                  format(sum(binned_rows[[i]]$n_pairs), big.mark = ","),
                  nrow(binned_rows[[i]])))
}
binned_df <- bind_rows(binned_rows) %>%
  arrange(chrom, bin_idx) %>%
  mutate(chrom_idx = chrom_labels[chrom],
         is_chr1   = chrom == chr1)

# --- Step 4: Per-chrom half-decay distance ----------------------------------

log_msg(sprintf("\n[4] Computing half-decay distance per chrom (r^2 threshold = %.2f)",
                HALF_DECAY_THRESHOLD))

# Half-decay = first bin (in increasing distance order) at which mean r^2
# drops below the threshold. If r^2 never drops below threshold within the
# window, report NA and the value at the largest bin instead.
half_decay_one <- function(d) {
  d <- d %>% arrange(bin_idx)
  below <- which(d$mean_r2 < HALF_DECAY_THRESHOLD)
  if (length(below) == 0) {
    return(tibble(half_decay_bp = NA_real_,
                  reached_threshold = FALSE,
                  min_observed_r2 = min(d$mean_r2),
                  min_observed_at_bp = d$bin_center_bp[which.min(d$mean_r2)]))
  }
  first_bin_below <- d[below[1], ]
  tibble(half_decay_bp = first_bin_below$bin_center_bp,
         reached_threshold = TRUE,
         min_observed_r2 = min(d$mean_r2),
         min_observed_at_bp = d$bin_center_bp[which.min(d$mean_r2)])
}

# Reference r^2 at four distance scales for the summary table.
reference_distances <- c(100, 5000, 25000, 100000)
mean_r2_at <- function(d, target_bp) {
  if (nrow(d) == 0) return(NA_real_)
  d$mean_r2[which.min(abs(d$bin_center_bp - target_bp))]
}

summary_rows <- vector("list", length(autosomes))
for (i in seq_along(autosomes)) {
  chrom <- autosomes[i]
  d <- binned_df %>% filter(chrom == !!chrom)
  hd <- half_decay_one(d)
  summary_rows[[i]] <- tibble(
    chrom         = chrom,
    chrom_idx     = chrom_labels[chrom],
    is_chr1       = chrom == chr1,
    n_pairs       = sum(d$n_pairs),
    mean_r2_at_100bp  = mean_r2_at(d, 100),
    mean_r2_at_5kb    = mean_r2_at(d, 5000),
    mean_r2_at_25kb   = mean_r2_at(d, 25000),
    mean_r2_at_100kb  = mean_r2_at(d, 100000),
    half_decay_bp     = hd$half_decay_bp,
    reached_threshold = hd$reached_threshold,
    min_observed_r2   = hd$min_observed_r2,
    min_observed_at_bp = hd$min_observed_at_bp
  )
}
summary_df <- bind_rows(summary_rows)
write_tsv(summary_df, out_summary, na = "NA")
log_msg(sprintf("  wrote: %s", out_summary))

log_msg("\n  per-chromosome summary:")
log_msg(sprintf("  %-12s  %12s  %8s  %8s  %8s  %8s  %12s",
                "chrom", "n_pairs", "r2@100bp", "r2@5kb", "r2@25kb", "r2@100kb",
                "half_decay"))
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  hd_str <- if (is.na(r$half_decay_bp)) "not reached"
            else sprintf("%6.0f bp", r$half_decay_bp)
  log_msg(sprintf("  chr%d %s  %12s  %8.3f  %8.3f  %8.3f  %8.3f  %12s",
                  r$chrom_idx, r$chrom,
                  format(r$n_pairs, big.mark = ","),
                  r$mean_r2_at_100bp, r$mean_r2_at_5kb,
                  r$mean_r2_at_25kb, r$mean_r2_at_100kb,
                  hd_str))
}

# --- Step 5: Plotting -------------------------------------------------------

log_msg("\n[5] Plotting decay curves")

# Discrete chromosome palette: highlight chr 1, others in greyscale.
chrom_color <- function(idx, is_c1) {
  if (is_c1) "#D55E00" else scales::alpha("#5B8AB7", 0.6)
}
binned_df$chrom_label <- sprintf("chr%d", binned_df$chrom_idx)
binned_df$chrom_label <- factor(binned_df$chrom_label,
                                  levels = sprintf("chr%d", 1:9))

# Panel A: per-chromosome decay curves.
pal <- setNames(rep("#5B8AB7", 9), sprintf("chr%d", 1:9))
pal["chr1"] <- "#D55E00"
linewidths <- setNames(rep(0.55, 9), sprintf("chr%d", 1:9))
linewidths["chr1"] <- 1.05

p_curves <- ggplot(binned_df,
                   aes(x = bin_center_bp, y = mean_r2,
                       color = chrom_label, group = chrom_label,
                       linewidth = chrom_label)) +
  geom_line(alpha = 0.95) +
  geom_point(size = 1.4) +
  geom_hline(yintercept = HALF_DECAY_THRESHOLD,
             linetype = "dashed", color = "grey20", linewidth = 0.4) +
  scale_x_log10(breaks = c(10, 100, 1e3, 1e4, 1e5),
                labels = c("10 bp", "100 bp", "1 kb", "10 kb", "100 kb")) +
  scale_color_manual(values = pal, name = NULL) +
  scale_linewidth_manual(values = linewidths, guide = "none") +
  labs(
    title    = sprintf("Per-chromosome LD decay (population = %s, n=%d)",
                       args$population, length(ids)),
    subtitle = sprintf("Mean r^2 within %d distance bins (log spaced); MAF >= %.2f; pairs within %d kb. Dashed line: r^2 = %.2f.",
                       n_bins, MAF_MIN, LD_WINDOW_KB, HALF_DECAY_THRESHOLD),
    x = "Pairwise SNP distance (bp; log scale)",
    y = expression(Mean ~ r^2)
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9),
        legend.position = "right")

ggsave(file.path(fig_dir,
                 sprintf("per_chrom_decay_curves_%s.png", args$population)),
       p_curves, width = 220, height = 130, units = "mm", dpi = 200)
log_msg(sprintf("  wrote: per_chrom_decay_curves_%s.png", args$population))

# Panel B: chr 1 vs pooled-others. Aggregate non-chr-1 chromosomes by
# weighted mean (weight = n_pairs in bin).
pooled_others <- binned_df %>%
  filter(!is_chr1) %>%
  group_by(bin_idx, bin_center_bp) %>%
  summarise(mean_r2 = weighted.mean(mean_r2, w = n_pairs, na.rm = TRUE),
            n_pairs = sum(n_pairs),
            .groups = "drop") %>%
  mutate(group = "other autosomes (pooled)")
chr1_only <- binned_df %>%
  filter(is_chr1) %>%
  transmute(bin_idx, bin_center_bp, mean_r2, n_pairs,
            group = "chr 1")
pooled_df <- bind_rows(chr1_only, pooled_others)

p_chr1 <- ggplot(pooled_df,
                 aes(x = bin_center_bp, y = mean_r2,
                     color = group, group = group)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_hline(yintercept = HALF_DECAY_THRESHOLD,
             linetype = "dashed", color = "grey20", linewidth = 0.4) +
  scale_x_log10(breaks = c(10, 100, 1e3, 1e4, 1e5),
                labels = c("10 bp", "100 bp", "1 kb", "10 kb", "100 kb")) +
  scale_color_manual(values = c("chr 1" = "#D55E00",
                                 "other autosomes (pooled)" = "#5B8AB7"),
                     name = NULL) +
  labs(
    title    = sprintf("LD decay: chr 1 vs other autosomes pooled (population = %s)",
                       args$population),
    subtitle = sprintf("Pooled others = weighted mean r^2 across chr 2-9. Dashed line: r^2 = %.2f.",
                       HALF_DECAY_THRESHOLD),
    x = "Pairwise SNP distance (bp; log scale)",
    y = expression(Mean ~ r^2)
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9),
        legend.position = "right")
ggsave(file.path(fig_dir,
                 sprintf("chr1_vs_others_decay_%s.png", args$population)),
       p_chr1, width = 220, height = 120, units = "mm", dpi = 200)
log_msg(sprintf("  wrote: chr1_vs_others_decay_%s.png", args$population))

# Panel C: half-decay distance per chromosome.
hd_plot_df <- summary_df %>%
  mutate(half_decay_kb = half_decay_bp / 1000,
         label_str = if_else(reached_threshold,
                             sprintf("%.1f kb", half_decay_kb),
                             "not reached"))
p_hd <- ggplot(hd_plot_df,
               aes(x = factor(chrom_idx),
                   y = ifelse(is.na(half_decay_kb),
                              max(half_decay_kb, na.rm = TRUE) * 1.1,
                              half_decay_kb),
                   fill = is_chr1)) +
  geom_col(color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = label_str), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("FALSE" = "#5B8AB7", "TRUE" = "#D55E00"),
                     guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = sprintf("LD half-decay distance per chromosome (r^2 < %.2f, population = %s)",
                    HALF_DECAY_THRESHOLD, args$population),
    subtitle = "Bars not reaching threshold within 100 kb shown at the maximum bin's location for visual reference.",
    x = "Chromosome (1-9; autosomes only)",
    y = sprintf("Distance at which mean r^2 first drops below %.2f (kb)",
                HALF_DECAY_THRESHOLD)
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8.5, color = "grey20"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank())
ggsave(file.path(fig_dir,
                 sprintf("half_decay_distances_%s.png", args$population)),
       p_hd, width = 200, height = 110, units = "mm", dpi = 200)
log_msg(sprintf("  wrote: half_decay_distances_%s.png", args$population))

# --- Run log ----------------------------------------------------------------

writeLines(log_lines, out_log)
log_msg(sprintf("\nRun log: %s", out_log))
