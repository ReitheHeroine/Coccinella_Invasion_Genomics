#!/usr/bin/env Rscript

# title: manhattan_mirrored_dataprep.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-04-28
# last modified: 2026-04-28
#
# purpose:
#   Data prep for the rebuilt mirrored two-panel Manhattan (Figure 1).
#   Builds per-SNP max-FST vectors per axis from the source per-SNP FST
#   tables, applies the CHI genotype-count filter to the native-vs-invasive
#   axis, computes per-axis top-1% thresholds, identifies the 2+ method
#   outlier SNPs (highlight layer), and resolves representative SNPs per
#   featured gene. Also adds a functional_theme column to the featured
#   candidate table. Plotting is deferred to a separate step pending the
#   author's check-in on the data-prep numbers.
#
# inputs:
#   - results/fst/{COMP}/fst_genomewide/{CHROM}.filtered.weir.fst
#       Per-SNP Weir & Cockerham FST. Columns: CHROM, POS,
#       WEIR_AND_COCKERHAM_FST. One file per (comparison, chromosome).
#       Comparisons: CHI_vs_EEU, CHI_vs_USA, CHI_vs_WEU,
#       EEU_vs_USA, EEU_vs_WEU, USA_vs_WEU.
#   - results/per_snp_selection_scan/chi_genotype_counts/{CHROM}.lmiss
#       Per-SNP CHI genotype counts (vcftools --missing-site). Columns:
#       CHR, POS, N_DATA (alleles, =10 for 5 CHI individuals),
#       N_GENOTYPE_FILTERED, N_MISS, F_MISS. CHI filter: keep if
#       N_DATA - N_MISS >= 6 (i.e., >= 3 of 5 CHI individuals genotyped).
#   - results/cross_method_concordance/master_candidate_table.tsv
#       Cross-method outlier table; n_methods column drives the 2+ method
#       outlier highlight layer.
#   - results/enrichment/featured_candidate_table.tsv
#       10 featured genes with chromosome, gene_start, gene_end,
#       functional_category. Used to find per-axis representative SNPs and
#       receives a new functional_theme column (in-place update).
#
# outputs:
#   - results/enrichment/plots/manhattan_axis_dataprep/
#       per_snp_axis_fst.tsv           : SNP-level table with native and
#                                        among-invasive max FST, cum_pos,
#                                        outlier-highlight flag.
#       featured_rep_snps.tsv          : Per-featured-gene representative
#                                        SNPs per axis (top and bottom panel).
#       axis_thresholds.tsv            : Per-axis top-1% threshold values.
#       chrom_offsets.tsv              : Cumulative chromosome offsets +
#                                        center positions for x-axis labels.
#       qc_summary.txt                 : Human-readable QC summary printed
#                                        in the check-in.
#   - results/enrichment/featured_candidate_table.tsv (updated in place,
#     adds functional_theme column).
#
# usage example:
#   Rscript scripts/manhattan_mirrored_dataprep.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# --- Paths ------------------------------------------------------------------

project_root <- if (file.exists("results/enrichment/foreground_genes.tsv")) {
  "."
} else if (file.exists("../results/enrichment/foreground_genes.tsv")) {
  ".."
} else {
  stop("Cannot locate results/enrichment/. Run from project root or scripts/.")
}

fst_base       <- file.path(project_root, "results/fst")
chi_lmiss_dir  <- file.path(project_root, "results/per_snp_selection_scan/chi_genotype_counts")
master_in      <- file.path(project_root, "results/cross_method_concordance/master_candidate_table.tsv")
featured_in    <- file.path(project_root, "results/enrichment/featured_candidate_table.tsv")

out_dir        <- file.path(project_root, "results/enrichment/plots/manhattan_axis_dataprep")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

per_snp_out    <- file.path(out_dir, "per_snp_axis_fst.tsv")
rep_snp_out    <- file.path(out_dir, "featured_rep_snps.tsv")
thresh_out     <- file.path(out_dir, "axis_thresholds.tsv")
offsets_out    <- file.path(out_dir, "chrom_offsets.tsv")
qc_out         <- file.path(out_dir, "qc_summary.txt")

# --- Configuration ----------------------------------------------------------

chi_comparisons      <- c("CHI_vs_EEU", "CHI_vs_USA", "CHI_vs_WEU")
invasive_comparisons <- c("EEU_vs_USA", "EEU_vs_WEU", "USA_vs_WEU")
all_comparisons      <- c(chi_comparisons, invasive_comparisons)

# 9 autosomes; sex chromosome NC_058198.1 excluded throughout (consistent
# with cross_method_concordance.R and the genome-wide FST treatment).
autosomes <- sprintf("NC_0581%02d.1", 89:97)

# CHI filter: at least 3 of 5 CHI individuals genotyped at the site. With
# N_DATA = 2 * 5 = 10 alleles, that translates to N_DATA - N_MISS >= 6.
chi_min_alleles_called <- 6

# --- Logging ----------------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

log_msg("=== Mirrored Manhattan data prep ===")
log_msg(sprintf("project root: %s", normalizePath(project_root)))
log_msg(sprintf("autosomes (n=%d): %s", length(autosomes),
                paste(autosomes, collapse = ", ")))

# --- Step 1: Load per-SNP FST for all 6 comparisons -------------------------

log_msg("\n[1/6] Loading per-SNP FST for 6 comparisons x 9 autosomes")

read_fst_one <- function(comp, chrom) {
  f <- file.path(fst_base, comp, "fst_genomewide",
                 sprintf("%s.filtered.weir.fst", chrom))
  if (!file.exists(f)) {
    stop(sprintf("missing per-SNP FST file: %s", f))
  }
  dt <- fread(f, sep = "\t", header = TRUE,
              colClasses = c(CHROM = "character",
                             POS = "integer",
                             WEIR_AND_COCKERHAM_FST = "double"))
  dt[, comparison := comp]
  dt[!is.na(WEIR_AND_COCKERHAM_FST)]
}

fst_list <- list()
for (comp in all_comparisons) {
  comp_chunks <- lapply(autosomes, read_fst_one, comp = comp)
  fst_list[[comp]] <- rbindlist(comp_chunks)
  log_msg(sprintf("  %-12s %d SNPs", comp, nrow(fst_list[[comp]])))
}

# --- Step 2: CHI genotype-count filter --------------------------------------

log_msg("\n[2/6] Loading CHI genotype counts and building filter set")

read_lmiss_one <- function(chrom) {
  f <- file.path(chi_lmiss_dir, sprintf("%s.lmiss", chrom))
  if (!file.exists(f)) {
    stop(sprintf("missing lmiss file: %s", f))
  }
  dt <- fread(f, sep = "\t", header = TRUE,
              colClasses = c(CHR = "character",
                             POS = "integer",
                             N_DATA = "integer",
                             N_GENOTYPE_FILTERED = "integer",
                             N_MISS = "integer",
                             F_MISS = "double"))
  setnames(dt, "CHR", "CHROM")
  dt
}

lmiss <- rbindlist(lapply(autosomes, read_lmiss_one))
lmiss[, n_alleles_called := N_DATA - N_MISS]
lmiss[, chi_pass := n_alleles_called >= chi_min_alleles_called]

chi_pass_count <- sum(lmiss$chi_pass)
log_msg(sprintf("  total autosomal SNPs in lmiss: %d", nrow(lmiss)))
log_msg(sprintf("  CHI filter pass (>=3/5 CHI genotyped): %d (%.2f%%)",
                chi_pass_count, 100 * chi_pass_count / nrow(lmiss)))

# Build a lookup keyed by (CHROM, POS) for fast filtering.
chi_pass_keys <- lmiss[chi_pass == TRUE, .(CHROM, POS)]
setkey(chi_pass_keys, CHROM, POS)

# --- Step 3: Per-axis max FST per SNP ---------------------------------------

log_msg("\n[3/6] Computing per-SNP max FST per axis")

# Native-vs-invasive axis: 3 CHI-vs-X comparisons, CHI-filtered.
native_long <- rbindlist(fst_list[chi_comparisons])
setkey(native_long, CHROM, POS)
# Inner join on CHI-pass keys.
native_long <- native_long[chi_pass_keys, nomatch = 0L]
native_axis <- native_long[, .(
  native_max_fst = max(WEIR_AND_COCKERHAM_FST),
  n_native_comps = .N
), by = .(CHROM, POS)]
setkey(native_axis, CHROM, POS)

# Among-invasive axis: 3 invasive-vs-invasive comparisons, no CHI filter.
inv_long <- rbindlist(fst_list[invasive_comparisons])
inv_axis <- inv_long[, .(
  invasive_max_fst = max(WEIR_AND_COCKERHAM_FST),
  n_invasive_comps = .N
), by = .(CHROM, POS)]
setkey(inv_axis, CHROM, POS)

log_msg(sprintf("  native-axis SNPs (post CHI-filter, max across 3 CHI-vs-X): %d",
                nrow(native_axis)))
log_msg(sprintf("  among-invasive axis SNPs (max across 3 inv-vs-inv): %d",
                nrow(inv_axis)))

# Outer-merge for the combined plotting table. SNPs may have one axis or both.
per_snp <- merge(native_axis, inv_axis, by = c("CHROM", "POS"), all = TRUE)
log_msg(sprintf("  union (either-axis) SNPs: %d", nrow(per_snp)))
log_msg(sprintf("  intersection (both-axis) SNPs: %d",
                sum(!is.na(per_snp$native_max_fst) &
                    !is.na(per_snp$invasive_max_fst))))

# Range and median per axis.
nat_vals <- per_snp$native_max_fst[!is.na(per_snp$native_max_fst)]
inv_vals <- per_snp$invasive_max_fst[!is.na(per_snp$invasive_max_fst)]

log_msg(sprintf("  native-axis max FST: range=[%.4f, %.4f]  median=%.4f  mean=%.4f",
                min(nat_vals), max(nat_vals), median(nat_vals), mean(nat_vals)))
log_msg(sprintf("  invasive-axis max FST: range=[%.4f, %.4f]  median=%.4f  mean=%.4f",
                min(inv_vals), max(inv_vals), median(inv_vals), mean(inv_vals)))

# --- Step 4: Per-axis top-1% thresholds -------------------------------------

log_msg("\n[4/6] Computing per-axis top-1% thresholds")

native_top1   <- as.numeric(quantile(nat_vals, 0.99))
invasive_top1 <- as.numeric(quantile(inv_vals, 0.99))

log_msg(sprintf("  native axis  top-1%% threshold: %.4f", native_top1))
log_msg(sprintf("  invasive axis top-1%% threshold: %.4f", invasive_top1))
log_msg(sprintf("  native axis SNPs >= threshold: %d (%.3f%%)",
                sum(nat_vals >= native_top1),
                100 * sum(nat_vals >= native_top1) / length(nat_vals)))
log_msg(sprintf("  invasive axis SNPs >= threshold: %d (%.3f%%)",
                sum(inv_vals >= invasive_top1),
                100 * sum(inv_vals >= invasive_top1) / length(inv_vals)))

thresh_df <- tibble(
  axis = c("native_vs_invasive", "among_invasive"),
  panel = c("top", "bottom"),
  top1pct_threshold = c(native_top1, invasive_top1),
  n_snps_axis = c(length(nat_vals), length(inv_vals)),
  n_at_or_above_threshold = c(sum(nat_vals >= native_top1),
                              sum(inv_vals >= invasive_top1))
)
write_tsv(thresh_df, thresh_out, na = "NA")

# --- Step 5: Cumulative position + 2+ method outlier flag -------------------

log_msg("\n[5/6] Building cumulative x-axis and outlier highlight flag")

# Chromosome lengths from per-SNP data (max POS observed). Add a small pad
# between chromosomes (matches existing plotdata convention -- offset is
# cumsum of preceding chrom lengths with no inter-chrom padding).
chrom_max_pos <- per_snp[, .(max_pos = max(POS, na.rm = TRUE)), by = CHROM]
chrom_max_pos <- chrom_max_pos[order(match(CHROM, autosomes))]
chrom_max_pos[, offset := cumsum(c(0, head(max_pos, -1)))]

setkey(chrom_max_pos, CHROM)
per_snp <- merge(per_snp, chrom_max_pos[, .(CHROM, offset)], by = "CHROM", all.x = TRUE)
per_snp[, cum_pos := POS + offset]

# Center positions for x-axis labels (per chromosome midpoint along cum_pos).
chrom_centers <- per_snp[, .(center = (min(cum_pos) + max(cum_pos)) / 2,
                             min_cum = min(cum_pos),
                             max_cum = max(cum_pos)), by = CHROM]
chrom_centers <- chrom_centers[order(match(CHROM, autosomes))]
chrom_centers[, chrom_idx := seq_len(.N)]
chrom_centers[, label := sprintf("%d", chrom_idx)]

# 2+ method outliers from the master candidate table (already excludes the
# sex chromosome and is the canonical foreground criterion).
master <- fread(master_in, sep = "\t", header = TRUE)
outlier_keys <- master[n_methods >= 2, .(CHROM = chrom, POS = pos)]
log_msg(sprintf("  2+ method outliers in master table: %d", nrow(outlier_keys)))
log_msg(sprintf("    n_methods distribution (>=2):"))
nm_tab <- table(master[n_methods >= 2]$n_methods)
for (lvl in names(nm_tab)) {
  log_msg(sprintf("      n_methods = %s : %d", lvl, nm_tab[[lvl]]))
}

setkey(outlier_keys, CHROM, POS)
setkey(per_snp, CHROM, POS)
per_snp[, is_outlier_2plus := FALSE]
per_snp[outlier_keys, is_outlier_2plus := TRUE]

# How many outliers landed on each axis (should be high but may not be 100%
# because some master-table SNPs are CHI-axis-only and may have been removed
# by the CHI filter, etc.).
n_out_on_native <- sum(per_snp$is_outlier_2plus & !is.na(per_snp$native_max_fst))
n_out_on_inv    <- sum(per_snp$is_outlier_2plus & !is.na(per_snp$invasive_max_fst))
n_out_total     <- sum(per_snp$is_outlier_2plus)
log_msg(sprintf("  outlier flag on per_snp: %d total; %d have native-axis FST; %d have invasive-axis FST",
                n_out_total, n_out_on_native, n_out_on_inv))

# --- Step 6: Featured-gene representative SNPs per axis ---------------------

log_msg("\n[6/6] Resolving featured-gene representative SNPs per axis")

featured <- read_tsv(featured_in, show_col_types = FALSE)
log_msg(sprintf("  featured genes loaded: %d", nrow(featured)))
stopifnot(nrow(featured) == 10)

# Functional theme mapping (handoff section 1.4).
theme_map <- c(
  "Fatty acid metabolism"             = "Metabolism",
  "Lipid metabolism / membrane"       = "Metabolism",
  "Metabolic transport"               = "Metabolism",
  "Ion channel / K+ transport"        = "Membrane / ion transport",
  "Ion homeostasis"                   = "Membrane / ion transport",
  "Inositol / phospholipid signaling" = "Cell signaling",
  "Immune / Toll signaling"           = "Cell signaling",
  "Growth / proliferation"            = "Growth & development",
  "Cytoskeletal / developmental"      = "Growth & development",
  "Sensory / chemoreception"          = "Sensory"
)

unmapped <- setdiff(featured$functional_category, names(theme_map))
if (length(unmapped) > 0) {
  stop(sprintf("featured table has unmapped functional_category values: %s",
               paste(unmapped, collapse = "; ")))
}

featured$functional_theme <- theme_map[featured$functional_category]
theme_counts <- featured %>% count(functional_theme, name = "n_genes") %>% arrange(desc(n_genes))
log_msg("  per-theme gene counts:")
for (i in seq_len(nrow(theme_counts))) {
  log_msg(sprintf("    %-25s %d", theme_counts$functional_theme[i],
                  theme_counts$n_genes[i]))
}

# Write featured table back in place with theme appended.
write_tsv(featured, featured_in, na = "NA")
log_msg(sprintf("  updated featured table written: %s", featured_in))

# Resolve per-axis rep SNPs by finding max-FST SNP within each gene's
# coordinate range. May be different SNPs on the two axes.
rep_rows <- vector("list", nrow(featured))
for (i in seq_len(nrow(featured))) {
  g <- featured[i, ]
  in_gene <- per_snp[CHROM == g$chromosome &
                       POS >= g$gene_start &
                       POS <= g$gene_end]

  if (nrow(in_gene) == 0) {
    log_msg(sprintf("  WARNING: no per-SNP rows in gene %s (%s:%d-%d)",
                    g$gene_id, g$chromosome, g$gene_start, g$gene_end))
    next
  }

  # Top-panel rep: max native-axis FST within the gene.
  has_nat <- !is.na(in_gene$native_max_fst)
  if (any(has_nat)) {
    nat_idx <- which(has_nat)[which.max(in_gene$native_max_fst[has_nat])]
    nat_pos       <- in_gene$POS[nat_idx]
    nat_cum       <- in_gene$cum_pos[nat_idx]
    nat_fst       <- in_gene$native_max_fst[nat_idx]
    nat_inv_atpos <- in_gene$invasive_max_fst[nat_idx]
  } else {
    nat_pos <- NA_integer_; nat_cum <- NA_real_
    nat_fst <- NA_real_;    nat_inv_atpos <- NA_real_
  }

  # Bottom-panel rep: max invasive-axis FST within the gene.
  has_inv <- !is.na(in_gene$invasive_max_fst)
  if (any(has_inv)) {
    inv_idx <- which(has_inv)[which.max(in_gene$invasive_max_fst[has_inv])]
    inv_pos       <- in_gene$POS[inv_idx]
    inv_cum       <- in_gene$cum_pos[inv_idx]
    inv_fst       <- in_gene$invasive_max_fst[inv_idx]
    inv_nat_atpos <- in_gene$native_max_fst[inv_idx]
  } else {
    inv_pos <- NA_integer_; inv_cum <- NA_real_
    inv_fst <- NA_real_;    inv_nat_atpos <- NA_real_
  }

  rep_rows[[i]] <- tibble(
    gene_id              = g$gene_id,
    chromosome           = g$chromosome,
    gene_start           = g$gene_start,
    gene_end             = g$gene_end,
    selection_block      = g$selection_block,
    functional_category  = g$functional_category,
    functional_theme     = g$functional_theme,
    functional_annotation = g$functional_annotation,
    n_snps_in_gene       = nrow(in_gene),
    top_rep_pos          = nat_pos,
    top_rep_cum_pos      = nat_cum,
    top_rep_native_fst   = nat_fst,
    top_rep_invasive_fst = nat_inv_atpos,
    bot_rep_pos          = inv_pos,
    bot_rep_cum_pos      = inv_cum,
    bot_rep_invasive_fst = inv_fst,
    bot_rep_native_fst   = inv_nat_atpos,
    same_snp_both_axes   = !is.na(nat_pos) & !is.na(inv_pos) & nat_pos == inv_pos
  )
}
rep_df <- bind_rows(rep_rows)

log_msg("  featured-gene representative SNPs:")
for (i in seq_len(nrow(rep_df))) {
  r <- rep_df[i, ]
  log_msg(sprintf(
    "    %-14s [%s] block=%-20s n_SNPs=%-4d  top@%d (nat=%.3f, inv@pos=%s)  bot@%d (inv=%.3f, nat@pos=%s)  same=%s",
    r$gene_id,
    r$functional_theme,
    as.character(r$selection_block),
    r$n_snps_in_gene,
    if (is.na(r$top_rep_pos)) 0L else r$top_rep_pos,
    if (is.na(r$top_rep_native_fst)) NA_real_ else r$top_rep_native_fst,
    if (is.na(r$top_rep_invasive_fst)) "NA" else format(round(r$top_rep_invasive_fst, 3), nsmall = 3),
    if (is.na(r$bot_rep_pos)) 0L else r$bot_rep_pos,
    if (is.na(r$bot_rep_invasive_fst)) NA_real_ else r$bot_rep_invasive_fst,
    if (is.na(r$bot_rep_native_fst)) "NA" else format(round(r$bot_rep_native_fst, 3), nsmall = 3),
    r$same_snp_both_axes
  ))
}

write_tsv(rep_df, rep_snp_out, na = "NA")
write_tsv(as_tibble(chrom_centers), offsets_out, na = "NA")

# Final per-SNP table: keep light columns to keep file size reasonable.
per_snp_out_df <- as_tibble(per_snp) %>%
  select(CHROM, POS, cum_pos,
         native_max_fst, invasive_max_fst,
         is_outlier_2plus)
write_tsv(per_snp_out_df, per_snp_out, na = "NA")

log_msg(sprintf("\nWrote: %s (%d rows)", per_snp_out, nrow(per_snp_out_df)))
log_msg(sprintf("Wrote: %s (%d rows)", rep_snp_out, nrow(rep_df)))
log_msg(sprintf("Wrote: %s",          thresh_out))
log_msg(sprintf("Wrote: %s",          offsets_out))

writeLines(log_lines, qc_out)
log_msg(sprintf("Wrote: %s", qc_out))
