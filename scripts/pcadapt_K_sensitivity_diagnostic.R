#!/usr/bin/env Rscript

# title: pcadapt_K_sensitivity_diagnostic.R
# project: BIOL624 Final Project -- Selection Detection in Lady Beetles
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-04-29
# last modified: 2026-04-29
#
# purpose:
#   Diagnostic battery to support a re-evaluation of the K=2 choice for the
#   pcadapt selection scan. Six sub-tasks:
#     (1) PC interpretation: per-pop mean/SD of PC1 and PC2 scores; identify
#         samples driving PC2; drop-one sensitivity.
#     (2) K=1 parallel run: re-run pcadapt with K=1, save outliers at q<0.05
#         and q<0.01 with _K1 suffix; report lambda and counts.
#     (3) K-sensitivity of cross-method concordance: intersection / set-diff
#         between K=1 and K=2 outlier sets; rebuild a SHADOW master
#         candidate table that substitutes K=1 for K=2 (without overwriting
#         the canonical master_candidate_table.tsv); quantify how SNPs
#         change n_methods status broken down by consensus_pattern.
#     (4) Fix scree annotation: drop the max(2, .) floor in the auto K
#         detection; report strict Cattell, Kaiser-Guttman, and Horn's
#         parallel analysis K side-by-side; redraw scree plot with all
#         three criteria annotated.
#     (5) PC2 sanity check: correlation between |PC2 loading| and per-SNP
#         among-invasive max FST; binned enrichment summary.
#     (6) Optional K=3 run: lambda + outlier counts + PC3 score plot.
#
#   IMPORTANT: this script does NOT modify any K=2 output, nor the canonical
#   master_candidate_table.tsv. All diagnostic outputs are written under
#   results/pcadapt/diagnostics/.
#
# inputs:
#   - data/plink/ladybug_snps_pcadapt.{bed,bim,fam} (81 individuals, 803,821 SNPs)
#   - data/plink/ladybug_snps_pcadapt.imiss
#   - metadata/pop_{CHI,EEU,WEU,USA}.txt
#   - results/pcadapt/pcadapt_outliers_q05.tsv      (canonical K=2 outliers)
#   - results/pcadapt/pcadapt_outliers_q01.tsv
#   - results/per_snp_selection_scan/snp_summary_1pct_chifiltered.tsv
#   - results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv
#   - results/outflank_ldpruned/{EEU_vs_WEU,EEU_vs_USA,USA_vs_WEU}/outflank_full_results.rds
#   - results/enrichment/plots/manhattan_axis_dataprep/per_snp_axis_fst.tsv
#       (per-SNP among-invasive max FST from the Figure 1 dataprep)
#   - results/cross_method_concordance/master_candidate_table.tsv
#       (canonical K=2-based master, used for diff comparison only; never
#       overwritten)
#
# outputs:
#   results/pcadapt/diagnostics/
#     # Task 1
#     pc_scores_overspecified.tsv       (per-sample PC1..PC10 from K=10 fit)
#     pc_per_pop_stats.tsv              (per-pop mean/SD on PC1..PC4)
#     pc2_extreme_samples.tsv           (|PC2 z-score| > 2 sample list)
#     pca_drop_one_sensitivity.tsv      (variance explained with each
#                                        candidate drop sample removed)
#     # Task 2
#     pcadapt_outliers_q05_K1.tsv
#     pcadapt_outliers_q01_K1.tsv
#     pcadapt_K1_summary.txt
#     # Task 3
#     K_set_comparison.tsv              (Venn-style intersection summary)
#     master_candidate_table_K1.tsv     (shadow master with K=1 substituted)
#     master_K_diff_summary.txt         (how n_methods changes per consensus_pattern)
#     # Task 4
#     scree_diagnostic.tsv              (per-PC variance + cumulative + parallel
#                                        analysis null + Kaiser threshold)
#     scree_diagnostic.png              (scree plot with all 3 K criteria)
#     scree_K_summary.txt               (strict Cattell, Kaiser, Horn K values)
#     # Task 5
#     pc2_loading_vs_fst.tsv            (per-SNP |PC2 loading| + among-inv FST)
#     pc2_loading_vs_fst_summary.txt    (Spearman correlation, binned enrichment)
#     pc2_loading_vs_fst.png            (scatter / hexbin)
#     # Task 6
#     pcadapt_K3_summary.txt
#     pca_K3_scores.png                 (PC1xPC2, PC1xPC3, PC2xPC3)
#     # Synthesis
#     diagnostic_report.md
#
# usage example:
#   Rscript scripts/pcadapt_K_sensitivity_diagnostic.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(pcadapt)
  library(qvalue)
  library(data.table)
})

set.seed(42)

# --- Paths ------------------------------------------------------------------

project_root <- if (file.exists("data/plink/ladybug_snps_pcadapt.bed")) {
  "."
} else if (file.exists("../data/plink/ladybug_snps_pcadapt.bed")) {
  ".."
} else {
  stop("Cannot locate data/plink/. Run from project root or scripts/.")
}

PLINK_PREFIX  <- file.path(project_root, "data/plink/ladybug_snps_pcadapt")
IMISS_FILE    <- file.path(project_root, "data/plink/ladybug_snps_pcadapt.imiss")
POP_DIR       <- file.path(project_root, "metadata")
OUT_DIR       <- file.path(project_root, "results/pcadapt/diagnostics")
PCADAPT_DIR   <- file.path(project_root, "results/pcadapt")
MASTER_TBL    <- file.path(project_root, "results/cross_method_concordance/master_candidate_table.tsv")
PCT_SUMMARY   <- file.path(project_root, "results/per_snp_selection_scan/snp_summary_1pct_chifiltered.tsv")
PCT_DETAIL    <- file.path(project_root, "results/per_snp_selection_scan/outliers_1pct_chifiltered.tsv")
OUTFLANK_RDS  <- list(
  EEU_vs_WEU = file.path(project_root, "results/outflank_ldpruned/EEU_vs_WEU/outflank_full_results.rds"),
  EEU_vs_USA = file.path(project_root, "results/outflank_ldpruned/EEU_vs_USA/outflank_full_results.rds"),
  USA_vs_WEU = file.path(project_root, "results/outflank_ldpruned/USA_vs_WEU/outflank_full_results.rds")
)
PER_SNP_FST   <- file.path(project_root, "results/enrichment/plots/manhattan_axis_dataprep/per_snp_axis_fst.tsv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEX_CHROM <- "NC_058198.1"

# --- Logging helper ---------------------------------------------------------

log_lines <- character()
log_msg <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n", sep = "")
}

# --- Load PLINK metadata ----------------------------------------------------

log_msg("\n=== pcadapt K-sensitivity diagnostic battery ===")
log_msg("Loading PLINK files and population assignments")

bim <- fread(paste0(PLINK_PREFIX, ".bim"), header = FALSE,
             col.names = c("chrom", "snp_id", "cm", "pos", "a1", "a2"))
fam <- fread(paste0(PLINK_PREFIX, ".fam"), header = FALSE,
             col.names = c("fid", "iid", "pid", "mid", "sex", "phen"))
imiss <- fread(IMISS_FILE)

pop_map <- bind_rows(
  tibble(iid = readLines(file.path(POP_DIR, "pop_CHI.txt")), pop = "CHI"),
  tibble(iid = readLines(file.path(POP_DIR, "pop_EEU.txt")), pop = "EEU"),
  tibble(iid = readLines(file.path(POP_DIR, "pop_WEU.txt")), pop = "WEU"),
  tibble(iid = readLines(file.path(POP_DIR, "pop_USA.txt")), pop = "USA")
)
fam_pop <- fam %>% left_join(pop_map, by = "iid")
stopifnot(!any(is.na(fam_pop$pop)))

log_msg(sprintf("  %d individuals (CHI=%d, EEU=%d, WEU=%d, USA=%d)",
                nrow(fam_pop),
                sum(fam_pop$pop == "CHI"),
                sum(fam_pop$pop == "EEU"),
                sum(fam_pop$pop == "WEU"),
                sum(fam_pop$pop == "USA")))
log_msg(sprintf("  %d SNPs in .bim", nrow(bim)))

# Pre-build snp_id <-> (chrom, pos) lookup.
bim_lookup <- bim %>%
  mutate(snp_id_cp = paste0(chrom, ":", pos)) %>%
  select(snp_idx = NULL, snp_id_bim = snp_id, chrom, pos, snp_id_cp)
bim_lookup$snp_idx <- seq_len(nrow(bim_lookup))

################################################################################
# TASK 1: PC interpretation check
################################################################################

log_msg("\n--- TASK 1: PC interpretation check ---")

geno <- read.pcadapt(paste0(PLINK_PREFIX, ".bed"), type = "bed")

log_msg("  Running pcadapt at K=10 (overspecified) for variance / scores")
fit_K10 <- pcadapt(geno, K = 10)

# Variance explained, scores, loadings.
sing <- fit_K10$singular.values
var_exp <- sing^2 / sum(sing^2)
log_msg(sprintf("  Variance explained by PC1..PC10: %s",
                paste(sprintf("%.3f", var_exp), collapse = " ")))

scores <- as.data.frame(fit_K10$scores)
colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
scores$iid <- fam_pop$iid
scores$pop <- fam_pop$pop
scores <- scores %>% select(iid, pop, everything())
write_tsv(scores, file.path(OUT_DIR, "pc_scores_overspecified.tsv"))
log_msg(sprintf("  wrote %s", file.path(OUT_DIR, "pc_scores_overspecified.tsv")))

per_pop_stats <- scores %>%
  pivot_longer(starts_with("PC"), names_to = "pc", values_to = "score") %>%
  group_by(pop, pc) %>%
  summarise(
    n         = n(),
    mean_pc   = mean(score),
    sd_pc     = sd(score),
    min_pc    = min(score),
    max_pc    = max(score),
    .groups   = "drop"
  ) %>%
  mutate(pc = factor(pc, levels = paste0("PC", 1:10))) %>%
  arrange(pc, pop)
write_tsv(per_pop_stats, file.path(OUT_DIR, "pc_per_pop_stats.tsv"))

log_msg("  Per-pop mean (SD) for PC1, PC2:")
for (pc_lvl in c("PC1", "PC2")) {
  for (p in c("CHI", "EEU", "WEU", "USA")) {
    r <- per_pop_stats %>% filter(pc == pc_lvl, pop == p)
    log_msg(sprintf("    %s  %s  n=%d  mean=%+0.4f  sd=%0.4f  range=[%+0.4f, %+0.4f]",
                    pc_lvl, p, r$n, r$mean_pc, r$sd_pc, r$min_pc, r$max_pc))
  }
}

# Identify samples with extreme PC2 scores (|z| > 2 within pop).
pc2_z <- scores %>%
  group_by(pop) %>%
  mutate(pc2_z = (PC2 - mean(PC2)) / sd(PC2)) %>%
  ungroup() %>%
  arrange(desc(abs(pc2_z))) %>%
  select(iid, pop, PC1, PC2, pc2_z) %>%
  mutate(extreme = abs(pc2_z) > 2)
write_tsv(pc2_z, file.path(OUT_DIR, "pc2_extreme_samples.tsv"))

extreme_samples <- pc2_z %>% filter(extreme) %>% pull(iid)
log_msg(sprintf("  Samples with |PC2 within-pop z| > 2 (n=%d): %s",
                length(extreme_samples), paste(extreme_samples, collapse = ", ")))

# Drop-one sensitivity. Test the 4 borderline samples + the most extreme PC2
# samples not already in the borderline list.
borderline <- c("DEU12", "NLD3", "SVK1", "WUS7")
drop_candidates <- unique(c(borderline,
                            head(pc2_z$iid[pc2_z$extreme], 6)))
log_msg(sprintf("  Drop-one candidates (n=%d): %s",
                length(drop_candidates),
                paste(drop_candidates, collapse = ", ")))

drop_results <- list()
drop_results[["__none__"]] <- tibble(
  dropped       = "(none)",
  n_indv        = nrow(fam_pop),
  pc1_var_pct   = 100 * var_exp[1],
  pc2_var_pct   = 100 * var_exp[2],
  pc3_var_pct   = 100 * var_exp[3],
  pc1_pc2_ratio = var_exp[1] / var_exp[2]
)

for (drop_iid in drop_candidates) {
  if (!drop_iid %in% fam_pop$iid) {
    log_msg(sprintf("    skip %s (not in dataset)", drop_iid))
    next
  }
  keep_idx <- which(fam_pop$iid != drop_iid)
  # pcadapt accepts a subset via the bed file; easiest path is to read the
  # bed into memory and subset. The pcadapt .bed reader doesn't support an
  # `indv_keep` argument, but we can subset post-hoc by re-running on a
  # genotype matrix. genotype matrix is samples x SNPs.
  geno_mat <- bed2matrix(paste0(PLINK_PREFIX, ".bed"))
  geno_sub <- geno_mat[keep_idx, , drop = FALSE]
  fit_sub <- pcadapt(read.pcadapt(geno_sub, type = "lfmm"), K = 10)
  sing_s <- fit_sub$singular.values
  var_s  <- sing_s^2 / sum(sing_s^2)
  drop_results[[drop_iid]] <- tibble(
    dropped       = drop_iid,
    n_indv        = length(keep_idx),
    pc1_var_pct   = 100 * var_s[1],
    pc2_var_pct   = 100 * var_s[2],
    pc3_var_pct   = 100 * var_s[3],
    pc1_pc2_ratio = var_s[1] / var_s[2]
  )
  log_msg(sprintf("    drop %-7s -> PC1=%5.2f%% PC2=%5.2f%% PC3=%5.2f%%  ratio=%.3f",
                  drop_iid, 100 * var_s[1], 100 * var_s[2], 100 * var_s[3],
                  var_s[1] / var_s[2]))
}
drop_df <- bind_rows(drop_results)
write_tsv(drop_df, file.path(OUT_DIR, "pca_drop_one_sensitivity.tsv"))
log_msg(sprintf("  wrote %s", file.path(OUT_DIR, "pca_drop_one_sensitivity.tsv")))

################################################################################
# TASK 2: K=1 parallel pcadapt run
################################################################################

log_msg("\n--- TASK 2: K=1 pcadapt run ---")

fit_K1 <- pcadapt(geno, K = 1)

# Genomic inflation lambda (median observed chi^2 / median expected chi^2).
chisq_obs_K1 <- qchisq(1 - fit_K1$pvalues, df = 1)
lambda_K1 <- median(chisq_obs_K1, na.rm = TRUE) / qchisq(0.5, df = 1)
log_msg(sprintf("  Lambda (K=1):  %.4f", lambda_K1))

# Q-values via Storey.
qv_K1 <- qvalue(fit_K1$pvalues)
qvals_K1 <- qv_K1$qvalues

snp_table_K1 <- bim %>%
  select(chrom, pos, snp_id) %>%
  mutate(
    snp_index   = seq_len(nrow(bim)),
    pvalue      = fit_K1$pvalues,
    qvalue      = qvals_K1,
    z_PC1       = fit_K1$loadings[, 1] /
                   sd(fit_K1$loadings[, 1], na.rm = TRUE)  # standardize
  )

out_q05_K1 <- snp_table_K1 %>% filter(!is.na(qvalue) & qvalue < 0.05)
out_q01_K1 <- snp_table_K1 %>% filter(!is.na(qvalue) & qvalue < 0.01)

log_msg(sprintf("  K=1 outliers q<0.05: %d  (vs K=2 q<0.05: %d)",
                nrow(out_q05_K1), 42288L))
log_msg(sprintf("  K=1 outliers q<0.01: %d  (vs K=2 q<0.01: %d)",
                nrow(out_q01_K1), 30157L))

write_tsv(out_q05_K1, file.path(OUT_DIR, "pcadapt_outliers_q05_K1.tsv"))
write_tsv(out_q01_K1, file.path(OUT_DIR, "pcadapt_outliers_q01_K1.tsv"))

K1_summary <- c(
  sprintf("K = 1 (diagnostic parallel run)"),
  sprintf("Date: %s", format(Sys.time())),
  sprintf("N SNPs: %d", nrow(snp_table_K1)),
  sprintf("Lambda (genomic inflation): %.4f", lambda_K1),
  sprintf("Outliers q<0.05: %d", nrow(out_q05_K1)),
  sprintf("Outliers q<0.01: %d", nrow(out_q01_K1)),
  sprintf("Outliers on sex chrom (NC_058198.1) at q<0.05: %d",
          sum(out_q05_K1$chrom == SEX_CHROM)),
  sprintf("Autosomal q<0.05: %d",
          sum(out_q05_K1$chrom != SEX_CHROM)),
  "",
  sprintf("K=2 comparison (canonical run, for reference):"),
  sprintf("  q<0.05 = 42,288   q<0.01 = 30,157"),
  ""
)
writeLines(K1_summary, file.path(OUT_DIR, "pcadapt_K1_summary.txt"))
log_msg(sprintf("  wrote %s", file.path(OUT_DIR, "pcadapt_K1_summary.txt")))

################################################################################
# TASK 3: K-sensitivity of cross-method concordance
################################################################################

log_msg("\n--- TASK 3: K-sensitivity of concordance ---")

pcadapt_q05_K2 <- fread(file.path(PCADAPT_DIR, "pcadapt_outliers_q05.tsv"))
pcadapt_q01_K2 <- fread(file.path(PCADAPT_DIR, "pcadapt_outliers_q01.tsv"))

# Set comparison
set_K1_q05 <- out_q05_K1$snp_id
set_K2_q05 <- pcadapt_q05_K2$snp_id
set_K1_q01 <- out_q01_K1$snp_id
set_K2_q01 <- pcadapt_q01_K2$snp_id

K_compare <- tibble(
  threshold        = c("q<0.05", "q<0.01"),
  n_K1             = c(length(set_K1_q05), length(set_K1_q01)),
  n_K2             = c(length(set_K2_q05), length(set_K2_q01)),
  intersection     = c(length(intersect(set_K1_q05, set_K2_q05)),
                       length(intersect(set_K1_q01, set_K2_q01))),
  K1_only          = c(length(setdiff(set_K1_q05, set_K2_q05)),
                       length(setdiff(set_K1_q01, set_K2_q01))),
  K2_only          = c(length(setdiff(set_K2_q05, set_K1_q05)),
                       length(setdiff(set_K2_q01, set_K1_q01))),
  union_total      = c(length(union(set_K1_q05, set_K2_q05)),
                       length(union(set_K1_q01, set_K2_q01))),
  jaccard          = c(length(intersect(set_K1_q05, set_K2_q05)) /
                         length(union(set_K1_q05, set_K2_q05)),
                       length(intersect(set_K1_q01, set_K2_q01)) /
                         length(union(set_K1_q01, set_K2_q01)))
)
write_tsv(K_compare, file.path(OUT_DIR, "K_set_comparison.tsv"))
log_msg("  K=1 vs K=2 outlier set comparison:")
print(K_compare)

# Build shadow master with K=1 substituted. Use the same pipeline logic as
# cross_method_concordance.R: percentile + outflank stay the same; pcadapt
# is replaced.

# Load percentile and outflank inputs
log_msg("  Loading percentile + outflank for shadow master rebuild")
pct_summary <- fread(PCT_SUMMARY)
pct_detail  <- fread(PCT_DETAIL)
outflank_list <- lapply(names(OUTFLANK_RDS), function(comp) {
  if (!file.exists(OUTFLANK_RDS[[comp]])) return(NULL)
  fr <- readRDS(OUTFLANK_RDS[[comp]])
  res <- fr$results
  outliers <- res[!is.na(res$qvalues) & res$qvalues < 0.05, ]
  if (nrow(outliers) == 0) return(NULL)
  data.frame(
    snp_id          = paste0(outliers$LocusName, ""),
    chrom           = sub(":.*$", "", outliers$LocusName),
    pos             = as.integer(sub("^.*:", "", outliers$LocusName)),
    comparison      = comp,
    qvalue          = outliers$qvalues,
    fst             = outliers$FST,
    stringsAsFactors = FALSE
  )
})
outflank_all <- bind_rows(outflank_list)
outflank_snps <- outflank_all %>%
  group_by(snp_id, chrom, pos) %>%
  summarise(
    outflank_comparisons = paste(sort(unique(comparison)), collapse = ";"),
    outflank_q_min       = min(qvalue, na.rm = TRUE),
    .groups = "drop"
  )
log_msg(sprintf("    percentile rows: %d  outflank SNPs: %d",
                nrow(pct_summary), nrow(outflank_snps)))

# Helper: build master with a given pcadapt outlier set.
build_master <- function(pcadapt_q05_df, pcadapt_K_label) {
  pct_for_join <- pct_summary %>%
    filter(chrom != SEX_CHROM) %>%
    select(snp_id, chrom, pos,
           percentile_comparisons = comparisons,
           percentile_pattern = pattern)
  of_for_join <- outflank_snps %>% filter(chrom != SEX_CHROM)
  pa_for_join <- pcadapt_q05_df %>%
    filter(chrom != SEX_CHROM) %>%
    transmute(snp_id, chrom, pos,
              pcadapt_q = qvalue)

  master <- pct_for_join %>% select(snp_id, chrom, pos) %>%
    full_join(of_for_join %>% select(snp_id, chrom, pos), by = c("snp_id","chrom","pos")) %>%
    full_join(pa_for_join %>% select(snp_id, chrom, pos), by = c("snp_id","chrom","pos")) %>%
    distinct() %>%
    left_join(pct_for_join, by = c("snp_id","chrom","pos")) %>%
    left_join(of_for_join,  by = c("snp_id","chrom","pos")) %>%
    left_join(pa_for_join,  by = c("snp_id","chrom","pos"))

  master %>% mutate(
    percentile_outlier = !is.na(percentile_comparisons),
    outflank_outlier   = !is.na(outflank_comparisons),
    pcadapt_outlier    = !is.na(pcadapt_q),
    n_methods = as.integer(percentile_outlier) +
                 as.integer(outflank_outlier) +
                 as.integer(pcadapt_outlier),
    pct_pat = ifelse(percentile_outlier, percentile_pattern, NA_character_),
    of_pat  = ifelse(outflank_outlier, "POST_INVASION", NA_character_),
    consensus_pattern = case_when(
      !is.na(pct_pat) & !is.na(of_pat) & pct_pat == "BOTH" ~ "BOTH",
      !is.na(pct_pat) & !is.na(of_pat) & pct_pat == "POST_INVASION" ~ "POST_INVASION",
      !is.na(pct_pat) & !is.na(of_pat) & pct_pat == "INVASION_SPECIFIC" ~ "DISCORDANT",
      !is.na(pct_pat) ~ pct_pat,
      !is.na(of_pat) ~ of_pat,
      TRUE ~ NA_character_
    ),
    K_label = pcadapt_K_label
  )
}

master_K2_shadow <- build_master(pcadapt_q05_K2, "K2")
master_K1_shadow <- build_master(out_q05_K1,    "K1")

log_msg(sprintf("  shadow master rows: K1=%d  K2=%d",
                nrow(master_K1_shadow), nrow(master_K2_shadow)))

# Diff: how many SNPs change n_methods status between K=1 and K=2 master.
diff_join <- master_K2_shadow %>%
  select(snp_id, chrom, pos,
         n_methods_K2 = n_methods,
         consensus_pattern_K2 = consensus_pattern) %>%
  full_join(
    master_K1_shadow %>% select(snp_id, chrom, pos,
                                n_methods_K1 = n_methods,
                                consensus_pattern_K1 = consensus_pattern),
    by = c("snp_id", "chrom", "pos")
  ) %>%
  mutate(
    n_methods_K1 = coalesce(n_methods_K1, 0L),
    n_methods_K2 = coalesce(n_methods_K2, 0L),
    diff_n_methods = n_methods_K1 - n_methods_K2,
    consensus_pattern = coalesce(consensus_pattern_K2, consensus_pattern_K1, "(none)")
  )

diff_summary <- diff_join %>%
  count(consensus_pattern, n_methods_K2, n_methods_K1, name = "n_snps") %>%
  arrange(consensus_pattern, n_methods_K2, n_methods_K1)
write_tsv(diff_summary, file.path(OUT_DIR, "master_K_diff_summary.txt"))

# Specifically the 2+ method outliers under each K.
n_2plus_K2 <- master_K2_shadow %>% filter(n_methods >= 2) %>% nrow()
n_2plus_K1 <- master_K1_shadow %>% filter(n_methods >= 2) %>% nrow()
n_2plus_intersect <- diff_join %>% filter(n_methods_K1 >= 2 & n_methods_K2 >= 2) %>% nrow()
n_2plus_K1only <- diff_join %>% filter(n_methods_K1 >= 2 & n_methods_K2 < 2) %>% nrow()
n_2plus_K2only <- diff_join %>% filter(n_methods_K2 >= 2 & n_methods_K1 < 2) %>% nrow()

log_msg(sprintf("  2+ method outlier counts: K2=%d  K1=%d  shared=%d  K1-only=%d  K2-only=%d",
                n_2plus_K2, n_2plus_K1, n_2plus_intersect, n_2plus_K1only, n_2plus_K2only))

# Among-invasive panel impact: how many of the K2-only or K1-only 2+ outliers
# would land on the invasive axis. Use per_snp_axis_fst.tsv to look up.
if (file.exists(PER_SNP_FST)) {
  axis_fst <- fread(PER_SNP_FST)
  axis_fst <- axis_fst %>%
    transmute(snp_id = paste0(CHROM, ":", POS),
              invasive_max_fst, native_max_fst, is_outlier_2plus)
  K1_only_snps <- diff_join %>% filter(n_methods_K1 >= 2 & n_methods_K2 < 2) %>% pull(snp_id)
  K2_only_snps <- diff_join %>% filter(n_methods_K2 >= 2 & n_methods_K1 < 2) %>% pull(snp_id)
  K1_only_inv <- axis_fst %>%
    filter(snp_id %in% K1_only_snps) %>%
    summarise(n = n(), median_inv = median(invasive_max_fst, na.rm = TRUE),
              n_above_thresh_0.19 = sum(invasive_max_fst >= 0.19, na.rm = TRUE))
  K2_only_inv <- axis_fst %>%
    filter(snp_id %in% K2_only_snps) %>%
    summarise(n = n(), median_inv = median(invasive_max_fst, na.rm = TRUE),
              n_above_thresh_0.19 = sum(invasive_max_fst >= 0.19, na.rm = TRUE))
  log_msg(sprintf("  K1-only 2+ outliers among-inv: n=%d  median FST=%.4f  n>=0.190=%d",
                  K1_only_inv$n, K1_only_inv$median_inv, K1_only_inv$n_above_thresh_0.19))
  log_msg(sprintf("  K2-only 2+ outliers among-inv: n=%d  median FST=%.4f  n>=0.190=%d",
                  K2_only_inv$n, K2_only_inv$median_inv, K2_only_inv$n_above_thresh_0.19))
}

write_tsv(master_K1_shadow %>% filter(n_methods >= 1) %>% arrange(chrom, pos),
          file.path(OUT_DIR, "master_candidate_table_K1.tsv"))
log_msg(sprintf("  wrote %s",
                file.path(OUT_DIR, "master_candidate_table_K1.tsv")))

################################################################################
# TASK 4: Fix scree annotation (strict Cattell + Kaiser + parallel analysis)
################################################################################

log_msg("\n--- TASK 4: scree annotation diagnostic ---")

# Strict Cattell: K = position of the elbow (largest drop, with the elbow
# *separating* high from low). Convention here: K_strict = number of PCs
# above the plateau. The plateau is identified as the first PC after which
# subsequent drops are < 1/3 of the PC1 -> PC2 drop. That is the strict
# scree-test reading without the max(2, .) floor that appeared in the
# original auto-detection.
diffs_var <- diff(var_exp)
elbow_threshold <- abs(diffs_var[1]) / 3
K_cattell_strict <- 1
for (k in 2:length(diffs_var)) {
  if (abs(diffs_var[k]) < elbow_threshold) {
    K_cattell_strict <- k - 1   # PCs above the plateau
    break
  }
}
log_msg(sprintf("  Cattell (strict, no min-floor): K = %d", K_cattell_strict))

# Kaiser-Guttman analog: keep PCs with var_explained > 1/N where N is number
# of PCs computed. This is equivalent to "eigenvalue > average" in the
# scaled-PCA sense. With K=10 PCs, the threshold is 0.1 (10%).
N_pcs <- length(var_exp)
kaiser_thresh <- 1 / N_pcs
K_kaiser <- sum(var_exp > kaiser_thresh)
log_msg(sprintf("  Kaiser-Guttman (var > 1/%d = %.3f): K = %d",
                N_pcs, kaiser_thresh, K_kaiser))

# Horn's parallel analysis. Permute genotypes within each SNP across samples
# (breaks any structure), re-run pcadapt at K=10, record eigenvalues. Repeat
# n_perm times. Compute the 95th percentile of permuted variance per PC; the
# observed PC k is "real" iff its var_exp exceeds the 95th-percentile null.
# K_horn = number of consecutive PCs from PC1 that exceed null.
n_perm <- 10
log_msg(sprintf("  Horn's parallel analysis (n_perm=%d, please wait)...", n_perm))

geno_mat_full <- bed2matrix(paste0(PLINK_PREFIX, ".bed"))   # 81 x 803,821
# Subsample SNPs for parallel-analysis speed: 50k random SNPs is plenty for
# stable eigenvalue estimates.
n_subsample_snps <- 50000
set.seed(123)
snp_subsample <- sample(seq_len(ncol(geno_mat_full)), n_subsample_snps)
geno_obs <- geno_mat_full[, snp_subsample, drop = FALSE]

# Observed eigenvalues on the subsample
fit_obs_sub <- pcadapt(read.pcadapt(geno_obs, type = "lfmm"), K = 10)
sing_obs <- fit_obs_sub$singular.values
var_obs  <- sing_obs^2 / sum(sing_obs^2)

permuted_var <- matrix(NA_real_, nrow = n_perm, ncol = 10)
for (p in seq_len(n_perm)) {
  geno_perm <- geno_obs
  # Permute each column (SNP) independently across samples.
  for (j in seq_len(ncol(geno_perm))) {
    geno_perm[, j] <- geno_perm[sample(nrow(geno_perm)), j]
  }
  fit_perm <- pcadapt(read.pcadapt(geno_perm, type = "lfmm"), K = 10)
  sp <- fit_perm$singular.values
  permuted_var[p, ] <- sp^2 / sum(sp^2)
  cat(sprintf("    perm %d / %d done\n", p, n_perm))
}

null_p95 <- apply(permuted_var, 2, quantile, probs = 0.95, na.rm = TRUE)
K_horn <- 0
for (k in seq_along(var_obs)) {
  if (var_obs[k] > null_p95[k]) {
    K_horn <- k
  } else {
    break
  }
}
log_msg(sprintf("  Horn's parallel analysis (50k SNPs subsample): K = %d", K_horn))

scree_dx <- tibble(
  PC                  = seq_along(var_exp),
  variance_explained  = var_exp,
  cumulative          = cumsum(var_exp),
  drop_to_next        = c(diffs_var, NA_real_),
  kaiser_threshold    = kaiser_thresh,
  exceeds_kaiser      = var_exp > kaiser_thresh,
  observed_var_50k    = var_obs,
  null_p95_50k        = null_p95,
  exceeds_horn_null   = var_obs > null_p95
)
write_tsv(scree_dx, file.path(OUT_DIR, "scree_diagnostic.tsv"))

# Redrawn scree plot with all 3 K criteria.
scree_long <- tibble(
  PC = rep(seq_along(var_exp), 2),
  series = c(rep("observed (full SNP set)", length(var_exp)),
             rep("Horn null 95% (50k SNPs)", length(null_p95))),
  variance_pct = 100 * c(var_exp, null_p95)
)
p_scree <- ggplot(scree_long,
                  aes(x = PC, y = variance_pct, color = series, linetype = series)) +
  geom_line() + geom_point(size = 2) +
  geom_hline(yintercept = 100 * kaiser_thresh, linetype = "dotted", color = "grey40") +
  annotate("text", x = max(scree_long$PC) - 0.5, y = 100 * kaiser_thresh,
           label = "Kaiser threshold", vjust = -0.4, hjust = 1, size = 3) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0, max(scree_long$variance_pct) * 1.1)) +
  scale_color_manual(values = c("observed (full SNP set)" = "black",
                                "Horn null 95% (50k SNPs)" = "tomato")) +
  scale_linetype_manual(values = c("observed (full SNP set)" = "solid",
                                   "Horn null 95% (50k SNPs)" = "dashed")) +
  labs(
    title = "Scree plot with multiple K criteria (K=2 vs K=1 audit)",
    subtitle = sprintf(
      "Cattell strict (no min-floor) = %d  |  Kaiser (var > 1/N) = %d  |  Horn (95%% null) = %d",
      K_cattell_strict, K_kaiser, K_horn),
    x = "Principal Component",
    y = "Variance explained (%)",
    color = NULL, linetype = NULL,
    caption = "Whichever K is used downstream is an analytical choice; this plot reports the data-driven candidates."
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "scree_diagnostic.png"),
       p_scree, width = 160, height = 110, units = "mm", dpi = 300)
ggsave(file.path(OUT_DIR, "scree_diagnostic.pdf"),
       p_scree, width = 160, height = 110, units = "mm")

writeLines(c(
  "Scree-plot K-selection diagnostic",
  sprintf("Date: %s", format(Sys.time())),
  "",
  "Three independent criteria applied to the same K=10 pcadapt fit:",
  sprintf("  Cattell strict (no min-floor): K = %d", K_cattell_strict),
  sprintf("    - elbow defined as first PC where drop < (PC1->PC2 drop) / 3"),
  sprintf("    - the original auto-detection had a max(2, .) floor that was"),
  sprintf("      preventing K=1 from ever being reported. Floor removed."),
  sprintf("  Kaiser-Guttman (var > 1/N=%.3f): K = %d", kaiser_thresh, K_kaiser),
  sprintf("    - keeps any PC explaining more than its average share."),
  sprintf("  Horn's parallel analysis (n_perm=%d, 50k SNPs): K = %d",
          n_perm, K_horn),
  sprintf("    - PCs whose observed variance exceeds the 95%% null"),
  sprintf("      from genotype-shuffled data. Most conservative."),
  "",
  "Whichever K is committed to downstream is an analytical choice on top",
  "of these diagnostics. This file does NOT pick a winner; it reports the",
  "facts so that the choice can be made and disclosed transparently."
), file.path(OUT_DIR, "scree_K_summary.txt"))
log_msg(sprintf("  wrote %s and scree_diagnostic.png/pdf",
                file.path(OUT_DIR, "scree_K_summary.txt")))

################################################################################
# TASK 5: PC2 sanity check -- |PC2 loading| vs among-invasive FST
################################################################################

log_msg("\n--- TASK 5: PC2 loading vs among-invasive FST ---")

# Refit at K=2 to get PC2 loadings. (We have fit_K10$loadings but those are
# for the K=10 fit and may differ; refit at K=2 to match the canonical run.)
fit_K2 <- pcadapt(geno, K = 2)
abs_pc2_load <- abs(fit_K2$loadings[, 2])
load_table <- bim %>%
  select(chrom, pos, snp_id) %>%
  mutate(abs_pc2_load = abs_pc2_load)

if (file.exists(PER_SNP_FST)) {
  axis_fst <- fread(PER_SNP_FST)
  axis_fst <- axis_fst %>%
    transmute(snp_id = paste0(CHROM, ":", POS), invasive_max_fst)
  load_fst <- load_table %>%
    inner_join(axis_fst, by = "snp_id") %>%
    filter(!is.na(invasive_max_fst), !is.na(abs_pc2_load))

  log_msg(sprintf("  joined SNPs (PC2 loading + invasive FST): %d",
                  nrow(load_fst)))

  spearman_rho <- cor(load_fst$abs_pc2_load, load_fst$invasive_max_fst,
                       method = "spearman", use = "complete.obs")
  pearson_rho  <- cor(load_fst$abs_pc2_load, load_fst$invasive_max_fst,
                       method = "pearson",  use = "complete.obs")
  log_msg(sprintf("  Spearman rho (|PC2 loading|, invasive_max_fst): %.4f",
                  spearman_rho))
  log_msg(sprintf("  Pearson  rho (|PC2 loading|, invasive_max_fst): %.4f",
                  pearson_rho))

  # Binned enrichment: deciles of |PC2 loading|, mean invasive FST per decile.
  load_fst <- load_fst %>%
    mutate(load_decile = ntile(abs_pc2_load, 10))
  bin_summary <- load_fst %>%
    group_by(load_decile) %>%
    summarise(
      n         = n(),
      load_mean = mean(abs_pc2_load),
      fst_mean  = mean(invasive_max_fst),
      fst_med   = median(invasive_max_fst),
      fst_q90   = quantile(invasive_max_fst, 0.9),
      .groups   = "drop"
    )
  log_msg("  Mean among-invasive max FST by |PC2 loading| decile:")
  for (i in seq_len(nrow(bin_summary))) {
    log_msg(sprintf("    decile %d  n=%6d  load_mean=%.4f  fst_mean=%.4f  fst_q90=%.4f",
                    bin_summary$load_decile[i], bin_summary$n[i],
                    bin_summary$load_mean[i], bin_summary$fst_mean[i],
                    bin_summary$fst_q90[i]))
  }

  write_tsv(load_fst, file.path(OUT_DIR, "pc2_loading_vs_fst.tsv"))
  writeLines(c(
    sprintf("PC2 sanity check"),
    sprintf("Date: %s", format(Sys.time())),
    sprintf("N SNPs: %d", nrow(load_fst)),
    sprintf("Spearman rho: %.4f", spearman_rho),
    sprintf("Pearson  rho: %.4f", pearson_rho),
    "",
    "Mean among-invasive max FST by |PC2 loading| decile:",
    paste(capture.output(print(bin_summary, n = 10)), collapse = "\n")
  ), file.path(OUT_DIR, "pc2_loading_vs_fst_summary.txt"))

  # Hexbin plot
  p_pc2 <- ggplot(load_fst, aes(abs_pc2_load, invasive_max_fst)) +
    geom_hex(bins = 80) +
    scale_fill_viridis_c(trans = "log10") +
    labs(
      title = "PC2 sanity check: |PC2 loading| vs per-SNP among-invasive max FST",
      subtitle = sprintf("Spearman rho = %.3f  |  N = %d  |  K=2 fit",
                         spearman_rho, nrow(load_fst)),
      x = "|PC2 loading|", y = "per-SNP among-invasive max FST",
      fill = "SNP count"
    ) +
    theme_bw(base_size = 9)
  ggsave(file.path(OUT_DIR, "pc2_loading_vs_fst.png"),
         p_pc2, width = 130, height = 100, units = "mm", dpi = 300)
} else {
  log_msg(sprintf("  per-SNP axis FST file not found at %s; skipping task 5",
                  PER_SNP_FST))
}

################################################################################
# TASK 6: Optional K=3 run
################################################################################

log_msg("\n--- TASK 6: Optional K=3 run ---")

fit_K3 <- pcadapt(geno, K = 3)
chisq_obs_K3 <- qchisq(1 - fit_K3$pvalues, df = 1)
lambda_K3 <- median(chisq_obs_K3, na.rm = TRUE) / qchisq(0.5, df = 1)
qv_K3 <- qvalue(fit_K3$pvalues)
n_q05_K3 <- sum(!is.na(qv_K3$qvalues) & qv_K3$qvalues < 0.05)
n_q01_K3 <- sum(!is.na(qv_K3$qvalues) & qv_K3$qvalues < 0.01)
log_msg(sprintf("  K=3 lambda=%.4f  q<0.05 outliers=%d  q<0.01 outliers=%d",
                lambda_K3, n_q05_K3, n_q01_K3))

# PC3 score plot to check structure.
scores_K3 <- as.data.frame(fit_K3$scores)
colnames(scores_K3) <- paste0("PC", 1:3)
scores_K3$iid <- fam_pop$iid
scores_K3$pop <- fam_pop$pop

p1 <- ggplot(scores_K3, aes(PC1, PC2, color = pop)) +
  geom_point(size = 2, alpha = 0.85) + theme_bw() + ggtitle("K=3 fit  -  PC1 vs PC2")
p2 <- ggplot(scores_K3, aes(PC1, PC3, color = pop)) +
  geom_point(size = 2, alpha = 0.85) + theme_bw() + ggtitle("K=3 fit  -  PC1 vs PC3")
p3 <- ggplot(scores_K3, aes(PC2, PC3, color = pop)) +
  geom_point(size = 2, alpha = 0.85) + theme_bw() + ggtitle("K=3 fit  -  PC2 vs PC3")

# Combine via patchwork if available, otherwise save individually.
if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- patchwork::wrap_plots(p1, p2, p3, ncol = 3) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(file.path(OUT_DIR, "pca_K3_scores.png"),
         combined, width = 220, height = 80, units = "mm", dpi = 300)
} else {
  ggsave(file.path(OUT_DIR, "pca_K3_scores_PC1xPC2.png"),
         p1, width = 90, height = 80, units = "mm", dpi = 300)
  ggsave(file.path(OUT_DIR, "pca_K3_scores_PC1xPC3.png"),
         p2, width = 90, height = 80, units = "mm", dpi = 300)
  ggsave(file.path(OUT_DIR, "pca_K3_scores_PC2xPC3.png"),
         p3, width = 90, height = 80, units = "mm", dpi = 300)
}

# PC3 per-pop summary
pc3_per_pop <- scores_K3 %>%
  group_by(pop) %>%
  summarise(mean = mean(PC3), sd = sd(PC3),
            min = min(PC3), max = max(PC3), .groups = "drop")

writeLines(c(
  sprintf("K=3 diagnostic"),
  sprintf("Date: %s", format(Sys.time())),
  sprintf("Lambda: %.4f", lambda_K3),
  sprintf("Outliers q<0.05: %d", n_q05_K3),
  sprintf("Outliers q<0.01: %d", n_q01_K3),
  "",
  "PC3 per-pop summary:",
  paste(capture.output(print(pc3_per_pop, n = 10)), collapse = "\n")
), file.path(OUT_DIR, "pcadapt_K3_summary.txt"))

################################################################################
# Synthesis
################################################################################

log_msg("\n=== Diagnostic battery complete ===")
log_msg(sprintf("Outputs in: %s", OUT_DIR))
writeLines(log_lines, file.path(OUT_DIR, "diagnostic_run.log"))
