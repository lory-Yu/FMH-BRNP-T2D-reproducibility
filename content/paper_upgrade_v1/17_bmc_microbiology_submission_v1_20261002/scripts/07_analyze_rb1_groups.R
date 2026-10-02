#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath("wei_2025_calibration")
results_dir <- file.path(project, "results")
figures_dir <- file.path(project, "figures")
provenance_dir <- file.path(project, "provenance")
logs_dir <- file.path(project, "logs")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
log_connection <- file(file.path(logs_dir, "07_analyze_rb1_groups.log"), open = "wt")
sink(log_connection, type = "output")
sink(log_connection, type = "message")

obj <- readRDS(file.path(results_dir, "wei2025_dada2_objects.rds"))
metadata <- obj$metadata
rownames(metadata) <- metadata$donor_id
counts <- obj$seqtab_clean
counts <- counts[metadata$donor_id, , drop = FALSE]
taxa_clean <- obj$taxonomy_clean
genus_labels <- taxa_clean[, "Genus"]
missing_genus <- is.na(genus_labels) | genus_labels == ""
family_labels <- taxa_clean[, "Family"]
usable_family <- missing_genus & !is.na(family_labels) & family_labels != ""
genus_labels[usable_family] <- paste0("norank_f_", family_labels[usable_family])
genus_labels[missing_genus & !usable_family] <- "Unclassified_genus"
genus_counts <- t(rowsum(t(counts), group = genus_labels, reorder = FALSE))
write.table(
  data.frame(donor_id = rownames(genus_counts), genus_counts, check.names = FALSE),
  file.path(results_dir, "genus_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(donor_id = rownames(genus_counts), genus_counts / rowSums(genus_counts), check.names = FALSE),
  file.path(results_dir, "genus_relative_abundance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

alpha <- data.frame(
  donor_id = rownames(counts),
  Observed_ASV = specnumber(counts),
  Chao1 = unname(estimateR(counts)["S.chao1", ]),
  Shannon = diversity(counts, index = "shannon"),
  Simpson = diversity(counts, index = "simpson"),
  library_size = rowSums(counts),
  stringsAsFactors = FALSE
)
alpha <- merge(alpha, metadata[, c("donor_id", "age", "gender", "rb1_group")], by = "donor_id")
alpha <- alpha[order(as.integer(sub("D", "", alpha$donor_id))), ]
write.table(alpha, file.path(results_dir, "alpha_diversity_by_donor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rank_biserial <- function(x, y) {
  comparisons <- outer(x, y, "-")
  (sum(comparisons > 0) - sum(comparisons < 0)) / length(comparisons)
}

hl <- alpha$rb1_group %in% c("High", "Low")
alpha_tests <- do.call(rbind, lapply(c("Observed_ASV", "Chao1", "Shannon", "Simpson"), function(metric) {
  high <- alpha[hl & alpha$rb1_group == "High", metric]
  low <- alpha[hl & alpha$rb1_group == "Low", metric]
  test <- wilcox.test(high, low, exact = FALSE, conf.int = TRUE)
  fit <- lm(alpha[[metric]][hl] ~ relevel(factor(alpha$rb1_group[hl]), ref = "Low") + alpha$age[hl] + factor(alpha$gender[hl]))
  coef_name <- grep("High", names(coef(fit)), value = TRUE)[1]
  data.frame(
    metric = metric,
    n_high = length(high), n_low = length(low),
    high_median = median(high), low_median = median(low),
    median_difference = median(high) - median(low),
    rank_biserial_high_vs_low = rank_biserial(high, low),
    wilcoxon_p = test$p.value,
    adjusted_mean_difference = unname(coef(fit)[coef_name]),
    adjusted_p = coef(summary(fit))[coef_name, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}))
alpha_tests$adjusted_q <- p.adjust(alpha_tests$adjusted_p, method = "BH")
write.table(alpha_tests, file.path(results_dir, "alpha_diversity_high_vs_low.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

rel <- counts / rowSums(counts)
bray <- vegdist(rel, method = "bray")
meta_hl <- metadata[rownames(rel), ]
meta_hl <- meta_hl[meta_hl$rb1_group %in% c("High", "Low"), ]
rel_hl <- rel[meta_hl$donor_id, , drop = FALSE]
bray_hl <- vegdist(rel_hl, method = "bray")
meta_hl$rb1_group <- relevel(factor(meta_hl$rb1_group), ref = "Low")
meta_hl$gender <- factor(meta_hl$gender)
set.seed(20250903)
permanova <- adonis2(bray_hl ~ rb1_group + age + gender, data = meta_hl, permutations = 999, by = "margin")
permanova_out <- data.frame(term = rownames(permanova), permanova, check.names = FALSE)
write.table(permanova_out, file.path(results_dir, "beta_bray_permanova_high_vs_low.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

dispersion <- betadisper(bray_hl, meta_hl$rb1_group)
set.seed(20250903)
dispersion_test <- permutest(dispersion, permutations = 999)
dispersion_out <- data.frame(
  statistic = unname(dispersion_test$tab[1, "F"]),
  p_value = unname(dispersion_test$tab[1, "Pr(>F)"]),
  permutations = 999
)
write.table(dispersion_out, file.path(results_dir, "beta_bray_dispersion_high_vs_low.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Method-matched sensitivity analysis: the original article used ANOSIM.
set.seed(20250903)
anosim_fit <- anosim(bray_hl, grouping = meta_hl$rb1_group, permutations = 999)
anosim_out <- data.frame(
  statistic_R = unname(anosim_fit$statistic),
  p_value = unname(anosim_fit$signif),
  permutations = 999
)
write.table(anosim_out, file.path(results_dir, "beta_bray_anosim_high_vs_low.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pcoa <- cmdscale(bray_hl, eig = TRUE, k = 2, add = TRUE)
pcoa_df <- data.frame(
  donor_id = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  rb1_group = meta_hl[rownames(pcoa$points), "rb1_group"],
  stringsAsFactors = FALSE
)
write.table(pcoa_df, file.path(results_dir, "beta_bray_pcoa_high_vs_low.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

theme_set(theme_classic(base_size = 12))
alpha_long <- reshape(
  alpha[hl, c("donor_id", "rb1_group", "Observed_ASV", "Chao1", "Shannon", "Simpson")],
  varying = c("Observed_ASV", "Chao1", "Shannon", "Simpson"),
  v.names = "value", timevar = "metric",
  times = c("Observed ASVs", "Chao1", "Shannon", "Simpson"), direction = "long"
)
p_alpha <- ggplot(alpha_long, aes(x = rb1_group, y = value, fill = rb1_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.08, size = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(Low = "#4E79A7", High = "#E15759")) +
  labs(x = "Rb1 conversion group", y = "Diversity value", title = "Wei 2025: alpha diversity in Rb1 high vs low converters") +
  guides(fill = "none")
ggsave(file.path(figures_dir, "Figure_Wei2025_alpha_diversity.png"), p_alpha, width = 8.2, height = 6.6, dpi = 300)
ggsave(file.path(figures_dir, "Figure_Wei2025_alpha_diversity.pdf"), p_alpha, width = 8.2, height = 6.6)

p_beta <- ggplot(pcoa_df, aes(PCoA1, PCoA2, color = rb1_group)) +
  geom_point(size = 3) +
  stat_ellipse(type = "norm", linewidth = 0.7, level = 0.8) +
  scale_color_manual(values = c(Low = "#4E79A7", High = "#E15759")) +
  labs(x = "PCoA axis 1", y = "PCoA axis 2", color = "Rb1 group", title = "Wei 2025: Bray-Curtis PCoA")
ggsave(file.path(figures_dir, "Figure_Wei2025_beta_Bray_PCoA.png"), p_beta, width = 6.4, height = 5.0, dpi = 300)
ggsave(file.path(figures_dir, "Figure_Wei2025_beta_Bray_PCoA.pdf"), p_beta, width = 6.4, height = 5.0)

prevalence <- colMeans(genus_counts[meta_hl$donor_id, , drop = FALSE] > 0)
mean_rel <- colMeans(genus_counts[meta_hl$donor_id, , drop = FALSE] / rowSums(genus_counts[meta_hl$donor_id, , drop = FALSE]))
keep_genera <- prevalence >= 0.20 & mean_rel >= 0.0001 & colnames(genus_counts) != "Unclassified_genus"
gc <- genus_counts[, keep_genera, drop = FALSE]
clr <- log(gc + 0.5) - rowMeans(log(gc + 0.5))

fit_genus <- function(genus) {
  dat <- data.frame(
    clr = clr[meta_hl$donor_id, genus],
    group = meta_hl$rb1_group,
    age = meta_hl$age,
    gender = meta_hl$gender
  )
  fit <- lm(clr ~ group + age + gender, data = dat)
  coef_name <- "groupHigh"
  rel_values <- genus_counts[meta_hl$donor_id, genus] / rowSums(genus_counts[meta_hl$donor_id, , drop = FALSE])
  wilcox <- wilcox.test(rel_values[meta_hl$rb1_group == "High"], rel_values[meta_hl$rb1_group == "Low"], exact = FALSE)
  data.frame(
    genus = genus,
    prevalence = prevalence[genus],
    mean_relative_abundance = mean_rel[genus],
    clr_difference_high_vs_low = unname(coef(fit)[coef_name]),
    clr_se = coef(summary(fit))[coef_name, "Std. Error"],
    adjusted_p = coef(summary(fit))[coef_name, "Pr(>|t|)"],
    wilcoxon_p = wilcox$p.value,
    high_median_relative_abundance = median(rel_values[meta_hl$rb1_group == "High"]),
    low_median_relative_abundance = median(rel_values[meta_hl$rb1_group == "Low"]),
    stringsAsFactors = FALSE
  )
}

genus_results <- do.call(rbind, lapply(colnames(gc), fit_genus))
genus_results$adjusted_q <- p.adjust(genus_results$adjusted_p, method = "BH")
genus_results$wilcoxon_q <- p.adjust(genus_results$wilcoxon_p, method = "BH")
genus_results <- genus_results[order(genus_results$adjusted_q, genus_results$adjusted_p), ]
write.table(genus_results, file.path(results_dir, "genus_high_vs_low_clr.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

normalize_taxon <- function(x) gsub("[^a-z0-9]+", "", tolower(x))
# Genera reported as higher in the original Wei 2025 Rb1 high-converter group
# (article Results/Fig. S2D); retained here as a prespecified replication set.
candidates <- c(
  "Eubacterium_hallii_group", "Anaerostipes",
  "norank_f_Eubacterium_coprostanoligenes_group", "Coprococcus", "Barnesiella",
  "norank_f_Oscillospiraceae", "UBA1819", "Family_XIII_AD3011_group"
)
candidate_patterns <- c(
  "^(eubacteriumhalliigroup|eubacteriumhalligroup)$",
  "^anaerostipes$",
  "^(norankf)?eubacteriumcoprostanoligenesgroup$",
  "^coprococcus$", "^barnesiella$",
  "^(norankf)?oscillospiraceae$",
  "^uba1819$", "^familyxiiiad3011group$"
)
candidate_rows <- lapply(candidates, function(candidate) {
  genus_norm <- normalize_taxon(genus_results$genus)
  pattern <- candidate_patterns[match(candidate, candidates)]
  idx <- which(grepl(pattern, genus_norm, perl = TRUE))
  if (!length(idx)) {
    missing_result <- as.data.frame(as.list(setNames(rep(NA, ncol(genus_results) - 1), setdiff(names(genus_results), "genus"))))
    cbind(
      data.frame(candidate = candidate, observed_genus = NA_character_, detected = FALSE),
      missing_result
    )
  } else {
    cbind(data.frame(candidate = candidate, observed_genus = genus_results$genus[idx], detected = TRUE), genus_results[idx, setdiff(names(genus_results), "genus"), drop = FALSE])
  }
})
candidate_results <- do.call(rbind, candidate_rows)
candidate_results$candidate_set_adjusted_q <- NA_real_
candidate_results$candidate_set_wilcoxon_q <- NA_real_
candidate_detected <- candidate_results$detected & !is.na(candidate_results$adjusted_p)
candidate_results$candidate_set_adjusted_q[candidate_detected] <- p.adjust(
  candidate_results$adjusted_p[candidate_detected], method = "BH"
)
candidate_results$candidate_set_wilcoxon_q[candidate_detected] <- p.adjust(
  candidate_results$wilcoxon_p[candidate_detected], method = "BH"
)
write.table(candidate_results, file.path(results_dir, "prespecified_candidate_genera_high_vs_low.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

candidate_plot_data <- candidate_results[candidate_results$detected, , drop = FALSE]
candidate_plot_data$display_name <- factor(
  candidate_plot_data$candidate,
  levels = rev(candidate_plot_data$candidate)
)
p_candidates <- ggplot(
  candidate_plot_data,
  aes(x = clr_difference_high_vs_low, y = display_name)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7A7A7A") +
  geom_errorbarh(
    aes(xmin = clr_difference_high_vs_low - 1.96 * clr_se,
        xmax = clr_difference_high_vs_low + 1.96 * clr_se),
    height = 0.16, color = "#4E79A7"
  ) +
  geom_point(size = 2.8, color = "#E15759") +
  labs(
    x = "Adjusted CLR difference (High - Low; 95% CI)", y = NULL,
    title = "Wei 2025: prespecified Rb1-associated taxa"
  )
ggsave(file.path(figures_dir, "Figure_Wei2025_candidate_taxa.png"), p_candidates, width = 8.2, height = 5.2, dpi = 300)
ggsave(file.path(figures_dir, "Figure_Wei2025_candidate_taxa.pdf"), p_candidates, width = 8.2, height = 5.2)

metadata$rb1_ordinal <- c(Low = 0, Middle = 1, High = 2)[metadata$rb1_group]
trend_results <- do.call(rbind, lapply(colnames(gc), function(genus) {
  dat <- data.frame(clr = clr[metadata$donor_id, genus], ordinal = metadata$rb1_ordinal, age = metadata$age, gender = factor(metadata$gender))
  fit <- lm(clr ~ ordinal + age + gender, data = dat)
  data.frame(
    genus = genus,
    ordinal_beta = coef(fit)["ordinal"],
    ordinal_se = coef(summary(fit))["ordinal", "Std. Error"],
    ordinal_p = coef(summary(fit))["ordinal", "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}))
trend_results$ordinal_q <- p.adjust(trend_results$ordinal_p, method = "BH")
trend_results <- trend_results[order(trend_results$ordinal_q, trend_results$ordinal_p), ]
write.table(trend_results, file.path(results_dir, "genus_low_middle_high_trend.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary <- list(
  primary_n = nrow(meta_hl),
  high_n = sum(meta_hl$rb1_group == "High"),
  low_n = sum(meta_hl$rb1_group == "Low"),
  total_n = nrow(metadata),
  tested_genera = nrow(genus_results),
  genus_fdr_significant = sum(genus_results$adjusted_q < 0.05),
  prespecified_taxa_reported = nrow(candidate_results),
  prespecified_taxa_detected = sum(candidate_results$detected),
  prespecified_taxa_positive_direction = sum(candidate_results$detected & candidate_results$clr_difference_high_vs_low > 0, na.rm = TRUE),
  prespecified_taxa_candidate_set_fdr_significant_adjusted = sum(candidate_results$candidate_set_adjusted_q < 0.05, na.rm = TRUE),
  prespecified_taxa_candidate_set_fdr_significant_wilcoxon = sum(candidate_results$candidate_set_wilcoxon_q < 0.05, na.rm = TRUE),
  alpha = alpha_tests,
  permanova_group_R2 = permanova["rb1_group", "R2"],
  permanova_group_p = permanova["rb1_group", "Pr(>F)"],
  anosim_R = anosim_out$statistic_R,
  anosim_p = anosim_out$p_value,
  dispersion_p = dispersion_out$p_value
)
write_json(summary, file.path(results_dir, "rb1_group_analysis_summary.json"), pretty = TRUE, auto_unbox = TRUE, dataframe = "rows")
writeLines(capture.output(sessionInfo()), file.path(provenance_dir, "statistics_sessionInfo.txt"))
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, dataframe = "rows"), "\n")
sink(type = "message")
sink(type = "output")
close(log_connection)
