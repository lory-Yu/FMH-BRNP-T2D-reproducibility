#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(vegan)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args) >= 1) normalizePath(args[[1]]) else normalizePath("wei_2025_calibration")
results_dir <- file.path(project, "results")
obj <- readRDS(file.path(results_dir, "wei2025_dada2_objects.rds"))
metadata <- obj$metadata
rownames(metadata) <- metadata$donor_id
counts <- obj$seqtab_clean[metadata$donor_id, , drop = FALSE]

rarefaction_depth <- 30000L
if (min(rowSums(counts)) < rarefaction_depth) stop("At least one sample is below rarefaction depth")
set.seed(20250903)
rare <- rrarefy(counts, sample = rarefaction_depth)

alpha <- data.frame(
  donor_id = rownames(rare),
  Observed_ASV = specnumber(rare),
  Chao1 = unname(estimateR(rare)["S.chao1", ]),
  Shannon = diversity(rare, index = "shannon"),
  Simpson = diversity(rare, index = "simpson"),
  rb1_group = metadata[rownames(rare), "rb1_group"],
  age = metadata[rownames(rare), "age"],
  gender = metadata[rownames(rare), "gender"],
  stringsAsFactors = FALSE
)

hl <- alpha$rb1_group %in% c("High", "Low")
alpha_tests <- do.call(rbind, lapply(c("Observed_ASV", "Chao1", "Shannon", "Simpson"), function(metric) {
  high <- alpha[hl & alpha$rb1_group == "High", metric]
  low <- alpha[hl & alpha$rb1_group == "Low", metric]
  wilcox <- wilcox.test(high, low, exact = FALSE)
  frame <- alpha[hl, c(metric, "rb1_group", "age", "gender")]
  names(frame)[1] <- "response"
  fit <- lm(response ~ relevel(factor(rb1_group), ref = "Low") + age + factor(gender), data = frame)
  coef_name <- grep("High", names(coef(fit)), value = TRUE)[1]
  data.frame(
    metric = metric, high_median = median(high), low_median = median(low),
    wilcoxon_p = wilcox$p.value,
    adjusted_difference_high_vs_low = unname(coef(fit)[coef_name]),
    adjusted_p = coef(summary(fit))[coef_name, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}))
alpha_tests$adjusted_q <- p.adjust(alpha_tests$adjusted_p, method = "BH")
write.table(alpha_tests, file.path(results_dir, "rarefied_30000_alpha_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

meta_hl <- metadata[rownames(rare), ]
meta_hl <- meta_hl[meta_hl$rb1_group %in% c("High", "Low"), ]
meta_hl$rb1_group <- relevel(factor(meta_hl$rb1_group), ref = "Low")
meta_hl$gender <- factor(meta_hl$gender)
bray <- vegdist(rare[meta_hl$donor_id, ] / rowSums(rare[meta_hl$donor_id, ]), method = "bray")
set.seed(20250903)
permanova <- adonis2(bray ~ rb1_group + age + gender, data = meta_hl, permutations = 999, by = "margin")
set.seed(20250903)
anosim_fit <- anosim(bray, grouping = meta_hl$rb1_group, permutations = 999)
beta <- data.frame(
  rarefaction_depth = rarefaction_depth,
  permanova_R2 = permanova["rb1_group", "R2"],
  permanova_p = permanova["rb1_group", "Pr(>F)"],
  anosim_R = unname(anosim_fit$statistic),
  anosim_p = unname(anosim_fit$signif)
)
write.table(beta, file.path(results_dir, "rarefied_30000_beta_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write_json(
  list(rarefaction_depth = rarefaction_depth, alpha = alpha_tests, beta = beta),
  file.path(results_dir, "rarefied_30000_sensitivity_summary.json"),
  pretty = TRUE, auto_unbox = TRUE, dataframe = "rows"
)

