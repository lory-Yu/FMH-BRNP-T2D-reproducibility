#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_path <- normalizePath(script_arg, mustWork = TRUE)
package <- dirname(dirname(script_path))
src <- file.path(package, "inputs", "three_cohort")
out <- file.path(package, "results", "three_cohort")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pathways <- c(
  "PWY-5676" = "acetyl-CoA fermentation to butanoate II",
  "CENTFERM-PWY" = "pyruvate fermentation to butanoate"
)

pathway_vector <- function(obj, pid) {
  ids <- sub(":.*$", "", rownames(obj$pathway_matrix))
  keep <- ids == pid & !grepl("\\|", rownames(obj$pathway_matrix))
  if (sum(keep) != 1L) stop(sprintf("Expected exactly one unstratified row for %s; found %d", pid, sum(keep)))
  y <- as.numeric(obj$pathway_matrix[keep, ])
  names(y) <- colnames(obj$pathway_matrix)
  if (!identical(names(y), obj$metadata$sample_id)) stop("Pathway matrix and metadata order differ")
  y
}

fit_one <- function(cohort, obj, pid, variant, subset_keep, formula_rhs) {
  y <- pathway_vector(obj, pid)
  dat <- obj$metadata[subset_keep, , drop = FALSE]
  y <- y[subset_keep]
  dat$detected <- as.integer(y > 0)
  dat$log2_positive <- ifelse(y > 0, log2(y), NA_real_)
  out_rows <- list()

  # Detection is retained for audit but never pooled when separation/no variation occurs.
  f_det <- as.formula(paste("detected ~", formula_rhs))
  complete_det <- complete.cases(model.frame(f_det, dat, na.action = na.pass))
  d_det <- dat[complete_det, , drop = FALSE]
  detection_status <- if (length(unique(d_det$detected)) < 2L) "not_estimable_no_detection_variation" else "estimated"
  if (detection_status == "estimated") {
    warns <- character()
    g <- withCallingHandlers(glm(f_det, data = d_det, family = binomial()), warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    z <- coef(summary(g))
    if (!("groupT2D" %in% rownames(z)) || !is.finite(z["groupT2D", "Std. Error"]) || z["groupT2D", "Std. Error"] > 100) {
      detection_status <- "not_pooled_separation_or_unstable"
      est <- se <- p <- NA_real_
    } else {
      est <- z["groupT2D", "Estimate"]
      se <- z["groupT2D", "Std. Error"]
      p <- z["groupT2D", "Pr(>|z|)"]
    }
  } else {
    warns <- character(); est <- se <- p <- NA_real_
  }
  out_rows[[1]] <- data.frame(
    cohort = cohort, pathway_id = pid, pathway_name = unname(pathways[pid]), variant = variant,
    metric = "detection", formula = paste("detected ~", formula_rhs), n = nrow(d_det),
    n_control = sum(d_det$group == "Control"), n_t2d = sum(d_det$group == "T2D"),
    estimate_log_odds = est, se = se, ci95_low = est - 1.96 * se, ci95_high = est + 1.96 * se,
    p_value = p, transformed_measure = "odds_ratio", transformed_estimate = exp(est),
    transformed_ci95_low = exp(est - 1.96 * se), transformed_ci95_high = exp(est + 1.96 * se),
    status = detection_status, warnings = paste(unique(warns), collapse = " | ")
  )

  f_pos <- as.formula(paste("log2_positive ~", formula_rhs))
  complete_pos <- complete.cases(model.frame(f_pos, dat, na.action = na.pass))
  d_pos <- dat[complete_pos, , drop = FALSE]
  lm_fit <- lm(f_pos, data = d_pos)
  z <- coef(summary(lm_fit))
  est <- z["groupT2D", "Estimate"]
  se <- z["groupT2D", "Std. Error"]
  p <- z["groupT2D", "Pr(>|t|)"]
  out_rows[[2]] <- data.frame(
    cohort = cohort, pathway_id = pid, pathway_name = unname(pathways[pid]), variant = variant,
    metric = "positive_abundance", formula = paste("log2_positive ~", formula_rhs), n = nrow(d_pos),
    n_control = sum(d_pos$group == "Control"), n_t2d = sum(d_pos$group == "T2D"),
    estimate_log_odds = est, se = se, ci95_low = est - 1.96 * se, ci95_high = est + 1.96 * se,
    p_value = p, transformed_measure = "positive_abundance_ratio",
    transformed_estimate = 2^est, transformed_ci95_low = 2^(est - 1.96 * se), transformed_ci95_high = 2^(est + 1.96 * se),
    status = "estimated", warnings = ""
  )
  do.call(rbind, out_rows)
}

q <- readRDS(file.path(src, "Qin_analysis_ready_raw.rds"))
k <- readRDS(file.path(src, "Karlsson_analysis_ready_raw.rds"))
m <- readRDS(file.path(src, "MetaCardis_analysis_ready_raw.rds"))

fits <- list()
for (pid in names(pathways)) {
  fits[[length(fits) + 1L]] <- fit_one("Qin", q, pid, "primary_demographic_adjusted", rep(TRUE, nrow(q$metadata)), "group + age_category + gender + BMI")
  fits[[length(fits) + 1L]] <- fit_one("Qin", q, pid, "known_no_metformin", q$metadata$treatment == "no" & !is.na(q$metadata$treatment), "group + age_category + gender + BMI")
  fits[[length(fits) + 1L]] <- fit_one("Karlsson", k, pid, "primary_demographic_adjusted", rep(TRUE, nrow(k$metadata)), "group + age + BMI")
  fits[[length(fits) + 1L]] <- fit_one("Karlsson", k, pid, "known_no_metformin", !grepl("metformin", k$metadata$treatment, ignore.case = TRUE) & !is.na(k$metadata$treatment), "group + age + BMI")
  fits[[length(fits) + 1L]] <- fit_one("MetaCardis", m, pid, "primary_demographic_adjusted", rep(TRUE, nrow(m$metadata)), "group + age_category + gender + BMI + country")
  fits[[length(fits) + 1L]] <- fit_one("MetaCardis", m, pid, "known_no_metformin", !grepl("metformin", m$metadata$treatment, ignore.case = TRUE) & !is.na(m$metadata$treatment), "group + age_category + gender + BMI + country")
}
cohort_effects <- do.call(rbind, fits)
write.table(cohort_effects, file.path(out, "cohort_two_part_effects.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

# REML intercept-only random-effects model, with ordinary and modified KH side by side.
reml_pool <- function(y, se) {
  v <- se^2; k <- length(y)
  nll <- function(tau2) {
    w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
    0.5 * (sum(log(v + tau2)) + log(sum(w)) + sum(w * (y - mu)^2))
  }
  upper <- max(2, var(y) * 10 + max(v))
  tau2 <- optimize(nll, interval = c(0, upper), tol = 1e-12)$minimum
  if (tau2 < 1e-8) tau2 <- 0
  w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
  q_hk <- sum(w * (y - mu)^2) / (k - 1)
  se_hk <- sqrt(q_hk / sum(w))
  se_mkh <- sqrt(max(1, q_hk) / sum(w))
  wf <- 1 / v; mu_f <- sum(wf * y) / sum(wf)
  Q <- sum(wf * (y - mu_f)^2)
  # metafor::rma.uni reports I2 from tau2 and the typical sampling variance.
  vt <- (k - 1) * sum(wf) / (sum(wf)^2 - sum(wf^2))
  I2 <- ifelse(tau2 > 0, 100 * tau2 / (tau2 + vt), 0)
  make <- function(kind, sx) data.frame(
    method = "REML", test = kind, k = k, estimate_log2_ratio = mu, se = sx,
    ci95_low = mu - qt(.975, k - 1) * sx, ci95_high = mu + qt(.975, k - 1) * sx,
    p_value = 2 * pt(-abs(mu / sx), df = k - 1), tau2 = tau2, I2 = I2, Q = Q,
    Q_p = pchisq(Q, df = k - 1, lower.tail = FALSE), q_scale = q_hk
  )
  rbind(make("ordinary_KH", se_hk), make("modified_KH", se_mkh))
}

meta_rows <- list()
for (variant in c("primary_demographic_adjusted", "known_no_metformin")) {
  for (pid in names(pathways)) {
    d <- cohort_effects[cohort_effects$variant == variant & cohort_effects$pathway_id == pid & cohort_effects$metric == "positive_abundance" & cohort_effects$status == "estimated", ]
    z <- reml_pool(d$estimate_log_odds, d$se)
    z$variant <- variant; z$pathway_id <- pid
    z$abundance_ratio <- 2^z$estimate_log2_ratio
    z$ratio_ci95_low <- 2^z$ci95_low; z$ratio_ci95_high <- 2^z$ci95_high
    meta_rows[[length(meta_rows) + 1L]] <- z
  }
}
meta <- do.call(rbind, meta_rows)
meta <- meta[, c("variant", "pathway_id", "method", "test", "k", "estimate_log2_ratio", "se", "ci95_low", "ci95_high", "p_value", "abundance_ratio", "ratio_ci95_low", "ratio_ci95_high", "tau2", "I2", "Q", "Q_p", "q_scale")]
write.table(meta, file.path(out, "three_cohort_REML_KH_rebuilt.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

loo_rows <- list()
for (variant in c("primary_demographic_adjusted", "known_no_metformin")) {
  for (pid in names(pathways)) {
    d <- cohort_effects[cohort_effects$variant == variant & cohort_effects$pathway_id == pid & cohort_effects$metric == "positive_abundance" & cohort_effects$status == "estimated", ]
    for (omitted in d$cohort) {
      z <- reml_pool(d$estimate_log_odds[d$cohort != omitted], d$se[d$cohort != omitted])
      z <- z[z$test == "modified_KH", ]
      z$variant <- variant; z$pathway_id <- pid; z$omitted_cohort <- omitted
      z$abundance_ratio <- 2^z$estimate_log2_ratio
      z$ratio_ci95_low <- 2^z$ci95_low; z$ratio_ci95_high <- 2^z$ci95_high
      loo_rows[[length(loo_rows) + 1L]] <- z
    }
  }
}
loo <- do.call(rbind, loo_rows)
loo <- loo[, c("variant", "pathway_id", "omitted_cohort", "method", "test", "k", "estimate_log2_ratio", "se", "ci95_low", "ci95_high", "p_value", "abundance_ratio", "ratio_ci95_low", "ratio_ci95_high", "tau2", "I2", "Q", "Q_p", "q_scale")]
write.table(loo, file.path(out, "leave_one_cohort_out.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Compare the reconstructed primary cohort effects and modified-KH results with the frozen release.
old_eff <- read.delim(file.path(src, "three_cohort_adjusted_meta_input_effects.tsv"), check.names = FALSE)
old_eff <- old_eff[old_eff$metric == "positive_abundance", ]
cmp <- merge(
  cohort_effects[cohort_effects$metric == "positive_abundance", c("cohort", "pathway_id", "variant", "estimate_log_odds", "se", "n")],
  transform(old_eff, variant = ifelse(meta_family == "no_metformin_sensitivity", "known_no_metformin", "primary_demographic_adjusted"))[, c("cohort", "pathway_id", "variant", "estimate", "SE", "n_analysis")],
  by = c("cohort", "pathway_id", "variant"), all.x = TRUE
)
cmp$delta_estimate <- cmp$estimate_log_odds - cmp$estimate
cmp$delta_se <- cmp$se - cmp$SE
cmp$delta_n <- cmp$n - cmp$n_analysis
cmp$reproduced_1e_8 <- abs(cmp$delta_estimate) < 1e-8 & abs(cmp$delta_se) < 1e-8 & cmp$delta_n == 0
write.table(cmp, file.path(out, "reproduction_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
cat(sprintf("Rebuilt %d cohort-model rows; %d/%d comparable positive-abundance rows reproduced at 1e-8.\n", nrow(cohort_effects), sum(cmp$reproduced_1e_8, na.rm = TRUE), sum(!is.na(cmp$estimate))))
