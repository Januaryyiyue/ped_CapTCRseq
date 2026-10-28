# Patient-aware associations computed by notebook 15 on a common PCA reference.
# All numerical quantities are recomputed; this file reads no previous results.

# Holm is retained for tests. Bonferroni gives simultaneous confidence intervals
# across the same family without assumptions about dependence between contrasts.
.pca_simultaneous_t_intervals <- function(table, lower, upper, level = .95) {
  n_tests <- nrow(table)
  stopifnot(n_tests > 0L, all(is.finite(table$estimate)),
            all(is.finite(table$SE)), all(table$SE >= 0), all(table$df > 0))
  table[[paste0(lower, "_pointwise")]] <- table[[lower]]
  table[[paste0(upper, "_pointwise")]] <- table[[upper]]
  half_width <- qt(1 - (1 - level) / (2 * n_tests), table$df) * table$SE
  table[[lower]] <- table$estimate - half_width
  table[[upper]] <- table$estimate + half_width
  table$ci_adjustment <- "Bonferroni"
  table$ci_family_level <- level
  table$ci_n_tests <- n_tests
  table$ci_individual_level <- 1 - (1 - level) / n_tests
  table
}

fit_pca_associations <- function(inputs, soc) {
  for (pkg in c("lme4", "lmerTest", "emmeans")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Install the required R package: ", pkg, call. = FALSE)
    }
  }
  datasets <- list(Original = inputs$real_metadata75,
                   missMDA = inputs$metadata235)

  # HC3 sandwich covariance with residual-t reference, retaining the manuscript
  # analysis specification. There is one X01 record per patient in each cohort.
  hc3_contrast <- function(fit, term) {
    X <- model.matrix(fit)
    bread <- solve(crossprod(X))
    adjusted_residuals <- residuals(fit) / (1 - hatvalues(fit))
    variance <- bread %*% crossprod(X, X * as.numeric(adjusted_residuals^2)) %*% bread
    estimate <- unname(coef(fit)[term])
    se <- unname(sqrt(diag(variance))[term])
    df <- df.residual(fit)
    half_width <- qt(.975, df) * se
    c(estimate = estimate, SE = se, df = df,
      lower = estimate - half_width, upper = estimate + half_width,
      p = 2 * pt(-abs(estimate / se), df))
  }

  baseline <- recorded_group <- list()
  for (cohort in names(datasets)) {
    z <- datasets[[cohort]]
    z <- z[z$cycle == "X01", ]
    z$Lymphoma <- as.integer(z$lymphoma == "Lymphoma")
    stopifnot(!anyDuplicated(z$Patient), length(unique(z$Lymphoma)) == 2L)
    fit <- lm(PC1 ~ Lymphoma + Age, data = z)
    result <- hc3_contrast(fit, "Lymphoma")
    lymph <- z$PC1[z$Lymphoma == 1L]
    other <- z$PC1[z$Lymphoma == 0L]
    rank_test <- suppressWarnings(wilcox.test(lymph, other, exact = TRUE))
    baseline[[cohort]] <- data.frame(
      dataset = cohort, n = nrow(z), n_lymphoma = length(lymph),
      n_non_lymphoma = length(other), t(result), p_rank = rank_test$p.value,
      rank_biserial = 2 * as.numeric(rank_test$statistic) /
        (length(lymph) * length(other)) - 1,
      standardized_shift = result[["estimate"]] / sigma(fit), row.names = NULL)

    # ALCL is clinically a lymphoma; this sensitivity retains the source's
    # recorded cancer classes, where ALCL belongs to the T-cell class.
    z$recordedLymphoma <- as.integer(z$cancergroup == "Lymphoma")
    sensitivity <- hc3_contrast(lm(PC1 ~ recordedLymphoma + Age, data = z),
                               "recordedLymphoma")
    recorded_group[[cohort]] <- data.frame(
      dataset = cohort, t(sensitivity), n_alcl_baseline = sum(z$Disease_type == "ALCL"),
      definition = "Recorded cancer group (ALCL retained in T-cell class)", row.names = NULL)
  }
  baseline <- do.call(rbind, baseline)
  rownames(baseline) <- NULL
  baseline$p_Holm <- p.adjust(baseline$p, "holm")
  baseline$p_rank_Holm <- p.adjust(baseline$p_rank, "holm")
  baseline <- .pca_simultaneous_t_intervals(baseline, "lower", "upper")
  recorded_group <- do.call(rbind, recorded_group)
  rownames(recorded_group) <- NULL
  recorded_group$p_Holm <- p.adjust(recorded_group$p, "holm")
  recorded_group <- .pca_simultaneous_t_intervals(recorded_group, "lower", "upper")

  # The cohorts overlap. Each is fitted separately; no statistical comparison
  # between cohorts is made. Include all estimable regimen-v-baseline contrasts
  # for both PCs in the multiplicity family, rather than only plotted hits.
  contrasts <- diagnostics <- list()
  model_warnings <- character()
  for (adjusted in c(FALSE, TRUE)) {
    for (cohort in names(datasets)) {
      for (pc in c("PC1", "PC2")) {
        z <- datasets[[cohort]]
        z$score <- z[[pc]]
        z$regimen <- droplevels(factor(z$regimen, levels = names(inputs$therapy_palette)))
        z$Patient <- factor(z$Patient)
        z$cancergroup <- factor(z$cancergroup)
        formula <- if (adjusted) {
          score ~ regimen + Age + cancergroup + (1 | Patient)
        } else score ~ regimen + (1 | Patient)
        fit <- withCallingHandlers(
          lmerTest::lmer(formula, data = z, REML = TRUE,
                        control = lme4::lmerControl(optimizer = "bobyqa",
                                                    optCtrl = list(maxfun = 100000))),
          warning = function(w) {
            model_warnings <<- c(model_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          })
        marginal <- emmeans::emmeans(fit, ~ regimen, lmer.df = "satterthwaite",
                                    weights = "proportional")
        result <- as.data.frame(summary(
          emmeans::contrast(marginal, "trt.vs.ctrl", ref = 1, adjust = "none"),
          infer = c(TRUE, TRUE), adjust = "none"))
        result$regimen <- levels(z$regimen)[-1]
        result$dataset <- cohort
        result$PC <- pc
        result$analysis <- if (adjusted) "Age_cancer_adjusted_sensitivity" else
          "Current_patient_intercept_primary"
        result$score_SD <- sd(z$score)
        result$standardized_shift <- result$estimate / result$score_SD
        result$standardized_lower <- result$lower.CL / result$score_SD
        result$standardized_upper <- result$upper.CL / result$score_SD
        contrasts[[length(contrasts) + 1L]] <- result
        diagnostics[[length(diagnostics) + 1L]] <- data.frame(
          dataset = cohort, PC = pc, analysis = result$analysis[1],
          model_scope = "Regimen categories",
          n_samples = nrow(z), n_patients = length(unique(z$Patient)),
          singular = lme4::isSingular(fit),
          convergence = paste(fit@optinfo$conv$lme4$messages, collapse = ";"))
      }
    }
  }
  # Append the additional pooled PC2 hypothesis before the shared correction.
  # The original regimen models retain all samples and their existing fits.
  pooled <- fit_pca_cyclophosphamide(inputs, soc)
  contrasts <- .pca_bind_rows(do.call(rbind, contrasts), pooled$regimen_contrasts)
  model_warnings <- c(model_warnings, pooled$warnings)
  contrasts$p_Holm <- NA_real_
  contrasts$standardized_lower_pointwise <- contrasts$standardized_lower
  contrasts$standardized_upper_pointwise <- contrasts$standardized_upper
  families <- list()
  for (analysis in unique(contrasts$analysis)) {
    rows <- which(contrasts$analysis == analysis)
    contrasts$p_Holm[rows] <- p.adjust(contrasts$p.value[rows], "holm")
    families[[analysis]] <- .pca_simultaneous_t_intervals(
      contrasts[rows, , drop = FALSE], "lower.CL", "upper.CL")
  }
  contrasts <- do.call(rbind, families)
  rownames(contrasts) <- NULL
  contrasts$standardized_lower <- contrasts$lower.CL / contrasts$score_SD
  contrasts$standardized_upper <- contrasts$upper.CL / contrasts$score_SD
  attr(contrasts, "mesg") <- c(
    "Degrees-of-freedom method: satterthwaite",
    "95% simultaneous confidence intervals: Bonferroni within each analysis family",
    "P values: Holm within each analysis family")

  # A paired sensitivity: one first eligible post-treatment visit per patient
  # and regimen, matched to that patient's own X01. Not the raincloud's full
  # sample population and not the primary mixed-model inference.
  pairs <- paired_rank <- list()
  for (cohort in names(datasets)) {
    for (regimen in c("Vincristine/Cyclophosphamide", "Vincristine")) {
      pc <- if (regimen == "Vincristine/Cyclophosphamide") "PC1" else "PC2"
      z <- datasets[[cohort]]
      baseline_rows <- z[z$cycle == "X01", ]
      post <- z[z$regimen == regimen & z$cycle != "X01", ]
      post <- post[order(post$Patient, post$cycle_num), ]
      post <- post[!duplicated(post$Patient), ]
      match_rows <- match(post$Patient, baseline_rows$Patient)
      post <- post[!is.na(match_rows), ]
      match_rows <- match_rows[!is.na(match_rows)]
      delta <- post[[pc]] - baseline_rows[[pc]][match_rows]
      stopifnot(length(delta) > 0L)
      pairs[[length(pairs) + 1L]] <- data.frame(
        dataset = cohort, regimen = regimen, PC = pc, Patient = post$Patient,
        baseline_sample = baseline_rows$sample_id[match_rows], post_sample = post$sample_id,
        cycle = post$cycle, baseline = baseline_rows[[pc]][match_rows],
        post = post[[pc]], delta = delta)
      nonzero <- delta[abs(delta) > 1e-12]
      rank_test <- if (length(nonzero)) suppressWarnings(
        wilcox.test(nonzero, mu = 0, exact = TRUE)) else NULL
      rank_effect <- function(values) {
        values <- values[abs(values) > 1e-12]
        ranks <- rank(abs(values))
        if (length(ranks)) sum(sign(values) * ranks) / sum(ranks) else 0
      }
      set.seed(20261022L + length(paired_rank))
      ci <- if (length(delta) >= 5L) quantile(
        replicate(2000, median(sample(delta, replace = TRUE))), c(.025, .975),
        names = FALSE) else c(NA_real_, NA_real_)
      rank_ci <- if (length(delta) >= 5L) quantile(
        replicate(2000, rank_effect(sample(delta, replace = TRUE))), c(.025, .975),
        names = FALSE) else c(NA_real_, NA_real_)
      paired_rank[[length(paired_rank) + 1L]] <- data.frame(
        dataset = cohort, regimen = regimen, PC = pc,
        n_paired_patients = length(delta), n_nonzero = length(nonzero),
        p = if (is.null(rank_test)) 1 else rank_test$p.value,
        median_delta = median(delta), median_lower = ci[1], median_upper = ci[2],
        rank_biserial = rank_effect(delta), rank_biserial_lower = rank_ci[1],
        rank_biserial_upper = rank_ci[2], n_decreased = sum(delta < 0),
        n_increased = sum(delta > 0),
        p_sign = if (length(nonzero)) binom.test(sum(nonzero > 0), length(nonzero), .5)$p.value else 1,
        test = if (is.null(rank_test)) "All paired changes zero" else rank_test$method)
    }
  }
  paired_rank <- .pca_bind_rows(do.call(rbind, paired_rank), pooled$paired_rank)
  paired_rank$p_Holm <- p.adjust(paired_rank$p, "holm")
  paired_rank$p_sign_Holm <- p.adjust(paired_rank$p_sign, "holm")

  list(baseline_lymphoma = baseline, recorded_group_sensitivity = recorded_group,
       regimen_contrasts = contrasts, paired_rank = paired_rank,
       paired_changes = .pca_bind_rows(do.call(rbind, pairs), pooled$paired_changes),
       diagnostics = .pca_bind_rows(do.call(rbind, diagnostics), pooled$diagnostics),
       cyclophosphamide = pooled[c("diagnostics", "analysis_source", "coverage",
         "composition", "exposure_audit", "paired_patient_changes", "warnings")],
       warnings = unique(model_warnings),
       families = data.frame(
         family = c("Baseline cancer contrasts", "Primary regimen contrasts",
                    "Adjusted regimen sensitivity", "Selected paired rank sensitivity"),
         n_tests = c(nrow(baseline), sum(contrasts$analysis == "Current_patient_intercept_primary"),
                     sum(contrasts$analysis == "Age_cancer_adjusted_sensitivity"), nrow(paired_rank)),
         correction = "Holm",
         ci_adjustment = c("Bonferroni", "Bonferroni", "Bonferroni",
                           "None (descriptive bootstrap median and rank effect)"),
         ci_family_level = c(.95, .95, .95, NA_real_)))
}

write_pca_association_tables <- function(statistics, tables_dir) {
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(baseline_lymphoma = "baseline_lymphoma_statistics.csv",
             recorded_group_sensitivity = "recorded_group_sensitivity.csv",
             regimen_contrasts = "all_regimen_contrasts_primary_and_adjusted.csv",
             paired_rank = "paired_rank_sensitivity_and_effects.csv",
             paired_changes = "paired_first_exposure_changes.csv",
             diagnostics = "mixed_model_diagnostics.csv",
             families = "multiplicity_families.csv")
  for (name in names(files)) {
    write.csv(statistics[[name]], file.path(tables_dir, files[[name]]), row.names = FALSE)
  }
  writeLines(statistics$warnings, file.path(tables_dir, "model_warnings.txt"))
  write_pca_cyclophosphamide_tables(statistics, file.path(tables_dir, "cyclophosphamide"))
  invisible(files)
}
