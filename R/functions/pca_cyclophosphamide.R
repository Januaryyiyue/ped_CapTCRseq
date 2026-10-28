# Pooled cyclophosphamide models and exposure audit used by notebook 15.
# All fits start from common PCA scores and the original treatment table.

.pca_bind_rows <- function(...) {
  tables <- list(...)
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(d) {
    class(d) <- "data.frame"
    for (column in setdiff(columns, names(d))) d[[column]] <- NA
    d[, columns, drop = FALSE]
  })
  result <- do.call(rbind, tables)
  rownames(result) <- NULL
  result
}

fit_pca_cyclophosphamide <- function(inputs, soc) {
  stopifnot(!anyDuplicated(soc$Patient))
  datasets <- list(Original = inputs$real_metadata75, missMDA = inputs$metadata235)
  exposed_label <- "Any cyclophosphamide"
  warnings <- character()
  capture_warning <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
  # Audit drug membership using the recorded previous-cycle drug list.
  sources <- coverage <- composition <- audits <- pairs <- list()
  for (cohort in names(datasets)) {
    d <- datasets[[cohort]]
    stopifnot(!anyDuplicated(d$sample_id), !anyDuplicated(paste(d$Patient, d$cycle)),
              all(is.finite(d$PC2)), !anyNA(d$Age), !anyNA(d$cancergroup))
    d$contains_cyclophosphamide <- grepl('cyclophosphamide', d$regimen,
                                       ignore.case = TRUE)
    d$prior_cycle_drugs <- vapply(seq_len(nrow(d)), function(i) {
      if (d$cycle_num[i] == 1L) return(NA_character_)
      patient <- match(d$Patient[i], soc$Patient)
      column <- paste0('Cycle ', d$cycle_num[i] - 1L, ' Drugs')
      if (is.na(patient) || !column %in% names(soc))
        stop('Cannot map prior-cycle drug list: ', d$sample_id[i])
      as.character(soc[[column]][patient])
    }, character(1))
    post <- d$cycle_num > 1L
    recorded <- post & !is.na(d$prior_cycle_drugs) & nzchar(d$prior_cycle_drugs)
    raw_cyclophosphamide <- grepl('cyclophosphamide', d$prior_cycle_drugs[recorded],
                                 ignore.case = TRUE)
    stopifnot(identical(d$contains_cyclophosphamide[recorded], raw_cyclophosphamide),
              all(recorded[d$contains_cyclophosphamide]))
    d$included_in_model <- d$cycle_num == 1L | recorded
    d$exclusion_reason <- ifelse(d$included_in_model, '', 'Missing prior-cycle drug list')
    d$pooled_regimen <- ifelse(d$contains_cyclophosphamide, exposed_label,
                               as.character(d$regimen))
    reference_levels <- unique(c('baseline', exposed_label,
      names(inputs$therapy_palette)[!grepl('cyclophosphamide',
        names(inputs$therapy_palette), ignore.case = TRUE)]))
    d$pooled_regimen <- droplevels(factor(d$pooled_regimen, levels = reference_levels))
    d$Patient <- factor(d$Patient)
    d$cancergroup <- factor(d$cancergroup)
    stopifnot(!anyNA(d$pooled_regimen), all(d$regimen[d$cycle_num == 1L] == 'baseline'))
    datasets[[cohort]] <- droplevels(d[d$included_in_model, ])
    b <- d[d$cycle_num == 1L, ]
    e <- d[d$contains_cyclophosphamide, ]
    first <- e[order(e$Patient, e$cycle_num), ]
    first <- first[!duplicated(first$Patient), ]
    first <- first[first$Patient %in% b$Patient, ]
    baseline_index <- match(first$Patient, b$Patient)
    pairs[[cohort]] <- data.frame(dataset = cohort, Patient = as.character(first$Patient),
      baseline_sample = b$sample_id[baseline_index], post_sample = first$sample_id,
      post_cycle = first$cycle, regimen = as.character(first$regimen),
      baseline_PC2 = b$PC2[baseline_index], post_PC2 = first$PC2,
      delta_PC2 = first$PC2 - b$PC2[baseline_index])
    coverage[[cohort]] <- data.frame(dataset = cohort,
      n_input_samples = nrow(d), n_excluded_unknown_drug_list = sum(!d$included_in_model),
      n_model_samples = sum(d$included_in_model),
      n_model_patients = length(unique(d$Patient[d$included_in_model])),
      n_baseline_samples = nrow(b), n_baseline_patients = length(unique(b$Patient)),
      n_exposed_samples = nrow(e), n_exposed_patients = length(unique(e$Patient)),
      n_paired_patients = nrow(first), baseline_PC2_median = median(b$PC2),
      exposed_PC2_median = median(e$PC2))
    composition[[cohort]] <- aggregate(sample_id ~ regimen + Disease_type,
      data = transform(e, regimen = as.character(regimen)), FUN = length)
    names(composition[[cohort]])[3] <- 'n_samples'
    composition[[cohort]]$dataset <- cohort
    audits[[cohort]] <- data.frame(dataset = cohort, n_followup_audited = sum(recorded),
                                   n_missing_followup_drug_lists = sum(post & !recorded),
                                   n_drug_label_mismatches = 0L)
    d$dataset <- cohort
    sources[[cohort]] <- d[, c('dataset','sample_id','Patient','cycle','cycle_num',
      'Age','cancergroup','Disease_type','regimen','pooled_regimen',
      'contains_cyclophosphamide','prior_cycle_drugs','included_in_model',
      'exclusion_reason','PC2','any_feature_imputed')]
  }

  # Keep other treatment categories separate so the model specification matches
  # the existing Figure 2 analysis except for pooling the three exposed categories.
  results <- diagnostics <- list()
  for (adjusted in c(FALSE, TRUE)) {
    analysis <- if (adjusted) 'Age_cancer_adjusted_sensitivity' else
      'Current_patient_intercept_primary'
    for (cohort in names(datasets)) {
      d <- datasets[[cohort]]
      formula <- if (adjusted) PC2 ~ pooled_regimen + Age + cancergroup + (1 | Patient)
        else PC2 ~ pooled_regimen + (1 | Patient)
      fit <- withCallingHandlers(lmerTest::lmer(formula, data = d, REML = TRUE,
        control = lme4::lmerControl(optimizer = 'bobyqa',
                                   optCtrl = list(maxfun = 100000))),
        warning = capture_warning)
      em <- emmeans::emmeans(fit, ~ pooled_regimen, lmer.df = 'satterthwaite',
                             weights = 'proportional')
      levels <- levels(d$pooled_regimen)
      contrast <- as.numeric(levels == exposed_label) - as.numeric(levels == 'baseline')
      r <- as.data.frame(summary(emmeans::contrast(em,
        method = list(any_cyclophosphamide_vs_baseline = contrast), adjust = 'none'),
        infer = c(TRUE, TRUE), adjust = 'none'))
      r$analysis <- analysis
      r$dataset <- cohort
      r$PC <- 'PC2'
      r$regimen <- exposed_label
      r$score_SD <- sd(d$PC2)
      r$standardized_shift <- r$estimate / r$score_SD
      r$standardized_lower <- r$lower.CL / r$score_SD
      r$standardized_upper <- r$upper.CL / r$score_SD
      results[[length(results) + 1L]] <- r
      vc <- as.data.frame(lme4::VarCorr(fit))
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        analysis = analysis, dataset = cohort, PC = "PC2",
        model_scope = "Pooled cyclophosphamide", n_samples = nrow(d),
        n_patients = length(unique(d$Patient)), singular = lme4::isSingular(fit),
        patient_intercept_variance = vc$vcov[vc$grp == 'Patient'][1],
        residual_SD = sigma(fit), convergence = paste(
          fit@optinfo$conv$lme4$messages, collapse = ';'))
    }
  }
  results <- do.call(rbind, results)


  paired <- paired_changes <- list()
  for (cohort in names(pairs)) {
    p <- pairs[[cohort]]
    delta <- p$delta_PC2
    nonzero <- delta[abs(delta) > 1e-12]
    rank_test <- if (length(nonzero)) withCallingHandlers(
      wilcox.test(nonzero, mu = 0, exact = TRUE), warning = capture_warning) else NULL
    ranks <- rank(abs(nonzero))
    paired[[cohort]] <- data.frame(dataset = cohort, regimen = exposed_label, PC = "PC2",
      n_paired_patients = nrow(p), n_nonzero = length(nonzero),
      p = if (is.null(rank_test)) 1 else rank_test$p.value,
      median_delta = median(delta), median_lower = NA_real_, median_upper = NA_real_,
      rank_biserial = if (length(ranks)) sum(sign(nonzero) * ranks) / sum(ranks) else 0,
      rank_biserial_lower = NA_real_, rank_biserial_upper = NA_real_,
      n_decreased = sum(delta < 0), n_increased = sum(delta > 0),
      p_sign = if (length(nonzero)) binom.test(sum(nonzero > 0), length(nonzero), .5)$p.value else 1,
      test = if (is.null(rank_test)) "All paired changes zero" else rank_test$method)
    paired_changes[[cohort]] <- data.frame(dataset = cohort, regimen = exposed_label,
      PC = "PC2", Patient = p$Patient, baseline_sample = p$baseline_sample,
      post_sample = p$post_sample, cycle = p$post_cycle,
      baseline = p$baseline_PC2, post = p$post_PC2, delta = p$delta_PC2)
  }
  list(regimen_contrasts = results, paired_rank = do.call(rbind, paired),
       paired_changes = do.call(rbind, paired_changes),
       diagnostics = do.call(rbind, diagnostics),
       analysis_source = do.call(rbind, sources), coverage = do.call(rbind, coverage),
       composition = do.call(rbind, composition), exposure_audit = do.call(rbind, audits),
       paired_patient_changes = do.call(rbind, pairs), warnings = unique(warnings))
}

write_pca_cyclophosphamide_tables <- function(statistics, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  csv <- function(x, filename) write.csv(x, file.path(output_dir, filename),
                                        row.names = FALSE, na = "")
  pooled <- statistics$cyclophosphamide
  selected <- statistics$regimen_contrasts$regimen == "Any cyclophosphamide"
  results <- statistics$regimen_contrasts[selected, ]
  results$p_Holm_extended <- results$p_Holm
  results$lower_pointwise <- results$lower.CL_pointwise
  results$upper_pointwise <- results$upper.CL_pointwise
  csv(results, "model_results.csv")
  csv(pooled$diagnostics, "model_diagnostics.csv")
  csv(pooled$coverage, "cohort_counts.csv")
  csv(pooled$composition, "exposed_regimen_cancer_composition.csv")
  csv(pooled$exposure_audit, "exposure_audit.csv")
  csv(pooled$analysis_source, "analysis_source.csv")
  csv(pooled$paired_patient_changes, "paired_patient_changes.csv")
  paired <- statistics$paired_rank[statistics$paired_rank$regimen == "Any cyclophosphamide", ]
  paired$p_Holm_extended <- paired$p_Holm
  paired$family_n_tests <- nrow(statistics$paired_rank)
  csv(paired, "paired_sensitivity.csv")
  family <- statistics$regimen_contrasts[, c("analysis", "dataset", "PC", "regimen", "p.value", "p_Holm", "ci_n_tests")]
  names(family)[5:7] <- c("p_raw", "p_Holm_extended", "n_tests")
  csv(family, "expanded_model_multiplicity_families.csv")
  family <- statistics$paired_rank[, c("dataset", "PC", "regimen", "p", "p_Holm")]
  names(family)[4:5] <- c("p_raw", "p_Holm_extended")
  family$n_tests <- nrow(family)
  csv(family, "expanded_paired_multiplicity_family.csv")
  writeLines(pooled$warnings, file.path(output_dir, "model_warnings.txt"))
  writeLines(c("# PC2 and cyclophosphamide-containing regimens", "",
    "Notebook 15 fits these models and writes the numerical data used by Figure 2G.",
    "The standalone analysis script uses the same helpers and correction families.",
    "Exposure pools Cyclophosphamide, Vincristine/Cyclophosphamide and",
    "Vincristine/Methotrexate/Cyclophosphamide. Prior-cycle drug lists confirm membership.",
    "Only the pooled models exclude unknown drug lists (2 Original; 4 missMDA).",
    "Primary and age/cancer-adjusted families each contain 32 tests; paired checks contain six.",
    "P values use Holm correction; model intervals are 95% simultaneous Bonferroni intervals.",
    "The pooled comparison was added after reviewing the earlier results and is exploratory.",
    "Cohorts overlap, co-administered drugs differ, and imputation uncertainty is not propagated.",
    "Paired median changes are descriptive, without a substituted confidence interval.",
    "Models combine within-patient and between-patient information; see the paired sensitivity.",
    "Reproduce: Rscript R/analysis_PC2_any_cyclophosphamide.R"), file.path(output_dir, "README.md"))
  invisible(results)
}
