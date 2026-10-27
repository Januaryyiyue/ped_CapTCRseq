# Prepared numerical outputs written by notebook 15 and read by figure scripts.
# Reading this bundle never runs imputation, PCA fitting or statistical models.

.pca_figure_schema <- 1L

.validate_pca_figure_inputs <- function(inputs) {
  required <- c("metadata235", "real_metadata75", "raw235", "completed235",
                "real_features75", "reference_pca", "scoredata310", "disease_order",
                "therapy_palette", "features", "variance_percent")
  if (!is.list(inputs) || !all(required %in% names(inputs)))
    stop("Incomplete notebook-15 PCA input bundle", call. = FALSE)
  raw <- inputs$raw235
  completed <- inputs$completed235
  full <- inputs$metadata235
  original <- inputs$real_metadata75
  scores <- inputs$scoredata310
  ref <- inputs$reference_pca
  stopifnot(is.matrix(raw), is.matrix(completed), identical(dimnames(raw), dimnames(completed)),
            all(is.finite(completed)), all(is.finite(raw) | is.na(raw)),
            identical(colnames(raw), inputs$features),
            !anyDuplicated(full$sample_id), !anyDuplicated(original$sample_id),
            identical(as.character(full$sample_id), rownames(raw)),
            identical(as.character(original$sample_id), rownames(inputs$real_features75)),
            !anyDuplicated(paste(full$Patient, full$cycle)),
            setequal(original$sample_id, full$sample_id[rowSums(is.na(raw)) == 0L]),
            identical(dimnames(inputs$real_features75),
                      dimnames(raw[match(original$sample_id, rownames(raw)), , drop = FALSE])),
            identical(rownames(ref$rotation), colnames(raw)),
            ncol(ref$rotation) >= 2L, ncol(ref$x) >= 2L,
            all(is.finite(ref$rotation)), all(is.finite(ref$x)),
            all(as.character(full$Disease_type) %in% inputs$disease_order),
            all(as.character(full$exposed_regimen) %in% names(inputs$therapy_palette)),
            all(as.character(scores$dataset) %in% c("Original", "missMDA")),
            !anyDuplicated(paste(scores$dataset, scores$sample_id)))
  observed <- !is.na(raw)
  stopifnot(max(abs(completed[observed] - raw[observed])) == 0,
            max(abs(inputs$real_features75 -
              raw[match(original$sample_id, rownames(raw)), , drop = FALSE])) == 0)
  for (cohort in c("Original", "missMDA")) {
    data <- if (cohort == "Original") original else full
    plotted <- scores[as.character(scores$dataset) == cohort, , drop = FALSE]
    stopifnot(setequal(plotted$sample_id, data$sample_id))
    plotted <- plotted[match(data$sample_id, plotted$sample_id), , drop = FALSE]
    stopifnot(max(abs(as.matrix(plotted[, c("PC1", "PC2")]) -
                       as.matrix(data[, c("PC1", "PC2")]))) < 1e-9)
  }
  full_scores <- ref$x[match(full$sample_id, rownames(ref$x)), 1:2, drop = FALSE]
  overlap <- full[match(original$sample_id, full$sample_id), c("PC1", "PC2")]
  expected_variance <- 100 * ref$sdev[1:2]^2 / sum(ref$sdev^2)
  stopifnot(all(is.finite(full_scores)),
            max(abs(full_scores - as.matrix(full[, c("PC1", "PC2")]))) < 1e-9,
            max(abs(as.matrix(original[, c("PC1", "PC2")]) - as.matrix(overlap))) < 1e-9,
            max(abs(inputs$variance_percent - expected_variance)) < 1e-9)
  invisible(TRUE)
}

.validate_pca_figure_statistics <- function(statistics) {
  required <- c("baseline_lymphoma", "regimen_contrasts", "paired_rank", "families")
  if (!is.list(statistics) || !all(required %in% names(statistics)))
    stop("Incomplete notebook-15 association result bundle", call. = FALSE)
  baseline <- statistics$baseline_lymphoma
  regimen <- statistics$regimen_contrasts
  paired <- statistics$paired_rank
  ci_columns <- c("SE", "df", "ci_adjustment", "ci_family_level", "ci_n_tests",
                  "ci_individual_level")
  if (!all(ci_columns %in% names(baseline)) || !all(ci_columns %in% names(regimen)))
    stop("Association results need simultaneous confidence intervals. ",
         "Rerun notebook 15 before plotting Figure 2.", call. = FALSE)
  stopifnot(setequal(baseline$dataset, c("Original", "missMDA")),
            !anyDuplicated(baseline$dataset),
            all(c("estimate", "lower", "upper", "p_Holm", "p_rank_Holm") %in% names(baseline)),
            all(c("analysis", "PC", "regimen", "estimate", "lower.CL", "upper.CL", "p_Holm") %in% names(regimen)),
            all(c("PC", "regimen", "p_Holm", "rank_biserial", "n_paired_patients") %in% names(paired)))
  for (data in list(baseline, regimen, paired)) {
    stopifnot(all(is.finite(data$p_Holm)), all(data$p_Holm >= 0 & data$p_Holm <= 1))
  }
  family_counts <- c(nrow(baseline),
    sum(regimen$analysis == "Current_patient_intercept_primary"),
    sum(regimen$analysis == "Age_cancer_adjusted_sensitivity"), nrow(paired))
  stopifnot(identical(as.integer(statistics$families$n_tests), as.integer(family_counts)),
            all(statistics$families$correction == "Holm"))
  for (data in list(baseline, regimen)) {
    lower <- if ("lower.CL" %in% names(data)) data$lower.CL else data$lower
    upper <- if ("upper.CL" %in% names(data)) data$upper.CL else data$upper
    half_width <- qt(1 - .05 / (2 * data$ci_n_tests), data$df) * data$SE
    stopifnot(all(data$ci_adjustment == "Bonferroni"),
              all(data$ci_family_level == .95),
              max(abs(data$ci_individual_level - (1 - .05 / data$ci_n_tests))) < 1e-12,
              max(abs(lower - (data$estimate - half_width))) < 1e-9,
              max(abs(upper - (data$estimate + half_width))) < 1e-9)
  }
  stopifnot(all(baseline$ci_n_tests == nrow(baseline)))
  for (analysis in unique(regimen$analysis)) {
    rows <- regimen$analysis == analysis
    stopifnot(all(regimen$ci_n_tests[rows] == sum(rows)))
  }
  invisible(TRUE)
}

.assert_numeric_figure_bundle <- function(object) {
  if (inherits(object, c("ggplot", "gg", "patchwork", "grob", "gTree", "gtable")) ||
      is.function(object) || is.environment(object) || isS4(object))
    stop("Notebook figure outputs must contain numerical results, not plots or fitted mixed models",
         call. = FALSE)
  if (is.list(object)) invisible(lapply(object, .assert_numeric_figure_bundle))
  invisible(TRUE)
}

write_pca_figure_outputs <- function(inputs, statistics, project_root = getwd(),
                                     output_dir = file.path(project_root, "data")) {
  .validate_pca_figure_inputs(inputs)
  .validate_pca_figure_statistics(statistics)
  # Store numerical results and the PCA reference; no saved plot objects.
  prepared <- inputs[setdiff(names(inputs), c("project", "old_real_pca"))]
  if (!is.null(prepared$source_paths)) {
    prepared$source_paths <- setNames(file.path("data", basename(prepared$source_paths)),
                                     names(prepared$source_paths))
  }
  if (!is.null(prepared$source_hashes)) {
    names(prepared$source_hashes) <- file.path("data", basename(names(prepared$source_hashes)))
  }
  .assert_numeric_figure_bundle(prepared)
  .assert_numeric_figure_bundle(statistics)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(inputs = "pca_inputs.rds", statistics = "association_results.rds")
  objects <- list(inputs = list(schema_version = .pca_figure_schema, inputs = prepared),
                  statistics = list(schema_version = .pca_figure_schema, statistics = statistics))
  for (name in names(files)) {
    temporary <- tempfile(pattern = ".notebook15_", tmpdir = output_dir)
    saveRDS(objects[[name]], temporary, version = 3)
    if (!file.rename(temporary, file.path(output_dir, files[[name]]))) {
      unlink(temporary)
      stop("Could not write notebook-15 output: ", files[[name]], call. = FALSE)
    }
  }
  invisible(file.path(output_dir, files))
}

read_pca_figure_outputs <- function(project_root = getwd(), include_statistics = TRUE) {
  filenames <- c("pca_inputs.rds", if (include_statistics) "association_results.rds")
  paths <- file.path(project_root, "data", filenames)
  instruction <- "Run all cells in notebooks/15_PCA_imputation_associations.ipynb first."
  if (!all(file.exists(paths))) {
    stop("Missing notebook-15 output: ",
         paste(file.path("data", filenames[!file.exists(paths)]), collapse = ", "),
         ". ", instruction, call. = FALSE)
  }
  bundle <- readRDS(paths[1])
  if (!identical(bundle$schema_version, .pca_figure_schema))
    stop("Unsupported notebook-15 PCA schema. ", instruction, call. = FALSE)
  .validate_pca_figure_inputs(bundle$inputs)
  statistics <- NULL
  if (include_statistics) {
    results <- readRDS(paths[2])
    if (!identical(results$schema_version, .pca_figure_schema))
      stop("Unsupported notebook-15 statistics schema. ", instruction, call. = FALSE)
    .validate_pca_figure_statistics(results$statistics)
    statistics <- results$statistics
  }
  list(inputs = bundle$inputs, statistics = statistics)
}
