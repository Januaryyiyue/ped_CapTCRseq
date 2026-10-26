# Numerical preparation run by notebook 15 for Figure 2, Extended 2 and Figure 3.
# Figure scripts read its prepared data/pca_inputs.rds output.
# Inputs: data/metadata.rds and data/soc_good_with_baseregimen.rds.
# All calculations start from these upstream data; no stored plot or analysis
# result is read. Outcome, treatment, cancer, age and patient labels are never
# predictors for imputation or PCA. This is a descriptive, all-timepoint PCA,
# not a prospective prediction model.

prepare_pca_inputs <- function(project_root = getwd(), ncp = 2L,
                               seed = 20261012L, rerun_rank_cv = FALSE) {
  if (!requireNamespace("missMDA", quietly = TRUE))
    stop("Install missMDA and its dependencies before rerunning notebook 15.")
  project_root <- normalizePath(project_root, mustWork = TRUE)
  source_paths <- c(
    metadata = file.path(project_root, "data", "metadata.rds"),
    soc = file.path(project_root, "data", "soc_good_with_baseregimen.rds")
  )
  if (!all(file.exists(source_paths)))
    stop("Missing required inputs: ", paste(source_paths[!file.exists(source_paths)], collapse = ", "))
  source_hashes <- tools::md5sum(source_paths)
  metadata <- readRDS(source_paths[["metadata"]])
  soc <- readRDS(source_paths[["soc"]])
  if (anyDuplicated(soc$Patient)) stop("Treatment table has duplicated patients")

  # PBMC anchors each patient/timepoint. The corresponding cfDNA row supplies
  # cfTCR Shannon diversity. Never impute between unrelated sample IDs.
  pbmc <- metadata[!is.na(metadata$sampletype) & metadata$sampletype == "PBMC", , drop = FALSE]
  cfdna <- metadata[!is.na(metadata$sampletype) & metadata$sampletype == "cfDNA", , drop = FALSE]
  if (anyDuplicated(pbmc$sample_id) || anyDuplicated(cfdna$sample_id))
    stop("Sample IDs must be unique within each sample type")
  pbmc$Patient <- as.character(pbmc$Patient)
  pbmc$cycle <- as.character(pbmc$cycle)
  pbmc$cfShannon <- cfdna$observed_Shannon[match(pbmc$sample_id, cfdna$sample_id)]
  pbmc$cycle_num <- as.integer(sub("^X0*", "", pbmc$cycle))
  if (anyNA(pbmc$cycle_num) || any(pbmc$cycle_num < 1L))
    stop("Expected cycle labels X01, X02, ...")

  therapy_palette <- c(
    "baseline" = "black", "Vincristine" = "#E41A1C",
    "Methotrexate" = "#377EB8", "Cyclophosphamide" = "#4DAF4A",
    "Vincristine/Methotrexate" = "#984EA3",
    "Vincristine/Cyclophosphamide" = "#fa9734",
    "Vincristine/Methotrexate/Cyclophosphamide" = "#f781bf",
    "No therapy" = "grey", "Others" = "#6c6c6c"
  )
  # At Xk the exposure label is the regimen administered in cycle k-1.
  pbmc$exposed_regimen <- vapply(seq_len(nrow(pbmc)), function(i) {
    if (pbmc$cycle_num[i] == 1L) return("baseline")
    patient_row <- match(pbmc$Patient[i], as.character(soc$Patient))
    regimen_column <- paste0("base_regimen_cycle", pbmc$cycle_num[i] - 1L)
    if (is.na(patient_row) || !regimen_column %in% names(soc)) return(NA_character_)
    as.character(soc[[regimen_column]][patient_row])
  }, character(1))

  feature_columns <- c("ATC", "Naïve%", "SCM%", "CM%", "EM%", "TE%",
                       "PD1%", "LAG3%", "TIM3%", "observed_Shannon", "cfShannon")
  if (!all(feature_columns %in% names(pbmc))) stop("Missing immune-feature columns")
  raw <- as.matrix(pbmc[, feature_columns, drop = FALSE])
  storage.mode(raw) <- "double"
  raw[!is.finite(raw)] <- NA_real_
  rownames(raw) <- as.character(pbmc$sample_id)
  colnames(raw)[10:11] <- c("TCR Diversity", "cfTCR Diversity")
  n_missing <- rowSums(is.na(raw))
  keep <- n_missing < 8L & rowSums(raw, na.rm = TRUE) != 0
  eligibility <- data.frame(
    sample_id = pbmc$sample_id, Patient = pbmc$Patient, cycle = pbmc$cycle,
    n_missing_features = n_missing, included = keep,
    reason = ifelse(keep, "included", ifelse(n_missing >= 8L,
      "eight_or_more_missing_features", "all_zero_observed_features")),
    stringsAsFactors = FALSE
  )
  raw <- raw[keep, , drop = FALSE]
  md <- pbmc[keep, , drop = FALSE]
  if (anyNA(md$exposed_regimen) || !all(md$exposed_regimen %in% names(therapy_palette)))
    stop("Unmapped treatment exposure in eligible samples")
  if (anyDuplicated(paste(md$Patient, md$cycle)))
    stop("Expected one eligible row per patient and timepoint")
  if (any(colSums(!is.na(raw)) < 2L)) stop("Too few observed values for an immune feature")

  # The manuscript fixes rank 2, selected in the original repeated masked-entry
  # cross-validation. The optional notebook step reruns that selection in full.
  rank_cv <- NULL
  if (rerun_rank_cv) {
    set.seed(seed)
    cv_input <- sweep(sweep(raw, 2L, colMeans(raw, na.rm = TRUE), "-"),
                      2L, apply(raw, 2L, stats::sd, na.rm = TRUE), "/")
    rank_cv <- missMDA::estim_ncpPCA(
      as.data.frame(cv_input), ncp.min = 0L, ncp.max = 8L, scale = TRUE,
      method = "Regularized", method.cv = "Kfold", nbsim = 20L,
      pNA = 0.10, threshold = 1e-6, verbose = FALSE
    )
    ncp <- as.integer(rank_cv$ncp)
  }
  if (length(ncp) != 1L || !is.finite(ncp) || ncp < 1L || ncp >= ncol(raw))
    stop("Invalid reconstruction rank")
  set.seed(seed)
  imputation <- missMDA::imputePCA(
    as.data.frame(raw), ncp = ncp, scale = TRUE, method = "Regularized",
    seed = seed, nb.init = 5L, threshold = 1e-6, maxiter = 1000L
  )
  completed <- as.matrix(imputation$completeObs)
  dimnames(completed) <- dimnames(raw)
  observed <- !is.na(raw)
  if (any(!is.finite(completed)) || !isTRUE(all.equal(
      completed[observed], raw[observed], tolerance = 0, check.attributes = FALSE)))
    stop("Imputation changed observed entries or left nonfinite values")

  complete_index <- which(rowSums(is.na(raw)) == 0L)
  real_features <- raw[complete_index, , drop = FALSE]
  original_pca <- stats::prcomp(real_features, center = TRUE, scale. = TRUE)
  reference <- stats::prcomp(completed, center = TRUE, scale. = TRUE)
  # PCA signs are arbitrary. Match the complete-case reference calculated above,
  # exactly as in the original analysis; this changes signs, never axis geometry.
  signs <- sign(colSums(reference$rotation * original_pca$rotation))
  signs[signs == 0] <- 1
  reference$rotation <- sweep(reference$rotation, 2L, signs, "*")
  reference$x <- sweep(reference$x, 2L, signs, "*")
  # Deterministic anchors keep the two manuscript axes stable across SVD builds.
  anchors <- c("TE%", "Naïve%")
  for (j in seq_len(2L)) if (reference$rotation[anchors[j], j] < 0) {
    reference$rotation[, j] <- -reference$rotation[, j]
    reference$x[, j] <- -reference$x[, j]
  }

  md$n_missing_original <- rowSums(is.na(raw))
  md$n_missing_features <- md$n_missing_original
  md$any_feature_imputed <- md$n_missing_original > 0L
  md$Age <- as.numeric(as.character(md$Age))
  md$regimen <- as.character(md$exposed_regimen)
  for (field in c("Remission", "Relapse", "vital_status")) {
    md[[field]] <- as.integer(as.character(md[[field]]))
    if (anyNA(md[[field]]) || !all(md[[field]] %in% c(0L, 1L)))
      stop("Expected binary metadata field: ", field)
  }
  if (any(!is.finite(md$Age))) stop("Missing age in eligible samples")
  md$lymphoma <- ifelse(as.character(md$Disease_type) %in%
    c("ALCL", "BL", "BLL", "DLBCL", "HD", "PMBCL"), "Lymphoma", "Non-lymphoma")
  md$recorded_lymphoma <- ifelse(as.character(md$cancergroup) == "Lymphoma",
                                "Lymphoma", "Other recorded groups")
  scores <- stats::predict(reference, newdata = completed)[, 1:2, drop = FALSE]
  md$PC1 <- scores[, 1L]
  md$PC2 <- scores[, 2L]
  real_md <- md[complete_index, , drop = FALSE]
  real_scores <- stats::predict(reference, newdata = real_features)[, 1:2, drop = FALSE]
  if (max(abs(real_scores - as.matrix(real_md[, c("PC1", "PC2")]))) > 1e-10)
    stop("Shared samples do not project to identical positions")
  for (patient in unique(md$Patient)) {
    rows <- md$Patient == patient
    for (field in c("Age", "Disease_type", "cancergroup", "Remission", "Relapse", "vital_status"))
      if (length(unique(md[rows, field])) != 1L)
        stop("Inconsistent patient-level metadata: ", patient, "/", field)
  }
  original_scores <- real_md
  original_scores$dataset <- "Original"
  completed_scores <- md
  completed_scores$dataset <- "missMDA"
  scoredata <- rbind(original_scores, completed_scores)

  disease_order <- c("ALL", "AML", "CML", "HR ALL", "BL", "DLBCL", "HD", "PMBCL",
    "ARMS", "ERMS", "EWS", "HB", "NB", "OS", "WILMS", "ALCL", "T-ALL", "RMS", "BLL", "MS NB")
  palette_env <- new.env(parent = baseenv())
  sys.source(file.path(project_root, "R_updated", "functions", "color_schemes.R"), envir = palette_env)
  if (!all(as.character(md$Disease_type) %in% disease_order)) stop("Unmapped cancer type")
  settings <- list(ncp = ncp, seed = seed, method = "Regularized", scale = TRUE,
    nb.init = 5L, threshold = 1e-6, maxiter = 1000L,
    exclusion_missing_features = 8L, rank_cv_rerun = rerun_rank_cv,
    rank_cv_method = "Kfold", rank_cv_nbsim = 20L, rank_cv_pNA = 0.10,
    rank_cv_range = 0:8, shared_axes = "missMDA-completed cohort PCA")
  versions <- data.frame(package = c("R", "missMDA", "FactoMineR"), version = c(
    as.character(getRversion()), as.character(utils::packageVersion("missMDA")),
    as.character(utils::packageVersion("FactoMineR"))), stringsAsFactors = FALSE)
  if (!identical(source_hashes, tools::md5sum(source_paths))) stop("Upstream inputs changed during analysis")
  list(project = project_root, metadata235 = md, real_metadata75 = real_md,
       completed235 = completed, raw235 = raw, real_features75 = real_features,
       reference_pca = reference, old_real_pca = original_pca,
       scoredata310 = scoredata, disease_order = disease_order,
       cancer_palette = palette_env$cohort_col[disease_order],
       group_palette = palette_env$group_col, therapy_palette = therapy_palette,
       features = colnames(raw), variance_percent =
         100 * reference$sdev[1:2]^2 / sum(reference$sdev^2),
       source_paths = source_paths, source_hashes = source_hashes,
       sample_eligibility = eligibility, rank_cv = rank_cv,
       settings = settings, versions = versions)
}

write_pca_source_data <- function(inputs, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(x, name) utils::write.csv(x, file.path(output_dir, name),
    row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write_matrix <- function(x, name) write_table(data.frame(
    sample_id = rownames(x), x, check.names = FALSE), name)
  write_matrix(inputs$raw235, "raw_immune_features.csv")
  write_matrix(inputs$completed235, "completed_immune_features_missMDA.csv")
  write_table(inputs$scoredata310, "shared_axis_scores_metadata.csv")
  write_table(inputs$sample_eligibility, "sample_eligibility.csv")
  write_table(as.data.frame(table(dataset = inputs$scoredata310$dataset,
    cycle = inputs$scoredata310$cycle)), "cohort_timepoint_counts.csv")
  reference <- inputs$reference_pca
  write_table(data.frame(feature = rownames(reference$rotation),
    reference$rotation, check.names = FALSE), "PCA_loadings.csv")
  variance <- reference$sdev^2 / sum(reference$sdev^2)
  write_table(data.frame(PC = seq_along(variance), variance_fraction = variance,
    cumulative_variance = cumsum(variance)), "PCA_variance.csv")
  if (!is.null(inputs$rank_cv)) write_table(data.frame(ncp = 0:8,
    prediction_error = as.numeric(inputs$rank_cv$criterion)), "missMDA_rank_crossvalidation.csv")
  write_table(data.frame(source = names(inputs$source_paths),
    file = file.path("data", basename(inputs$source_paths)),
    md5 = unname(inputs$source_hashes)), "input_checksums.csv")
  write_table(inputs$versions, "software_versions.csv")
  writeLines(capture.output(dput(inputs$settings)), file.path(output_dir, "PCA_settings.R"))
  observed <- !is.na(inputs$raw235)
  audit <- data.frame(
    measure = c("eligible_samples", "eligible_patients", "complete_samples",
      "complete_patients", "X01_completed_samples", "X01_complete_samples",
      "imputed_cells", "observed_entry_max_abs_change", "PC1_variance_percent", "PC2_variance_percent"),
    value = c(nrow(inputs$raw235), length(unique(inputs$metadata235$Patient)),
      nrow(inputs$real_features75), length(unique(inputs$real_metadata75$Patient)),
      sum(inputs$metadata235$cycle == "X01"), sum(inputs$real_metadata75$cycle == "X01"),
      sum(!observed), max(abs(inputs$completed235[observed] - inputs$raw235[observed])),
      inputs$variance_percent))
  write_table(audit, "PCA_preparation_audit.csv")
  invisible(audit)
}
