#!/usr/bin/env Rscript

# Title: Figure 4 - GLIPH TCR specificity group analysis
#
# Inputs: 
#   - `data/metadata.rds`
#   - `data/soc_good_with_baseregimen.rds`
#   - `data/gliph_tigerdb_annot_stat.rds`
#   - `data/int_gliph_tigerdb_cleaned.rds`
#   - `data/int_samples_gliph_tigerdb_annot_freq_stat.rds`
#   - Color and theme helpers from `R/functions/color_schemes.R` and `R/functions/ggplot2_theme.R`
# Outputs:
#   - `plots/Fig4A.pdf` - Volcano-style plot of TCR specificity-group clonal change vs Fisher score
#   - `plots/Fig4B.pdf` - Paired antigen fraction boxplots per treatment regimen (2 x 2 grid)
#   - PNG/SVG copies of each panel and a combined Figure4 overview
# Usage:
#   Rscript R/Figure4.R

# --- Setup --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggsignif)
  library(ggrepel)
  library(cowplot)
  library(tidyselect)
})

Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

# Load helper palettes and themes
source(file.path(getwd(), "R/functions", "color_schemes.R"))
source(file.path(getwd(), "R/functions", "ggplot2_theme.R"))
source(file.path(getwd(), "R/functions", "Misc_functions.R"))
source(file.path(getwd(), "R/functions", "figure_exports.R"))

# Project-rooted paths
project_root <- getwd()
plots_dir <- file.path(project_root, "plots")
datapath <- file.path(project_root, "data/")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# --- Palettes -----------------------------------------------------------------

colpal <- c(
  "Vincristine" = "#E41A1C",
  "Methotrexate" = "#377EB8",
  "Cyclophosphamide" = "#4DAF4A",
  "Vincristine/Methotrexate" = "#984EA3",
  "Vincristine/Cyclophosphamide" = "#fa9734",
  "Vincristine/Methotrexate/Cyclophosphamide" = "#f781bf",
  "No therapy" = "grey",
  "Others" = "#6c6c6c"
)

# --- Helpers ------------------------------------------------------------------

to_num <- function(x) as.integer(gsub("\\D", "", x))

# Standardize a cycle-specific pattern data frame for row-binding.
standardize_cycle_df <- function(df, cycle_label) {
  out <- df %>% select(Sample, pattern, Patient, MeanFreq, exposed_regimen)
  out$cycle <- cycle_label
  out
}

# Build paired-antigen plot data + wilcoxon for a target regimen
make_regimen_paired_antigen_plot <- function(data, target_regimen,
                                             antigens = c("Human", "CMV", "EBV", "HCV", "HPV",
                                                          "HSV", "MCPyV", "Influenza", "YFV", "CEF"),
                                             p_alpha = 0.15, adjust_method = "BH") {

  # Create _frac columns if missing
  if (!all(paste0(antigens, "_frac") %in% names(data)) && "cdr3_reads" %in% names(data)) {
    for (ag in antigens) {
      if (ag %in% names(data) && !(paste0(ag, "_frac") %in% names(data))) {
        data[[paste0(ag, "_frac")]] <- data[[ag]] / data$cdr3_reads
      }
    }
  }

  df_reg <- data %>%
    dplyr::filter(.data$exposed_regimen == target_regimen) %>%
    dplyr::select(Sample = dplyr::all_of(if ("sample_id" %in% names(.)) "sample_id" else "Sample"),
                  Patient, cycle, exposed_regimen,
                  dplyr::all_of(paste0(antigens, "_frac"))) %>%
    dplyr::mutate(cycle_num = to_num(.data$cycle))

  if (nrow(df_reg) == 0) {
    return(list(pairs_long = tibble::tibble(), wilcox_results = tibble::tibble(), plot = ggplot2::ggplot()))
  }

  from_tbl <- df_reg %>%
    dplyr::transmute(Patient, exposed_regimen,
                     cycle_from = .data$cycle, cycle_num_from = .data$cycle_num,
                     dplyr::across(tidyselect::ends_with("_frac"), identity, .names = "from.{.col}"))

  to_tbl <- df_reg %>%
    dplyr::transmute(Patient, exposed_regimen,
                     cycle_to = .data$cycle, cycle_num_to = .data$cycle_num,
                     dplyr::across(tidyselect::ends_with("_frac"), identity, .names = "to.{.col}"))

  pairs <- try(
    from_tbl %>%
      dplyr::inner_join(to_tbl, by = c("Patient", "exposed_regimen"), relationship = "many-to-many") %>%
      dplyr::filter(.data$cycle_num_to == .data$cycle_num_from + 1L),
    silent = TRUE
  )
  if (inherits(pairs, "try-error")) {
    pairs <- from_tbl %>%
      dplyr::inner_join(to_tbl, by = c("Patient", "exposed_regimen")) %>%
      dplyr::filter(.data$cycle_num_to == .data$cycle_num_from + 1L)
  }
  if (nrow(pairs) == 0) {
    return(list(pairs_long = tibble::tibble(), wilcox_results = tibble::tibble(), plot = ggplot2::ggplot()))
  }

  pairs_long_from <- pairs %>%
    dplyr::select(Patient, exposed_regimen, cycle_from, cycle_to, tidyselect::starts_with("from.")) %>%
    tidyr::pivot_longer(tidyselect::starts_with("from."), names_to = "var", values_to = "fraction") %>%
    dplyr::mutate(antigen = sub("^from\\.(.+)_frac$", "\\1", .data$var),
                  timepoint = "from",
                  cycle = .data$cycle_from) %>%
    dplyr::select(Patient, exposed_regimen, interval = .data$cycle_to, cycle, antigen, timepoint, fraction)

  pairs_long_to <- pairs %>%
    dplyr::select(Patient, exposed_regimen, cycle_from, cycle_to, tidyselect::starts_with("to.")) %>%
    tidyr::pivot_longer(tidyselect::starts_with("to."), names_to = "var", values_to = "fraction") %>%
    dplyr::mutate(antigen = sub("^to\\.(.+)_frac$", "\\1", .data$var),
                  timepoint = "to",
                  cycle = .data$cycle_to) %>%
    dplyr::select(Patient, exposed_regimen, interval = .data$cycle_to, cycle, antigen, timepoint, fraction)

  pairs_long <- dplyr::bind_rows(pairs_long_from, pairs_long_to) %>%
    dplyr::filter(.data$antigen %in% antigens) %>%
    dplyr::mutate(
      antigen = factor(.data$antigen, levels = antigens),
      pair_id = paste(.data$Patient, .data$interval, sep = "_"),
      cycle_order = to_num(.data$cycle)
    )

  wilcox_results <- pairs_long %>%
    dplyr::select(Patient, antigen, interval, timepoint, fraction) %>%
    dplyr::filter(!is.na(.data$fraction)) %>%
    tidyr::pivot_wider(names_from = timepoint, values_from = fraction) %>%
    dplyr::group_by(antigen) %>%
    dplyr::summarise(
      n_pairs = sum(!is.na(.data$from) & !is.na(.data$to)),
      p = if (n_pairs >= 2)
        tryCatch(stats::wilcox.test(.data$to, .data$from, paired = TRUE)$p.value,
                 error = function(e) NA_real_) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_adj = p.adjust(.data$p, method = adjust_method)) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(pairs_long,
                       ggplot2::aes(x = factor(.data$timepoint, levels = c("from", "to")),
                                    y = .data$fraction)) +
    ggplot2::geom_line(ggplot2::aes(color = .data$exposed_regimen, group = .data$pair_id), alpha = 0.7) +
    ggplot2::geom_point(alpha = 0.7, shape = 16) +
    ggplot2::geom_boxplot(ggplot2::aes(group = .data$timepoint), inherit.aes = TRUE,
                          outlier.shape = NA, alpha = 0.25) +
    ggplot2::facet_grid(. ~ antigen, scales = "free_y",
                        labeller = ggplot2::labeller(antigen = function(x) ifelse(x == "Influenza", "Flu", x))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = "", y = "Fraction of TCR reads", color = "", title = target_regimen) +
    ggplot2::scale_color_manual(values = colpal) +
    ggplot2::scale_x_discrete(labels = c("from" = "pre", "to" = "post"))

  annot_sig <- wilcox_results %>%
    dplyr::filter(!is.na(.data$p_adj) & .data$p_adj < p_alpha) %>%
    dplyr::mutate(label = dplyr::if_else(.data$p_adj < 0.01,
                                         sprintf("p = %.2e", .data$p_adj),
                                         sprintf("p = %.3f", .data$p_adj)))

  if (nrow(annot_sig) > 0) {
    y_map <- sapply(levels(pairs_long$antigen), function(a) {
      vals <- pairs_long$fraction[pairs_long$antigen == a]
      if (all(is.na(vals))) return(NA_real_)
      max(vals, na.rm = TRUE) * 1.08
    })
    signif_df <- data.frame(
      antigen = annot_sig$antigen,
      y_position = y_map[as.character(annot_sig$antigen)],
      label = annot_sig$label
    )
    signif_df <- signif_df[is.finite(signif_df$y_position), , drop = FALSE]
    if (nrow(signif_df) > 0) {
      p <- p + ggsignif::geom_signif(
        data = signif_df,
        ggplot2::aes(xmin = 1, xmax = 2, annotations = .data$label, y_position = .data$y_position),
        manual = TRUE, textsize = 3.2, tip_length = 0.015
      )
    }
  }

  list(pairs_long = pairs_long, wilcox_results = wilcox_results, plot = p)
}

# --- Load inputs --------------------------------------------------------------

meta_div <- readr::read_rds(paste0(datapath, "metadata.rds"))
soc_good <- readr::read_rds(paste0(datapath, "soc_good_with_baseregimen.rds"))

soc_good_long <- soc_good %>%
  pivot_longer(
    cols = starts_with("base_regimen_"),
    names_to = "cycle",
    values_to = "regimen",
    names_prefix = "base_regimen_"
  ) %>%
  select(Patient, cycle, regimen)

meta_div$cycle_mapped <- paste0("cycle", gsub("X0?", "", meta_div$cycle))
meta_div <- meta_div %>%
  left_join(soc_good_long, by = c("Patient" = "Patient", "cycle_mapped" = "cycle"))
meta_div$cycle_mapped <- NULL

# --- Derive exposed_regimen (previous cycle's regimen) ------------------------

meta_div$exposed_regimen <- NA
for (i in seq_len(nrow(meta_div))) {
  current_patient <- meta_div$Patient[i]
  current_cycle <- meta_div$cycle[i]
  if (current_cycle == "X01") {
    meta_div$exposed_regimen[i] <- "baseline"
  } else {
    cycle_num <- as.numeric(gsub("X0?", "", current_cycle))
    prev_cycle_col <- paste0("cycle", cycle_num - 1)
    prev_regimen <- as.character(soc_good[soc_good$Patient == current_patient,
                                          paste0("base_regimen_", prev_cycle_col), drop = TRUE])
    if (length(prev_regimen) > 0 && !is.na(prev_regimen)) {
      meta_div$exposed_regimen[i] <- prev_regimen
    }
  }
}

# Restrict to samples with TCRseq data
meta_div <- meta_div[!is.na(meta_div$TRB), ]

# --- Load GLIPH outputs -------------------------------------------------------

gliphdf   <- readr::read_rds(paste0(datapath, "gliph_tigerdb_annot_stat.rds"))
gliph_out <- readr::read_rds(paste0(datapath, "int_gliph_tigerdb_cleaned.rds"))

chp_gliph <- gliph_out[grepl("CHP", gliph_out$Sample), ]
chp_gliph_pbmc <- chp_gliph[grepl("PBMC", chp_gliph$Sample), ]

chp_gliph_pbmc$Patient         <- meta_div$Patient[match(chp_gliph_pbmc$Sample, meta_div$basename)]
chp_gliph_pbmc$cycle           <- meta_div$cycle[match(chp_gliph_pbmc$Sample, meta_div$basename)]
chp_gliph_pbmc$sampletype      <- meta_div$sampletype[match(chp_gliph_pbmc$Sample, meta_div$basename)]
chp_gliph_pbmc$cancergroup     <- meta_div$cancergroup[match(chp_gliph_pbmc$Sample, meta_div$basename)]
chp_gliph_pbmc$regimen         <- meta_div$regimen[match(chp_gliph_pbmc$Sample, meta_div$basename)]
chp_gliph_pbmc$exposed_regimen <- meta_div$exposed_regimen[match(chp_gliph_pbmc$Sample, meta_div$basename)]

# Patients with at least one pair of consecutive timepoints
mytab <- table(chp_gliph_pbmc$Patient, chp_gliph_pbmc$cycle)
consecutive_matrix <- mytab > 0
consecutive_matrix <- consecutive_matrix[apply(consecutive_matrix, 1,
                                               function(x) any(x[-length(x)] & x[-1])), , drop = FALSE]
baseline_patients <- rownames(consecutive_matrix)

df1 <- chp_gliph_pbmc[chp_gliph_pbmc$Patient %in% baseline_patients, ]

# --- Fig 4A: TCR specificity-group clonal change vs Fisher score -------------

# Average per (Sample, pattern) at each cycle
cycles <- c("X01", "X02", "X03", "X04", "X05")
per_cycle <- lapply(cycles, function(cc) {
  df1 %>%
    dplyr::filter(cycle == cc) %>%
    dplyr::group_by(Sample, pattern) %>%
    dplyr::summarise(MeanFreq = mean(Freq), .groups = "drop") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      exposed_regimen = meta_div$exposed_regimen[match(Sample, meta_div$basename)],
      cancergroup     = meta_div$cancergroup[match(Sample, meta_div$basename)],
      Disease_type    = meta_div$Disease_type[match(Sample, meta_div$basename)],
      Patient         = meta_div$Patient[match(Sample, meta_div$basename)]
    )
})
names(per_cycle) <- cycles

df_all <- do.call(rbind, lapply(names(per_cycle),
                                function(nm) standardize_cycle_df(per_cycle[[nm]], nm)))

sample_pat_map <- unique(df_all[, c("Sample", "Patient", "cycle", "exposed_regimen")])

df_pc2 <- df_all %>% mutate(cycle_num = to_num(cycle))

from <- df_pc2 %>%
  transmute(Patient, pattern,
            cycle_from = cycle, cycle_num_from = cycle_num,
            MeanFreq_from = MeanFreq,
            next_cycle_num = cycle_num + 1L)

to <- df_pc2 %>%
  transmute(Patient, pattern,
            cycle_to = cycle, cycle_num_to = cycle_num,
            MeanFreq_to = MeanFreq)

delta_long <- from %>%
  inner_join(to, by = c("Patient", "pattern", "next_cycle_num" = "cycle_num_to")) %>%
  transmute(Patient, pattern,
            cycle_from, cycle_to, next_cycle_num,
            MeanFreq_from, MeanFreq_to,
            Delta = MeanFreq_to - MeanFreq_from)

# Fisher score per pattern
pattern_fisher <- df1 %>%
  distinct(pattern, Fisher_score) %>%
  group_by(pattern) %>%
  summarise(fisher = suppressWarnings(min(Fisher_score, na.rm = TRUE)), .groups = "drop") %>%
  mutate(fisher = ifelse(is.infinite(fisher), NA_real_, fisher))

# Regimen per patient x cycle (for the "to" cycle)
regimen_pc <- sample_pat_map %>%
  select(Patient, cycle, exposed_regimen) %>%
  group_by(Patient, cycle) %>%
  summarise(exposed_regimen = dplyr::first(na.omit(exposed_regimen)), .groups = "drop")

# Antigen flags per pattern
ant_lookup <- gliphdf %>%
  select(gliph, CMV, Influenza, HCV, EBV, HPV, Human, CEF, HSV, TB, SP, MCPyV, YFV) %>%
  rename(Flu = Influenza, HS = Human)

delta_long <- delta_long %>%
  left_join(pattern_fisher, by = "pattern") %>%
  left_join(regimen_pc, by = c("Patient", "cycle_to" = "cycle")) %>%
  left_join(ant_lookup, by = c("pattern" = "gliph"))

viral_cols     <- c("CMV", "Flu", "HCV", "EBV", "HPV", "MCPyV", "HSV", "YFV")
bacterial_cols <- c("TB", "SP")

delta_long$antigen <- "Unknown"
delta_long$antigen[rowSums(delta_long[, viral_cols], na.rm = TRUE) > 0]     <- "Viral"
delta_long$antigen[rowSums(delta_long[, bacterial_cols], na.rm = TRUE) > 0] <- "Bacterial"
delta_long$antigen[delta_long$HS > 0] <- "HomoSapiens"
delta_long$antigen[
  (delta_long$HS > 0) &
    (rowSums(delta_long[, c(viral_cols, bacterial_cols)], na.rm = TRUE) > 0)
] <- "Cross-reactive"

# Unique patient count per pattern
pattern_patient_count <- delta_long %>%
  dplyr::distinct(Patient, pattern) %>%
  dplyr::count(pattern, name = "n_patients_with_pattern")
delta_long <- delta_long %>%
  dplyr::left_join(pattern_patient_count, by = "pattern")

delta_long$exposed_regimen <- factor(delta_long$exposed_regimen,
                                     levels = c("No therapy", setdiff(names(colpal), "No therapy")))

quantiles <- quantile(delta_long$Delta, c(0.1, 0.9), na.rm = TRUE)

p_fig4a <- ggplot() +
  geom_point(data = delta_long,
             aes(x = Delta, y = -log10(fisher + 1e-25),
                 color = exposed_regimen, shape = antigen, size = n_patients_with_pattern),
             alpha = 0.8) +
  myplot + myaxis +
  geom_vline(xintercept = c(quantiles[1], quantiles[2]), linetype = 3) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_text(size = 10)) +
  labs(y = "-log10(Fisher score)",
       x = "TCR specificity group clonal change \n[post-therapy - pre-therapy]") +
  guides(color = guide_legend(title = "exposed_regimen"),
         shape = guide_legend(title = "Antigen"),
         size  = guide_legend(title = "Unique CDR3s (count)")) +
  scale_color_manual(values = colpal) +
  scale_size_continuous(range = c(3, 6))

# Retain the labeled panel for the shared export block below.
panel_fig4a <- cowplot::plot_grid(p_fig4a, labels = "A")

# --- Fig 4B: Paired antigen fraction boxplots per treatment regimen ----------

samples_tcr_freq <- readr::read_rds(paste0(datapath, "int_samples_gliph_tigerdb_annot_freq_stat.rds"))
samples_tcr_freq$sample_id       <- rownames(samples_tcr_freq)
samples_tcr_freq$Relapse         <- meta_div$Relapse[match(samples_tcr_freq$sample_id, meta_div$basename)]
samples_tcr_freq$exposed_regimen <- meta_div$exposed_regimen[match(samples_tcr_freq$sample_id, meta_div$basename)]

df_pbmc <- samples_tcr_freq[samples_tcr_freq$sampletype == "PBMC", ]

antigens_used <- c("Human", "CMV", "EBV", "HCV", "HPV", "HSV", "MCPyV", "Influenza")

vin       <- make_regimen_paired_antigen_plot(df_pbmc, target_regimen = "Vincristine",                  antigens = antigens_used)
mtx       <- make_regimen_paired_antigen_plot(df_pbmc, target_regimen = "Methotrexate",                 antigens = antigens_used)
vincyclo  <- make_regimen_paired_antigen_plot(df_pbmc, target_regimen = "Vincristine/Cyclophosphamide", antigens = antigens_used)
others    <- make_regimen_paired_antigen_plot(df_pbmc, target_regimen = "Others",                       antigens = antigens_used)

combined_boxplots <- cowplot::plot_grid(
  vin$plot, mtx$plot, vincyclo$plot, others$plot,
  ncol = 2, nrow = 2, align = "v", axis = "tblr", labels = "B"
)

# Retain panel A and each treatment subplot's native dimensions and typography.
# Panel B wraps four treatment subplots into two rows of two.
export_figure_panels(
  panels = list(panel_fig4a, combined_boxplots),
  dimensions = list(c(10, 5.6), c(8.5, 9.6)),
  plots_dir = plots_dir, panel_stems = c("Fig4A", "Fig4B"),
  figure_stem = "Figure4")
message("Figure 4 panels and overview saved in ", plots_dir)

