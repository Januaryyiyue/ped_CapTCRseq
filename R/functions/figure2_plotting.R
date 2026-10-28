# Plot construction for Figure 2 and Extended Figure 2.
# Inputs are numeric results returned by prepare_pca_inputs() and
# fit_pca_associations(). No previously rendered plot is read by this file.

figure2_theme <- function() {
  # The same theme components are used by Figure 1 and Figure 4. Horizontal
  # axis labels are appropriate for the continuous PCA and raincloud axes.
  theme_minimal(base_size = 13, base_family = "sans") + myplot + myaxis +
    theme(axis.text.x = element_text(angle = 0, hjust = .5),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 13, colour = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 13, hjust = .5),
          plot.tag = element_text(size = 16, face = "bold"),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin = margin(7, 7, 7, 7))
}

figure2_palettes <- function(inputs) {
  # Darker categorical colours retain readability of the 20 cancer types.
  disease_order <- c("ALL", "AML", "CML", "HR ALL", "BL", "DLBCL", "HD", "PMBCL",
    "ARMS", "ERMS", "EWS", "HB", "NB", "OS", "WILMS", "ALCL", "T-ALL", "RMS", "BLL", "MS NB")
  cancer <- setNames(c("#1A6297", "#479506", "#CC681C", "#DF45B0", "#423B26",
    "#632643", "#3C918B", "#868194", "#3D219A", "#928355", "#6B240C", "#D95A61",
    "#8572E9", "#326941", "#1C4144", "#393A4E", "#705956", "#845405", "#A403B5", "#258EB5"),
    disease_order)
  list(cancer = cancer, therapy = inputs$therapy_palette,
       shape = c(Original = 16, missMDA = 1))
}

pca_loading_arrows <- function(reference_pca) {
  # Correlation-vector coordinates, expanded to 70% of the shared score range.
  # This is the biplot scaling used for the corresponding individual PCA plots.
  xy <- sweep(reference_pca$rotation[, 1:2, drop = FALSE], 2,
              reference_pca$sdev[1:2], "*")
  score_range <- apply(reference_pca$x[, 1:2, drop = FALSE], 2, function(v) diff(range(v)))
  vector_range <- apply(xy, 2, function(v) diff(range(v)))
  multiplier <- .7 * min(score_range / vector_range)
  data.frame(name = rownames(xy), x = xy[, 1] * multiplier,
             y = xy[, 2] * multiplier, row.names = NULL, check.names = FALSE)
}

figure2_score_data <- function(inputs) {
  d <- inputs$scoredata310
  stopifnot(all(c("sample_id", "Patient", "cycle", "PC1", "PC2", "Disease_type",
    "exposed_regimen", "dataset") %in% names(d)))
  d$cohort <- factor(d$dataset, levels = c("Original", "missMDA"))
  if (!"cycle_num" %in% names(d)) d$cycle_num <- as.integer(sub("^X0*", "", d$cycle))
  stopifnot(!anyNA(d$cohort), !anyNA(d$cycle_num), !anyNA(d$PC1), !anyNA(d$PC2),
    !anyDuplicated(paste(d$dataset, d$sample_id)))
  d
}

make_pca_biplot <- function(data, inputs, colour_column = "Disease_type",
                           y_limit = c(-4, 5), show_cohort_legend = TRUE) {
  palettes <- figure2_palettes(inputs)
  palette <- if (colour_column == "Disease_type") palettes$cancer else palettes$therapy
  title <- if (colour_column == "Disease_type") "Cancer type" else "Regimen"
  arrows <- pca_loading_arrows(inputs$reference_pca)
  variance <- 100 * inputs$reference_pca$sdev^2 / sum(inputs$reference_pca$sdev^2)
  # The cohorts overlap: draw the filled Original marks after the open marks.
  data <- data[order(data$cohort, decreasing = TRUE), , drop = FALSE]
  stopifnot(all(data[[colour_column]] %in% names(palette)),
    all(data$PC1 > -7 & data$PC1 < 4), all(data$PC2 > y_limit[1] & data$PC2 < y_limit[2]))
  p <- ggplot(data, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = .5, colour = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = .5, colour = "black") +
    geom_point(aes(colour = .data[[colour_column]], shape = cohort), size = 3.2,
      stroke = .8, alpha = if (colour_column == "Disease_type") .85 else .6,
      show.legend = TRUE) +
    geom_segment(data = arrows, aes(x = 0, y = 0, xend = x, yend = y),
      inherit.aes = FALSE, linewidth = .5, colour = "black",
      arrow = grid::arrow(length = grid::unit(.2, "cm"), type = "open")) +
    ggrepel::geom_text_repel(data = arrows, aes(x, y, label = name),
      inherit.aes = FALSE, colour = "black", size = 13 / .pt, family = "sans",
      seed = 123, max.overlaps = Inf, box.padding = .25,
      point.padding = .000001, segment.size = .2) +
    scale_colour_manual(values = palette, limits = names(palette), breaks = names(palette),
      drop = FALSE, name = title) +
    scale_shape_manual(values = palettes$shape, limits = names(palettes$shape),
      drop = FALSE, name = "Cohort") +
    scale_x_continuous(limits = c(-7, 4), breaks = c(-6, -4, -2, 0, 2, 4)) +
    scale_y_continuous(limits = y_limit,
      breaks = if (y_limit[2] == 5) c(-4, -2, 0, 2, 4, 5) else c(-4, -2, 0, 2, 4, 6),
      expand = if (y_limit[2] == 5) c(0, 0) else waiver()) +
    labs(x = sprintf("Inverted Immune Checkpoint axis\n(PC1 = %.1f%%)", variance[1]),
         y = sprintf("Diversity axis\n(PC2 = %.1f%%)", variance[2])) +
    figure2_theme() +
    guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, alpha = 1)),
           shape = guide_legend(order = 2, override.aes = list(colour = "black", alpha = 1)))
  if (colour_column == "Disease_type") {
    p <- p + geom_point(data = data[data$cohort == "Original", , drop = FALSE],
      aes(PC1, PC2), inherit.aes = FALSE, shape = 21, size = 3.2,
      stroke = .35, colour = "grey40", fill = NA, alpha = .65, show.legend = FALSE)
  }
  if (!show_cohort_legend) p <- p + guides(shape = "none")
  p
}

polygon_hatch_segments <- function(polygon, spacing) {
  # Vertical vector segments clipped to a polygon (PC2 rotates them horizontally).
  ans <- list(); n <- nrow(polygon); j <- c(2:n, 1L)
  levels <- seq(ceiling(min(polygon$x) / spacing) * spacing, max(polygon$x), by = spacing)
  for (level in levels) {
    hits <- which((polygon$x <= level & polygon$x[j] > level) |
                    (polygon$x[j] <= level & polygon$x > level))
    if (length(hits) < 2L) next
    values <- sort(polygon$y[hits] + (level - polygon$x[hits]) *
      (polygon$y[j[hits]] - polygon$y[hits]) / (polygon$x[j[hits]] - polygon$x[hits]))
    stopifnot(length(values) %% 2L == 0L)
    for (k in seq.int(1L, length(values), by = 2L)) {
      if (values[k + 1L] - values[k] < 1e-8) next
      ans[[length(ans) + 1L]] <- data.frame(x = level, y = values[k],
        xend = level, yend = values[k + 1L])
    }
  }
  if (length(ans)) do.call(rbind, ans) else
    data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
}

make_pca_marginal_data <- function(data, therapy_palette) {
  curves <- polygons <- hatches <- status <- list()
  for (pc in c("PC1", "PC2")) for (regimen in names(therapy_palette)) {
    position <- match(regimen, names(therapy_palette))
    values <- lapply(c("Original", "missMDA"), function(co)
      data[data$cohort == co & data$exposed_regimen == regimen, pc])
    names(values) <- c("Original", "missMDA")
    full <- values$missMDA
    bandwidth <- if (length(full) > 1L && sd(full) > 0) bw.nrd0(full) else NA_real_
    densities <- lapply(values, function(v) {
      if (length(v) > 1L && sd(v) > 0 && is.finite(bandwidth))
        density(v, bw = bandwidth, from = if (pc == "PC1") -7 else -4,
                to = if (pc == "PC1") 4 else 6, n = 512) else NULL
    })
    peak <- max(c(0, unlist(lapply(densities, function(v) if (!is.null(v)) v$y else 0))))
    for (cohort in names(values)) {
      values_cohort <- values[[cohort]]; de <- densities[[cohort]]
      status[[length(status) + 1L]] <- data.frame(PC = pc, regimen = regimen,
        cohort = cohort, n_samples = length(values_cohort), bandwidth = bandwidth,
        density_drawn = !is.null(de))
      if (is.null(de)) next
      direction <- if (cohort == "Original") 1 else -1
      y <- position + direction * .43 * de$y / peak
      key <- paste(pc, position, cohort, sep = ":")
      curves[[key]] <- data.frame(PC = pc, regimen = regimen, cohort = cohort,
        coordinate = de$x, density = de$y, ridge_coordinate = y,
        regimen_position = position, bandwidth = bandwidth)
      poly <- data.frame(x = c(de$x[1], de$x, tail(de$x, 1)), y = c(position, y, position))
      polygons[[key]] <- data.frame(poly, PC = pc, regimen = regimen, cohort = cohort, group = key)
      if (cohort == "missMDA") {
        hs <- polygon_hatch_segments(poly, .24)
        if (nrow(hs)) hatches[[key]] <- data.frame(hs, PC = pc, regimen = regimen, group = key)
      }
    }
  }
  list(curves = do.call(rbind, curves), polygons = do.call(rbind, polygons),
       hatches = do.call(rbind, hatches), status = do.call(rbind, status))
}

make_pca_marginal <- function(pc, data, marginal_data, palette) {
  polygons <- marginal_data$polygons[marginal_data$polygons$PC == pc, ]
  data$value <- data[[pc]]
  data$regimen_position <- match(data$exposed_regimen, names(palette)) +
    ifelse(data$cohort == "Original", .035, -.035)
  p <- ggplot() +
    geom_hline(yintercept = seq_along(palette), colour = "grey85", linewidth = .3) +
    # PC2 score coordinates rotate into horizontal tick-aligned background lines.
    (if (pc == "PC2") geom_vline(xintercept = c(-4, -2, 0, 2, 4, 6),
      colour = "grey85", linewidth = .3)) +
    geom_polygon(data = polygons[polygons$cohort == "Original", ],
      aes(x, y, group = group, fill = regimen), colour = "black", linewidth = .35, alpha = .4) +
    geom_polygon(data = polygons[polygons$cohort == "missMDA", ],
      aes(x, y, group = group, colour = regimen), fill = NA, linewidth = .65) +
    geom_point(data = data, aes(value, regimen_position), shape = "|", size = 2.3,
      colour = "black", alpha = .6) +
    scale_fill_manual(values = palette) + scale_colour_manual(values = palette) +
    scale_x_continuous(limits = if (pc == "PC1") c(-7, 4) else c(-4, 6),
      breaks = if (pc == "PC1") c(-6, -4, -2, 0, 2, 4) else c(-4, -2, 0, 2, 4, 6)) +
    scale_y_continuous(limits = c(.5, length(palette) + .5), breaks = seq_along(palette)) +
    labs(x = NULL, y = NULL) + figure2_theme() +
    theme(panel.grid.major.x = element_line(colour = "grey85", linewidth = .3),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line = element_blank(), legend.position = "none", plot.margin = margin(0, 0, 0, 0)) +
    guides(fill = "none", colour = "none")
  if (pc == "PC2") p <- p + coord_flip() +
    theme(panel.grid.major.x = element_blank(), axis.text.x = element_blank(),
      axis.text.y = element_text(size = 13, colour = "black"))
  p
}

make_pca_marginal_legend <- function(palette) {
  data <- data.frame(regimen = names(palette), y = seq(15.5, 7.5, by = -1))
  ggplot() +
    annotate("text", x = .045, y = 16.65, label = "Regimen", hjust = 0, size = 10 / .pt) +
    geom_point(data = data, aes(x = .09, y = y, colour = regimen), size = 3.2, show.legend = FALSE) +
    geom_text(data = data, aes(x = .19, y = y, label = regimen), hjust = 0, size = 10 / .pt) +
    scale_colour_manual(values = palette, limits = names(palette), guide = "none") +
    annotate("text", x = .045, y = 5.8, label = "Cohort", hjust = 0, size = 10 / .pt) +
    annotate("point", x = .09, y = 4.7, shape = 16, size = 3.2) +
    annotate("text", x = .19, y = 4.7, label = "Original", hjust = 0, size = 10 / .pt) +
    annotate("point", x = .09, y = 3.7, shape = 1, size = 3.2, stroke = .8) +
    annotate("text", x = .19, y = 3.7, label = "missMDA", hjust = 0, size = 10 / .pt) +
    annotate("rect", xmin = .045, xmax = .14, ymin = 1.55, ymax = 2.05,
      fill = "grey75", colour = "black", linewidth = .3) +
    annotate("text", x = .19, y = 1.8, label = "Original: solid density", hjust = 0, size = 9.5 / .pt) +
    annotate("rect", xmin = .045, xmax = .14, ymin = .55, ymax = 1.05,
      fill = NA, colour = "grey40", linewidth = .3) +
    annotate("text", x = .19, y = .8, label = "missMDA: outlined density", hjust = 0, size = 9.5 / .pt) +
    coord_cartesian(xlim = c(0, 1.2), ylim = c(0, 17.2), expand = FALSE, clip = "off") +
    theme_void(base_family = "sans") + theme(plot.margin = margin(0, 0, 0, 0))
}

make_raincloud_data <- function(data, category, pc, levels, tag) {
  data$category <- factor(data[[category]], levels = levels)
  data$score <- data[[pc]]
  data$dataset <- factor(data$dataset, levels = c("Original", "missMDA"))
  set.seed(20261022L + match(tag, letters))
  data$rain_x <- as.integer(data$category) + ifelse(data$dataset == "Original", -.045, .045) +
    runif(nrow(data), -.022, .022)
  clouds <- summary <- list()
  for (k in seq_along(levels)) for (cohort in c("Original", "missMDA")) {
    values <- data$score[data$category == levels[k] & data$dataset == cohort]
    if (!length(values)) next
    side <- if (cohort == "Original") -1 else 1
    anchor <- k + side * .10
    if (length(values) >= 2 && sd(values) > 1e-10) {
      de <- density(values, n = 256, from = min(values), to = max(values))
      width <- .18 * de$y / max(de$y)
      clouds[[length(clouds) + 1L]] <- data.frame(
        x = c(rep(anchor, length(de$x)), rev(anchor + side * width)),
        y = c(de$x, rev(de$x)), category = levels[k], cohort = cohort,
        id = paste(k, cohort))
    }
    q <- quantile(values, c(.25, .5, .75), names = FALSE)
    summary[[length(summary) + 1L]] <- data.frame(x = k + side * .082,
      q25 = q[1], median = q[2], q75 = q[3], category = levels[k], cohort = cohort, n = length(values))
  }
  list(source = data, cloud_polygons = do.call(rbind, clouds), median_IQR = do.call(rbind, summary))
}

figure2_statistics_strip <- function(table, primary_note = NULL, paired_note = NULL) {
  format_p <- function(p) ifelse(is.na(p), "NE", format.pval(p, digits = 2, eps = .0001))
  x <- c(.13, .36, .58)
  headers <- c("Cohort", "Holm P", "Shift [95% Bonferroni CI]")
  text <- function(label, x, y, bold = FALSE, size = 10.5)
    grid::textGrob(label, x = x, y = y, hjust = 0,
      gp = grid::gpar(fontfamily = "sans", fontsize = size,
        fontface = if (bold) "bold" else "plain"))
  grobs <- list(grid::segmentsGrob(x0 = .125, x1 = .99, y0 = .98, y1 = .98,
    gp = grid::gpar(col = "grey65", lwd = .7)))
  grobs <- c(grobs, lapply(seq_along(x), function(j) text(headers[j], x[j], .86, TRUE)))
  for (i in seq_len(nrow(table))) {
    values <- c(table$dataset[i], format_p(table$primary_p_Holm[i]),
      sprintf("%.2f [%.2f, %.2f]", table$estimate[i], table$lower[i], table$upper[i]))
    grobs <- c(grobs, lapply(seq_along(x), function(j) text(values[j], x[j], .65 - (i - 1) * .21)))
  }
  notes <- c(primary_note, paired_note)
  for (i in seq_along(notes)) grobs <- c(grobs,
    list(text(notes[i], .13, .22 - (i - 1) * .17, size = 9.7)))
  patchwork::wrap_elements(full = do.call(grid::grobTree, grobs), clip = FALSE)
}

make_association_raincloud <- function(tag, data, inputs, statistics) {
  cohort_order <- c("Original", "missMDA")
  regimen <- switch(tag, d = NA_character_, e = "Vincristine/Cyclophosphamide",
                    f = "Vincristine", g = "Any cyclophosphamide")
  levels <- if (tag == "d") c("Non-lymphoma", "Lymphoma") else
    c("Baseline (X01)", switch(tag, e = "Vincristine /\ncyclophosphamide",
                              f = "Vincristine", g = "Cyclophosphamide-\ncontaining"))
  pc <- if (tag %in% c("f", "g")) "PC2" else "PC1"
  if (tag == "d") {
    d <- data[data$cycle_num == 1, , drop = FALSE]
    d$category_show <- d$lymphoma
    colours <- c("Non-lymphoma" = "#7570B3", "Lymphoma" = "#D95F02")
    r <- statistics$baseline_lymphoma[match(cohort_order, statistics$baseline_lymphoma$dataset), ]
    table <- data.frame(panel = "Fig2d", dataset = r$dataset, estimate = r$estimate,
      lower = r$lower, upper = r$upper, primary_p_Holm = r$p_Holm,
      ci_adjustment = r$ci_adjustment, ci_family_level = r$ci_family_level,
      ci_n_tests = r$ci_n_tests)
    strip <- figure2_statistics_strip(table)
  } else {
    if (tag == "g") {
      audited <- statistics$cyclophosphamide$analysis_source
      keep <- audited$included_in_model &
        (audited$cycle_num == 1L | audited$contains_cyclophosphamide)
      keys <- paste(audited$dataset[keep], audited$sample_id[keep])
      d <- data[paste(data$dataset, data$sample_id) %in% keys, , drop = FALSE]
      stopifnot(nrow(d) == sum(keep))
    } else d <- data[data$regimen %in% c("baseline", regimen), , drop = FALSE]
    d$category_show <- ifelse(d$regimen == "baseline", levels[1], levels[2])
    therapy_colour <- inputs$therapy_palette[if (tag == "g") "Cyclophosphamide" else regimen]
    colours <- setNames(c("grey50", therapy_colour), levels)
    all_contrasts <- statistics$regimen_contrasts
    r <- all_contrasts[all_contrasts$analysis == "Current_patient_intercept_primary" &
      all_contrasts$regimen == regimen & all_contrasts$PC == pc, ]
    r <- r[match(cohort_order, r$dataset), ]
    table <- data.frame(panel = paste0("Fig2", tag), dataset = r$dataset, estimate = r$estimate,
      lower = r$lower.CL, upper = r$upper.CL, primary_p_Holm = r$p_Holm,
      ci_adjustment = r$ci_adjustment, ci_family_level = r$ci_family_level,
      ci_n_tests = r$ci_n_tests)
    strip <- figure2_statistics_strip(table)
  }
  pd <- make_raincloud_data(d, "category_show", pc, levels, tag)
  labels <- vapply(levels, function(category) paste0(category, "\nOriginal ",
    sum(pd$source$category == category & pd$source$dataset == "Original"), "\nmissMDA ",
    sum(pd$source$category == category & pd$source$dataset == "missMDA")), character(1))
  main <- ggplot() +
    geom_polygon(data = pd$cloud_polygons, aes(x, y, group = id, fill = category),
      alpha = .40, colour = "grey40", linewidth = .3) +
    geom_segment(data = pd$median_IQR, aes(x = x, xend = x, y = q25, yend = q75),
      linewidth = .85, colour = "grey25") +
    geom_segment(data = pd$median_IQR, aes(x = x - .03, xend = x + .03, y = median, yend = median),
      linewidth = .8, colour = "grey15") +
    geom_point(data = pd$source, aes(rain_x, score, shape = dataset),
      colour = "black", size = 1.65, stroke = .6, alpha = .80) +
    scale_shape_manual(values = c(Original = 16, missMDA = 1), limits = cohort_order,
      labels = c("Original", "missMDA"), name = NULL) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_x_continuous(breaks = seq_along(levels), labels = labels,
      limits = c(.5, length(levels) + .5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(.06, .07))) +
    labs(title = switch(tag, d = "Baseline PC1: lymphoma",
      e = "PC1: vincristine /\ncyclophosphamide", f = "PC2: vincristine",
      g = "PC2: cyclophosphamide-\ncontaining regimens"),
      x = NULL, y = paste0(pc, " score"), tag = tag) +
    figure2_theme() + theme(legend.position = "top", legend.key.width = grid::unit(12, "pt"),
      legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(3, "pt"),
      plot.margin = margin(7, 8, 3, 8))
  pd$plot <- (main / strip) + patchwork::plot_layout(heights = c(4.35, 1.35))
  pd$statistics <- table
  pd
}

build_figure2_panels <- function(inputs, statistics) {
  data <- figure2_score_data(inputs)
  palettes <- figure2_palettes(inputs)
  subsets <- list(Fig2a = data[data$cycle_num == 1, ],
    Fig2b = data[data$cycle_num > 1, ], Fig2c = data)
  plots <- list(
    Fig2a = make_pca_biplot(subsets$Fig2a, inputs) + labs(tag = "a") +
      theme(legend.key.height = grid::unit(.3, "cm")),
    Fig2b = make_pca_biplot(subsets$Fig2b, inputs) + labs(tag = "b") +
      theme(legend.key.height = grid::unit(.3, "cm")))
  central <- make_pca_biplot(data, inputs, "exposed_regimen", y_limit = c(-4, 6)) +
    guides(colour = "none", shape = "none") + theme(plot.margin = margin(0, 0, 0, 0))
  marginal <- make_pca_marginal_data(data, palettes$therapy)
  mx <- make_pca_marginal("PC1", data, marginal, palettes$therapy)
  my <- make_pca_marginal("PC2", data, marginal, palettes$therapy)
  p_c <- cowplot::plot_grid(mx, make_pca_marginal_legend(palettes$therapy), central, my,
    ncol = 2, nrow = 2, rel_widths = c(1.65, 1), rel_heights = c(1, 1.65), align = "hv", axis = "tblr")
  plots$Fig2c <- cowplot::plot_grid(p_c, labels = "c", label_size = 16,
    label_fontfamily = "sans")
  association_data <- setNames(lapply(c("d", "e", "f", "g"), function(tag)
    make_association_raincloud(tag, data, inputs, statistics)), paste0("Fig2", c("d", "e", "f", "g")))
  plots <- c(plots, lapply(association_data, `[[`, "plot"))
  counts <- do.call(rbind, lapply(names(subsets), function(panel)
    do.call(rbind, lapply(split(subsets[[panel]], subsets[[panel]]$dataset), function(d)
      data.frame(panel = panel, dataset = d$dataset[1], n_samples = nrow(d),
        n_patients = length(unique(d$Patient)))))))
  list(plots = plots, dimensions = list(Fig2a = c(6, 4.8), Fig2b = c(6, 4.8),
    Fig2c = c(12, 9.6), Fig2d = c(4.65, 5.7), Fig2e = c(4.65, 5.7),
    Fig2f = c(4.65, 5.7), Fig2g = c(4.65, 5.7)),
    source = subsets, loading_arrows = pca_loading_arrows(inputs$reference_pca),
    marginal_data = marginal, association_data = association_data, counts = counts,
    cancer_palette = palettes$cancer)
}

build_extended2_panels <- function(inputs) {
  data <- figure2_score_data(inputs)
  original <- data[data$dataset == "Original", , drop = FALSE]
  sources <- list(Ext2a = original[original$cycle_num == 1, , drop = FALSE],
                  Ext2b = original[original$cycle_num > 1, , drop = FALSE])
  stopifnot(all(vapply(sources, function(x) all(x$n_missing_original == 0L), logical(1))),
            nrow(sources$Ext2a) + nrow(sources$Ext2b) == nrow(original))
  plots <- lapply(names(sources), function(name) {
    make_pca_biplot(sources[[name]], inputs, show_cohort_legend = FALSE) +
      labs(tag = sub("Ext2", "", name))
  })
  names(plots) <- names(sources)
  counts <- do.call(rbind, lapply(names(sources), function(name) {
    data.frame(panel = name, dataset = "Original", n_samples = nrow(sources[[name]]),
               n_patients = length(unique(sources[[name]]$Patient)))
  }))
  list(plots = plots, source = sources,
       dimensions = list(Ext2a = c(10, 7), Ext2b = c(10, 7)), counts = counts)
}

# Backward-compatible single baseline panel for interactive inspection.
build_extended2_panel <- function(inputs) {
  result <- build_extended2_panels(inputs)
  list(plot = result$plots$Ext2a, source = result$source$Ext2a,
       dimensions = result$dimensions$Ext2a)
}

export_extended2 <- function(result, plots_dir) {
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(result$plots)) {
    dimensions <- result$dimensions[[name]]
    stem <- paste0("Ext2", toupper(sub("Ext2", "", name)))
    save_figure_plot(result$plots[[name]], file.path(plots_dir, stem),
                    dimensions[1], dimensions[2])
  }
  grDevices::cairo_pdf(file.path(plots_dir, "Extended2.pdf"),
                      width = 10, height = 7, onefile = TRUE, family = "sans")
  for (plot in result$plots) print(plot)
  grDevices::dev.off()
  overview <- cowplot::plot_grid(plotlist = result$plots, ncol = 2, align = "hv")
  save_figure_plot(overview, file.path(plots_dir, "Extended2_overview"), 20, 7, dpi = 180)
  for (extension in c("png", "svg")) file.copy(
    file.path(plots_dir, paste0("Extended2_overview.", extension)),
    file.path(plots_dir, paste0("Extended2.", extension)), overwrite = TRUE)
  invisible(overview)
}

save_figure_plot <- function(plot, stem, width, height, dpi = 300) {
  for (extension in c("pdf", "png", "svg")) {
    set.seed(123)
    ggsave(paste0(stem, ".", extension), plot, width = width, height = height,
      dpi = dpi, bg = "white", limitsize = FALSE,
      device = switch(extension, pdf = grDevices::cairo_pdf, png = ragg::agg_png,
        svg = svglite::svglite))
  }
}

write_figure2_source_tables <- function(result, tables_dir) {
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  csv <- function(data, name) write.csv(data, file.path(tables_dir, name), row.names = FALSE, na = "")
  for (name in names(result$source)) csv(result$source[[name]], paste0(name, "_source.csv"))
  csv(result$loading_arrows, "shared_PCA_loading_arrows.csv")
  csv(result$counts, "panel_cohort_counts.csv")
  for (name in names(result$marginal_data))
    csv(result$marginal_data[[name]], paste0("Fig2c_marginal_", name, ".csv"))
  for (name in names(result$association_data)) {
    data <- result$association_data[[name]]
    for (kind in c("source", "cloud_polygons", "median_IQR", "statistics"))
      csv(data[[kind]], paste0(name, "_", kind, ".csv"))
  }
  csv(data.frame(cancer_type = names(result$cancer_palette), colour = unname(result$cancer_palette)),
      "cancer_palette.csv")
}

export_figure2 <- function(result, plots_dir) {
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(result$plots)) {
    stem <- paste0("Fig2", toupper(sub("Fig2", "", name)))
    dimensions <- result$dimensions[[name]]
    save_figure_plot(result$plots[[name]], file.path(plots_dir, stem), dimensions[1], dimensions[2])
  }
  # Keep native panel dimensions on each PDF page, including A/B at half C's width.
  grDevices::cairo_pdf(file.path(plots_dir, "Figure2.pdf"), width = 12, height = 9.6,
    onefile = TRUE, family = "sans")
  for (name in names(result$plots)) {
    dimensions <- result$dimensions[[name]]
    page <- cowplot::ggdraw() + cowplot::draw_plot(result$plots[[name]],
      x = (12 - dimensions[1]) / 24, y = (9.6 - dimensions[2]) / 19.2,
      width = dimensions[1] / 12, height = dimensions[2] / 9.6)
    print(page)
  }
  grDevices::dev.off()
  # Stack A/B at left with exactly the same combined height as C at right.
  # D-G share the next row, at one quarter of the existing overview width each.
  dimensions <- result$dimensions
  bottom_keys <- c("Fig2d", "Fig2e", "Fig2f", "Fig2g")
  upper_height <- dimensions$Fig2c[2]
  lower_height <- max(vapply(dimensions[bottom_keys],
    function(d) d[2], numeric(1)))
  canvas_width <- sum(vapply(dimensions[bottom_keys],
    function(d) d[1], numeric(1)))
  row_gap <- .3
  canvas <- c(canvas_width, upper_height + row_gap + lower_height)
  column_gap <- canvas_width - dimensions$Fig2a[1] - dimensions$Fig2c[1]
  stopifnot(abs(dimensions$Fig2a[2] + dimensions$Fig2b[2] - upper_height) < 1e-9,
    dimensions$Fig2a[1] == dimensions$Fig2b[1], column_gap >= 0)
  upper_bottom <- lower_height + row_gap
  cells <- data.frame(key = names(result$plots),
    x = c(0, 0, dimensions$Fig2a[1] + column_gap,
      head(cumsum(c(0, vapply(dimensions[bottom_keys], function(d) d[1], numeric(1)))), -1)),
    y = c(upper_bottom + dimensions$Fig2b[2], upper_bottom, upper_bottom,
          rep(0, length(bottom_keys))))
  overview <- cowplot::ggdraw()
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]; panel_size <- dimensions[[cell$key]]
    overview <- overview + cowplot::draw_plot(result$plots[[cell$key]],
      x = cell$x / canvas[1], y = cell$y / canvas[2],
      width = panel_size[1] / canvas[1], height = panel_size[2] / canvas[2])
  }
  save_figure_plot(overview, file.path(plots_dir, "Figure2_overview"), canvas[1], canvas[2], dpi = 150)
  file.copy(file.path(plots_dir, "Figure2_overview.png"),
    file.path(plots_dir, "Figure2.png"), overwrite = TRUE)
  file.copy(file.path(plots_dir, "Figure2_overview.svg"),
    file.path(plots_dir, "Figure2.svg"), overwrite = TRUE)
  invisible(overview)
}
