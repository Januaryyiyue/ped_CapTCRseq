# Export already constructed figure panels at their native dimensions, plus a
# vertically arranged overview that preserves each panel's width and height.
# This helper performs no numerical analysis and reads no cached plot objects.

export_figure_panels <- function(panels, dimensions, plots_dir, panel_stems,
                                 figure_stem, overview_labels = NULL,
                                 dpi = 220) {
  stopifnot(length(panels) == length(dimensions),
            length(panels) == length(panel_stems),
            all(vapply(dimensions, function(x) length(x) == 2L && all(x > 0), logical(1))))
  for (package in c("cowplot", "ragg", "svglite")) {
    if (!requireNamespace(package, quietly = TRUE))
      stop("Install the required plotting package: ", package, call. = FALSE)
  }
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  save_plot <- function(plot, stem, width, height) {
    for (extension in c("pdf", "png", "svg")) {
      device <- switch(extension, pdf = grDevices::cairo_pdf,
                       png = ragg::agg_png, svg = svglite::svglite)
      ggplot2::ggsave(file.path(plots_dir, paste0(stem, ".", extension)),
                     plot, width = width, height = height, units = "in",
                     dpi = dpi, bg = "white", device = device, limitsize = FALSE)
    }
  }
  for (i in seq_along(panels)) {
    save_plot(panels[[i]], panel_stems[i], dimensions[[i]][1], dimensions[[i]][2])
  }
  widths <- vapply(dimensions, `[[`, numeric(1), 1L)
  heights <- vapply(dimensions, `[[`, numeric(1), 2L)
  width <- max(widths)
  height <- sum(heights)
  top <- height
  overview <- cowplot::ggdraw()
  for (i in seq_along(panels)) {
    top <- top - heights[i]
    x <- (width - widths[i]) / 2
    overview <- overview + cowplot::draw_plot(panels[[i]],
      x = x / width, y = top / height,
      width = widths[i] / width, height = heights[i] / height)
    if (!is.null(overview_labels)) {
      overview <- overview + cowplot::draw_plot_label(overview_labels[i],
        x = (x + 0.08) / width, y = (top + heights[i] - 0.08) / height,
        hjust = 0, vjust = 1, size = 13, fontface = "bold")
    }
  }
  save_plot(overview, figure_stem, width, height)
  invisible(list(plot = overview, width = width, height = height))
}
