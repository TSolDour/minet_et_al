# 05_figures.R ----------------------------------------------------------------
# Panels specific to each variable type, building and saving.
# -----------------------------------------------------------------------------

#' Brightness panel
make_brightness_panel <- function(data, y_var, year = NULL, day = NULL,
                                  y_label = "", use_pattern = FALSE) {
  make_panel(
    data, y_var,
    ylim        = c(0.6, 1.15),
    y_label     = y_label,
    annotations = condition_annotations(labels_brightness, y = 1.02, size = 6),
    year = year, year_y = 1.15,
    day  = day,  day_y  = 1.1,
    ref_line    = 1,
    use_pattern = use_pattern
  )
}

#' Disruptive panel
make_disruptive_panel <- function(data, y_var, labels, year = NULL, day = NULL,
                                  y_label = "", use_pattern = FALSE) {
  make_panel(
    data, y_var,
    ylim        = c(0, 43),
    y_label     = y_label,
    annotations = condition_annotations(labels, y = 35, size = 4),
    year = year, year_y = 43,
    day  = day,  day_y  = 40,
    use_pattern = use_pattern
  )
}

#' Laterality index (LI) panel
make_li_panel <- function(data, y_var, year = NULL, day = NULL, y_label = "",
                          use_pattern = FALSE,
                          stars = NULL, group_letters = NULL) {
  annotations <- dplyr::bind_rows(
    if (!is.null(stars))         condition_annotations(stars, y = 0.8, size = 6),
    if (!is.null(group_letters)) condition_annotations(group_letters, y = 1.2, size = 4)
  )
  make_panel(
    data, y_var,
    ylim        = c(-1, 1.8),
    y_label     = y_label,
    annotations = annotations,
    year = year, year_y = 1.8,
    day  = day,  day_y  = 1.5, day_size = 4,
    ref_line    = 0,
    use_pattern = use_pattern
  )
}

#' Panels for 0-1 variables (AbsLI, vigilance monoculaire)
make_abs_panel <- function(data, y_var, year = NULL, day = NULL, y_label = "",
                           use_pattern = FALSE, group_letters = NULL) {
  annotations <- if (!is.null(group_letters)) {
    condition_annotations(group_letters, y = 1.1, size = 4)
  }
  make_panel(
    data, y_var,
    ylim        = c(0, 1.3),
    y_label     = y_label,
    annotations = annotations,
    year = year, year_y = 1.3,
    day  = day,  day_y  = 1.2, day_size = 4,
    use_pattern = use_pattern
  )
}

#' assembling several panels (patchwork)
combine_panels <- function(..., nrow = NULL, ncol = NULL) {
  if (is.null(nrow) && is.null(ncol)) nrow <- 1
  patchwork::wrap_plots(list(...), nrow = nrow, ncol = ncol)
}

#' save a figure and return path
save_figure <- function(plot, path, width = 12, height = 4, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  path
}

#'Locomotion panel (distance, mobility)
make_locomotion_panel <- function(data, y_var, ylim, y_label, letter_y,
                                  group_letters = labels_locomotion) {
  make_panel(
    data, y_var,
    ylim        = ylim,
    y_label     = y_label,
    annotations = condition_annotations(group_letters, y = letter_y, size = 4)
  )
}

#' Thigmotaxy panel
make_thigmotaxis_panel <- function(data) {
  annotations <- tibble::tibble(x = 1.5, y = 107,
                                label = thigmotaxis_p_label, size = 4)
  make_panel(
    data, "thigmotaxis",
    ylim        = c(50, 110),
    y_label     = "Thigmotaxis (%)",
    annotations = annotations
  ) +
    ggplot2::geom_segment(x = 1, xend = 2, y = 103, linewidth = 1)
}
