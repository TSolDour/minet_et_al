# 04_plot_helpers.R -----------------------------------------------------------
# Fonction générique de construction d'un panneau (boxplot par condition).
# Nécessite ggplot2 >= 3.5.0 (argument `staplewidth`).
# -----------------------------------------------------------------------------

#' Annotations (une valeur par condition) sous forme de tibble x / y / label / size
#' Les NA dans `labels` sont ignorés (pas d'annotation pour cette condition).
condition_annotations <- function(labels, y, size) {
  tibble::tibble(x = seq_along(labels), y = y, label = labels, size = size) |>
    dplyr::filter(!is.na(label))
}

#' Ajoute à un graphique des annotations textuelles (tibble x / y / label / size)
add_annotations <- function(p, annotations) {
  if (is.null(annotations)) return(p)
  for (i in seq_len(nrow(annotations))) {
    p <- p + ggplot2::annotate(
      "text",
      x = annotations$x[i], y = annotations$y[i],
      label = annotations$label[i], size = annotations$size[i]
    )
  }
  p
}

#' Boxplot par condition avec annotations
#'
#' @param data        data.frame
#' @param y_var       nom de la variable à tracer (chaîne)
#' @param ylim        limites de l'axe y (coord_cartesian)
#' @param y_label     titre de l'axe y ("" pour le masquer)
#' @param annotations tibble x / y / label / size (cf. condition_annotations()),
#'                    ou NULL
#' @param year,day    textes d'en-tête (NULL pour ne pas les afficher)
#' @param year_y,day_y positions y des textes d'en-tête
#' @param day_size    taille du texte "day"
#' @param ref_line    ordonnée de la ligne pointillée de référence (NULL = aucune)
#' @param use_pattern boxplots avec motifs (Hg / CO2_Hg rayés), pour 2021
make_panel <- function(data, y_var, ylim, y_label = "", annotations = NULL,
                       year = NULL, day = NULL,
                       year_y = NULL, day_y = NULL, day_size = 3,
                       ref_line = NULL, use_pattern = FALSE) {

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = Condition, y = .data[[y_var]], fill = Condition)
  )

  if (use_pattern) {
    p <- p +
      ggpattern::geom_boxplot_pattern(
        ggplot2::aes(pattern = Condition),
        pattern_fill = "white", pattern_angle = 45, staplewidth = 0.5
      ) +
      ggpattern::scale_pattern_manual(values = condition_patterns)
  } else {
    p <- p + ggplot2::geom_boxplot(staplewidth = 0.5)
  }

  p <- p +
    ggplot2::ylab(y_label) +
    ggplot2::xlab("Conditions") +
    ggplot2::scale_fill_manual(values = condition_colors) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))
    ) +
    ggplot2::coord_cartesian(ylim = ylim)

  if (!is.null(ref_line)) {
    p <- p + ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed",
                                 linewidth = 0.5)
  }

  if (!is.null(year)) {
    p <- p + ggplot2::annotate("text", x = 2.5, y = year_y, label = year,
                               size = 4, fontface = "bold")
  }
  if (!is.null(day)) {
    p <- p + ggplot2::annotate("text", x = 2.5, y = day_y, label = day,
                               size = day_size)
  }

  add_annotations(p, annotations)
}
