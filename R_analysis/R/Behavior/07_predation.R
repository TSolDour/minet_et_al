# 07_predation.R --------------------------------------------------------------
# Succès de prédation (crevettes aux jours 3, 5, 7, 12 ; crabe) :
# tables longues, tests et figures.
# -----------------------------------------------------------------------------

#' Table longue : une ligne par individu et par jour
#' @param data   data.frame 2020 ou 2021
#' @param id_col nom de la colonne d'identifiant ("Ind" en 2020, "Individu" en 2021)
#' @param days   jours testés (ex. c(3, 5, 7, 12) -> colonnes PredaCre3, ...)
#' @return tibble : Individuals, Condition, Day (facteur Day3, Day5, ...), success
make_predation_long <- function(data, id_col, days) {
  cols <- paste0("PredaCre", days)
  data |>
    dplyr::select(Individuals = dplyr::all_of(id_col), Condition,
                  dplyr::all_of(cols)) |>
    tidyr::pivot_longer(
      cols         = dplyr::all_of(cols),
      names_to     = "Day",
      names_prefix = "PredaCre",
      values_to    = "success"
    ) |>
    dplyr::mutate(Day = factor(paste0("Day", Day),
                               levels = paste0("Day", days))) |>
    dplyr::arrange(Day)   # même ordre que le script d'origine (jour par jour)
}

#' Q de Cochran (comparaison des jours) par condition, et McNemar par paires
#' (permutation, FDR) pour les conditions de `post_hoc`
#' @return liste nommée par condition (cochran, mcnemar)
test_predation_by_condition <- function(long, post_hoc = condition_levels) {
  conds <- levels(long$Condition)
  names(conds) <- conds
  lapply(conds, function(cond) {
    d <- as.data.frame(long[long$Condition == cond, ])
    out <- list(
      cochran = RVAideMemoire::cochran.qtest(success ~ Day | Individuals, data = d)
    )
    if (cond %in% post_hoc) {
      out$mcnemar <- rcompanion::pairwiseMcnemar(
        success ~ Day | Individuals, data = d,
        test = "permutation", method = "fdr", digits = 3
      )
    }
    out
  })
}

#' Test binomial (p = 0.5) du succès, pour chaque condition et chaque jour
binom_by_day <- function(long) {
  long |>
    dplyr::group_by(Condition, Day) |>
    dplyr::summarise(
      n          = dplyr::n(),
      success    = sum(ifelse(success == "Yes", 1, 0)),
      proportion = success / n,
      p_value    = binom.test(success, n, p = 0.5)$p.value,
      .groups    = "drop"
    )
}

#' Chi-deux entre conditions, séparément pour chaque jour
chisq_by_day <- function(long) {
  long |>
    dplyr::group_by(Day) |>
    dplyr::group_modify(~ broom::tidy(
      chisq.test(table(.x$Condition, .x$success))
    )) |>
    dplyr::ungroup()
}

#' Prédation sur le crabe (2020) : chi-deux entre conditions et binomial par condition
test_crab <- function(data) {
  tab <- table(data$Condition, data$PredaCrab_success)
  list(
    table = tab,
    chisq = broom::tidy(chisq.test(tab)),
    binom = tibble::tibble(
      Condition = rownames(tab),
      p_value   = unname(apply(tab, 1, function(row) {
        binom.test(x = row[["Yes"]], n = sum(row), p = 0.5)$p.value
      }))
    )
  )
}

#' Un panneau (une condition) : barres empilées de succès (%) par jour
make_predation_panel <- function(pct, x_label, y_label = "",
                                 annotations = NULL, show_legend = FALSE) {
  p <- ggplot2::ggplot(
    pct,
    ggplot2::aes(x = Day, y = percentage, fill = success)
  ) +
    ggplot2::geom_bar(stat = "identity", color = "black") +
    ggplot2::coord_cartesian(ylim = c(0, 120)) +
    ggplot2::scale_fill_manual(values = predation_colors) +
    ggplot2::labs(y = y_label, x = x_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))
    )
  add_annotations(p, annotations)
}

#' Figure complète : un panneau par condition (grille 2 x 2)
#' @param annotations liste par condition (group_letters, stars), cf. 00_config.R
make_predation_figure <- function(long, annotations) {
  conds  <- condition_levels
  panels <- lapply(seq_along(conds), function(i) {
    cond <- conds[i]
    pct <- long |>
      dplyr::filter(Condition == cond) |>
      dplyr::group_by(Day) |>
      dplyr::count(success) |>
      dplyr::mutate(percentage = n / sum(n) * 100) |>
      dplyr::ungroup()

    ann <- annotations[[cond]]
    notes <- dplyr::bind_rows(
      if (!is.null(ann$group_letters))
        condition_annotations(ann$group_letters, y = 115, size = 4),
      if (!is.null(ann$stars))
        condition_annotations(ann$stars, y = 105, size = 6)
    )

    make_predation_panel(
      pct,
      x_label     = predation_condition_labels[[cond]],
      y_label     = if (i %% 2 == 1) "Success (%)" else "",
      annotations = notes,
      show_legend = (i == length(conds))
    )
  })
  do.call(combine_panels, c(panels, list(ncol = 2)))
}
