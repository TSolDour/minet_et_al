# 06_categories.R -------------------------------------------------------------
# Catégories de préférence d'œil (gauche / droite / pas de préférence) :
# construction des tables, tests et figure en barres empilées.
# -----------------------------------------------------------------------------

#' Construit la table longue condition / Category / numbers
#' @param counts liste nommée par condition, chaque élément étant un vecteur
#'               nommé (left, right, no_preference) -- cf. 00_config.R
make_categories_table <- function(counts) {
  tibble::tibble(
    condition = factor(rep(names(counts), each = 3), levels = categories_levels),
    Category  = rep(c("left", "right", "no preference"), times = length(counts)),
    numbers   = unlist(
      lapply(counts, function(x) x[c("left", "right", "no_preference")]),
      use.names = FALSE
    )
  )
}

#' Ajoute les pourcentages par condition et ordonne les catégories (pour le graphique)
categories_to_percent <- function(tbl) {
  tbl |>
    dplyr::group_by(condition) |>
    dplyr::mutate(percentage = numbers / sum(numbers) * 100) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      Category = factor(Category, levels = c("no preference", "right", "left"))
    )
}

#' Tests sur la table de catégories
#'  - fisher : différence de répartition entre conditions
#'  - chisq  : écart à 1/3 par catégorie, pour chaque condition
#'  - binom  : parmi les individus latéralisés (gauche + droite), écart à 50 %
test_categories <- function(tbl) {
  xt <- xtabs(numbers ~ condition + Category, data = tbl)

  chisq <- dplyr::bind_rows(lapply(rownames(xt), function(cond) {
    broom::tidy(chisq.test(xt[cond, ], p = rep(1 / 3, 3))) |>
      dplyr::mutate(condition = cond, .before = 1)
  }))

  binom <- dplyr::bind_rows(lapply(rownames(xt), function(cond) {
    left  <- unname(xt[cond, "left"])
    right <- unname(xt[cond, "right"])
    broom::tidy(binom.test(left, left + right, p = 0.5)) |>
      dplyr::mutate(condition = cond, left = left, right = right, .before = 1)
  }))

  list(
    table   = xt,
    fisher  = broom::tidy(fisher.test(xt)),
    chisq   = chisq,
    binom   = binom
  )
}

#' Barres empilées (%) de préférence d'œil par condition
#' @param tbl_percent sortie de categories_to_percent()
#' @param star_x      positions x des étoiles (NULL = aucune)
make_category_plot <- function(tbl_percent, year, star_x = NULL,
                               y_label = "", show_legend = FALSE) {
  p <- ggplot2::ggplot(
    tbl_percent,
    ggplot2::aes(x = condition, y = percentage, fill = Category)
  ) +
    ggplot2::geom_bar(stat = "identity", color = "black") +
    ggplot2::coord_cartesian(ylim = c(0, 120)) +
    ggplot2::labs(y = y_label, x = "Conditions") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10))
    ) +
    ggplot2::annotate("text", x = 2.5, y = 118, label = year,
                      size = 4, fontface = "bold") +
    ggplot2::scale_fill_manual(
      values = as.vector(wesanderson::wes_palette("GrandBudapest1", n = 3))
    )

  for (x in star_x) {
    p <- p + ggplot2::annotate("text", x = x, y = 105, label = "*", size = 6)
  }

  p
}
