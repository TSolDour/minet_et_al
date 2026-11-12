# 02_stats_helpers.R ----------------------------------------------------------
# Briques statistiques élémentaires. Chaque fonction prend un data.frame et le
# NOM (chaîne) de la variable réponse, et retourne un objet / tibble.
# -----------------------------------------------------------------------------

#' Wilcoxon par condition contre une valeur de référence (ex. mu = 1)
test_against_value <- function(data, var, mu = 1) {
  data |>
    dplyr::group_by(Condition) |>
    dplyr::summarise(
      res     = list(wilcox.test(.data[[var]], mu = mu)),
      W       = unname(res[[1]]$statistic),
      p_value = res[[1]]$p.value,
      .groups = "drop"
    ) |>
    dplyr::select(-res)
}

#' Shapiro-Wilk par condition
check_normality_by_group <- function(data, var) {
  data |>
    dplyr::group_by(Condition) |>
    dplyr::summarise(
      res     = list(shapiro.test(.data[[var]])),
      W       = unname(res[[1]]$statistic),
      p_value = res[[1]]$p.value,
      .groups = "drop"
    ) |>
    dplyr::select(-res)
}

#' ANOVA à un facteur (var ~ Condition)
run_aov <- function(data, var) {
  aov(reformulate("Condition", response = var), data = data)
}

#' Shapiro-Wilk sur les résidus d'un modèle
check_normality_residuals <- function(model) {
  broom::tidy(shapiro.test(residuals(model)))
}

#' ANOVA par permutation (permuco::aovperm)
run_perm_anova <- function(data, var, np = 5000) {
  permuco::aovperm(reformulate("Condition", response = var),
                   data = data, np = np)
}

#' Tests t par permutation par paires, p ajustées (FDR)
#' Retourne un tibble : group1, group2, p_adj
#' `group_var` : facteur de comparaison ("Condition" par défaut, "Days" possible)
run_pairwise_perm <- function(data, var, nperm = 10000, group_var = "Condition") {
  res <- RVAideMemoire::pairwise.perm.t.test(
    resp     = data[[var]],
    fact     = as.factor(data[[group_var]]),
    p.method = "fdr",
    nperm    = nperm,
    progress = FALSE
  )
  as.data.frame(as.table(res$p.value)) |>
    dplyr::filter(!is.na(Freq)) |>
    dplyr::transmute(group1 = Var1, group2 = Var2, p_adj = Freq) |>
    tibble::as_tibble()
}

#' Passe deux colonnes (jours) au format long pour les analyses à 2 facteurs
pivot_two_way <- function(data, cols, values_to) {
  tidyr::pivot_longer(
    data,
    cols      = dplyr::all_of(cols),
    names_to  = "Days",
    values_to = values_to
  )
}
