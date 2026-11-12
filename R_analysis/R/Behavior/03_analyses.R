# 03_analyses.R ---------------------------------------------------------------
# Analyses complètes, assemblées à partir des briques de 02_stats_helpers.R.
# -----------------------------------------------------------------------------

#' Analyse complète d'une variable en fonction de Condition
#'
#' @param data     data.frame
#' @param var      nom de la variable réponse (chaîne)
#' @param mu       si non NULL : Wilcoxon par condition contre cette valeur
#' @param post_hoc si TRUE : Tukey HSD + tests par permutation par paires (FDR)
#' @param np       nb de permutations pour aovperm
#' @param nperm    nb de permutations pour les tests par paires
#' @return liste nommée de résultats (tibbles / objets de modèles)
analyse_variable <- function(data, var, mu = NULL, post_hoc = TRUE,
                             np = 5000, nperm = 10000) {
  model <- run_aov(data, var)

  out <- list(
    variable            = var,
    normality_groups    = check_normality_by_group(data, var),
    perm_anova          = run_perm_anova(data, var, np = np),
    anova               = broom::tidy(model),
    normality_residuals = check_normality_residuals(model)
  )

  if (!is.null(mu)) {
    out$vs_value <- test_against_value(data, var, mu = mu)
  }

  if (post_hoc) {
    out$tukey         <- broom::tidy(TukeyHSD(model))
    out$pairwise_perm <- run_pairwise_perm(data, var, nperm = nperm)
  }

  out
}

#' ANOVA à mesures répétées Condition x Days (formule Error(1/id) d'origine)
#'
#' @param data_long format long issu de pivot_two_way()
#' @param response  nom de la variable réponse (chaîne)
#' @param id        nom de la colonne identifiant l'individu ("Ind" en 2020,
#'                  "Individu" en 2021)
twoway_aov <- function(data_long, response, id) {
  fml <- as.formula(
    paste0(response, " ~ Condition * Days + Error(1 / ", id, ")")
  )
  mod <- aov(fml, data = data_long)
  list(model = mod, summary = summary(mod))
}

#' ANOVA par permutation à 2 facteurs, mesures répétées : Error(id/Days)
twoway_perm <- function(data_long, response, id) {
  fml <- as.formula(
    paste0(response, " ~ Condition * Days + Error(", id, " / Days)")
  )
  permuco::aovperm(fml, data = data_long)
}
