# 03_analyses.R ---------------------------------------------------------------
# Full analyses built from 02_stats_helpers.R basic functions
# -----------------------------------------------------------------------------

#' Full analysis of a variable using Condition as explaining variable
#' @param data     data.frame
#' @param var      response variable's name (string)
#' @param mu       if != NULL : Wilcoxon per condition against this value
#' @param post_hoc if TRUE : Tukey HSD + permutation pairwise t tests (FDR)
#' @param np       nb of permutation for aovperm
#' @param nperm    nb of permutations for pairwise tests
#' @return named list of results (tibbles / models objects)
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

#' Repeated measures ANOVA Condition x Days (formula Error(1/id))
#'
#' @param data_long long format from pivot_two_way()
#' @param response  response variable's name (string)
#' @param id        Individual variable ("Ind" in 2020, "Individu" in 2021)
twoway_aov <- function(data_long, response, id) {
  fml <- as.formula(
    paste0(response, " ~ Condition * Days + Error(1 / ", id, ")")
  )
  mod <- aov(fml, data = data_long)
  list(model = mod, summary = summary(mod))
}

#' 2-way permutation ANOVA, repeated : Error(id/Days)
twoway_perm <- function(data_long, response, id) {
  fml <- as.formula(
    paste0(response, " ~ Condition * Days + Error(", id, " / Days)")
  )
  permuco::aovperm(fml, data = data_long)
}
