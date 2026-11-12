# 01_import.R -----------------------------------------------------------------
# Import, typage et vérification des données.
# -----------------------------------------------------------------------------

#' Vérifie la présence des colonnes attendues et la validité de Condition
check_data <- function(data, required_cols) {
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Colonnes manquantes : ", paste(missing_cols, collapse = ", "))
  }
  if (anyNA(data$Condition)) {
    warning("Valeurs de 'Condition' inattendues (NA après conversion en facteur).")
  }
  data
}

#' Importe et type les données 2020
read_data_2020 <- function(path, sheet = sheet_2020) {
  readxl::read_excel(path, sheet = sheet) |>
    dplyr::mutate(
      Ind       = as.factor(Ind),
      Condition = factor(as.character(Condition), levels = condition_levels),
      Replicate = as.factor(Replicate),
      dplyr::across(dplyr::any_of(factor_cols_2020), as.factor),
      dplyr::across(dplyr::any_of(numeric_cols_2020), as.numeric)
    ) |>
    check_data(c("Ind", "Condition", "Replicate",
                 factor_cols_2020, numeric_cols_2020))
}

#' Importe et type les données 2021
read_data_2021 <- function(path, sheet = sheet_2021) {
  readxl::read_excel(path, sheet = sheet) |>
    dplyr::mutate(
      Individu  = as.factor(Individu),
      Condition = factor(as.character(Condition), levels = condition_levels),
      Replicate = as.factor(Replicate),
      dplyr::across(dplyr::any_of(factor_cols_2021), as.factor),
      dplyr::across(dplyr::any_of(numeric_cols_2021), as.numeric)
    ) |>
    check_data(c("Individu", "Condition", "Replicate",
                 factor_cols_2021, numeric_cols_2021))
}

#' Importe la feuille de l'ACP combinée 2020-2021
#' Condition y contient 6 modalités (Hglow, Hghigh, ...) : on ne force donc pas
#' les niveaux de condition_levels (ordre alphabétique, comme dans l'original).
read_data_acp <- function(path, sheet = sheet_acp) {
  readxl::read_excel(path, sheet = sheet) |>
    dplyr::mutate(
      Condition = as.factor(Condition),
      Ind       = as.factor(Ind),
      dplyr::across(dplyr::any_of(pca_vars_acp), as.numeric)
    ) |>
    check_data(c("Ind", "Condition", pca_vars_acp))
}
