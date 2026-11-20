# 01_import.R -----------------------------------------------------------------
# Import, typing & data check.
# -----------------------------------------------------------------------------

#' Check if expected cols are present and Condition validy
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

#' Import and format 2020 data
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

#' Import and format 2021 data
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

#' Import 202-2021 PCA sheet
read_data_acp <- function(path, sheet = sheet_acp) {
  readxl::read_excel(path, sheet = sheet) |>
    dplyr::mutate(
      Condition = as.factor(Condition),
      Ind       = as.factor(Ind),
      dplyr::across(dplyr::any_of(pca_vars_acp), as.numeric)
    ) |>
    check_data(c("Annee", "Ind", "Condition", pca_vars_acp))
}
