# 00_config.R -----------------------------------------------------------------
# Constantes partagées par tout le pipeline.
# -----------------------------------------------------------------------------

# Feuilles du fichier Excel
sheet_2020 <- "data2020"
sheet_2021 <- "data2021"

# Ordre, couleurs et motifs des conditions
condition_levels   <- c("Ctrl", "CO2", "Hg", "CO2_Hg")
condition_colors   <- c("green4", "khaki2", "cadetblue", "black")
condition_patterns <- c("none", "none", "stripe", "stripe")

# Annotations des figures.
# ATTENTION : ces valeurs sont reportées à la main depuis les résultats
# statistiques du script d'origine. Elles ne se mettent PAS à jour si les
# données ou les analyses changent.
labels_brightness <- rep("*", 4)   # différence significative vs 1 (Wilcoxon)

labels_disruptive <- list(         # lettres de groupes homogènes (post-hoc)
  d15 = c("a", "b", "b",  "b"),
  d30 = c("a", "b", "bc", "c"),
  d35 = c("a", "c", "b",  "d")
)

# Colonnes numériques à typer à l'import
numeric_cols_2020 <- c(
  "Disr15", "Disr30", "Brightness_score15", "Brightness_score30",
  "LI5", "LI14", "LI25", "AbsLI5", "AbsLI14", "AbsLI25",
  "vigilance_mono5", "vigilance_mono14", "vigilance_mono25",
  "thigmotaxis", "distance", "nb_frame_int", "nb_frame_ext",
  "Mobility_time", "Mobility_percent"
)
numeric_cols_2021 <- c(
  "Disr35", "Brightness_score", "LI14", "AbsLI14", "LI31", "AbsLI31"
)

# --- Laterality : annotations manuelles (une valeur par condition, NA = rien)
stars_li <- list(
  d25 = c("*", NA, NA, NA),   # 2020, jour 25 : Ctrl différent de 0 (Wilcoxon)
  d31 = c("*", NA, NA, NA)    # 2021, jour 31 : idem
)
letters_laterality_d31 <- c("a", "b", "ab", "b")  # LI et AbsLI, 2021, jour 31

# --- Catégories de préférence d'œil -------------------------------------------
# Effectifs saisis à la main (repris du script d'origine).
categories_levels <- c("Ctrl", "CO2", "Hg", "CO2Hg")

categories_counts_2020 <- list(
  Ctrl  = c(left = 16, right = 3,  no_preference = 8),
  CO2   = c(left = 9,  right = 9,  no_preference = 9),
  Hg    = c(left = 8,  right = 8,  no_preference = 11),
  CO2Hg = c(left = 6,  right = 10, no_preference = 9)
)
categories_counts_2021 <- list(
  Ctrl  = c(left = 16, right = 3, no_preference = 7),
  CO2   = c(left = 10, right = 9, no_preference = 8),
  Hg    = c(left = 16, right = 3, no_preference = 8),
  CO2Hg = c(left = 7,  right = 9, no_preference = 9)
)

stars_categories <- list(y2020 = 1, y2021 = c(1, 3))  # positions x des étoiles

# =============================================================================
# PARTIE 3 : locomotion, prédation, PCA
# =============================================================================

sheet_acp <- "acp2020_2021"

# --- Colonnes de prédation (facteurs Yes / No) -------------------------------
factor_cols_2020 <- c("PredaCre3", "PredaCre5", "PredaCre7", "PredaCre12",
                      "PredaCrab_success")
factor_cols_2021 <- c("PredaCre3", "PredaCre5", "PredaCre7")

predation_days_2020 <- c(3, 5, 7, 12)
predation_days_2021 <- c(3, 5, 7)

predation_colors <- c(No = "paleturquoise4", Yes = "paleturquoise3")

# Libellés d'axe x des panneaux de prédation
predation_condition_labels <- c(Ctrl = "Ctrl", CO2 = "CO2", Hg = "Hg",
                                CO2_Hg = "CO2Hg")

# Conditions pour lesquelles les comparaisons par paires entre jours sont
# calculées (2021 : pas pour Hg, Cochran non significatif)
predation_posthoc_2021 <- c("Ctrl", "CO2", "CO2_Hg")

# Annotations manuelles par condition et par jour (NA = rien).
# ATTENTION : reportées à la main depuis les résultats du script d'origine.
predation_annotations <- list(
  y2020 = list(
    Ctrl   = list(group_letters = c("a", "ab", "bc", "c"),
                  stars         = c(NA, NA, "*", "*")),
    CO2    = list(group_letters = c("a", "b", "b", "b"),
                  stars         = c(NA, "*", "*", "*")),
    Hg     = list(group_letters = c("a", "ab", "b", "b"),
                  stars         = c(NA, "*", "*", "*")),
    CO2_Hg = list(group_letters = c("a", "b", "b", "c"),
                  stars         = c(NA, "*", "*", "*"))
  ),
  y2021 = list(
    Ctrl   = list(group_letters = c("a", "a", "b"),
                  stars         = c(NA, NA, "*")),
    CO2    = list(group_letters = c("a", "b", "b"),
                  stars         = c(NA, "*", NA)),
    Hg     = list(group_letters = NULL,
                  stars         = c(NA, "*", "*")),
    CO2_Hg = list(group_letters = c("a", "b", "b"),
                  stars         = c(NA, "*", "*"))
  )
)

# --- Locomotion ---------------------------------------------------------------
labels_locomotion   <- c("a", "b", "a", "ab")  # distance et mobilité
thigmotaxis_p_label <- "p = 0.056"

# --- PCA ----------------------------------------------------------------------
# Les colonnes étaient sélectionnées par position dans le script d'origine
# (ex. data[, c(5, 8, 9, 13, 21)]) ; elles sont ici nommées explicitement.
# À VÉRIFIER : ce sont les variables que vous aviez en tête.
pca_vars_2020 <- c("thigmotaxis", "distance", "Disr30", "Brightness_score30",
                   "LI25")
pca_vars_2021 <- c("Disr35", "Brightness_score", "LI31")
pca_vars_acp  <- c("Disr", "Brightness", "LI")   # feuille "acp2020_2021"

# Libellés d'affichage des conditions (noms = valeurs dans les données,
# ordre = ordre d'affichage)
cluster_condition_labels_2020 <- c(Ctrl = "Ctrl", CO2 = "CO2", Hg = "Hg",
                                   CO2_Hg = "CO2Hg")
cluster_condition_labels_acp <- c(
  Ctrl = "Ctrl", CO2 = "CO2", Hglow = "Hglow", Hghigh = "Hghigh",
  CO2_Hglow = "CO2Hglow", CO2_Hghigh = "CO2Hghigh"
)

# Annotations manuelles des figures de clusters
cluster_annotations_2020 <- tibble::tibble(
  x     = c(1:4, 1:4, 1:4),
  y     = rep(c(115, 105, 50), each = 4),
  label = c("a", "b", "b", "b", rep("*", 4), rep("#", 4)),
  size  = rep(c(4, 7, 4), each = 4)
)
cluster_annotations_acp <- tibble::tibble(
  x     = c(1:6, 1:6, 1, 3, 6),
  y     = c(rep(115, 6), rep(105, 6), 25, 25, 75),
  label = c("a", "c", "b", "c", "c", "c", rep("*", 6), "#", "#", "#"),
  size  = c(rep(4, 6), rep(7, 6), 4, 4, 4)
)
