# 00_config.R -----------------------------------------------------------------
# Constantes partagées par tout le pipeline.
# -----------------------------------------------------------------------------

# Excel sheets
sheet_2020 <- "data2020"
sheet_2021 <- "data2021"

# Condition aesthetics
condition_levels   <- c("Ctrl", "CO2", "Hg", "CO2_Hg")
condition_colors   <- c("green4", "khaki2", "cadetblue", "black")
condition_patterns <- c("none", "none", "stripe", "stripe")

# Figure annotation
# ATTENTION : Theses values are reported manually, they do not update if analyses change.
labels_brightness <- rep("*", 4)   # significant difference vs 1 (Wilcoxon)

labels_disruptive <- list(         
  d15 = c("a", "b", "b",  "b"),
  d30 = c("a", "b", "bc", "c"),
  d35 = c("a", "c", "b",  "d")
)

# numerical columns to be formated when imported
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

# --- Laterality : manual annotations
stars_li <- list(
  d25 = c("*", NA, NA, NA),   # 2020, day 25 : Ctrl different from 0 (Wilcoxon)
  d31 = c("*", NA, NA, NA)    # 2021, day 31 : idem
)
letters_laterality_d31 <- c("a", "b", "ab", "b")  # LI et AbsLI, 2021, day 31

# --- Eye use preference category -------------------------------------------
# Manually written.
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

stars_categories <- list(y2020 = 1, y2021 = c(1, 3))  # stars x coordinates

# =============================================================================
# PART 3 : locomotion, predation, PCA
# =============================================================================

sheet_acp <- "acp2020_2021"

# --- Predation cols (facteurs Yes / No) -------------------------------
factor_cols_2020 <- c("PredaCre3", "PredaCre5", "PredaCre7", "PredaCre12",
                      "PredaCrab_success")
factor_cols_2021 <- c("PredaCre3", "PredaCre5", "PredaCre7")

predation_days_2020 <- c(3, 5, 7, 12)
predation_days_2021 <- c(3, 5, 7)

predation_colors <- c(No = "paleturquoise4", Yes = "paleturquoise3")

# X axes titles
predation_condition_labels <- c(Ctrl = "Ctrl", CO2 = "CO2", Hg = "Hg",
                                CO2_Hg = "CO2Hg")

# Conditions for which daily pairwise comparisons are comupted
predation_posthoc_2021 <- c("Ctrl", "CO2", "CO2_Hg")

# Per condition/day manual annotations
# ATTENTION : manually reported
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
labels_locomotion   <- c("a", "b", "a", "ab")  # distance et mobility
thigmotaxis_p_label <- "p = 0.056"

# --- PCA ----------------------------------------------------------------------
pca_vars_2020 <- c("thigmotaxis", "distance", "Disr30", "Brightness_score30",
                   "LI25")
pca_vars_2021 <- c("Disr35", "Brightness_score", "LI31")
pca_vars_acp  <- c("Disr", "Brightness", "LI")   # feuille "acp2020_2021"

# 2020 PCA cluster figure
cluster_condition_labels_2020 <- c(Ctrl = "Ctrl", CO2 = "CO2", Hg = "Hg",
                                   CO2_Hg = "CO2Hg")
# manual annotation (a / b, étoiles, #)
cluster_annotations_2020 <- tibble::tibble(
  x     = c(1:4, 1:4, 1:4),
  y     = rep(c(115, 105, 50), each = 4),
  label = c("a", "b", "b", "b", rep("*", 4), rep("#", 4)),
  size  = rep(c(4, 7, 4), each = 4)
)

# --- 2020-2021 PCA : Figure 4 (panels a-e) ------------------------------------
condition_levels_acp <- c("Ctrl", "CO2", "Hglow", "Hghigh",
                          "CO2_Hglow", "CO2_Hghigh")
condition_labels_acp <- c("Ctrl", "CO\u2082", "HgLow", "HgHigh",
                          "CO\u2082+HgLow", "CO\u2082+HgHigh")
condition_axis_labels_acp <- c("Control", "CO\u2082", "HgLow", "HgHigh",
                               "CO\u2082+\nHgLow", "CO\u2082+\nHgHigh")
condition_cols_acp <- c(Ctrl = "#45A844", CO2 = "#DDD24D",
                        Hglow = "#19BEE0", Hghigh = "#14829C",
                        CO2_Hglow = "#F59A4A", CO2_Hghigh = "#E15D23")
cluster_cols_acp <- c("1" = "#E5A84D", "2" = "#58AABD", "3" = "#245C87")
vector_cols_acp  <- c(Brightness = "#6943C1", LI = "#1779B9", Disr = "#E46C30")

ellipse_level  <- 0.95   # 95 % confidence interval
acp_n_clusters <- 3      # number of clusters from HCPC

# Original effectif : 
reference_counts_acp <- matrix(
  c(3, 6, 45,   29, 23, 2,   5, 7, 15,
    16, 10, 1,  16, 11, 0,   19, 7, 1),
  nrow = 6, byrow = TRUE,
  dimnames = list(condition_levels_acp, c("1", "2", "3"))
)
letters_clusters_acp  <- c("a", "c", "b", "c", "c", "c")
expected_variance_acp <- c(46.3, 33.1)   # Variance % for PC1 & PC2

# Panel c (vectors) : 
pca_vector_labels <- c(Brightness = "Brightness score",
                       LI         = "Laterality index",
                       Disr       = "Disruptive score")
vector_label_offsets <- tibble::tibble(
  variable = c("Brightness", "LI", "Disr"),
  dx       = c(0,   -0.08,  0.18),
  dy       = c(0.17, 0.14, -0.08),
  hjust    = c(0.5,  1,     0.5)
)

# Figure 4 output folder
fig4_out_dir_fig <- "results/Behavior/figures"
fig4_out_dir_tab <- "results/Behavior/tables"
