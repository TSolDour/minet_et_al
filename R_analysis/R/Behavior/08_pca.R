# 08_pca.R --------------------------------------------------------------------
# ACP (FactoMineR) avec imputation des valeurs manquantes (missMDA),
# classification hiérarchique sur composantes principales (HCPC),
# comparaison des proportions de clusters entre conditions, figures.
# -----------------------------------------------------------------------------

#' ACP après imputation des valeurs manquantes
#' @param data  data.frame contenant `vars` et `group`
#' @param vars  noms des variables actives (chaînes)
#' @param group nom de la variable de groupe (supplémentaire qualitative)
#' @param ncp   nombre de dimensions conservées
#' @return liste : table (données imputées + Condition), pca, dimdesc
#' Remarque : Condition est re-factorisée dans l'ordre alphabétique, comme dans
#' le script d'origine (ce qui fixe l'attribution des couleurs de la palette).
run_pca <- function(data, vars, group = "Condition", ncp = 5) {
  imputed <- missMDA::imputePCA(as.data.frame(data[, vars]))$completeObs

  tbl <- data.frame(imputed, Condition = factor(as.character(data[[group]])))

  pca <- FactoMineR::PCA(
    tbl, scale.unit = TRUE, ncp = ncp, graph = FALSE,
    quali.sup = ncol(tbl)
  )

  list(table = tbl, pca = pca, dimdesc = FactoMineR::dimdesc(pca))
}

#' Classification hiérarchique sur composantes principales (nb de clusters automatique)
run_hcpc <- function(pca) {
  FactoMineR::HCPC(pca, graph = FALSE)
}

#' Table de contingence cluster x condition
cluster_table <- function(hcpc, pca_table) {
  table(hcpc$data.clust$clust, pca_table$Condition)
}

#' Tests sur la table cluster x condition
#'  - chisq / pairwise_nominal : proportions de clusters différentes entre conditions
#'  - within_condition (optionnel) : par condition, répartition des clusters
#'    vs équiprobabilité (chi-deux) et comparaisons par paires de clusters
test_cluster_table <- function(tab, within_condition = FALSE) {
  out <- list(
    table            = tab,
    chisq            = broom::tidy(chisq.test(tab)),
    pairwise_nominal = rcompanion::pairwiseNominalIndependence(t(tab),
                                                               method = "fdr")
  )

  if (within_condition) {
    conds <- colnames(tab)
    names(conds) <- conds
    out$within_condition <- lapply(conds, function(cond) {
      x <- tab[, cond]
      list(
        chisq         = broom::tidy(chisq.test(x)),
        pairwise_prop = pairwise.prop.test(x, rep(sum(x), length(x)),
                                           p.adjust.method = "fdr")
      )
    })
  }

  out
}

#' Barres empilées (%) de la répartition des clusters par condition
#' @param tab              table cluster x condition
#' @param annotations      tibble x / y / label / size (cf. 00_config.R)
#' @param condition_labels vecteur nommé : noms dans les données -> libellés
#'                         (l'ordre du vecteur fixe l'ordre d'affichage)
#' @param palette          nom d'une palette wesanderson
make_cluster_plot <- function(tab, annotations, condition_labels, palette) {
  tbl <- as.data.frame(tab) |>
    dplyr::rename(Clusters = Var1, condition = Var2, numbers = Freq) |>
    dplyr::mutate(
      condition = factor(condition_labels[as.character(condition)],
                         levels = unname(condition_labels))
    ) |>
    dplyr::group_by(condition) |>
    dplyr::mutate(percentage = numbers / sum(numbers) * 100) |>
    dplyr::ungroup()

  p <- ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = condition, y = percentage, fill = Clusters)
  ) +
    ggplot2::geom_bar(stat = "identity", color = "black") +
    ggplot2::coord_cartesian(ylim = c(0, 120)) +
    ggplot2::labs(y = "Proportions (%)", x = "Conditions") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10))
    ) +
    ggplot2::scale_fill_manual(
      values = as.vector(wesanderson::wes_palette(palette, n = 3))
    )

  add_annotations(p, annotations)
}

#' Nuage des individus (ellipses de confiance à 95 %) sur deux axes
pca_ind_plot <- function(pca, axes = c(1, 2)) {
  factoextra::fviz_pca_ind(
    pca, axes = axes, geom = "point",
    habillage = "Condition",
    addEllipses = TRUE, ellipse.level = 0.95,
    palette = "jco", pointshape = 19, pointsize = 2, repel = TRUE
  ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right",
                   text = ggplot2::element_text(size = 14))
}

#' Sauvegarde les 4 figures d'une ACP (cercles des variables axes 1-2 et 1-3,
#' nuages d'individus axes 1-2 et 1-3) et retourne les 4 chemins
save_pca_figures <- function(pca_res, prefix) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  axes_list <- list(c(1, 2), c(1, 3))
  paths <- paste0(prefix, c("_var_axes12.png", "_var_axes13.png",
                            "_ind_axes12.png", "_ind_axes13.png"))

  # Cercles des variables (graphiques de base de FactoMineR)
  for (k in 1:2) {
    grDevices::png(paths[k], width = 2000, height = 2000, res = 300)
    plot(pca_res$pca, choix = "var", axes = axes_list[[k]])
    grDevices::dev.off()
  }

  # Nuages d'individus (ggplot)
  for (k in 1:2) {
    ggplot2::ggsave(paths[2 + k],
                    plot = pca_ind_plot(pca_res$pca, axes = axes_list[[k]]),
                    width = 7, height = 5, dpi = 300)
  }

  paths
}
