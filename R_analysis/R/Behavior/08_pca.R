# 08_pca.R --------------------------------------------------------------------
# PCA (FactoMineR) with imputation of missing values (missMDA) + HCPC + figure
# Part 1 : ACP 2020 & 2021.
# Part 2 : ACP 2020-2021 (Figure 4, panels a-e).
# -----------------------------------------------------------------------------

# =============================================================================
# PART 1 : ACP 2020, 2021 & 2020-2021
# =============================================================================

#' PCA after imputation of missing values
#' @param data            data.frame containing `vars` et `group`
#' @param vars            variables nuames
#' @param group           grouping variable
#' @param ncp             number of PCs
#' @param condition_levels levels of the condition variable. 
run_pca <- function(data, vars, group = "Condition", ncp = 5,
                    condition_levels = NULL) {
  imputed <- missMDA::imputePCA(
    as.data.frame(data[, vars]),
    ncp = 2, scale = TRUE, method = "Regularized"
  )$completeObs

  if (is.null(condition_levels)) {
    cond <- factor(as.character(data[[group]]))
  } else {
    cond <- factor(as.character(data[[group]]), levels = condition_levels)
    if (anyNA(cond)) stop("Unrecognized experimental condition in spreadsheet.")
  }

  tbl <- data.frame(imputed, Condition = cond)

  pca <- FactoMineR::PCA(
    tbl, scale.unit = TRUE, ncp = ncp, graph = FALSE,
    quali.sup = ncol(tbl)
  )

  list(table = tbl, imputed = as.data.frame(imputed), pca = pca,
       dimdesc = FactoMineR::dimdesc(pca))
}

#' HCPC
#' @param nb_clust Nb of clusters
run_hcpc <- function(pca, nb_clust = -1) {
  FactoMineR::HCPC(pca, nb.clust = nb_clust, consol = TRUE, graph = FALSE)
}

#' Contingency table cluster x condition (ACP 2020 & 2021)
cluster_table <- function(hcpc, pca_table) {
  table(hcpc$data.clust$clust, pca_table$Condition)
}

#' Tests on the cluster x condition table (ACP 2020 & 2021)
#'  - chisq / pairwise_nominal : test for difference of cluster proportions between conditions
#'  - within_condition : per condition, clusters distribution
#'    vs equiprobability (chi-square) and ccluster pairwise comparisons
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

#' barplot (%) of cluster distribution between conditions (ACP 2020)
#' @param tab              table cluster x condition
#' @param annotations      tibble x / y / label / size (cf. 00_config.R)
#' @param condition_labels labelled vector
#' @param palette          wesanderson color palette
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

#' dot of individuals (95 % CI)
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

#' Saving the 4 figures from PCA
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

# =============================================================================
# PART 2 : ACP 2020-2021 -- clusters, tests, tables (Figure 4)
# =============================================================================

#' Relabel clusters from HCPC based on biological profil
#' @return factor Cluster (levels 1, 2, 3) ordered based on PCA lines
relabel_clusters_by_profile <- function(hcpc, pca_res,
                                        high_var = "Disr",
                                        order_var = "Brightness") {
  pca_ids  <- rownames(pca_res$pca$ind$coord)
  hcpc_ids <- rownames(hcpc$data.clust)
  idx <- match(pca_ids, hcpc_ids)
  if (anyNA(idx)) stop("HCPC and PCA individual identifiers do not align.")

  raw_cluster <- as.character(hcpc$data.clust$clust[idx])
  if (length(raw_cluster) != nrow(pca_res$table)) {
    stop("HCPC result has unexpected length.")
  }

  raw_means <- aggregate(pca_res$imputed,
                         by = list(Cluster = raw_cluster), FUN = mean)

  cluster3_old <- as.character(raw_means$Cluster[which.max(raw_means[[high_var]])])
  other <- raw_means[as.character(raw_means$Cluster) != cluster3_old, ,
                     drop = FALSE]
  other <- other[order(other[[order_var]]), , drop = FALSE]
  old_to_new <- setNames(c("1", "2", "3"),
                         c(as.character(other$Cluster), cluster3_old))

  cluster <- factor(unname(old_to_new[raw_cluster]), levels = c("1", "2", "3"))
  stopifnot(!anyNA(cluster))
  cluster
}

#' Contingency table Condition x Cluster
make_cluster_counts <- function(pca_res, cluster) {
  table(Condition = pca_res$table$Condition, Cluster = cluster)
}

#' Compare computed effectif with original ones
check_counts_against_reference <- function(counts, reference) {
  exact_match <- identical(dim(counts), dim(reference)) &&
    isTRUE(all(unname(counts) == unname(reference)))
  message("Exact match with historical PCA.R counts: ", exact_match)
  list(
    exact_match = exact_match,
    recomputed  = as.data.frame.matrix(counts),
    reference   = as.data.frame.matrix(reference)
  )
}

#' Find individuals with missing values
audit_missing <- function(data, vars, condition) {
  data <- as.data.frame(data)
  audit <- data.frame(
    Annee     = data$Annee,
    Ind       = data$Ind,
    Condition = condition,
    N_missing = rowSums(is.na(data[, vars, drop = FALSE]))
  )
  audit[audit$N_missing > 0, , drop = FALSE]
}

#' Orienting PCA scores
orient_pca_scores <- function(pca_res, data, cluster,
                              expected_variance = expected_variance_acp) {
  pca <- pca_res$pca
  data <- as.data.frame(data)

  scores    <- pca$ind$coord[, 1:3, drop = FALSE]
  var_coord <- pca$var$coord[, 1:3, drop = FALSE]

  signs <- rep(1, 3)
  if (var_coord["Disr", 1] < 0)       signs[1] <- -1
  if (var_coord["Brightness", 2] < 0) signs[2] <- -1
  if (var_coord["LI", 3] < 0)         signs[3] <- -1
  scores    <- sweep(scores, 2, signs, "*")
  var_coord <- sweep(var_coord, 2, signs, "*")
  colnames(scores) <- colnames(var_coord) <- c("PC1", "PC2", "PC3")

  pct <- round(pca$eig[1:3, 2], 1)
  message("Explained variance PC1/PC2/PC3: ", paste(pct, collapse = " / "), "%")
  if (abs(pct[1] - expected_variance[1]) > 1 ||
      abs(pct[2] - expected_variance[2]) > 1) {
    warning("PCA explained variance differs from supplied Figure 4: ",
            "check imputation settings.")
  }

  inds <- data.frame(
    Annee = data$Annee, Ind = data$Ind,
    Condition = pca_res$table$Condition, Cluster = cluster,
    PC1 = scores[, 1], PC2 = scores[, 2], PC3 = scores[, 3],
    row.names = NULL
  )
  vars <- data.frame(
    variable = rownames(var_coord),
    PC1 = var_coord[, 1], PC2 = var_coord[, 2], PC3 = var_coord[, 3],
    row.names = NULL
  )

  list(inds = inds, vars = vars, pct = pct, signs = signs)
}

#' Chi-square global Condition x Cluster, & pairwise chi-square (p.adjusted BH)
test_condition_clusters <- function(counts) {
  chi_global <- suppressWarnings(chisq.test(counts))
  if (any(chi_global$expected < 5)) {
    message("Low expected counts: verify with simulated/exact tests.")
  }

  pairs <- combn(rownames(counts), 2, simplify = FALSE)
  pairwise <- do.call(rbind, lapply(pairs, function(pair) {
    tab <- counts[pair, , drop = FALSE]
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    p_raw <- suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
    data.frame(Group1 = pair[1], Group2 = pair[2], p_value = p_raw)
  }))
  pairwise$p_FDR_BH <- p.adjust(pairwise$p_value, method = "BH")

  list(global = chi_global, pairwise = pairwise)
}

# =============================================================================
# PARTIE 3 : Figure 4 (panels a-e)
# =============================================================================

#' Theme common to all panels
theme_fig4 <- function() {
  ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.line  = ggplot2::element_line(linewidth = .95, colour = "#202020"),
      axis.ticks = ggplot2::element_line(linewidth = .8),
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12, colour = "black"),
      axis.text  = ggplot2::element_text(size = 10, colour = "black"),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text  = ggplot2::element_text(size = 10),
      plot.tag = ggplot2::element_text(size = 15, face = "bold"),
      plot.tag.position = c(.005, .995),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}

#' Dot clound of individuals (panels a, b et d). `group` : "Condition" or "Cluster".
make_pca_score_panel <- function(inds, x, y, group, colors, tag,
                                 x_label, y_label,
                                 legend_title = group,
                                 legend_breaks = ggplot2::waiver(),
                                 legend_labels = ggplot2::waiver(),
                                 legend_position = "right",
                                 legend_key_cm = NULL,
                                 ellipse_alpha = .11, ellipse_lw = .68,
                                 point_size = 1.8, point_alpha = .78,
                                 mean_shape = 16,
                                 ref_linetype = "dashed", ref_colour = "black") {
  p <- ggplot2::ggplot(
    inds,
    ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[group]])
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(fill = .data[[group]]), geom = "polygon",
      level = ellipse_level, type = "norm", alpha = ellipse_alpha,
      colour = NA, show.legend = FALSE
    ) +
    ggplot2::stat_ellipse(level = ellipse_level, type = "norm",
                          linewidth = ellipse_lw) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 4.0,
                          shape = mean_shape, show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 0, linetype = ref_linetype,
                        colour = ref_colour, linewidth = .45) +
    ggplot2::geom_hline(yintercept = 0, linetype = ref_linetype,
                        colour = ref_colour, linewidth = .45) +
    ggplot2::scale_colour_manual(values = colors, breaks = legend_breaks,
                                 labels = legend_labels) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(tag = tag, x = x_label, y = y_label, colour = legend_title) +
    theme_fig4() +
    ggplot2::theme(legend.position = legend_position)

  if (!is.null(legend_key_cm)) {
    p <- p + ggplot2::theme(legend.key.size = grid::unit(legend_key_cm, "cm"))
  }
  p
}

#' Labels "C1" / "C2" / "C3" in the middle of each cluster (panel d)
add_cluster_labels <- function(p, inds) {
  centers <- aggregate(inds[, c("PC1", "PC2")],
                       by = list(Cluster = inds$Cluster), FUN = mean)
  p + ggplot2::geom_text(
    data = centers,
    ggplot2::aes(x = PC1 + .16, y = PC2 + .12,
                 label = paste0("C", Cluster), colour = Cluster),
    inherit.aes = FALSE, size = 4.2, fontface = "bold", show.legend = FALSE
  )
}

#' Panel c : correlation circles
make_pca_vector_panel <- function(vars, colors = vector_cols_acp,
                                  labels = pca_vector_labels,
                                  offsets = vector_label_offsets) {
  varplot <- vars |>
    dplyr::left_join(offsets, by = "variable") |>
    dplyr::mutate(
      lab   = unname(labels[variable]),
      lx    = PC1 + dplyr::coalesce(dx, 0),
      ly    = PC2 + dplyr::coalesce(dy, 0),
      hjust = dplyr::coalesce(hjust, 0.5)
    )

  circle <- data.frame(t = seq(0, 2 * pi, length.out = 361))
  circle$x <- cos(circle$t)
  circle$y <- sin(circle$t)

  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = .45) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey75", linewidth = .45) +
    ggplot2::geom_path(data = circle, ggplot2::aes(x, y), linetype = "dotted",
                       linewidth = .6, colour = "grey55") +
    ggplot2::geom_segment(
      data = varplot,
      ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2, colour = variable),
      linewidth = 1.05, arrow = grid::arrow(length = grid::unit(.20, "cm")),
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = varplot,
      ggplot2::aes(x = lx, y = ly, label = lab, colour = variable,
                   hjust = hjust),
      size = 3.45, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = colors) +
    ggplot2::coord_equal(xlim = c(-1.37, 1.37), ylim = c(-1.34, 1.34),
                         expand = FALSE, clip = "off") +
    ggplot2::scale_x_continuous(breaks = c(-1, -.5, 0, .5, 1)) +
    ggplot2::scale_y_continuous(breaks = c(-1, -.5, 0, .5, 1)) +
    ggplot2::labs(tag = "c", x = "PC1 loading (correlation)",
                  y = "PC2 loading (correlation)") +
    theme_fig4()
}

#' Panel e : barplot (%) per condition clusters
make_cluster_stack_panel <- function(counts, exact_match,
                                     conditions = condition_levels_acp,
                                     colors = cluster_cols_acp,
                                     axis_labels = condition_axis_labels_acp,
                                     group_letters = letters_clusters_acp) {
  bar_df <- as.data.frame(counts, stringsAsFactors = FALSE)
  names(bar_df) <- c("Condition", "Cluster", "n")
  bar_df$Condition <- factor(bar_df$Condition, levels = conditions)
  # ggplot empile les niveaux dans l'ordre inverse : cluster 3 en bas, 1 en haut
  bar_df$Cluster <- factor(bar_df$Cluster, levels = c("1", "2", "3"))
  bar_df$group_n <- as.integer(rowSums(counts))[match(as.character(bar_df$Condition),
                                                      rownames(counts))]
  bar_df$pct <- 100 * bar_df$n / bar_df$group_n
  bar_df$inner_label <- ifelse(bar_df$pct >= 10,
                               paste0(round(bar_df$pct), "%"), "")

  annotation <- data.frame(
    Condition = factor(conditions, levels = conditions),
    n         = as.integer(rowSums(counts)),
    letter    = if (exact_match) group_letters else rep("", length(conditions))
  )

  ggplot2::ggplot(bar_df, ggplot2::aes(x = Condition, y = pct, fill = Cluster)) +
    ggplot2::geom_col(width = .78, colour = "white", linewidth = .65,
                      position = "stack") +
    ggplot2::geom_text(ggplot2::aes(label = inner_label),
                       position = ggplot2::position_stack(vjust = .5),
                       fontface = "bold", size = 3.3, colour = "#0C2936") +
    ggplot2::geom_text(data = annotation,
                       ggplot2::aes(x = Condition, y = 104,
                                    label = paste0("n=", n)),
                       inherit.aes = FALSE, size = 3.3) +
    ggplot2::geom_text(data = annotation,
                       ggplot2::aes(x = Condition, y = 115, label = letter),
                       inherit.aes = FALSE, size = 4, fontface = "bold") +
    ggplot2::scale_fill_manual(values = colors, breaks = c("1", "2", "3"),
                               drop = FALSE) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::scale_y_continuous(limits = c(0, 120),
                                breaks = c(0, 20, 40, 60, 80, 100),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(tag = "e", x = NULL, y = "Individuals (%)", fill = "Cluster") +
    theme_fig4() +
    ggplot2::theme(legend.position = "right",
                   axis.text.x = ggplot2::element_text(size = 9),
                   legend.key.size = grid::unit(.37, "cm"))
}

#' Building the 5 panels a-e of Figure 4
#' @param oriented   orient_pca_scores() output
#' @param counts     table Condition x Cluster
make_figure4_panels <- function(oriented, counts, exact_match) {
  inds <- oriented$inds
  pct  <- oriented$pct
  axis1 <- sprintf("Dim1 (%.1f%%)", pct[1])
  axis2 <- sprintf("Dim2 (%.1f%%)", pct[2])
  axis3 <- sprintf("Dim3 (%.1f%%)", pct[3])

  a <- make_pca_score_panel(
    inds, "PC1", "PC2", "Condition", condition_cols_acp, tag = "a",
    x_label = axis1, y_label = axis2,
    legend_breaks = condition_levels_acp, legend_labels = condition_labels_acp,
    legend_position = "none"
  )
  b <- make_pca_score_panel(
    inds, "PC1", "PC3", "Condition", condition_cols_acp, tag = "b",
    x_label = axis1, y_label = axis3, legend_title = "Conditions",
    legend_breaks = condition_levels_acp, legend_labels = condition_labels_acp,
    legend_key_cm = .35
  )
  c <- make_pca_vector_panel(oriented$vars)
  d <- make_pca_score_panel(
    inds, "PC1", "PC2", "Cluster", cluster_cols_acp, tag = "d",
    x_label = sprintf("PC1 (%.1f%%)", pct[1]),
    y_label = sprintf("PC2 (%.1f%%)", pct[2]),
    ellipse_alpha = .12, ellipse_lw = .9, point_size = 1.5, point_alpha = .63,
    mean_shape = 19, ref_linetype = "dotted", ref_colour = "grey70"
  ) |>
    add_cluster_labels(inds)
  e <- make_cluster_stack_panel(counts, exact_match)

  list(a = a, b = b, c = c, d = d, e = e)
}

#' Assembling panels
combine_figure4 <- function(panels) {
  top <- (panels$a | panels$b | panels$c) +
    patchwork::plot_layout(widths = c(1.25, 1.60, 1.20))
  bottom <- (panels$d | panels$e) +
    patchwork::plot_layout(widths = c(1.15, 1.45))
  (top / bottom) + patchwork::plot_layout(heights = c(1.02, 1.0))
}

#' Saving Figure 4
save_figure4 <- function(figure4, panels, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, "Figure4_complete_a_to_e")

  ggplot2::ggsave(paste0(base, ".png"), figure4, width = 18.5, height = 10.2,
                  dpi = 400, units = "in", bg = "white")
  ggplot2::ggsave(paste0(base, ".pdf"), figure4, width = 18.5, height = 10.2,
                  units = "in", device = "pdf", bg = "white",
                  useDingbats = FALSE)
  paths <- paste0(base, c(".png", ".pdf"))

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(paste0(base, ".svg"), figure4, width = 18.5, height = 10.2,
                    units = "in", device = "svg", bg = "white")
    paths <- c(paths, paste0(base, ".svg"))
  }

  for (name in names(panels)) {
    png_path <- file.path(out_dir, paste0("Figure4_panel_", name, ".png"))
    pdf_path <- file.path(out_dir, paste0("Figure4_panel_", name, ".pdf"))
    ggplot2::ggsave(png_path, panels[[name]], width = 7.5, height = 5.5,
                    dpi = 400, bg = "white")
    ggplot2::ggsave(pdf_path, panels[[name]], width = 7.5, height = 5.5,
                    bg = "white", useDingbats = FALSE)
    paths <- c(paths, png_path, pdf_path)
  }

  paths
}

#' Exporting tables and résults text associated with Figure 4.
export_figure4_tables <- function(out_dir, check_counts, audit, oriented,
                                  pca_res, hcpc, tests) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  f <- function(name) file.path(out_dir, name)

  write.csv(check_counts$recomputed, f("HCPC_counts_recomputed.csv"))
  write.csv(check_counts$reference,  f("HCPC_counts_original_PCA_R.csv"))
  write.csv(audit, f("missing_data_audit.csv"), row.names = FALSE)
  write.csv(oriented$inds, f("Figure4_individual_PCA_HCPC.csv"), row.names = FALSE)
  write.csv(oriented$vars, f("Table_S1_variable_correlations.csv"),
            row.names = FALSE)
  write.csv(as.data.frame(pca_res$pca$eig), f("PCA_eigenvalues.csv"))
  writeLines(capture.output(print(hcpc$desc.var)),
             f("HCPC_description_vtests.txt"))
  writeLines(capture.output(print(hcpc$desc.axes)), f("HCPC_PC_vtests.txt"))
  writeLines(capture.output(print(tests$global)),
             f("HCPC_condition_chisquare.txt"))
  write.csv(tests$pairwise, f("Pairwise_condition_chisquare_FDR.csv"),
            row.names = FALSE)

  c(f("HCPC_counts_recomputed.csv"), f("HCPC_counts_original_PCA_R.csv"),
    f("missing_data_audit.csv"), f("Figure4_individual_PCA_HCPC.csv"),
    f("Table_S1_variable_correlations.csv"), f("PCA_eigenvalues.csv"),
    f("HCPC_description_vtests.txt"), f("HCPC_PC_vtests.txt"),
    f("HCPC_condition_chisquare.txt"), f("Pairwise_condition_chisquare_FDR.csv"))
}
