###################################################

# Plotting results from the DE analysis

###################################################

################################
# Plotting heatmap of DEGs
HeatDEG <- function(dds, X, DEG, Annot, ncutrow, ncutcol, ngap, rows = TRUE) {
  
  
  vsd <- vst(dds)
  a <- assay(vsd) %>%
    as.data.frame() %>%
    select(all_of(X)) %>%
    tibble::rownames_to_column(var = "ID") %>%
    filter(ID %in% DEG$ID)

  b <- a %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(ID = paste0(ID, "_", dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    as.data.frame()

  rownames(b) <- b$ID

  Matrix <- b %>%
    select(-ID)
  

  pheatmap::pheatmap((Matrix), # Modify parameters as convenience
    cluster_cols = FALSE,
    gaps_col = ngap,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    cluster_rows = TRUE,
    scale = "row",
    border_color = "grey",
    # color=colorRampPalette(c("#7294D4","#D5D5D3","#E6A0C4")) (100),
    color = colorRampPalette(c("#046C9A", "#D5D5D3", "#FD6467"))(100),
    cutree_rows = ncutrow,
    cutree_cols = ncutcol,
    drop_levels = TRUE,
    legend = TRUE,
    fontsize = 12,
    cellwidth = 20,
    angle_col = 0,
    show_rownames = rows
  )
}

################################
# Summarise number of DEG for each contrast
SumDEG <- function(A, B, C) {
  A <- A %>%
    as.data.frame()

  B <- B %>%
    as.data.frame()

  C <- C %>%
    as.data.frame()

  A$Trend <- ifelse(A$log2FoldChange > 0, "Induced", "Repressed")
  A$Condition <- "pCO2+MeHg"
  B$Trend <- ifelse(B$log2FoldChange > 0, "Induced", "Repressed")
  B$Condition <- "MeHg"
  C$Trend <- ifelse(C$log2FoldChange > 0, "Induced", "Repressed")
  C$Condition <- "pCO2"

  dat <- rbind(A, B, C) %>%
   dplyr:: mutate(nb = 1) %>%
    dplyr::select(Condition, Trend, nb) %>%
    dplyr::group_by(Condition, Trend) %>%
    dplyr::summarise(nobs = sum(nb)) %>%
    dplyr::mutate(N = ifelse(Trend == "Repressed", -nobs, nobs))

  dat %>%
    ggplot(aes(x = Condition, y = N, fill = Trend)) +
    geom_col(color = "white") +
    theme_bw(base_size=14) +
    scale_y_continuous(labels = abs) +
    scale_x_discrete(labels=c(expression(Hg[Low]), expression(CO[2]), expression(CO[2]+Hg[Low]))) +
    scale_fill_manual(values = c("#FD6467", "#046C9A")) +
    annotate("text", x = 1, y = 35, label = "3", color = "black", size = 6, fontface = "bold") +
    annotate("text", x = 1, y = -35, label = "5", color = "black", size = 6, fontface = "bold") +
    annotate("text", x = 2, y = 35, label = "5", color = "black", size = 6, fontface = "bold") +
    annotate("text", x = 2, y = -35, label = "25", color = "black", size = 6, fontface = "bold") +
    annotate("text", x = 3, y = 35, label = "114", color = "black", size = 6, fontface = "bold") +
    annotate("text", x = 3, y = -35, label = "361", color = "black", size = 6, fontface = "bold") +
    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=14)) +
    labs(
      x = "",
      y = "Differentially expressed genes"
    )
}

################################
### Interaction Plot - finaly not used here
IntPlot <- function(dds, NAME){
  
  vsd <- vst(dds, blind = FALSE)
  expr <- assay(vsd)
  
  res_lrt <- results(dds, name = NAME)
  sig_int <- res_lrt[!is.na(res_lrt$padj) & res_lrt$padj < 0.05, ]
  
  genes <- rownames(sig_int)
  
    df <- expr[genes, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
    left_join(
      as.data.frame(colData(dds)) %>% tibble::rownames_to_column("sample"),
      by = "sample" )
    
    profile_means <- df %>%
      group_by(gene, Mercury, pCO2) %>%
      summarise(mean_expr = mean(expression), .groups = "drop") %>%
      tidyr::unite(condition, Mercury, pCO2) %>%
      tidyr::pivot_wider(names_from = condition, values_from = mean_expr)
    
    profile_mat <- profile_means %>%
      tibble::column_to_rownames("gene") %>%
      scale(center = TRUE, scale = TRUE)
    
    dist_mat <- dist(profile_mat)
    hc <- hclust(dist_mat, method = "ward.D2")
    
    plot(hc, labels = FALSE, hang = -1, main = "Dendrogramme des gènes à interaction")
    rect.hclust(hc, k = 3, border = "red")
    
    k <- 3
    clusters <- cutree(hc, k = k)
    
    res_df <- as.data.frame(sig_int) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(cluster = clusters[gene])
    
    top_genes_by_cluster <- res_df %>%
      filter(!is.na(cluster)) %>%
      group_by(cluster) %>%
      arrange(padj) %>%
      slice_head(n = 1) %>%
      pull(gene)
    
    genes_to_plot <- c("sepoff.g111096", "sepoff.g108796", "sepoff.g110782")
    
    df_plot <- expr[genes_to_plot, ] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      tidyr::pivot_longer(
        -gene,
        names_to = "sample",
        values_to = "expression"
      ) %>%
      dplyr::left_join(
        as.data.frame(colData(dds)) %>% tibble::rownames_to_column("sample"),
        by = "sample"
      )
    
    # Plot interaction
    ggplot(df_plot, aes(x = Mercury, y = expression, color = pCO2, group = pCO2)) +
      stat_summary(fun = mean, geom = "line", linewidth = 1, linetype = "dotted") +
      stat_summary(fun = mean, geom = "point", size=5) +
      facet_wrap(~gene, scales = "free_y", nrow=1) +  # 1 gène par cluster
      theme_bw(base_size = 14) +
      scale_color_manual(values = c("#46ACC8", "#E58601")) +
      labs(
        y = "Expression (VST)",
        x = "Mercury treatment")
    
    
    

}