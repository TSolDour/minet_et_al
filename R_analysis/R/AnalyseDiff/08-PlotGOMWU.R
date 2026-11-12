###################################################

# Plotting results from the GOMWU test

###################################################

################################
# Loading GOMWU results file
loadGOMWU <- function(file) {
  res <- read.delim(file, sep = "\t") %>%
    data.frame()

  res$GOterms <- factor(res$GOterms, levels = res$GOterms)

  return(res)
}

################################
# Customized plot of GOMWU results
PlotGOMWU <- function(data) {
  
  
  data %>%
    mutate(
      Trend = factor(Trend, levels = c("Induced", "Repressed")),
      GOterms = factor(
        GOterms,
        levels = data %>%
          mutate(Trend = factor(Trend, levels = c("Induced", "Repressed"))) %>%
          arrange(Trend) %>%
          pull(GOterms) %>%
          unique() %>%
          rev())) %>%
    ggplot(aes(Trend, GOterms, color = pval)) +
    geom_point() +
    scale_color_gradient(low = "#E6A0C4", high = "#1E1E1E") +
    theme_bw(base_size=14) +
    theme(axis.text.x = element_text(size=10, angle = 45, vjust=0.7),
          axis.text.y = element_text(size=12)) +
    labs(
      x = "",
      y = "Biological Process"
    )
}


################################
## Retrieve list of DEG involved in the GOMWU signal
#DEG <- tar_read(DEG1)
#annot <- tar_read(AnnotationFile)
#gene2go <- tar_read(GOforMWU1)
#GO <- tar_read(HgCO2_GOMWU_BP)
gomwu2deg <- function(DEG, annot, gene2go, GO){
  
  GOsmall <- GO %>%
    select(GOID, GOterms, Trend)
  
  gene2go.filt <- gene2go %>%
    separate_rows(GOs, sep=";") %>%
    filter(GOs %in% GO$GOID)
  
  DEG.GO <- DEG %>%
    filter(ID %in% gene2go.filt$gene)
  
  gene2go.DEG <- gene2go.filt %>%
    filter(gene %in% DEG.GO$ID) %>%
    left_join(GOsmall, by = c("GOs" = "GOID"))
  
  A <- annot %>%
    dplyr::select(gene, name, product, Description) %>%
    group_by(gene) %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    )) %>%
    filter(gene %in% gene2go.DEG$gene)
  
  A. <- A %>%
    filter(name != "<unknown>")
  
  gene2go.DEG$gene2 <- gene2go.DEG$gene
  
  idx <- match(gene2go.DEG$gene, A.$gene)
  gene2go.DEG$gene2 <- ifelse(is.na(idx), gene2go.DEG$gene, A.$name[idx])
  
  X <- gene2go.DEG %>%
    dplyr::select(gene2, GOs, GOterms, Trend)
  
  return(X)
    
  
}