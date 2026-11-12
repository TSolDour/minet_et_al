#####################################################
#
#   Retrieving Hub genes for each trait             #
#
#####################################################


## After having defined the significant modules for each trait based on the correlation matrix

################################
## Calculate Module Membership
ModuleMembership <- function(ModuleEigen, Exp.Matrix, cols, R) {
  modNames <- substring(names(ModuleEigen), 3) # extract module names

  nSamples <- nrow(Exp.Matrix)

  geneModuleMembership <- as.data.frame(cor(Exp.Matrix, ModuleEigen, use = "p"))

  names(geneModuleMembership) <- paste("MM", modNames, sep = "")

  return(geneModuleMembership)

}

################################
# Retrieve hub genes based on KME:
HubGenes <- function(Colors, Mod, MM, Module, Variable, Percent, Annot) {
  Col <- Colors %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene") %>%
    filter(. == Mod)

  Hub <- MM %>%
    as.data.frame() %>%
    select({{ Module }}) %>%
    rownames_to_column(var = "Gene") %>%
    filter(Gene %in% Col$Gene) %>%
    filter({{ Variable }} >= quantile({{ Variable }}, Percent, na.rm = TRUE))
  
  a <- Annot
  A <- a %>%
    select(gene, name, product, Description) %>%
    group_by(gene) %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    )) %>%
    filter(gene %in% Hub$Gene)
  
  X <- merge(Hub, A, by.x="Gene", by.y="gene", all.x=TRUE)
  
#  Hub$ID2 <- Hub$Gene
  
#  idx <- match(Hub$Gene, A.$gene)
#  Hub$ID2 <- ifelse(is.na(idx), Hub$Gene, A.$name[idx])
  
#  Hub$Gene <- Hub$ID2
  
#  Hub$ID2 <- NULL
  
  return(X)
}

################################
## Create reduced adjacency and TOM matrices
TOMsub <- function(Col, Matrix, Module, softPower) {
  # Retrieving genes names
  a <- Col
  b <- Matrix

  moduleGenes <- names(a)[a == Module]

  datMod <- Matrix[, moduleGenes]
  adj <- adjacency(datMod, power = softPower, type = "signed")

  TOM <- TOMsimilarity(adj)
  colnames(TOM) <- moduleGenes
  rownames(TOM) <- moduleGenes

  return(TOM)
}



################################
## Export reduced TOM toward Cytoscape
ExportTOMtoCyto <- function(TOMSUB, Annot, Mod, PATH1, PATH2, Colors, THRLD) {
  Col <- Colors %>%
    as.data.frame() %>%
    rownames_to_column(var = "Gene") %>%
    dplyr::filter(. == Mod)

  a <- Annot
  A <- a %>%
    select(gene, name, product, Description) %>%
    group_by(gene) %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    )) %>%
    filter(gene %in% Col$Gene)
  
  A. <- A %>%
    filter(name != "<unknown>")
  
  Col$ID2 <- Col$Gene
  
  idx <- match(Col$Gene, A.$gene)
  Col$ID2 <- ifelse(is.na(idx), Col$Gene, A.$name[idx])
  
  Col$Gene <- Col$ID2
  Col$ID2 <- NULL

  selected_genes <- Col[, 1]

  colnames(Col) <- c("Gene", "Module")

  ColData <- setNames(Col$Module, Col$Gene)

  exportNetworkToCytoscape(TOMSUB,
    edgeFile = PATH1,
    nodeFile = PATH2,
    weighted = TRUE,
    threshold = THRLD,
    nodeNames = selected_genes,
    nodeAttr = ColData[selected_genes]
  )
}

################################
# Node attribute table for Cytoscape
NodAttr <- function(Colors, Annot, GMM, Mod, Vars, PATH) {
  Col <- Colors %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene") %>%
    dplyr::filter(. == Mod)

  selected_genes <- Col[, 1]

  nodAttr <- GMM %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene") %>%
    dplyr::filter(Gene %in% selected_genes) %>%
    dplyr::select(all_of(Vars)) %>%
    as.data.frame() %>%
    mutate(Size = 70)

  colnames(nodAttr) <- c("name", "kME", "Size")

  a <- Annot
  A <- a %>%
    select(gene, name, product, Description) %>%
    group_by(gene) %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    )) %>%
    filter(gene %in% nodAttr$name)
  
  A. <- A %>%
    filter(name != "<unknown>")
  
  nodAttr$ID2 <- nodAttr$name
  
  idx <- match(nodAttr$name, A.$gene)
  nodAttr$ID2 <- ifelse(is.na(idx), nodAttr$name, A.$name[idx])
  
  nodAttr$name <- nodAttr$ID2
  nodAttr$ID2 <- NULL



  write.table(nodAttr, PATH, sep = "\t", row.names = FALSE, quote = FALSE)
}

################################
#Colors <- tar_read(mergeColors)
#Mod <- "royalblue"
#Annot <- tar_read(AnnotationFile)
## Retrieve annotated list of modules' genes:
AnnotModules <- function(Colors, Mod, Annot){
  
  Col <- Colors %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Gene") %>%
    dplyr::filter(. == Mod) %>%
    select(Gene)
  
  a <- Annot
  A <- a %>%
    select(gene, name, product, Description) %>%
    group_by(gene) %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    )) %>%
    filter(gene %in% Col$Gene)
  
  X <- merge(Col, A, by.x = "Gene", by.y = "gene", all.x = TRUE)
  
}

################################
## Calculate Gene-Trait Significance
TraitSignificance <- function(TraitData, Exp.Matrix, R) {
  TraitVar <- as.data.frame(TraitData)

  geneTraitSignificance <- as.data.frame(cor(Exp.Matrix, TraitVar, use = "p")) %>%
    filter(abs(TraitData) >= R)

  names(geneTraitSignificance) <- paste("GS.", names(TraitVar), sep = "")

  return(geneTraitSignificance)
}
