####################################################
#
#   Constructing module eigengenes                 #
#
####################################################

################################
# Build dissimilarity matrix
TOMdissim <- function(ADJACENCY) {
  TOM <- TOMsimilarity(ADJACENCY) # Similarity
  NegTOM <- 1 - TOM # Dissimilarity

  rownames(NegTOM) <- rownames(ADJACENCY)
  colnames(NegTOM) <- colnames(ADJACENCY)

  return(NegTOM)
}

################################
# Clustering of genes
TreeGene <- function(matrix) {
  geneTree <- hclust(as.dist(matrix), method = "average")
}

################################
# Building modules
BuildModules <- function(DENDRO, matrix) {
  cutreeDynamic(dendro = DENDRO, distM = matrix, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
}

################################
MEColors <- function(Modules, Matrix) {
  NET <- Modules
  lab <- labels2colors(NET$colors)

  names(lab) <- colnames(Matrix)

  return(lab)
}

################################
MEdissim <- function(Exp.Matrix, ModColors) {
  MElist <- moduleEigengenes(Exp.Matrix, colors = ModColors)
  ME.dissimilarity <- 1 - cor(MElist$eigengenes, use = "complete") # Calculate eigengene dissimilarity

  return(ME.dissimilarity)
}

################################
mergeModules <- function(Exp.Matrix, ModColors) {
  mergeCloseModules(Exp.Matrix, ModColors, cutHeight = .25)
}

################################
retrieveMergedME <- function(obj) {
  a <- obj$newMEs

  return(a)
}

################################
retrieveMergedColors <- function(obj, ModColors) {
  a <- obj$colors

  names(a) <- names(ModColors)
  return(a)
}


################################
# Plot dendogram & modules
plotDendroColors <- function(GENETREE, ModCol, MergeCol, path) {
  Tree <- GENETREE$dendrograms[[1]]
  Col1 <- ModCol[GENETREE$blockGenes[[1]]]
  Col2 <- MergeCol[GENETREE$blockGenes[[1]]]

  png(path)

  plotDendroAndColors(Tree, cbind(Col1, Col2),
    c("Original Module", "Merged Module"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and module colors for original and merged modules"
  )
  dev.off()
  path
}
