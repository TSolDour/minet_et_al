#####################################################
#
#   External Traits Matching                        #
#
#####################################################


################################
# Load external traits
loadTraits <- function(path, SEP, DEC, Exp.Matrix) {
  a <- read.delim(path, sep = SEP, dec = DEC)

  b <- a %>%
    filter(SampleName %in% rownames(Exp.Matrix))
  rownames(b) <- b$SampleName
  b[, -1]
}

################################
# Computing Module-Traits correlation
ModTraitCor <- function(Exp.Matrix, MEmerged, Traits) {
  nGenes <- ncol(Exp.Matrix)
  nSamples <- nrow(Exp.Matrix)
  module.trait.correlation <- cor(MEmerged, Traits, use = "p") # p for pearson correlation coefficient
  module.trait.Pvalue <- corPvalueStudent(module.trait.correlation, nSamples) # calculate the p-value associated with the correlation

  # Will display correlations and their p-values
  textMatrix <- paste(signif(module.trait.correlation, 2), "\n(",
    signif(module.trait.Pvalue, 1), ")",
    sep = ""
  )
  dim(textMatrix) <- dim(module.trait.correlation)

  return(textMatrix)
}

################################
# Plotting the matrix
plotcor <- function(Traits, MEmerged, textMatrix, path, WIDTH, HEIGHT, UNITS, RES) {
  module.trait.correlation <- cor(MEmerged, Traits, use = "p") # p for pearson correlation coefficient

  png(path, width = WIDTH, height = HEIGHT, units = UNITS, res = RES)
  par(mar = c(6, 8.5, 3, 1))
  # Display the correlation values within a heatmap plot
  labeledHeatmap(
    Matrix = module.trait.correlation,
    xLabels = names(Traits),
    yLabels = names(MEmerged),
    ySymbols = names(MEmerged),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.6,
    zlim = c(-1, 1),
    main = paste("Module-trait relationships")
  )
  dev.off()
  path
}
