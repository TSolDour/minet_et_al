###########################################################################################
#
# Script dedicated to the display of PCA as preliminary visualization
#
#
###########################################################################################

makePCA <- function(STAR, top){
  vsd <- vst(STAR)
  pcaData <- plotPCA(vsd, intgroup=c("Treatment"), returnData=TRUE, top) # Performing PCA using the 1000 more variable genes
  pcaData$Condition <- ifelse(pcaData$Treatment == "CT8.1", "Control",
                              ifelse(pcaData$Treatment == "CT7.7", "pCO2",
                                     ifelse(pcaData$Treatment == "Hg8.1", "Hg",
                                            ifelse(pcaData$Treatment == "Hg7.7", "Hg+pCO2", NA))))
  percentVar <- round(100 * attr(pcaData, "percentVar")) # Retrieving the per component variance information
  p <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) + # Plotting the PCA
    geom_point(size=5) +
    stat_ellipse(level=0.95) +
    scale_color_manual(values=c("#50aa3a", "#30bcd7", "#e97132", "#ece87c"),
                       labels=c("Ctrl", expression(Hg[Low]), expression(CO[2]+Hg[Low]), expression(CO[2]))) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_bw() +
    theme(axis.title = element_text(size=18),
          legend.title = element_text(size=18),
          legend.position = "top",
          axis.text = element_text(size=12),
          legend.text = element_text(size=12)) +
    labs(color="Condition") +
    coord_fixed()
  
  return(p)
}