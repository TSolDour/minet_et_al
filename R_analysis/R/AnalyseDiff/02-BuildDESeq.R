###################################################

# Building DESeq file

###################################################

################################
# Build the DESeq operational file:
BuildDESeq <- function(ST, DIR, DES, relevel_treatment=FALSE, factor_treatment=FALSE, relevel_cluster=FALSE, factor_cluster=FALSE ) {
 X <- DESeqDataSetFromHTSeqCount(
    sampleTable = ST,
    directory = DIR,
    design = DES
  )
 
 X <- estimateSizeFactors(X)
 X <- estimateDispersions(X)
 
 idx <- rowSums(counts(X,normalized=TRUE) >= 10 ) >= 2  # check use normalize
 X <- X[idx,]
 
 if(factor_treatment == TRUE){
   X$Treatment <- as.factor(X$Treatment)
 }
 
 if(relevel_treatment == TRUE){
   X$Treatment <- relevel(X$Treatment, ref="CT8.1")

 }
  
  if(factor_cluster == TRUE) {
    X$Cluster <- as.factor(X$Cluster)
  }

  if(relevel_cluster == TRUE) {
    X$Cluster <- relevel(X$Cluster, ref="Cluster3")
  }
 
 return(X)
}

################################
# From deseq object to normalized counts matrix. Used for visualization.
StarToExp <- function(STAR) {
  star <- STAR
   star <- estimateSizeFactors(star)
   star <- estimateDispersions(star)

   idx <- rowSums(counts(star,normalized=TRUE) >= 10 ) >= 2
   star <- star[idx,]

  vsd <- vst(star)

  assay(vsd) %>%
    t()
}
