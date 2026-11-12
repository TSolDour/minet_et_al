###################################################

# Running DESeq

###################################################

################################
# Preforming differential expression analysis
DEanalysis <- function(X, TEST, reduced = NULL) {
  
  if(TEST == "Wald"){
    res <- DESeq(X, test = "Wald")
  } else if(TEST == "LRT"){
    if (is.null(reduced)) {
      stop("Pour le test LRT, vous devez spécifier l'argument 'reduced'.")
    }
    res <- DESeq(X, test = "LRT", reduced = reduced)
  } else {
    stop("L'argument TEST doit avoir la valeur 'Wald' ou 'LRT'")
  }
  
  return(res)
  
}

################################
# Retrieving list of differentialy expressed genes
RetrieveDEG <- function(dds, NAME, THRESHOLD) {
  Res <- results(dds, name = NAME)

  Results <- Res %>% # This objet contains every genes influenced by  Treatment
    data.frame() %>%
    rownames_to_column(var = "ID") %>%
    arrange(padj) %>%
    filter(padj < 0.05 &
      abs(log2FoldChange) > THRESHOLD) %>%
    select(ID, log2FoldChange)


  return(Results)
}

################################
# Retrieving list of interacting genes - finaly not used in the project.
RetrieveINT <- function(dds, NAME, THRESHOLD) {
  Res <- results(dds, name = NAME)
  
  Results <- Res %>% # This objet contains every genes influences Treatment
    data.frame() %>%
    rownames_to_column(var = "ID") %>%
    arrange(padj) %>%
    filter(padj < 0.05) %>%
    select(ID, padj)
  
  
  return(Results)
}

################################
# To add annotation to the differentialy expressed genes
AnnotDE <- function(DEG, Annot, fromWald=TRUE) {
  a <- Annot
  b <- DEG
  
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
    filter(gene %in% b$ID)
  
  A. <- A %>%
    filter(name != "<unknown>")
  
  b$ID2 <- b$ID
  
  idx <- match(b$ID, A.$gene)
  b$ID2 <- ifelse(is.na(idx), b$ID, A.$name[idx])
  
  if(fromWald == TRUE){
    B <- b %>%
      dplyr::select(ID2, log2FoldChange) %>%
      dplyr::rename(gene = ID2)
  } else {
    B <- b %>%
    dplyr::select(ID2, padj) %>%
    dplyr::rename(gene = ID2)
  }
  
  
  return(B)
#  valid_mapping <- Annot[!is.na(Annot$ProteinCode) & Annot$ProteinCode != "-", ]
#  transcript_to_preferred <- setNames(valid_mapping$ProteinCode, valid_mapping$Transcript)

#  Res <- DEG
#  Res$ID <- ifelse(Res$ID %in% names(transcript_to_preferred),
#    transcript_to_preferred[Res$ID],
#    Res$ID
#  )

#  return(Res)
}