#####################################

# Preparing files for GOMWU analysis

####################################

################################
# Gene to GO list
goMWU <- function(dds, NAME, GO) {
  Res <- results(dds, name = NAME) %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "ID") %>%
    dplyr::filter(pvalue != "NA")

GO %>%
  data.frame() %>%
  dplyr::filter(GOs != "NA" & GOs != "-" & gene %in% Res$ID) %>%
  dplyr::select(gene, GOs) %>%
  tidyr::separate_rows(GOs, sep=",") %>%
  unique() %>%
  #dplyr::mutate(GOs = stringr::str_replace_all(GOs, ",", ";")) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(across(
    .cols = everything(),
    .fns = ~ {
      vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
      if (length(vals) == 0) "-" else paste(vals, collapse = ";")
    }
  ))
}

################################
# Complete list of analyzed genes
geneMWU <- function(dds, NAME, THRESHOLD, Ref) {
  Res <- results(dds, name = NAME) %>% # This objet contains all the genes influences by Treatment
    data.frame() %>%
    tibble::rownames_to_column(var = "ID") %>%
    plyr::arrange(padj) %>%
    dplyr::filter(pvalue != "NA") %>%
    dplyr::select(ID, log2FoldChange) %>%
    dplyr::filter(ID %in% Ref$gene)
  
  return(Res)
}

################################
#ALL <- tar_read(MMall)
# Complete list of WGCNA analysed genes with value of belonging (1, 0)
moduleMWU <- function(ALL, module, value) {
  Gene <- rownames(ALL) %>%
    as.data.frame() %>%
    dplyr::rename(ID = ".") %>%
    dplyr::mutate(Value = ifelse(ID %in% rownames(module),
      ALL[[value]],
      0))

  colnames(Gene) <- c("ID", "Value")

  return(Gene)
}

################################
# Module to GO list
ModuleGO <- function(GO, ALL) {
  
  Gene <- rownames(ALL) %>%
    as.data.frame() %>%
    dplyr::rename(ID = ".")
  
  GO %>%
    data.frame() %>%
    dplyr::filter(GOs != "NA" & GOs != "-" & gene %in% Gene$ID) %>%
    dplyr::select(gene, GOs) %>%
    tidyr::separate_rows(GOs, sep=",") %>%
    unique() %>%
    #dplyr::mutate(GOs = stringr::str_replace_all(GOs, ",", ";")) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = ";")
      }
    ))
  
}


################################
# Export Go table
ExportGO <- function(X, FILE) {
  write.table(X, file = FILE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

################################
# Export Gene table
ExporteGENE <- function(X, FILE) {
  write.csv(X, file = FILE, sep = ",", quote = FALSE, row.names = FALSE)
}

################################
# Gene to KOG class list - don't used
KOG <- function(KO, CLASS) {
  class <- read.delim(CLASS) %>%
    as.data.frame()

  kog <- KO %>%
    data.frame() %>%
    dplyr::filter(COG_category != "NA" & COG_category != "-") %>%
    dplyr::select(Transcript, COG_category) %>%
    dplyr::mutate(COG_category = sapply(COG_category, function(x) paste(strsplit(COG_category, "")[[1]], collapse = ","))) %>%
    tidyr::separate_rows(COG_category, sep = ",")

  kog$COG_category <- ifelse(kog$COG_category %in% class$KOG_letter,
    class$KOG_class,
    kog$COG_category
  )
  return(kog)
}

################################
# Running KOGMWU - don't used.
KOGMWU <- function(GeneToKOG, ModuleGeneToKOG) {
  kog.mwu(ModuleGeneToKOG, GeneToKOG)
}