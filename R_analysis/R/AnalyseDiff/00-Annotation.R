
#### Build Annotation File from raw outputs ####

BuildAnnot <- function(path1, path2) {
  
  GOfile <- as.data.frame(read_excel(path1))
  annot <- read.delim(path2)
  
  GOfile$query <- sub("^([^.]+\\.[^.]+\\.[^.]+)\\..*$", "\\1", GOfile$query)
  
  GOfile <- GOfile %>%
    select(query, COG_category, Description, Preferred_name, GOs, KEGG_ko, KEGG_Pathway, KEGG_Module, PFAMs) %>%
    unique() %>%
    group_by(query) %>%
    summarise(.groups="drop",
              across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    )) %>%
    dplyr::rename(mrna = query)
  
  GOfile <- GOfile[-c(1:3), ]
  
  annot <- annot %>%
    select(gene, mrna, name, product) %>%
    group_by(gene, mrna) %>%
    summarise(.groups="drop",
              across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    ))
  
  Annotation <- base::merge.data.frame(annot, GOfile, by="mrna", all=TRUE)
  
  return(Annotation)
}