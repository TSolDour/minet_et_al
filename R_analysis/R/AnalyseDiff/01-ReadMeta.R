###################################################

# Loading the data and metadata files

###################################################

################################
# Function to read the metadata file
ReadMeta <- function(file) {
  read.table(file, header = T, sep = "\t")
}

# To read the emapper file
ReadMapper <- function(file) {
  read_excel(file)
}

################################
# To build the annotation file. Obsolete, see script 00-Annotation.R instead.
BuildAnnotation <- function(fileA, fileB) {
  
  GOfile <- fileA
  
  GOfile$query <- sub("^([^.]+\\.[^.]+)\\..*$", "\\1", GOfile$query)
  
  GOfile <- GOfile %>%
    dplyr::select(query, COG_category, Description, Preferred_name, GOs, KEGG_ko, KEGG_Pathway, KEGG_Module, PFAMs) %>%
    unique() %>%
    group_by(query) %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ {
        vals <- unique(.x[.x != "-" & .x != "" & !is.na(.x)])  # retire NA, "-", vide
        if (length(vals) == 0) "-" else paste(vals, collapse = "/")
      }
    )) %>%
    dplyr::rename(gene = query)
  
  GOfile <- GOfile[-c(1:3), ]
  
  annot <- fileB
  
  Annot <- annot %>%
    dplyr::select(gene, product) %>%
    group_by(gene) %>%
    summarise(product = paste(sort(unique(product)), collapse = "/"))
  
  Annotation <- base::merge(Annot, GOfile, by="gene", all.x=TRUE)
  
  return(Annotation)
}

################################
# Function to tidy metadata file - used for WGCNA
tidyMeta <- function(meta) {
  Meta <- meta %>%
    filter(SampleName != "80") %>%
    filter(SampleName != "39")
    #mutate(SampleName = paste(SampleName, "-", Treatment))
    #filter(Treatment == A | Treatment == B | Treatment == C)
  rownames(Meta) <- Meta$SampleName

  return(Meta)
}

################################
# Function to tidy metadata file - used for DESeq2
tidyMeta2 <- function(meta) {
  Meta <- meta %>%
    mutate(SampleName = paste(SampleName, "-", Treatment))
  rownames(meta) <- meta$SampleName

  return(meta)
}
