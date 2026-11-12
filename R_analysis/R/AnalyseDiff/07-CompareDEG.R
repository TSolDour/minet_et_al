###########################################################################################
#
# Script dedicated to the display of Venn diagram to compare DEGs between three conditions
#
#
###########################################################################################

################################
# Build Venn data from the DESeqDEG outputw
BuildVennData <- function(A, B, C) {
  x <- list(
    A$ID,
    B$ID,
    C$ID
  )

  return(x)
}


################################
# Displaying the VennDiagram
displayVenn <- function(Data) {
  
  
  names(Data) <- c("CO2+HgLow", "HgLow", "CO2")
  
  
  
a <-  ggvenn::ggvenn(
    Data,
    columns = c("CO2+HgLow", "HgLow", "CO2"),
    show_stats = "c",
    fill_color = c("#e97132", "#30bcd7", "#ece87c"),
    fill_alpha=0.7,
    auto_scale = FALSE,
    stroke_size = 0.9,
    stroke_color = "white",
    set_name_size = 10,
    text_size = 10,
    padding=0.2
  )

 return(a)
  
}