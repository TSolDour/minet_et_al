#####################################

# Export graphics

####################################

################################
# Exporting plot
PlotExport <- function(file, plot, W=NA, H=NA, U) {
  ggsave(file, plot, dpi = 300, width=W, height=H, units=U)
}