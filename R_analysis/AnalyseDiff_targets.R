# ==============================================================================
# Minet et al. - Unmasked in a changing ocean: Acidification and mercury impair cuttlefish defense mechanisms
# Main script of the targets pipeline managing transcriptomic analyses.
# Created by Thomas Sol Dourdin
# All figures used for the publication are exported with the prefix "Figure"
# ==============================================================================

rm(list = ls())
# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(visNetwork)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("readxl", "tidyverse", "DESeq2", "ggplot2", "pheatmap",
               "pcaExplorer", "RColorBrewer", "VennDiagram", "wesanderson", "ggvenn",
               "WGCNA", "pbapply", "circlize", "dplyr") # Packages that your targets need for their tasks.
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source(here::here("R", "AnalyseDiff", "00-Annotation.R"))
tar_source(here::here("R", "AnalyseDiff", "01-ReadMeta.R"))
tar_source(here::here("R", "AnalyseDiff", "02-BuildDESeq.R"))
tar_source(here::here("R", "AnalyseDiff", "02bis-PCA.R"))
tar_source(here::here("R", "AnalyseDiff", "03-DEanalysis.R"))
tar_source(here::here("R", "AnalyseDiff", "04-PlotDE.R"))
tar_source(here::here("R", "AnalyseDiff", "05-PlotExport.R"))
tar_source(here::here("R", "AnalyseDiff", "06-GeneMWU.R"))
tar_source(here::here("R", "AnalyseDiff", "07-CompareDEG.R"))
tar_source(here::here("R", "AnalyseDiff", "08-PlotGOMWU.R"))
tar_source(here::here("R", "AnalyseDiff", "09-CoExpNetwork.R"))
tar_source(here::here("R", "AnalyseDiff", "10-ModuleConstruction.R"))
tar_source(here::here("R", "AnalyseDiff", "11-MatchExtTraits.R"))
tar_source(here::here("R", "AnalyseDiff", "12-RetrieveHubGenes.R"))
tar_source(here::here("R", "AnalyseDiff", "13-RNAediting.R"))

list(
  # ============================================================================
  # Diff. analysis: Preparing data and DESeq2 object HgCO2 vs CT 
  # ============================================================================
  tar_target(AnnotationFile, BuildAnnot(here::here("configFiles", "mpibr_sepoff_v1.emapper.annotations.xlsx"), here::here("configFiles", "mpibr_sepoff_v1_annotations.tsv"))), # Import & read Annotation file 
  ################################################################################################################################################
  tar_target(file, here::here("configFiles", "Design_Deseq.tsv")),                                                                                # Script 1
  tar_target(meta, ReadMeta(file)),                                                                                                              # Read and prepare SampleTable
  tar_target(sampleTable1, tidyMeta2(meta)),                                                                                                     #
  ################################################################################################################################################
  tar_target(directory, here::here("data", "analyseDiff")),                                                                                      # Script 2
  tar_target(star1, BuildDESeq(ST = sampleTable1, DIR = directory, DES=~Treatment, relevel_treatment = TRUE, factor_treatment = TRUE)),          # Build the DESeq object from counts files
  
  # ============================================================================
  # Diff. analysis: PCA plot - 5000 most variable genes 
  # ============================================================================
  tar_target(ACP, makePCA(star1, 5000)),                                                                                                         # Script 2bis
  tar_target(ExportACP, PlotExport(here::here("results", "Analysediff", "figures", "Figure5a_pcaPlot.png"), ACP, H = 5, W = 6, U="in")),

  # ============================================================================
  # Diff. analysis: HgCO2 vs CT - Wald test
  # ============================================================================
  tar_target(dds1, DEanalysis(star1, "Wald")),                                                                                                   # Script3
  tar_target(DEG1, RetrieveDEG(dds1, "Treatment_Hg7.7_vs_CT8.1", 0.5)),                                                                          # Running DE analysis and retrieving DEGs
  tar_target(DEG1annot, AnnotDE(DEG1, AnnotationFile)),                                                                                          # 475 DEGs
  ################################################################################################################################################
  tar_target(Heat1, HeatDEG(dds1, X=c("5","6","8","9","37","38","80","28","30","31","32","68"), DEG1, AnnotationFile, 2, 2, 7, rows = FALSE)),   # Scripts 4 & 5
  tar_target(ExportHeat1, PlotExport(here::here("results", "Analysediff", "figures", "Figure5d_HgCO2_heatmap.png"), Heat1, W=5, H=7, U="in")),   #
  ################################################################################################################################################
  tar_target(GOforMWU1, goMWU(dds1, "Treatment_Hg7.7_vs_CT8.1", AnnotationFile)),                                                                #
  tar_target(GENEforMWU1, geneMWU(dds1, "Treatment_Hg7.7_vs_CT8.1", 0.5, GOforMWU1)),                                                            # Script 6
  tar_target(GOTable, ExportGO(GOforMWU1, "R/AnalyseDiff/GOMWU/GeneToGO.tab")),                                                                  # Build & Export Go & Gene files to be used by GOMWU (external scripts, cf. https://github.com/z0on/GO_MWU)
  tar_target(GENEtable1, ExporteGENE(GENEforMWU1, here::here("R", "AnalyseDiff", "GOMWU", "GeneToValue.csv"))),                                  #
  
  # ============================================================================
  # Diff. analysis: Hg vs CT - Wald test
  # ============================================================================
  tar_target(DEG2, RetrieveDEG(dds1, "Treatment_Hg8.1_vs_CT8.1", 0.5)),                                                                          # Script3
  tar_target(DEG2annot, AnnotDE(DEG2, AnnotationFile)),                                                                                          # Retrieving DEGs - 8 DEGs
  ################################################################################################################################################
  tar_target(Heat2, HeatDEG(dds1, X=c("5","6","8","9","37","38","80","19","20","21","24","25","55","56"), DEG2, AnnotationFile, 2, 2, 7, rows=FALSE)), # Scripts 4 & 5
  tar_target(ExportHeat2, PlotExport(here::here("results", "Analysediff", "figures", "Hg_heatmap.png"), Heat2, W=5, H=3, U="in")),                     # Build & export heatmap of DEGs
  # Very few DEGs so don't need to go further on functional enrichement analysis.
  
  # ============================================================================
  # Diff. analysis: CO2 vs CT - Wald test
  # ============================================================================
  tar_target(DEG3, RetrieveDEG(dds1, "Treatment_CT7.7_vs_CT8.1", 0.5)),                                                                          # Script3
  tar_target(DEG3annot, AnnotDE(DEG3, AnnotationFile)),                                                                                          # Retrieving DEGs - 30 DEGs
  ################################################################################################################################################
  tar_target(Heat3, HeatDEG(dds1, X=c("5","6","8","9","37","38","80","10","11","12","13","46","48","85","86"), DEG3, AnnotationFile, 2, 2, 7, rows=FALSE)), # Scripts 4 & 5
  tar_target(ExportHeat3, PlotExport(here::here("results", "Analysediff", "figures", "CO2_heatmap.png"), Heat3, W=5, H=3, U="in")),                         # Build & export heatmap of DEGs
  # Too few DEGs so don't need to go further on functional enrichement analysis.
  
  # ============================================================================
  # Diff. analysis: Comparing DEGs from all contrasts
  # ============================================================================
  tar_target(SumDEGplot, SumDEG(DEG1, DEG2, DEG3)),                                                                                              # Script 4
  tar_target(SumDEGExport, PlotExport(here::here("results", "Analysediff", "figures", "Figure5c_SumDEG.png"), SumDEGplot, W=5, H=5, U="in")),    #
  ################################################################################################################################################
  tar_target(VennData, BuildVennData(DEG1, DEG2, DEG3)),                                                                                         #
  tar_target(VennDiagram, displayVenn(VennData)),                                                                                                # Script 7
  tar_target(ExportVenn, PlotExport(here::here("results", "AnalyseDiff", "figures", "Figure5b_VennDEGS.png"), VennDiagram, W=5, H=5, U="in")),   # 
  ################################################################################################################################################
  tar_target(file2, here::here("results", "AnalyseDiff", "files", "HgCO2_BP_best_results.txt")),                                                             #
  tar_target(HgCO2_GOMWU_BP, loadGOMWU(file2)),                                                                                                              # Script 8
  tar_target(HgCO2_BP_Plot, PlotGOMWU(HgCO2_GOMWU_BP)),                                                                                                      # Computing/Plotting GOMWU results
  tar_target(ExportHgCO2_BP_Plot, PlotExport(here::here("results", "AnalyseDiff", "figures", "Figure5e_HgCO2_BP.png"), HgCO2_BP_Plot, W=7, H=6.5, U="in")),  #
  tar_target(file2a, here::here("results", "AnalyseDiff", "files", "HgCO2_BP_all_results.txt")),                                                             #
  tar_target(HgCO2_GOMWU_BP.all, loadGOMWU(file2a)),                                                                                                         #
  tar_target(gene2go.DEG, gomwu2deg(DEG1, AnnotationFile, GOforMWU1,HgCO2_GOMWU_BP.all)),                                                                    # find the DEGs involved in GOMWU signal (n=52).
  
  # ============================================================================
  # Diff. analysis: Comparing behavioral clusters - Wald test
  # ============================================================================
  tar_target(star2, BuildDESeq(ST = sampleTable1, DIR = directory, DES=~Cluster, relevel_cluster = TRUE, factor_cluster = TRUE)), # Build the DESeq object from counts files
  ################################################################################################################################################
  tar_target(dds2, DEanalysis(star2, "Wald")),                                                                                                   # Script3
  tar_target(DEG4, RetrieveDEG(dds2, "Cluster_Cluster1_vs_Cluster3", 0.3)),                                                                      # Running DE analysis and retrieving DEGs
  tar_target(DEG4annot, AnnotDE(DEG4, AnnotationFile)),                                                                                          #
  ################################################################################################################################################
  tar_target(DEG5, RetrieveDEG(dds2, "Cluster_Cluster2_vs_Cluster3", 0.3)),
  ## No DEGs
  
  # ============================================================================
  # Weighted Gene Coexpression Network Analysis (WGCNA)
  # ============================================================================
  tar_target(meta2, tidyMeta(meta)),                                                                                               #
  #tar_target(star, BuildDESeq(ST = meta2, DIR = directory, DES=~Treatment, relevel_treatment = FALSE)),                           # Script 09
  tar_target(exp, StarToExp(star1)),                                                                                               #
  tar_target(ExprFilt, ExpFilt(exp)),                                                                                              # To keep top 5000 more variable genes                                                         #
  tar_target(check, checkSamples(ExprFilt)),                                                                                       # If TRUE -> Don't need to remove outliers (then forget the three next targets).                                                    #
  #tar_target(sampleTree, findOut(ExprFilt)),                                                                                      # 
  #tar_target(PlotTree, plotTree("results/AnalyseDiff/figures/SampleTree.png", sampleTree, H = 64), format = "file"),              # Building coexpression network
  #tar_target(ExpFiltOut, rmOut(sampleTree, ExprFilt, H = 59)),                                                                    #
  tar_target(spt, softThres(ExprFilt)),                                                                                            #
  tar_target(R2, plotR2(here::here("results", "AnalyseDiff", "figures", "R2.png"), spt), format = "file"),                         #
  tar_target(meanConnect, plotConnect(here::here("results", "AnalyseDiff", "figures", "meanConnect.png"), spt), format = "file"),  #
  # tar_target(Adj, Adjacency(ExprFilt, 8)),                                                                                       # Adj is the weighted gene co-expression network and contains 23,106 nodes (genes).
  tar_target(net, moduleConstruct(ExprFilt, 7)),                                                                                   # Automatic construction of modules - can avoid the three next targets, go line 154.
  #################################################################################################################################################################
  # tar_target(TOM.dissim,TOMdissim(Adj)),                                                                                                                          #
  # tar_target(geneTree, TreeGene(TOM.dissim)),                                                                                                                     #
  # tar_target(ME, BuildModules(geneTree, TOM.dissim)),                                                                                                             #
  tar_target(ModuleColors, MEColors(net, ExprFilt)),                                                                                                                #
  tar_target(ME.dissim, MEdissim(ExprFilt, ModuleColors)),                                                                                                          # Script 10
  tar_target(merged, mergeModules(ExprFilt, ModuleColors)),                                                                                                         # Constructing Eigengene modules
  tar_target(mergedMEs, retrieveMergedME(merged)),                                                                                                                  #
  tar_target(mergeColors, retrieveMergedColors(obj = merged, ModuleColors)),                                                                                        #
  tar_target(DendroColors, plotDendroColors(net, ModuleColors, mergeColors, here::here("results", "AnalyseDiff", "figures", "dendroColors.png")), format = "file"), #
  ########################################################################################################################################################################
  tar_target(file3, here::here("data", "analyseDiff", "Other", "ExternalTraits.csv")),                                                                                                        #
  tar_target(ExtTraits, loadTraits(path = file3, SEP = ";", DEC = ",", Exp.Matrix = ExprFilt)),                                                                                               # Script 11
  tar_target(CorrelationMatrix, ModTraitCor(ExprFilt, mergedMEs, ExtTraits)),                                                                                                                 # Matching eigengene modules and traits
  tar_target(Matrix, plotcor(ExtTraits, mergedMEs, CorrelationMatrix, here::here("results", "AnalyseDiff", "figures", "Module-Traits-Correlation.png"), 12, 10, "in", 300), format = "file"), #
  ########################################################################################################################################################################
  tar_target(MMall, ModuleMembership(mergedMEs, ExprFilt)),                                                                                                              #
  tar_target(HubGeneGreenyellow, HubGenes(mergeColors, "greenyellow", MMall, c("MMgreenyellow"), Variable = MMgreenyellow, Percent = 0.90, AnnotationFile)),             # Script 12 - Retrieving Hub Genes for each module eigengene of interest
  tar_target(HubGeneRoyalblue, HubGenes(mergeColors, "royalblue", MMall, c("MMroyalblue"), Variable = MMroyalblue, Percent = 0.95, AnnotationFile)),                     #
  #################################################################################################################################################################################################################################################
  tar_target(TOM_subGreenyellow, TOMsub(mergeColors, ExprFilt, "greenyellow", 5)),                                                                                                                                                                #
  tar_target(exportGreenyellowCytoscape, ExportTOMtoCyto(TOM_subGreenyellow, AnnotationFile, "greenyellow",                                                                                                                                       #
    PATH1 = here::here("results", "analyseDiff", "files", "WGCNA", "GreenyellowCytoscapeEdges.txt"),                                                                                                                                              # Export TOM for blue module
    PATH2 = here::here("results","analyseDiff", "files", "WGCNA", "GreenyellowCytoscapeNodes.txt"), Colors = mergeColors, THRLD = 0.15), format = "rds"),                                                                                         #
  tar_target(GreenyellowNodeAttr, NodAttr(mergeColors, AnnotationFile, MMall, "greenyellow", Vars = c("Gene", "MMgreenyellow"), PATH = here::here("results", "analyseDiff", "files", "WGCNA", "GreenyellowNodeAttributes.txt")), format = "rds"), #
  tar_target(GreenyellowAnnot, AnnotModules(mergeColors, "greenyellow", AnnotationFile)),                                                                                                                                                         #
  #################################################################################################################################################################################################################################################
  tar_target(TOM_subRoyalblue, TOMsub(mergeColors, ExprFilt, "royalblue", 5)),                                                                                                                                                                    #
  tar_target(exportRoyalblueCytoscape, ExportTOMtoCyto(TOM_subRoyalblue, AnnotationFile, "royalblue",                                                                                                                                             #
                                                         PATH1 = here::here("results", "analyseDiff", "files", "WGCNA", "RoyalblueCytoscapeEdges.txt"),                                                                                           # Export TOM for blue module
                                                         PATH2 = here::here("results", "analyseDiff", "files", "WGCNA", "RoyalblueCytoscapeNodes.txt"), Colors = mergeColors, THRLD = 0.15), format = "rds"),                                     #
  tar_target(RoyalblueNodeAttr, NodAttr(mergeColors, AnnotationFile, MMall, "royalblue", Vars = c("Gene", "MMroyalblue"), PATH = here::here("results", "analyseDiff", "files", "WGCNA", "RoyalblueNodeAttributes.txt")), format = "rds"),         #
  tar_target(RoyalblueAnnot, AnnotModules(mergeColors, "royalblue", AnnotationFile)),                                                                                                                                                             #
  #################################################################################################################################################################################################################################################
  tar_target(moduleGO, ModuleGO(AnnotationFile, MMall)),                                                                                    #
  tar_target(moduleGOTable, ExportGO(moduleGO, here::here("R", "AnalyseDiff", "GOMWU", "moduleGeneToGO.tab"))),                             #
  tar_target(GreenyellowForMWU, moduleMWU(MMall, TOM_subGreenyellow, value = "MMgreenyellow")),                                             # Export genes & GOs for GOMWU - Blue module   
  tar_target(GreenyellowGENEtable, ExporteGENE(GreenyellowForMWU, here::here("R", "AnalyseDiff", "GOMWU", "GreenyellowGeneToValue.csv"))),  #
  tar_target(RoyalblueForMWU, moduleMWU(MMall, TOM_subRoyalblue, value = "MMroyalblue")),                                                   # Export genes & GOs for GOMWU - Blue module   
  tar_target(RoyalblueGENEtable, ExporteGENE(RoyalblueForMWU, here::here("R", "AnalyseDiff", "GOMWU", "RoyalblueGeneToValue.csv"))),        #
  
  # ============================================================================
  # mRNA editing analysis
  # ============================================================================
  tar_target(file6, here::here("data", "analyseDiff", "Other", "Soff_head_ES.csv")),                                                          #
  tar_target(file7, here::here("configFiles", "SampleTable.tsv")),                                                                            #
  tar_target(d_all, readES(file6)),                                                                                                           #
  tar_target(d_sample, ReadMeta(file7)),                                                                                                      #
  tar_target(d_sampleOK, tidySample(d_sample)),                                                                                               #
  tar_target(sampleTreatment, SampTreat(d_sampleOK)),                                                                                         #
  tar_target(knownES, retrieveKnownES(here::here("configFiles","editing_sites_Shoshan2021.txt"))),                                            #
  tar_target(CovHist, displayCovHist(d_all)),                                                                                                 #
  tar_target(dCov, rmLowCov(d_all, 2, 20)),                                                                                                   # Script 13
  tar_target(final.mat, buildMat(dCov, d_sampleOK)),                                                                                          #
  tar_target(dEL, rundEL(final.mat, dCov, 0.1, 0.05)),                                                                                        #
  tar_target(dELtornado, displayTornado(dEL)),                                                                                                #
  tar_target(Export_dELtornado, PlotExport(here::here("results", "AnalyseDiff", "figures", "dELtornado.png"), dELtornado, W=4, H=5, U="in")), #
  tar_target(dELcircos, displayCircos(dEL, "HgOA_vs_CT", "results/AnalyseDiff/figures/RNAeditingChordHgCO2rep.png")),                         #
  tar_target(EI, editingIndex(d_all, d_sampleOK)),                                                                                            # mean(tar_read(EI)$EI)=0.24; sd(tar_read(EI)$EI)=0.004
  tar_target(testEI, anovaEI(EI)),                                                                                                            #
  tar_target(EIplot, plotEI(EI)),                                                                                                             #
  tar_target(EIplotExport, PlotExport(here::here("results", "Analysediff", "figures", "Figure5f_EI.png"), EIplot, W=4, H=2.5, U="in"))        #             #
)

# Sys.setenv(TAR_PROJECT="AnalyseDiff")

# tar_manifest(fields=command)

# tar_visnetwork(physics=TRUE)

# tar_make(ExportTan_BP_Plot)

# tar_read(Tan_BP_Plot)