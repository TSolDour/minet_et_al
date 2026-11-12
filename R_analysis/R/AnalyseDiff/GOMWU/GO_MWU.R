
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.


# Edit these to match your data file names: 
input="GeneToValue.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="GeneToGO.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previous runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
	perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
	largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
	smallest=5,   # a GO category should contain at least this many genes to be considered
	clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
#	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
#	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
#	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# ----------- Plotting results

#png("HgCO2_CC.png", width=9, height=10, units="in", res=300)
results=gomwuPlot(input,goAnnotations,goDivision,
          #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
          #absValue=0.001,
          absValue=1, # un-remark this if you are using log2-fold changes
          level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
          level2=0.01, # FDR cutoff to print in regular (not italic) font.
          level3=0.005, # FDR cutoff to print in large bold font.
          txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
          treeHeight=0.5, # height of the hierarchical clustering tree
          colors=c("#046C9A","#FD6467","#3B9AB2","#E6A0C4") # these are default colors, un-remark and change if needed
)
#devname=dev.off()

dim(results[[1]])

######## There are too many enriched GO terms
### Extracting representative GOs: 
pcut=1e-2
hcut=0.1 # To be adjusted

plot(results[[2]], cex=0.6)
abline(h=hcut, col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
View(bestGOs)
nrow(bestGOs)

### To produce a customized results dataframe.
library(tidyverse)
library(plyr)
library(ggplot2)

 # Retrieving and preparing dataframe results table for custome plot
res <- results[[1]] %>%
  data.frame() %>%
  tibble::rownames_to_column("GOterms") %>%
  mutate(Trend = ifelse(direction == 0, "Repressed", "Induced"))

# Function to split a col into 2 cols without using `separate`
split_into_two <- function(x, SEP) {
  parts <- strsplit(as.character(x), SEP, fixed = TRUE)[[1]]
  data.frame(A = parts[1], B = paste(parts[-1], collapse = " "), stringsAsFactors = FALSE)
}

separated_GOcol <- do.call(rbind, lapply(res$GOterms, SEP=" ", split_into_two))

res$GOterms = factor(separated_GOcol$B, levels = separated_GOcol$B)
res$Ratio = separated_GOcol$A

separated_ratio <- do.call(rbind, lapply(res$Ratio, SEP="/", split_into_two))

separated_ratio$A <- as.numeric(separated_ratio$A)
separated_ratio$B <- as.numeric(separated_ratio$B)
separated_ratio$Ratio = separated_ratio$A/separated_ratio$B

res$Ratio = separated_ratio$A/separated_ratio$B
res$GeneNumber = separated_ratio$A

res.best <- res %>%
  filter(GOterms %in% bestGOs$name)

# Add IDs to GOs:
library(GO.db)
go_map <- AnnotationDbi::select(
  GO.db,
  keys = keys(GO.db, keytype = "GOID"),
  columns = c("GOID", "TERM", "ONTOLOGY"),
  keytype = "GOID"
)

res.ok <- res %>%
  left_join(go_map, by = c("GOterms" = "TERM"))
res.ok.best <- res.best %>%
  left_join(go_map, by = c("GOterms" = "TERM"))

write.table(res.ok.best, "~/Documents/Recherche/Current/rnasep-LIENSs/rnasep/R_analyses/results/AnalyseDiff/files/HgCO2_BP_best_results.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")



