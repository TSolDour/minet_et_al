# **Minet et al. — Unmasked in a changing ocean Acidification and mercury impair cuttlefish defence mechanisms**

**authors:** Antoine Minet,*, Thomas Sol Dourdin,*, Matthew A. Birk, Aude Andouche, Julien Bacque-Cazenave, Sébastien Baratte, Marina Fernandez-Declerck, Eric Pante, Paco Bustamante, Christelle Jozet-Alves, Thomas Lacoue-Labarthe

## Introduction:

This repository brings together the scripts and configuration files used to produce the analyses presented in the manuscript *"Unmasked in a changing ocean: Acidification and mercury impair cuttlefish defense mechanisms"*.

The data used for this project are available on **XXXX**


The project is structured in two main parts:

- **Bioinformatics processing** of transcriptomics data - scripts written for *SLURM* environment.

- **Post-processing analysis of transcriptomic data** and **analysis of behavioral data** carried out in R.

## Repository structure

```text
minet_et_al-main/
├── README.Rmd
├── bioinformatics/
│   ├── 1-QC/                 # Quality check FastQC / MultiQC
│   ├── 2-Trimming/           # Read trimming
│   ├── 3-Alignment/          # Genome indexation and reads alignment
│   └── 4-mRNAediting/        # mRNA editing estimation
└── R_analysis/
    ├── AnalyseDiff_targets.R # Main `targets` script of transcriptomic pipeline
    ├── Behavior_targets.R    # Main `targets` script for behavior pipeline
    ├── _targets.yaml         # Configuration file for `targets`
    ├── renv.lock             # R package versions and metadata
    ├── renv/                 # Scripts to initialize `renv`
    ├── configFiles/          # Some configuration files required to execute the pipeline
    ├── data/
    │   ├── analyseDiff/      # Counts files and mRNA editing data - to pick up from **XXXXX**
    │   └── behavior/         # Behavioral data - to pick up from **XXXXX**
    ├── R/
    │   ├── AnalyseDiff/      # Source scripts for transcriptomic pipeline
    │   └── Behavior/         # Source scripts for behavior pipeline
    ├── results/
    │   ├── AnalyseDiff/      # `targets` storage / transcriptomic outputs
    │   └── Behavior/         # `targets` storage / behavior outputs
    └── docs/                 # Some reference documentation (DESeq2, WGCNA, targets)
```

## Downloading the project:

Download the ZIP archive directly from GitHub, or from your terminal execute the following:

```bash
git clone <https://github.com/TSolDour/minet_et_al.git>
cd minet_et_al-main
```

## Bioinformatics data preprocessing

The `bioinformatics/` folder contains shell scripts describing a RNA-seq workflow.

| Step 					| Script 				| Function 								|
|---------------------------------------|---------------------------------------|-----------------------------------------------------------------------|
| Quality check 			| `1-QC/1-fastqc-0.12.1.sh`		| Execute FastQC on FASTQ files then MultiQC to aggregate reports.	|
| Trimming 				| `2-Trimming/1-Fastp_0.23.1_array.sh`	| Clean paired-end reads using fast, using a configuration file.	|
| Reference indexing			| `3-Alignment/1-STAR-Index-2.7.11a.sh`	| Build STAR index of the reference genome used to map the reads.	|
| Alignment 				| `3-Alignment/2-STAR-Align-2.7.11a.sh`	| Reads alignment using STAR.						|
| Retrieving coding genes counts	| `3-Alignment/filter.sh`		| Filter counts files with a coding genes ID list.			|
| mRNA editing 				| `4-mRNAediting/MASTER.sh`		| Quantify mRNA editing of previously known edited sites.		|

RNA-seq raw data used in this study are available from the European nucleotide archive (ENA) under the **accession number: PRJEB88225**.

## Preparing data for R pipeline

Paths are relative to the `R_analysis/` root folder.

### 1. Transcriptomic analysis (`AnalyseDiff`)

| Required file					| Retrieved from			| Location			| Usage										|
|-----------------------------------------------|---------------------------------------|-------------------------------|-------------------------------------------------------------------------------|
| `Design_Deseq.tsv`				| Available in the GitHub repository	| `configFiles/` 		| Design table for DESeq2 containing sample metadata.				|
| `*.counts`					| Available from XXXXX			| `data/analyseDiff		| Raw counts retrieved from STAR alignment and used as input for DESeq2.	|
| `mpibr_sepoff_v1.emapper.annotations.xlsx`	| Available in the GitHub repository	| `configFiles/`		| Functional annotation computed by EggNOG-mapper.				|
| `mpibr_sepoff_v1_annotations.tsv`		| Available in the GitHub repository	| `configFiles/`		| Functional annotation provided with the reference genome.			|
| `SampleTable.tsv`				| Available in the GitHub repository 	| `configFiles/`		| Metadata used in the mRNA editing part.					|
| `Soff_head_ES.csv`				| Available from XXXXX			| `data/analyseDiff/Other/`	| Editing levels used as input for the mRNA editing part.			|
| `editing_sites_Shoshan2021.txt`		| Available in the GitHub repository	| `configFiles/`		| Known mRNA editing sites.							|

The transcriptomic pipeline produces files to be used by the GO_MWU utils and process GO_MWU outputs. But GO_MWU itself is not included in the `targets` pipeline and the GO_MWU.R script should be independently executed. It is located in `R/AnalyseDiff/GOMWU` with its dependencies. 

### 2. Behavior analysis (`Behavior`)


| Required file		| Retrieved from	| Location		| Usage								|
|-----------------------|-----------------------|-----------------------|---------------------------------------------------------------|
| `Minet_et_al.xlsx`	| Available from XXXXX	| `data/behavior/`	| File containing all data used in the behavioral analysis	|

`Minet_et_al.xlsx`should contain 3 sheets: 

- `data2020`
- `data2021`
- `acp2020_2021`

## Setting up R environment

Open R and RStudio - or Positron - from the project's root folder `R_analysis`. If your work with RStudio, you can open the project by opening the `R_analysis.Rproj` file.

```r
## Install renv if necessary
install.packages("renv")

## Restore the project library stored in renv.lock
renv::restore()
```

The restoration can take some time. Accept package installation whenever suggested. If the installation fails, follow the instructions provided by renv to fix it.  

## Executing `{targets}` pipelines

Both transcriptomics and behavior pipelines are independents. They are managed by `{targets}`.

First, execute the *_targets.R script: `CMD + A` then `Ctrl + Enter`

Second, specify the configuration: `Sys.setenv(TAR_PROJECT="AnalyseDiff")` for the transcriptomic pipeline or `Sys.setenv(TAR_PROJECT="Behavior")` for the behavior pipeline.

Third, execute `targets::tar_make(script="AnalyseDiff_Targets.R")` to run the transcriptomic pipeline or `targets::tar_make(script="Behavior_targets.R") to ru the behavior pipeline.

Once the pipeline successfully ran, you can find the outputs in the `results/`folder. Figures used for the paper are labelled with the prefixe "Figure". 

To inspect an intermediate output, use the `targets::tar_read()` function with the name of a specific targets as argument.