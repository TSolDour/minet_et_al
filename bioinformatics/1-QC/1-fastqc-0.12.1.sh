#!/bin/bash
#SBATCH -p fast
#SBATCH --time=04:00:00
#SBATCH	-J FrstQCReport
#SBATCH -e error_out
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH	--mail-user=t.soldourdin@gmail.com
#SBATCH --mail-type=END

# Purge any previous module
module purge

# Load fastqc and multiqc environments
module load fastqc/0.12.1
module load multiqc/1.26

cd /shared/projects/rnasep/1-QC

# perform FASTQC Quality Check on all fastq.gz files
fastqc -o /shared/projects/rnasep/1-QC/ /shared/projects/rnasep/RAW/*.fastq.gz

# multiqc report
multiqc .
