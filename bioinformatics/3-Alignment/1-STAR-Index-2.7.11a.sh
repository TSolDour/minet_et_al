#!/bin/bash
#SBATCH -p long
#SBATCH -J Index-mpibr_sepoff_v1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=20G
#SBATCH -t 03:00:00
#SBATCH --error=index_%j.err
#SBATCH --output=index_%j.out

#Purge unused modules
module purge

#Load modules
module load star/2.7.11a

cd /shared/projects/rnasep/6-Alignement/rnasep2/GenomeBased/

#Defining variables
REFERENCE=/shared/projects/rnasep/Ref/mpibr_sepoff_v1.fa
GENOME_DIR=/shared/projects/rnasep/Ref/IndexSTAR
GTF=/shared/projects/rnasep/Ref/mpibr_sepoff_v1.gtf

STAR --runThreadN 10 --genomeSAindexNbases 12  --runMode genomeGenerate --genomeDir ${GENOME_DIR} --genomeFastaFiles ${REFERENCE} --sjdbGTFfile ${GTF} --sjdbOverhang 149
