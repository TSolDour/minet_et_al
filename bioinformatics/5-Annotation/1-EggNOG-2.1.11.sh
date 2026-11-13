#!/bin/bash
#SBATCH -p long
#SBATCH -J EggNOG
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=20G
#SBATCH -t 20:00:00
#SBATCH --error=eggnog_%j.err
#SBATCH --output=eggnog_%j.out

#Purge unused modules
module purge

#Load modules
module load eggnog-mapper/2.1.11

cd /shared/projects/rnasep/5-Annotation/

export EGGNOG_DATA_DIR=/shared/projects/rnasep/5-Annotation/EggNOGdata/

#Defining variables
QUERY=mpibr_sepoff_v1_proteins.faa
OUTDIR=/shared/projects/rnasep/5-Annotation/output

emapper.py -i ${QUERY} -m diamond --evalue 0.001 --tax_scope 33213 --target_orthologs one2one --go_evidence all --output_dir ${OUTDIR} -o mpibr_sepoff_v1 --excel --cpu 20 
