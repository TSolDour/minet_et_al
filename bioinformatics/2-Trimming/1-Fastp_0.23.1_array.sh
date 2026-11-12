#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=fastpTestArray_2
#SBATCH -t 05:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=15
#SBATCH --array=1-28

#Load modules
module load fastp/0.23.1

config=Fastp-Config_2.txt
READ1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print$2}' $config)
READ2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print$3}' $config)
OUTPUT=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print$4}' $config)
ADAPTER=TruSeq3-PE-2.fa

cd /shared/projects/rnasep/2-Trimming/rnasep2/fastp/

fastp --thread 15 -i ${READ1} -o ${OUTPUT}_R1_pe.fastq.gz -I ${READ2} -O ${OUTPUT}_R2_pe.fastq.gz --qualified_quality_phred 30 --unqualified_percent_limit 10 --detect_adapter_for_pe --adapter_fasta ${ADAPTER} -g -x --length_required 120
