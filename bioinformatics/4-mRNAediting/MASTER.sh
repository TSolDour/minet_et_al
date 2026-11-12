#!/bin/bash

#SBATCH -p long
#SBATCH --job-name=EditCalc
#SBATCH -t 20:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=40

n_CPUs=$SLURM_CPUS_PER_TASK

module load bowtie2/2.5.4
module load samtools/1.21
module load r/4.5.1

bowtie2-build Sepia_officinalis_ORFs.fa S_officinalis_Liscovitch2017	  # create a bowtie2 index from the transcriptome published by Liscovitch-Brauer et al. 2017 (Cell)

# loop through every sample

cat samples.txt | while read sample
do
	find /shared/projects/rnasep/2-Trimming/rnasep2/fastp/ -name "${sample}_*_pe.fastq.gz" | sort > tmp	       	# this will need changed probably to match the cleaned FASTQ file naming convention you have
	r1=$(head tmp -n 1)			# get the name of read 1
	r2=$(tail tmp -n 1)			# get the name of read 2
	rm tmp

	echo "========================${sample}==========================" >> bowtie_alignment_stats.log

	bowtie2 -x S_officinalis_Liscovitch2017 -p $n_CPUs \
		-1 $r1 \
		-2 $r2 \
		-S ${sample}.sam &>> bowtie_alignment_stats.log			# align reads to transcriptome with bowtie2

	samtools view -b -f 2 -@ $n_CPUs ${sample}.sam > ${sample}_mapped.bam		# retain only the reads mapped in proper pairs (those with flag 2)

	rm ${sample}.sam			# get rid of this big SAM file because it takes up lots of space

	samtools sort ${sample}_mapped.bam -o ${sample}_sorted.bam -@ $n_CPUs		# sort the mapped reads by position

	rm ${sample}_mapped.bam		# get rid of the unsorted BAM file
	
	samtools index ${sample}_sorted.bam -@ $n_CPUs

	samtools idxstats ${sample}_sorted.bam > ${sample}_sorted_stats.log

	samtools mpileup -f Sepia_officinalis_ORFs.fa -s ${sample}_sorted.bam -d 0 -Q 20 --no-output-ins --no-output-del --no-output-ends > ${sample}_pile.txt	# for every position in the transcriptome, count how many reads (with base call phred score >20) had what base call
done

Rscript calc_EL_from_mpileup.R			# take the (sample)_pile.txt output from samtools mpileup and calculate editing levels for every known editing site in the editing site catalog from Shoshan et al. 2021 (MBE)
