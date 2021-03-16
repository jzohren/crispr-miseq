#!/bin/bash

# example execution: ./crispr.sh A15 info_file.csv true genome.fasta fastqDir outDir
# where:
# 'A15' is the sample name
# 'info_file.csv' contains information about the guide sequence
# 'true/false' indicates whether the mapping step should be executed
# 'genome.fasta' is the full path to the genome file (bwa and fai index files need to be in the same directory)
# 'fastqDir' is the directory with the fastq files
# 'outDir' is the directory for the results


# variable definitions from command line parameters

sample=$1
infoFile=$2
map=$3
genome=$4
fastqDir=$5
outDir=$6

# only run mapping stage if 'map' variable is set to true
# explanation: usually, the read mapping only has to be done once
# however, you might want to re-run the analysis in R several times

if [ $map = true ]
then
	echo "running read mapping with bwa and samtools sort-index"
	
	# map reads with bwa and postprocess output with samtools
	bwa mem -t 4 $genome $fastqDir/${sample}_R1.fastq.gz $fastqDir/${sample}_R2.fastq.gz | samtools view -b - | samtools sort - -o $outDir/${sample}.bam
	samtools index $outDir/${sample}.bam
fi

# continue analysis in R
# relevant parameters are provided to the R script 

echo "running CrispRVariants analysis including creation of plots"
Rscript crisprvar.r $sample $infoFile $genome $outDir
