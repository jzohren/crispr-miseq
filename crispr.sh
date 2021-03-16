#!/bin/bash

# example execution: ./crispr.sh A15 info_file.csv true
# where 'A15' is the sample name
# and 'info_file.csv' contains information about the guide sequence
# and 'true/false' indicates whether the mapping step should be executed

# variable definitions from command line parameters

sample=$1
infoFile=$2
map=$3

# user-defined variable definitions

# provide full path to genome file
# bwa and fai index files need to be in the same directory
genome=<path to genome.fa>
# provide full path to directory with fastq files
fastqDir=<path to fastq files>
# provide an output directory for the results
outDir=<path to output directory>

# only run mapping stage if variable set to true
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
# provide relevant parameters to R script 

echo "running CrispRVariants analysis including creation of plots"
Rscript crisprvar.r $sample $infoFile $genome $outDir
