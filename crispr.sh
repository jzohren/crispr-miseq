#!/bin/bash

# example execution: ./crispr.sh A15 info_file.csv
# where 'A15' is the sample name
# and 'info_file.csv' contains information about the guide sequence

# variable definitions

sample=$1
infoFile=$2
map=true
genome=genome.fa

# directory definitions

fastqDir=/user/home/working_dir/fastq
outDir=/user/home/working_dir/results

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
