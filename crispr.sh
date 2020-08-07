#!/bin/bash

# example execution: ./crispr.sh A15 true info_file_example

# variable definitions

sample=$1
map=$2
infoFile=$3
genome=/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa

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
