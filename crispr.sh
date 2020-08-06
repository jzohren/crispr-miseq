#!/bin/bash

# example execution:
# ./crispr.sh A15 true info_file_example

sample=$1
map=$2
infoFile=$3

genome=/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa

fastqDir=/user/home/working_dir/fastq
outDir=/user/home/working_dir/results

if [ $map = true ]
then
	echo "running read mapping with bwa and samtools sort-index"

	bwa mem -t 4 $genome $fastqDir/${sample}_R1.fastq.gz $fastqDir/${sample}_R2.fastq.gz | samtools view -b - | samtools sort - -o $outDir/${sample}.bam
	samtools index $outDir/${sample}.bam
fi

echo "running CrispRVariants analysis including creation of plots"
Rscript crisprvar.r $sample $infoFile $genome $outDir
