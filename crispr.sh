#!/bin/bash

# example execution: ./crispr.sh -s A15 -i info_file.csv -m -g genome.fasta -f fastqDir -o outDir
# where:
# -s is the sample name
# -i contains information about the guide sequence
# -m indicates whether the mapping step should be executed (omit if you don't want to map)
# -g is the full path to the genome file (bwa and fai index files need to be in the same directory)
# -f is the directory with the fastq files
# -o is the directory for the results


# variable definitions from command line parameters

# default for mapping option unless -m flag is set
map=false

while getopts ":s:i:mg:f:o:" opt; do
	case $opt in
		s) sample="$OPTARG"
		;;
		i) infoFile="$OPTARG"
		;;
		m) map=true
		;;
		g) genome="$OPTARG"
		;;
		f) fastqDir="$OPTARG"
		;;
		o) outDir="$OPTARG"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
		;;
		:) echo "Option -$OPTARG requires an argument." >&2
	esac
done

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
