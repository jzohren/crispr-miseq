#!/bin/bash

# example execution: ./crispr.sh -n A15 -i info_file.csv -m -s species -f fastqDir -o outDir
# where:
# -n is the sample name [string]
# -i contains information about the guide sequence [string]
# -m indicates whether the mapping step should be executed (omit if you don't want to map) [boolean]
# -s is the name of the species being analysed [mouse/human]
# -f is the directory with the fastq files [string]
# -o is the directory for the results [string]


# variable definitions from command line parameters

# default for mapping option unless -m flag is set
map=false

while getopts ":n:i:ms:f:o:" opt; do
	case $opt in
		n) sample="$OPTARG"
		;;
		i) infoFile="$OPTARG"
		;;
		m) map=true
		;;
		s) species_tmp="$OPTARG"
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

species=`echo $species_tmp | tr '[:upper:]' '[:lower:]'`

if [ $species = mouse ]
then
	ref_1='/camp/svc/reference/Genomics/aws-igenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa'
elif [ $species = human ]
then
	ref_1='/camp/svc/reference/Genomics/aws-igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'
else
	echo "ERROR: unknown species provided; only "mouse" or "human" allowed"
fi

split_by='Sequence/'
prefix=${ref_1%$split_by*}
suffix='WholeGenomeFasta/genome.fa'
ref_2=$prefix$split_by$suffix

# validate
if [[ ! -f "$ref_2" ]]; 
then
	echo "ERROR: wrong ref provided"
	echo "e.g.: /camp/svc/reference/Genomics/aws-igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
	exit 1;
fi


# only run mapping stage if 'map' variable is set to true
# explanation: usually, the read mapping only has to be done once
# however, you might want to re-run the analysis in R several times

if [ $map = true ]
then
	echo "running read mapping with bwa and samtools sort-index"
	
	# map reads with bwa and postprocess output with samtools
	bwa mem -t 4 $ref_1 $fastqDir/${sample}_R1.fastq.gz $fastqDir/${sample}_R2.fastq.gz | samtools view -b - | samtools sort - -o $outDir/${sample}.bam
	samtools index $outDir/${sample}.bam
fi

# continue analysis in R
# relevant parameters are provided to the R script 

echo "running CrispRVariants analysis including creation of plots"
Rscript crisprvar.r $sample $infoFile $ref_2 $outDir
