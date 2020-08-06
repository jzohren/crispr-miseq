#!/bin/bash
# crispr analysis of MiSeq data
#SBATCH --job-name=crispr
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --partition=cpu
#SBATCH --array=1-10
#SBATCH --output=crispr_%A_%a.out

# action required:
# change array parameter (above) to the number of files to analyse (i.e. number of lines in readFiles)

echo 'JOB STARTED'
date


# action required:
# select your genome file here and modify path to go through your lab
GENOME=/camp/svc/reference/Genomics/aws-igenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa

SEQRUN=200630
MAP=true

if [ $MAP = true ]
then

	echo "running read mapping"
	
	# action required:
	# change name and potentially path to readFiles
	READ1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" readFiles_$SEQRUN)
	READ2="${READ1/_R1_/_R2_}"
	LIB=$(echo $READ1 | rev | cut -d "/" -f 1 | rev | cut -d "_" -f 1)


	# action required:
	# specify path where to store bam files
	OUT=/camp/lab/turnerj/working/Charlotte/bams/$SEQRUN/
	mkdir -p $OUT
	cd $OUT
	echo "cd $OUT"

	ml BWA
	bwa mem -t 4 $GENOME $READ1 $READ2 > "$LIB".sam

	ml purge    
	ml SAMtools

	samtools view -b "$LIB".sam -o "$LIB".bam
	samtools sort "$LIB".bam -o "$LIB".sorted.bam
	samtools index "$LIB".sorted.bam
	rm "$LIB".sam "$LIB".bam 

fi


ml Anaconda2 SAMtools
source activate r_3.6

DIR=/camp/lab/turnerj/working/Charlotte/scripts
cd $DIR
echo "cd $DIR"

Rscript crisprvar.r $SLURM_ARRAY_TASK_ID $SEQRUN

date
echo 'JOB ENDED'
