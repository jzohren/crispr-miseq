# crispr-miseq
Analysis of MiSeq amplicon data from CRISPR experiments. This branch is specific to using it within the Nautilus Suite of The Francis Crick Institute. 

## Prerequisites 

This pipeline utilises [BWA](http://bio-bwa.sourceforge.net/) to map short paired-end amplicon data to a genome sequence. After postprocessing with [samtools](http://samtools.sourceforge.net/), the data is further analysed using the [CrispRVariants](https://bioconductor.org/packages/release/bioc/html/CrispRVariants.html) R package. In addition to CrispRVariants, the following R packages need to be installed:

* rtracklayer
* Biostrings
* seqinr
* GenomicFeatures

## info_file

All sample information needs to be provided in a comma-separted `info_file.csv`. It has to include the sample name, gene name, guide sequence excluding PAM, chromosome location of the guide sequence, start position of the guide sequence, end position of the guide sequence, strand of the guide sequence, and how many base pairs up- and downstream of the guide sequence should be analysed. An example of such a file is in `info_file.csv`. Crucially, the sample name needs to be the whole first part of its corresponding fastq file. E.g. if the fastq file is `sample_name_1_R1_001.fastq.gz`, then the sample name in the `info_file.csv` should be `sample_name_1` (the remaining file name is appended in the `crispr.sh` script). It is also important that the chromosome names in the `info_file.csv` match those in the genome reference (e.g. _chrY_ and _>chrY_, rather than _Y_ and _>chrY_ or vice-versa). 

## Execution

The pipeline can be run using the following command: 
```sh
./crispr.sh -n sample_name -i info_file.csv -s species -f fastqDir -o outDir -m
```

## Parameter definitions
```
-n  is the sample name, i.e. the sample to be analysed [string]
-i  is the name of the `info_file.csv` that contains metainformation about the samples (see above) [string]
-s  is the name of the species being analysed [mouse/human]
-f  is the directory with the fastq files [string]
-o  is the directory for the results [string]
-m  is a flag indicating whether the mapping step should be executed (omit if you don't want to map)  [boolean]
```
