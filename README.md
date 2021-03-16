# crispr-miseq
Analysis of MiSeq amplicon data from CRISPR experiments.

## Prerequisites 

This pipeline utilises BWA (http://bio-bwa.sourceforge.net/) to map short paired-end amplicon data to a genome sequence. After postprocessing with samtools (http://samtools.sourceforge.net/), the data is further analysed using the CrispRVariants R package (https://bioconductor.org/packages/release/bioc/html/CrispRVariants.html). In addition to CrispRVariants, the following R packages need to be installed:

* rtracklayer
* Biostrings
* seqinr
* GenomicFeatures
* glue

## info_file

All sample information needs to be provided in a comma-separted `info_file.csv`. It has to include the sample name, gene name, guide sequence excluding PAM, chromosome location of the guide sequence, start position of the guide sequence, end position of the guide sequence, strand of the guide sequence, and how many base pairs up- and downstream of the guide sequence should be analysed. An example of such a file is in `info_file.csv`. 

The directory where `genome.fasta` is located also needs to contain the BWA index, as well as a FAI index.
It is crucial that the chromosome names in the `info_file.csv` match those in `genome.fasta` (e.g. _chrY_ and _>chrY_, rather than _Y_ and _>chrY_ or vice-versa). 

## Execution

The pipeline can be run using the following command: 
```
./crispr.sh sample_name info_file.csv true/false genome.fasta fastqDir outDir
```

## Parameter definitions
* **sample_name** is the sample name
* **info_file.csv** contains information about the guide sequence (see above)
* **true/false** indicates whether the mapping step should be executed
* **genome.fa** is the full path to the genome file (bwa and fai index files need to be in the same directory)
* **fastqDir** is the directory with the fastq files
* **outDir** is the directory for the results
