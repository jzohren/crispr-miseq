# crispr-miseq
Analysis of MiSeq amplicon data from CRISPR experiments.

This pipeline utilises BWA (http://bio-bwa.sourceforge.net/) to map short paired-end amplicon data to a genome sequence. After postprocessing with samtools (http://samtools.sourceforge.net/), the data is further analysed using the CrispRVariants R package (https://bioconductor.org/packages/release/bioc/html/CrispRVariants.html). 

In addition to the two read files in FASTA format, a tab-separted "info_file" is needed, which should include the sample name, gene name, guide sequence including PAM, chromosome location of guide sequence, start position of guide sequence, end location of guide sequence, strand of guide sequence, and how many base pairs up- and downstream of the guide sequence should be analysed. An example of such a file is in `info_file_example`. 

The directory where `genome.fa` is located also needs to contain the BWA index, as well as a FAI index.

The pipeline can be iniitated by running `./crispr.sh <sample_name> <info_file>`. 
