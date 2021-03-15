# crispr-miseq
Analysis of MiSeq amplicon data from CRISPR experiments.

This pipeline utilises BWA (http://bio-bwa.sourceforge.net/) to map short paired-end amplicon data to a genome sequence. After postprocessing with samtools (http://samtools.sourceforge.net/), the data is further analysed using the CrispRVariants R package (https://bioconductor.org/packages/release/bioc/html/CrispRVariants.html). 

All sample information needs to be provided in a comma-separted "info_file". It has to include the sample name, gene name, guide sequence including PAM, chromosome location of the guide sequence, start position of the guide sequence, end position of the guide sequence, strand of the guide sequence, and how many base pairs up- and downstream of the guide sequence should be analysed. An example of such a file is in `info_file.csv`. 

The directory where `genome.fa` is located also needs to contain the BWA index, as well as a FAI index.
It is crucial that the chromosome names in the `info _file.csv` match those in `genome.fa` (e.g. _chrY_ and _>chrY_, rather than _Y_ and _>chrY_ or vice-versa). 

The pipeline can be run using the following command: 
```
./crispr.sh <sample_name> <info_file.csv>
```
