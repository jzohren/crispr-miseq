# library loads

library(CrispRVariants)
library(rtracklayer)
library(Biostrings)
library(seqinr)
library(GenomicFeatures)
library(glue)

# variable definitions (provided to/from crispr.sh script)

args 		<- commandArgs(trailingOnly = TRUE)
sample 		<- args[1]
info_file 	<- read.table(args[2], sep = "\t")
genome 		<- args[3]
out_dir 	<- args[4]

setwd(out_dir)
colnames(info_file) <- c("sample", "gene", "seq", "chr", "start", "end", "strand", "plus", "minus")
gene 				<- as.character(info_file$gene)
sample_name 		<- glue("{sample}_{gene}")
bam 				<- glue("{out_dir}/{sample}.bam")

# create GRanges object of area around guide sequence

region <- GRanges(seqnames = info_file$chr, ranges = IRanges(info_file$start - info_file$minus, info_file$end + info_file$plus), strand = info_file$strand)

# the following steps are processed outside R, in your unix environment
# the 'region' object is used to extract a short FASTA sequence from the genome file

system(sprintf(glue("samtools faidx {genome} %s:%s-%s > %s.seq"), seqnames(region), start(region), end(region), sample_name))
fasta 		<- read.fasta(glue("{sample_name}.seq"), as.string = T)
reference 	<- DNAString(unlist(fasta))
system(sprintf("rm %s.seq", sample_name))

if (info_file$strand == "-") {
  reference <- reverseComplement(reference)
}

# the main CrispRVariants object is being created

crispr_set <- readsToTarget(bam, region, reference = reference, target.loc = 17 + info_file$minus, collapse.pairs = T, chimera.to.target = 200, names = sample_name)

var_counts 	<- variantCounts(crispr_set)
eff_raw 	<- mutationEfficiency(crispr_set)

# information about mutation efficancy and allele frequencies are written to simple text files

write("mutation efficiency:", file = glue("{sample_name}_mutEffic.txt"))
write(rbind(names(eff_raw), eff_raw), file = glue("{sample_name}_mutEffic.txt"), append = T)

allele_freqs <- barplotAlleleFreqs(var_counts, include.table = F)
write.table(allele_freqs$data, file = glue("{sample_name}_allele_freqs.txt"), sep = "\t", quote = F, row.names = F)

# plots are saved in a PDF file

pdf(file = glue("{sample_name}.pdf"), width = 25, height = 30)
  plotVariants(crispr_set, plotAlignments.args = list(tile.height = 1, ins.size = 3, legend.cols = 5))
  barplotAlleleFreqs(var_counts, include.table = F)
  plotAlignments(crispr_set, top.n = nrow(crispr_set$cigar_freqs), add.other = F, ins.size = 3, legend.cols = 13, tile.height = 1)
dev.off()

# a list of all possible variants is written to a text file

write.table(var_counts, file = glue("{sample_name}_variants.txt"), sep = "\t", quote = F)

# the variants occurring most frequently are plotted in a pie chart
# both with and without the "no variant" category

selection 		<- which(var_counts > 99)
top_mutations 	<- var_counts[selection]
mut_percentages <- round(100 * (var_counts / sum(var_counts)), 2)

pie_labels <- rownames(var_counts)[selection]
pie_labels <- glue("{mut_percentages[selection]} % {pie_labels}")
pie_colors <- sample(colors(), length(top_mutations))

pdf(file = glue("{sample_name}_pie.pdf"), width = 9.5, height = 7)
 pie(top_mutations, labels = NA, col = pie_colors, radius = 0.9, main = sample_name)
 legend("right", pie_labels, col = pie_colors, pch = 15)
 pie(top_mutations[-1], labels = NA, col = pie_colors, radius = 0.9)
 legend("right", pie_labels[-1], col = pie_colors, pch = 15)
dev.off()
