library(CrispRVariants)
library(rtracklayer)
library(Biostrings)
library(seqinr)
library(GenomicFeatures)


args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
out_dir <- args[2]
info_file <- args[3]
genome <- args[4]

setwd(out_dir)
colnames(info_file) <- c("sample", "gene", "seq", "chr", "start", "end", "strand", "plus", "minus")

sample_name <- paste0(sample, "_", as.character(info_file$gene))

bam <- paste0(out_dir, sample, ".bam")   

region <- GRanges(seqnames = info_file$chr, ranges = IRanges(info_file$start - info_file$minus, info_file$end + info_file$plus), strand = info_file$strand)

system(sprintf(paste("samtools faidx", genome, "%s:%s-%s > %s.seq"), seqnames(region), start(region), end(region), sample_name))
fasta <- read.fasta(paste0(sample_name, ".seq"), as.string = T)
reference <- DNAString(unlist(fasta))
system(sprintf("rm %s.seq", sample_name))

if (row$strand == "-") {
  reference <- reverseComplement(reference)
}

crispr_set <- readsToTarget(bam, region, reference = reference, target.loc = 17 + info_file$minus, collapse.pairs = T, chimera.to.target = 200, names = sample_name)

var_counts <- variantCounts(crispr_set)
eff_raw <- mutationEfficiency(crispr_set)

write("mutation efficiency:", file = paste0(sample_name, "_mutEffic.txt"))
write(rbind(names(eff_raw), eff_raw), file = paste0(sample_name, "_mutEffic.txt"), append = T)

write.table(var_counts, file = paste0(sample_name, "_variants.txt"), sep = "\t", quote = F)


selection <- which(var_counts > 99)
top_mutations <- var_counts[selection]
mut_percentages <- round(100 * (var_counts / sum(var_counts)), 2)
pie_labels <- rownames(var_counts)[selection]
pie_labels <- paste0(mut_percentages[selection], "% ", pie_labels)
pie_colors <- sample(colors(), length(top_mutations))

pdf(file = paste0(sample_name, "_pie.pdf"), width = 9.5, height = 7)
 pie(top_mutations, labels = NA, col = pie_colors, radius = 0.9, main = sample_name)
 legend("right", pie_labels, col = pie_colors, pch = 15)
 pie(top_mutations[-1], labels = NA, col = pie_colors, radius = 0.9)
 legend("right", pie_labels[-1], col = pie_colors, pch = 15)
dev.off()


pdf(file = paste0(sample_name, ".pdf"), width = 25, height = 30)
  plotVariants(crispr_set, plotAlignments.args = list(tile.height = 1, ins.size = 3, legend.cols = 5))
  barplotAlleleFreqs(var_counts, include.table = F)
  plotAlignments(crispr_set, top.n = nrow(crispr_set$cigar_freqs), add.other = F, ins.size = 3, legend.cols = 13, tile.height = 1)
dev.off()

allele_freqs <- barplotAlleleFreqs(var_counts, include.table = F)

write.table(allele_freqs$data, file = paste0(sample_name, "_allele_freqs.txt"), sep = "\t", quote = F, row.names = F)
