library(data.table)
library(ggplot2)

aln_len_ab <- fread("../jun2018_ab_pcr_vc2010_op50_all_alignments/mapped.alignmentset.merged.high_mapq.bam.alignment_length", header = FALSE, col.names = "alignment_length")
aln_len_cd <- fread("../jun2018_cd_pcr_vc2010_all_alignments/mapped.alignmentset.merged.high_mapq.bam.alignment_length", header = FALSE, col.names = "alignment_length")
aln_len_k <- fread("../jun2018_k_nopcr_vc2010_op50_all_alignments/mapped.alignmentset.merged.high_mapq.bam.alignment_length", header = FALSE, col.names = "alignment_length")
aln_len_l <- fread("../jun2018_l_nopcr_vc2010_all_alignments/mapped.alignmentset.merged.high_mapq.bam.alignment_length", header = FALSE, col.names = "alignment_length")
aln_len_PD2182sequel <- fread("../PD2182-sequel_all_alignments/mapped.alignmentset.merged.high_mapq.bam.alignment_length", header = FALSE, col.names = "alignment_length")

aln_len_ab[, "sample" := .("VC2010+OP50/WGA")]
aln_len_cd[, "sample" := .("VC2010/WGA")]
aln_len_k[, "sample" := .("VC2010+OP50/native")]
aln_len_l[, "sample" := .("VC2010/native")]
aln_len_PD2182sequel[, "sample" := .("PD2182/native")]
aln_len <- rbindlist(list(aln_len_ab, aln_len_cd, aln_len_k, aln_len_l, aln_len_PD2182sequel))

aln_len_means <- aln_len[, .(mean = mean(alignment_length)), by = sample]
fwrite(aln_len_means, file = "alignment_length.means.csv")

aln_len$sample <- factor(aln_len$sample, levels = c("VC2010+OP50/WGA", "VC2010/WGA", "VC2010+OP50/native", "VC2010/native", "PD2182/native"))
g1 <- ggplot(aln_len) + geom_freqpoly(aes(alignment_length, color = sample), binwidth = 500) + theme_classic() + xlab("Alignment length")
pdf("alignment_length.pdf", height = 5, width = 5)
print(g1)
print(g1 + scale_y_continuous(trans = "log10"))
invisible(dev.off())
