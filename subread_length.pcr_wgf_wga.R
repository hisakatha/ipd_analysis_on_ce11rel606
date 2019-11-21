library(ggplot2)
library(grid)
library(cowplot)

# Subread length
sc_len <- function(x){sc <- gsub("[0-9]+[=IDXH]","",x); ifelse(sc == "", 0, eval(parse(text=gsub("([0-9]+)S", "+\\1", sc))))}
pcr_alignment <- Rsamtools::scanBam("1_pcr/mapped.alignmentset.merged.bam")
cig <- pcr_alignment[[1]]$cigar
softclip_len <- sapply(cig, sc_len, USE.NAMES=FALSE)
pcr_subread_len <- pcr_alignment[[1]]$qwidth
#pcr_aligned_len <- pcr_subread_len - softclip_len

wgf_alignment <- Rsamtools::scanBam("2_whole_genome_fragments/mapped.alignmentset.merged.bam")
cig <- wgf_alignment[[1]]$cigar
softclip_len <- sapply(cig, sc_len, USE.NAMES=FALSE)
wgf_subread_len <- wgf_alignment[[1]]$qwidth
#wgf_aligned_len <- wgf_subread_len - softclip_len

wga_alignment <- Rsamtools::scanBam("3_whole_genome_pcr/mapped.alignmentset.merged.bam")
cig <- wga_alignment[[1]]$cigar
softclip_len <- sapply(cig, sc_len, USE.NAMES=FALSE)
wga_subread_len <- wga_alignment[[1]]$qwidth
#wga_aligned_len <- wga_subread_len - softclip_len

#subread_len_pval <- wilcox.test(pcr_subread_len, wgf_subread_len, alternative = "less")$p.value
#aligned_len_pval <- wilcox.test(pcr_aligned_len, wgf_aligned_len, alternative = "less")$p.value

#box_pl <- ggplot() + geom_boxplot(data = data.frame(PCR=pcr_aligned_len), aes("PCR", PCR, col="PCR")) + geom_boxplot(data = data.frame(WGF=wgf_aligned_len), aes("WGF", WGF, col="WGF")) + theme_bw() + xlab("Sample") + ylab("Aligned region length") + ggtitle(paste0("Box plots of aligned region lengths (p = ", sprintf("%.3G", aligned_len_pval), "; one-sided Wilcoxon rank sum test)")) + labs(col = "Sample")
box_pl2 <- ggplot() + geom_boxplot(data = data.frame(PCR=pcr_subread_len), aes("PCR", PCR, col="PCR")) + geom_boxplot(data = data.frame(WGF=wgf_subread_len), aes("WGF", WGF, col="WGF")) + geom_boxplot(data = data.frame(WGA=wga_subread_len), aes("WGA", WGA, col="WGA"))+ theme_bw() + xlab("Sample") + ylab("Subread length") + ggtitle(paste0("Box plots of subread lengths")) + labs(col = "Sample")
#ecdf_pl <- ggplot() + stat_ecdf(data = data.frame(PCR=pcr_aligned_len), aes(PCR, col="PCR")) + stat_ecdf(data = data.frame(WGF=wgf_aligned_len), aes(WGF, col="WGF")) + theme_bw() + xlab("Aligned region length") + ylab("F(length)") + ggtitle("ECDF of aligned region lengths") + labs(col = "Sample") + coord_cartesian(xlim = c(0, 15000))#+ theme(plot.margin = unit(c(6,3,6,3), "inches"))
ecdf_pl2 <- ggplot() + stat_ecdf(data = data.frame(PCR=pcr_subread_len), aes(PCR, col="PCR")) + stat_ecdf(data = data.frame(WGF=wgf_subread_len), aes(WGF, col="WGF")) + stat_ecdf(data = data.frame(WGA=wga_subread_len), aes(WGA, col="WGA")) + theme_bw() + xlab("Subread length") + ylab("F(length)") + ggtitle("ECDF of subread lengths") + labs(col = "Sample") + coord_cartesian(xlim = c(0, 15000))#+ theme(plot.margin = unit(c(6,3,6,3), "inches"))

alignment_pl <- plot_grid(box_pl2, ecdf_pl2, ncol = 1, labels = c("A", "B"), scale = c(0.8, 0.8))
#alignment_pl <- plot_grid(box_pl, box_pl2, ecdf_pl, ecdf_pl2, ncol = 1, labels = paste0("8", c("A", "B", "C", "D")))

pdf("subread_length.pcr_wgf_wga.pdf", width = 14, height = 18, onefile = TRUE)
print(alignment_pl)
dev.off()
