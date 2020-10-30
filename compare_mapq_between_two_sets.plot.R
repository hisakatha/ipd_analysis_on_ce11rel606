library(data.table)
library(ggplot2)

# geom_hex does not reflect binwidth well when the data is sparse
compare_mapq_hex <- function(data, title1) {
    g1 <- ggplot(data, aes(top_mapq_celegans, top_mapq_ecoli)) +
        geom_hex(binwidth = c(5, 5)) +
        xlab("Top MapQ per read in C. elegans") +
        ylab("Top MapQ per read in E. coli") +
        ggtitle(title1) +
        scale_fill_continuous(trans = "log10") +
        theme(plot.title = element_text(size = 8))
    return(g1)
}

compare_mapq_bin2d <- function(data, title1) {
    g1 <- ggplot(data, aes(top_mapq_celegans, top_mapq_ecoli)) +
        geom_bin2d(binwidth = c(5, 5)) +
        xlab("Top MapQ per read in C. elegans") +
        ylab("Top MapQ per read in E. coli") +
        ggtitle(title1) +
        scale_fill_continuous(trans = "log10") +
        theme(plot.title = element_text(size = 8))
    return(g1)
}

compare_mapq_hist <- function(data, title1) {
    g1 <- ggplot(data, aes(abs_diff_mapq)) +
        geom_histogram(binwidth = 5) +
        xlab("Absolute difference of top MapQ per read between C. elegans and E. coli") +
        ggtitle(title1) +
        theme(plot.title = element_text(size = 8))
    return(g1)
}

compare_mapq_ecdf <- function(data, title1) {
    g1 <- ggplot(data, aes(abs_diff_mapq)) +
        stat_ecdf(n = 1000) +
        xlab("Absolute difference of top MapQ per read between C. elegans and E. coli") +
        ylab("Empirical cumulative distribution") +
        ggtitle(title1) +
        theme(plot.title = element_text(size = 8))
    return(g1)
}

args <- commandArgs(trailingOnly = TRUE)
sample_name <- args[1]

title_full <- sprintf("All the reads in %s (MapQ for an unmapped read is regarded as 0)", sample_name)
title_subset <- sprintf("Reads aligned to both C. elegans and E. coli in %s", sample_name)

path <- "mapped.alignmentset.merged.sorted_by_name.bam.compare_celegans_ecoli.csv"
d1 <- fread(path)

d1_subset <- d1[count_celegans > 0 & count_ecoli > 0]

pdf_height <- 6
pdf_width <- 6

pdf("compare_mapq_between_two_sets.plot.pdf", height = pdf_height, width = pdf_width)
print(compare_mapq_bin2d(d1, title_full))
print(compare_mapq_hist(d1, title_full))
print(compare_mapq_ecdf(d1, title_full))
print(compare_mapq_bin2d(d1_subset, title_subset))
print(compare_mapq_hist(d1_subset, title_subset))
print(compare_mapq_ecdf(d1_subset, title_subset))
invisible(dev.off())
