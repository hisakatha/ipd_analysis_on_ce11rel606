library(data.table)
library(ggplot2)
plot.enrichment <- function(data, sample, ylab_text) {
    # Avoid a name conflict, because "data" has a column "sample".
    tmp_sample <- sample
    sub_data <- data[sample == tmp_sample]
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sample <- ifelse(sample == "k_normBy_ab", "VC2010+OP50", ifelse(sample == "l_normBy_cd", "VC2010", sample))
    sample <- ifelse(sample == "ab", "VC2010+OP50 (WGA)", ifelse(sample == "cd", "VC2010 (WGA)", sample))
    sample <- ifelse(sample == "PD2182", "PD2182 (PacBio RS II)", ifelse(sample == "PD2182sequel", "PD2182 (PacBio Sequel v1.2)", sample))
    ggplot(sub_data, aes(base, enrichment)) + geom_point() + facet_grid(. ~ region) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab("Base") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        ggtitle(sample) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5))
}

plot.enrichment.dual_facet <- function(data, sample, ylab_text) {
    # Avoid a name conflict, because "data" has a column "sample".
    tmp_sample <- sample
    sub_data <- data[sample == tmp_sample]
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sample <- ifelse(sample == "k_normBy_ab", "VC2010+OP50", ifelse(sample == "l_normBy_cd", "VC2010", sample))
    sample <- ifelse(sample == "ab", "VC2010+OP50 (WGA)", ifelse(sample == "cd", "VC2010 (WGA)", sample))
    sample <- ifelse(sample == "PD2182", "PD2182 (PacBio RS II)", ifelse(sample == "PD2182sequel", "PD2182 (PacBio Sequel v1.2)", sample))
    ggplot(sub_data, aes("", enrichment)) + geom_point() + facet_grid(base ~ region) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab(NULL) + scale_x_discrete(breaks = NULL) + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        ggtitle(sample) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5))
}

plot.enrichment.dual_facet.all_sample <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", "VC2010+OP50", ifelse(sample == "l_normBy_cd", "VC2010",
        ifelse(sample == "ab", "VC2010+OP50 (WGA)", ifelse(sample == "cd", "VC2010 (WGA)",
        ifelse(sample == "PD2182", "PD2182 (PacBio RS II)", ifelse(sample == "PD2182sequel", "PD2182 (PacBio Sequel v1.2)", sample))))))]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5))
}

plot.intron_details <- function(data, sample) {
    # Avoid a name conflict, because "data" has a column "sample".
    tmp_sample <- sample
    sub_data <- data[sample == tmp_sample]
    intron_callable <- sub_data[region == "intron", callableNum]
    intron_modified <- sub_data[region == "intron", modifiedNum]
    intron_outside_exon_callable <- sub_data[region == "intron outside exon", callableNum]
    intron_outside_exon_modified <- sub_data[region == "intron outside exon", modifiedNum]
    intron_overlapped_with_exon_callable <- intron_callable - intron_outside_exon_callable
    intron_overlapped_with_exon_modified <- intron_modified - intron_outside_exon_modified
    plot_data <- data.table(region = c("intron outside exon", "intron overlapped with exon"), fraction = c(intron_outside_exon_modified / intron_outside_exon_callable, intron_overlapped_with_exon_modified / intron_overlapped_with_exon_callable))
    pvalue <- fisher.test(rbind(c(intron_outside_exon_callable - intron_outside_exon_modified, intron_outside_exon_modified), c(intron_overlapped_with_exon_callable - intron_overlapped_with_exon_modified, intron_overlapped_with_exon_modified)))$p.value
    ggplot(plot_data, aes(region, fraction)) + geom_point() + ggtitle(sprintf("%s: Modification fraction per valid adenine within features on both the strands", sample), subtitle = sprintf("Fisher's exact test p-value: %.3g", pvalue)) +
        xlab("Region") + ylab("Modification fraction")
}

rename_region <- function(d) {
    d[region == "five prime utr", "region" := .("5' UTR")]
    d[region == "three prime utr", "region" := .("3' UTR")]
    d[region == "tandem repeat", "region" := .("tandem\nrepeat")]
    d[region == "pseudogene exon", "region" := .("pseudogene\nexon")]
    d[region == "TSS [-500, 0]", "region" := .("TSS\n[-500, 0]")]
    d[region == "TSS [0, 500]", "region" := .("TSS\n[0, 500]")]
    d[region == "TTS [-500, 0]", "region" := .("TTS\n[-500, 0]")]
    d[region == "TTS [0, 500]", "region" := .("TTS\n[0, 500]")]
    d[, "region" := factor(region, levels = c("enhancer", "promoter", "TSS\n[-500, 0]", "TSS\n[0, 500]", "5' UTR", "exon", "intron", "3' UTR", "TTS\n[-500, 0]", "TTS\n[0, 500]", "tandem\nrepeat", "ncRNA", "pseudogene", "pseudogene\nexon"))]
}
