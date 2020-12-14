library(data.table)
library(ggplot2)

name_ab <- "replicate 1/WGA"
name_cd <- "replicate 2/WGA"
name_k <- "replicate 1/native"
name_l <- "replicate 2/native"
name_k_normBy_ab <- "replicate 1"
name_l_normBy_cd <- "replicate 2"
name_PD2182 <- "PD2182 (RS II)/native"
name_PD2182sequel <- "PD2182/native"

plot.enrichment <- function(data, target_sample, ylab_text) {
    # You have to avoid a name conflict, because "data" has a column "sample".
    sub_data <- data[sample == target_sample]
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    target_sample <- ifelse(target_sample == "k_normBy_ab", name_k_normBy_ab, ifelse(target_sample == "l_normBy_cd", name_l_normBy_cd, target_sample))
    target_sample <- ifelse(target_sample == "ab", name_ab, ifelse(target_sample == "cd", name_cd, target_sample))
    target_sample <- ifelse(target_sample == "PD2182", name_PD2182, ifelse(target_sample == "PD2182sequel", name_PD2182sequel, target_sample))
    ggplot(sub_data, aes(base, enrichment)) + geom_point() + facet_grid(. ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab("Base") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        ggtitle(target_sample) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5))
}

plot.enrichment.dual_facet <- function(data, target_sample, ylab_text) {
    sub_data <- data[sample == target_sample]
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    target_sample <- ifelse(target_sample == "k_normBy_ab", name_k_normBy_ab, ifelse(target_sample == "l_normBy_cd", name_l_normBy_cd, target_sample))
    target_sample <- ifelse(target_sample == "ab", name_ab, ifelse(target_sample == "cd", name_cd, target_sample))
    target_sample <- ifelse(target_sample == "PD2182", name_PD2182, ifelse(target_sample == "PD2182sequel", name_PD2182sequel, target_sample))
    ggplot(sub_data, aes("", enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab(NULL) + scale_x_discrete(breaks = NULL) + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        ggtitle(target_sample) +
        theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5))
}

plot.enrichment.dual_facet.all_sample <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

plot.enrichment.dual_facet.all_sample.v2 <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "k", name_k, ifelse(sample == "l", name_l,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))))]
    sub_data[, "sample" := factor(sample, levels = c(name_ab, name_cd, name_k, name_l, name_PD2182, name_PD2182sequel, name_k_normBy_ab, name_l_normBy_cd))]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5, size = 8 / ggplot2::.pt) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

# For WGA samples
plot.enrichment.dual_facet.all_sample.subset1 <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data <- sub_data[sample == "ab" | sample == "cd"]
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))]
    sub_data <- sub_data[region == "enhancer" | region == "five prime utr" | region == "exon" | region == "intron" | region == "three prime utr" | region == "tandem repeat"]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

# For WGA samples
plot.enrichment.dual_facet.all_sample.subset2 <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data <- sub_data[sample == "ab" | sample == "cd"]
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))]
    sub_data <- sub_data[region == "enhancer" | region == "promoter" | region == "five prime utr" | region == "exon" | region == "intron" | region == "three prime utr" | region == "tandem repeat"]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

# For native samples
plot.enrichment.dual_facet.all_sample.subset3 <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data <- sub_data[sample == "k_normBy_ab" | sample == "l_normBy_cd"]
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))]
    sub_data <- sub_data[region == "promoter" | region == "exon" | region == "intron" | region == "three prime utr" | region == "TTS [0, 500]" | region == "tandem repeat"]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

# For WGA samples
plot.enrichment.dual_facet.all_sample.subset1.v2 <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data <- sub_data[sample == "ab" | sample == "cd" | sample == "k" | sample == "l"]
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "k", name_k, ifelse(sample == "l", name_l,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))))]
    sub_data[, "sample" := factor(sample, levels = c(name_ab, name_cd, name_k, name_l, name_PD2182, name_PD2182sequel, name_k_normBy_ab, name_l_normBy_cd))]
    sub_data <- sub_data[region == "enhancer" | region == "five prime utr" | region == "exon" | region == "intron" | region == "three prime utr" | region == "tandem repeat"]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5, size = 8 / ggplot2::.pt) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

# For WGA samples
plot.enrichment.dual_facet.all_sample.subset2.v2 <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data <- sub_data[sample == "ab" | sample == "cd" | sample == "k" | sample == "l"]
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "k", name_k, ifelse(sample == "l", name_l,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))))]
    sub_data[, "sample" := factor(sample, levels = c(name_ab, name_cd, name_k, name_l, name_PD2182, name_PD2182sequel, name_k_normBy_ab, name_l_normBy_cd))]
    sub_data <- sub_data[region == "enhancer" | region == "promoter" | region == "five prime utr" | region == "exon" | region == "intron" | region == "three prime utr" | region == "tandem repeat"]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5, size = 8 / ggplot2::.pt) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

# For native samples
plot.enrichment.dual_facet.all_sample.subset3.v2 <- function(data, ylab_text) {
    sub_data <- data
    #sub_data[region != "All", label := .(sprintf("%s\np = %.3g", region, pvalue))]
    #sub_data[region == "All", label := .(sprintf("%s\n(Control)", region))]
    sub_data[, qvalue := .(p.adjust(pvalue, method = "BH"))][, qvalue_mark := .(ifelse(qvalue > 0.05, "", ifelse(qvalue > 0.01, "*", ifelse(qvalue > 0.001, "**", "***"))))]
    sub_data$base <- factor(sub_data$base, levels = c("A", "C", "G", "T", "ALL"))
    sub_data <- sub_data[sample == "k_normBy_ab" | sample == "l_normBy_cd"]
    sub_data[, "sample" := ifelse(sample == "k_normBy_ab", name_k_normBy_ab, ifelse(sample == "l_normBy_cd", name_l_normBy_cd,
        ifelse(sample == "ab", name_ab, ifelse(sample == "cd", name_cd,
        ifelse(sample == "k", name_k, ifelse(sample == "l", name_l,
        ifelse(sample == "PD2182", name_PD2182, ifelse(sample == "PD2182sequel", name_PD2182sequel, sample))))))))]
    sub_data[, "sample" := factor(sample, levels = c(name_ab, name_cd, name_k, name_l, name_PD2182, name_PD2182sequel, name_k_normBy_ab, name_l_normBy_cd))]
    sub_data <- sub_data[region == "promoter" | region == "exon" | region == "intron" | region == "three prime utr" | region == "TTS [0, 500]" | region == "tandem repeat"]
    ggplot(sub_data, aes(sample, enrichment)) + geom_point() + facet_grid(base ~ region_label) +
        geom_text(aes(label = qvalue_mark), vjust = 0, hjust = 0.5, size = 10 / ggplot2::.pt) +
        xlab("Sample") + ylab(ylab_text) + geom_hline(yintercept = 1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), plot.title = element_text(hjust = 0.5)) +
        scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.2)))
}

plot.intron_details <- function(data, target_sample) {
    sub_data <- data[sample == target_sample]
    intron_callable <- sub_data[region == "intron", callableNum]
    intron_modified <- sub_data[region == "intron", modifiedNum]
    intron_outside_exon_callable <- sub_data[region == "intron outside exon", callableNum]
    intron_outside_exon_modified <- sub_data[region == "intron outside exon", modifiedNum]
    intron_overlapped_with_exon_callable <- intron_callable - intron_outside_exon_callable
    intron_overlapped_with_exon_modified <- intron_modified - intron_outside_exon_modified
    plot_data <- data.table(region = c("intron outside exon", "intron overlapped with exon"), fraction = c(intron_outside_exon_modified / intron_outside_exon_callable, intron_overlapped_with_exon_modified / intron_overlapped_with_exon_callable))
    pvalue <- fisher.test(rbind(c(intron_outside_exon_callable - intron_outside_exon_modified, intron_outside_exon_modified), c(intron_overlapped_with_exon_callable - intron_overlapped_with_exon_modified, intron_overlapped_with_exon_modified)))$p.value
    ggplot(plot_data, aes(region, fraction)) + geom_point() + ggtitle(sprintf("%s: Modification fraction per valid adenine within features on both the strands", target_sample), subtitle = sprintf("Fisher's exact test p-value: %.3g", pvalue)) +
        xlab("Region") + ylab("Modification fraction")
}

set_region_label <- function(d) {
    d[, "region_label" := .(region)]
    d[region == "five prime utr", "region_label" := .("5' UTR")]
    d[region == "three prime utr", "region_label" := .("3' UTR")]
    d[region == "tandem repeat", "region_label" := .("tandem\nrepeat")]
    d[region == "repeat region", "region_label" := .("repeat\nregion")]
    d[region == "pseudogene exon", "region_label" := .("pseudogene\nexon")]
    d[region == "TSS [-500, 0]", "region_label" := .("TSS\n[-500, 0]")]
    d[region == "TSS [0, 500]", "region_label" := .("TSS\n[0, 500]")]
    d[region == "TTS [-500, 0]", "region_label" := .("TTS\n[-500, 0]")]
    d[region == "TTS [0, 500]", "region_label" := .("TTS\n[0, 500]")]
    d[, "region_label" := factor(region_label, levels = c("enhancer", "promoter", "TSS\n[-500, 0]", "TSS\n[0, 500]", "5' UTR", "exon", "intron", "3' UTR", "TTS\n[-500, 0]", "TTS\n[0, 500]", "tandem\nrepeat", "repeat\nregion", "ncRNA", "pseudogene", "pseudogene\nexon"))]
}
