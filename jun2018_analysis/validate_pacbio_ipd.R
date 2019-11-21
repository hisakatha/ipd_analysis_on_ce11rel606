source("load_jun2018_data.saved.fst.R")

library(data.table)

score_thres <- 20
ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
coverage_thres <- 25
bases <- c("A", "C", "G", "T")

#k_normBy_ab_data[, deviated_flags := list(score >= score_thres)]
#k_normBy_ab_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]

#l_normBy_cd_data[, deviated_flags := list(score >= score_thres)]
#l_normBy_cd_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]

#kl_normBy_abcd_data[, deviated_flags := list(score >= score_thres)]
#kl_normBy_abcd_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]


#union_data <- cbind(a = k_normBy_ab_data, b = l_normBy_cd_data)

#deviated_flags_intersection <- k_normBy_ab_data$deviated_flags & l_normBy_cd_data$deviated_flags

#modified_flags_intersection <- k_normBy_ab_data$modified_flags & l_normBy_cd_data$modified_flags

library(ggplot2)

plot_ipdratio_with_or_without_pacbio_estimate <- function(data, name) {
    x.label <- sprintf("%s log2 (native IPD / control IPD)", name)
    y.label <- sprintf("%s log2 (native IPD / IPD predicted using the PacBio software)", name)
    subdata <- data[native_tMean > 0 & control_tMean > 0]
    subdata[, c("log2_ipdRatio_without_pacbio", "log2_ipdRatio_with_pacbio") := .(log2(native_tMean / control_tMean), log2(native_tMean / modelPrediction))]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2_ipdRatio_without_pacbio], base_data[,log2_ipdRatio_with_pacbio], method = "pearson")
            subdata_stat[base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
    }
    p <- ggplot(subdata, aes(log2_ipdRatio_without_pacbio, log2_ipdRatio_with_pacbio)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d, cor = %.3g (p = %.3g)", N, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s: IPD ratio comparison with or without PacBio estimates per base", name), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    subdata <- subdata[native_coverage >= coverage_thres & control_coverage >= coverage_thres]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2_ipdRatio_without_pacbio], base_data[,log2_ipdRatio_with_pacbio], method = "pearson")
            subdata_stat[base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
    }
    p <- ggplot(subdata, aes(log2_ipdRatio_without_pacbio, log2_ipdRatio_with_pacbio)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d, cor = %.3g (p = %.3g)", N, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s (coverage >= %g): IPD ratio comparison with or without PacBio estimates per base", name, coverage_thres), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    cat(sprintf("Printed plots for %s\n", name))
}

plot_ipd_vs_pacbio_estimate <- function(data, name) {
    x.label <- sprintf("%s log2 (IPD)", name)
    y.label <- sprintf("log2 (IPD predicted using the PacBio software)")
    subdata <- data[tMean > 0]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2(tMean)], base_data[,log2(modelPrediction)], method = "pearson")
            subdata_stat[base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
    }
    p <- ggplot(subdata, aes(log2(tMean), log2(modelPrediction))) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d, cor = %.3g (p = %.3g)", N, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s: log_2 (IPD) vs. log_2 (IPD predicted using the PacBio software) per base", name), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    subdata <- subdata[coverage >= coverage_thres]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2(tMean)], base_data[,log2(modelPrediction)], method = "pearson")
            subdata_stat[base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
    }
    p <- ggplot(subdata, aes(log2(tMean), log2(modelPrediction))) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d, cor = %.3g (p = %.3g)", N, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s (coverage >= %g): log_2 (IPD) vs. log_2 (IPD predicted using the PacBio software) per base", name, coverage_thres), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    cat(sprintf("Printed plots for %s\n", name))
}

pdf(sprintf("validate_pacbio_ipd.coverage%g.pdf", coverage_thres), width = 14, height = 14)
plot_ipd_vs_pacbio_estimate(ab_data, "AB_WGA_VC2010_OP50")
plot_ipd_vs_pacbio_estimate(ab_data[refName != ecoli_chr], "AB_WGA_VC2010_OP50 in C. elegans")
plot_ipd_vs_pacbio_estimate(ab_data[refName == ecoli_chr], "AB_WGA_VC2010_OP50 in E. coli")
plot_ipd_vs_pacbio_estimate(cd_data, "CD_WGA_VC2010")
plot_ipd_vs_pacbio_estimate(cd_data[refName != ecoli_chr], "CD_WGA_VC2010 in C. elegans")
plot_ipd_vs_pacbio_estimate(cd_data[refName == ecoli_chr], "CD_WGA_VC2010 in E. coli")
plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data, "k_normBy_ab (VC2010+OP50)")
plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data[refName != ecoli_chr], "k_normBy_ab (VC2010+OP50) in C. elegans")
plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data[refName == ecoli_chr], "k_normBy_ab (VC2010+OP50) in E. coli")
plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data, "l_normBy_cd (VC2010)")
plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data[refName != ecoli_chr], "l_normBy_cd (VC2010) in C. elegans")
plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data[refName == ecoli_chr], "l_normBy_cd (VC2010) in E. coli")
invisible(dev.off())

