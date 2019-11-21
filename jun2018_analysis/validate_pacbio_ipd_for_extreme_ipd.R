library(data.table)

#score_thres <- 20
#ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
coverage_thres <- 25
collect_frac <- 0.01
bases <- c("A", "C", "G", "T")

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

plot_extreme_ipd <- function(data, name){
    data_high <- data[min(ceiling((1 - collect_frac) * .N) + 1, .N):.N]
    plot_ipd_vs_pacbio_estimate(data_high, sprintf("%s (Top %g IPD)", name, collect_frac))
    data_low <- data[1:max(floor(collect_frac * .N), 1)]
    plot_ipd_vs_pacbio_estimate(data_low, sprintf("%s (Bottom %g IPD)", name, collect_frac))
}

source("load_jun2018_data.saved.fst.wga.R")

#cd_high_cov <- fread(cmd = "cat extreme_ipd.high.cd.coverage25.gff | grep -v '^#' | sed -E 's/.*tMean=([\\.0-9]+);.*modelPrediction=([\\.0-9]+).*/\\1,\\2/; 1i\\tMean,modelPrediction'")
#cd_low_cov <- fread(cmd = "cat extreme_ipd.low.cd.coverage25.gff | grep -v '^#' | sed -E 's/.*tMean=([\\.0-9]+);.*modelPrediction=([\\.0-9]+).*/\\1,\\2/; 1i\\tMean,modelPrediction'")

ab_data <- ab_data[tMean > 0]
cd_data <- cd_data[tMean > 0]
abcd_data <- abcd_data[tMean > 0]

setkey(ab_data, tMean)
setkey(cd_data, tMean)
setkey(abcd_data, tMean)

pdf(sprintf("validate_pacbio_ipd_for_extreme_ipd.coverage%g.pdf", coverage_thres), width = 14, height = 14)
plot_extreme_ipd(ab_data, "AB_WGA_VC2010_OP50")
plot_extreme_ipd(ab_data[refName != ecoli_chr], "AB_WGA_VC2010_OP50 in C. elegans")
plot_extreme_ipd(ab_data[refName == ecoli_chr], "AB_WGA_VC2010_OP50 in E. coli")
plot_extreme_ipd(cd_data, "CD_WGA_VC2010")
plot_extreme_ipd(cd_data[refName != ecoli_chr], "CD_WGA_VC2010 in C. elegans")
plot_extreme_ipd(cd_data[refName == ecoli_chr], "CD_WGA_VC2010 in E. coli")
plot_extreme_ipd(abcd_data, "ABCD_WGA")
plot_extreme_ipd(abcd_data[refName != ecoli_chr], "ABCD_WGA in C. elegans")
plot_extreme_ipd(abcd_data[refName == ecoli_chr], "ABCD_WGA in E. coli")
invisible(dev.off())

