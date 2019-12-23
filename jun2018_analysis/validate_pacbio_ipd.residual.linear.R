library(data.table)

score_thres <- 20
ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
coverage_thres <- 25
bases <- c("A", "C", "G", "T")

#csv_name <- sprintf("validate_pacbio_ipd.coverage%g.csv", coverage_thres)

library(ggplot2)

plot_ipdratio_with_or_without_pacbio_estimate <- function(data, name) {
    x.label <- bquote(log[2] ~ "(IPD observed in" ~ .(name) ~ ")")
    y.label <- bquote(log[2] ~ "(IPD predicted using the PacBio software)")
    subdata <- data[native_tMean > 0 & control_tMean > 0]
    subdata[, c("log2_ipdRatio_without_pacbio", "log2_ipdRatio_with_pacbio") := .(log2(native_tMean / control_tMean), log2(native_tMean / modelPrediction))]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, r_square = Inf, rmse = Inf, sample = name, coverage_thres = 0), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2_ipdRatio_without_pacbio], base_data[,log2_ipdRatio_with_pacbio], method = "pearson")
            base_r_square <- base_data[, 1 - sum((log2_ipdRatio_without_pacbio - log2_ipdRatio_with_pacbio)^2) / sum((log2_ipdRatio_without_pacbio - mean(log2_ipdRatio_without_pacbio))^2)]
            base_rmse <- base_data[, sqrt(sum((log2_ipdRatio_without_pacbio - log2_ipdRatio_with_pacbio)^2) / .N)]
            subdata_stat[base == b, c("cor.value", "p.value", "r_square", "rmse") := .(cor.ret$estimate, cor.ret$p.value, base_r_square, base_rmse)]
        }
    }
    fwrite(subdata_stat, file = csv_name, append = TRUE)
    p <- ggplot(subdata, aes(log2_ipdRatio_without_pacbio, log2_ipdRatio_with_pacbio)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nRMSE = %.3g\nR2 = %.3g\ncor = %.3g (p = %.3g)", N, rmse, r_square, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    #p <- p + ggtitle(sprintf("%s: IPD ratio comparison with or without PacBio estimates per base", name), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    p <- p + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    subdata <- subdata[native_coverage >= coverage_thres & control_coverage >= coverage_thres]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, r_square = Inf, rmse = Inf, sample = name, coverage_thres = coverage_thres), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2_ipdRatio_without_pacbio], base_data[,log2_ipdRatio_with_pacbio], method = "pearson")
            base_r_square <- base_data[, 1 - sum((log2_ipdRatio_without_pacbio - log2_ipdRatio_with_pacbio)^2) / sum((log2_ipdRatio_without_pacbio - mean(log2_ipdRatio_without_pacbio))^2)]
            base_rmse <- base_data[, sqrt(sum((log2_ipdRatio_without_pacbio - log2_ipdRatio_with_pacbio)^2) / .N)]
            subdata_stat[base == b, c("cor.value", "p.value", "r_square", "rmse") := .(cor.ret$estimate, cor.ret$p.value, base_r_square, base_rmse)]
        }
    }
    fwrite(subdata_stat, file = csv_name, append = TRUE)
    p <- ggplot(subdata, aes(log2_ipdRatio_without_pacbio, log2_ipdRatio_with_pacbio)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nRMSE = %.3g\nR2 = %.3g\ncor = %.3g (p = %.3g)", N, rmse, r_square, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    #p <- p + ggtitle(sprintf("%s (coverage >= %g): IPD ratio comparison with or without PacBio estimates per base", name, coverage_thres), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    x.label <- bquote(log[2] ~ "(IPD observed in" ~ .(name) ~ "with its depth" >= ~ .(coverage_thres) ~ ")")
    p <- p + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    cat(sprintf("Printed plots for %s\n", name))
}

plot_ipd_vs_pacbio_estimate <- function(data, name) {
    x.label <- bquote("IPD observed in" ~ .(name))
    y.label <- bquote("(observed IPD) - (predicted IPD)")
    subdata <- data[tMean > 0]
    #subdata[, c("log2tMean", "log2modelPrediction") := .(log2(tMean), log2(modelPrediction))]
    #subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, r_square = Inf, rmse = Inf, sample = name, coverage_thres = 0), by = base]
    #for (b in bases) {
    #    base_data <- subdata[base == b]
    #    if (base_data[,.N] >= 3L) {
    #        cor.ret <- cor.test(base_data[,log2tMean], base_data[,log2modelPrediction], method = "pearson")
    #        base_r_square <- base_data[, 1 - sum((log2tMean - log2modelPrediction)^2) / sum((log2tMean - mean(log2tMean))^2)]
    #        base_rmse <- base_data[, sqrt(sum((log2tMean - log2modelPrediction)^2) / .N)]
    #        subdata_stat[base == b, c("cor.value", "p.value", "r_square", "rmse") := .(cor.ret$estimate, cor.ret$p.value, base_r_square, base_rmse)]
    #    }
    #}
    #fwrite(subdata_stat, file = csv_name, append = TRUE)
    p <- ggplot(subdata, aes(tMean, tMean - modelPrediction)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both)
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    #p <- p + ggtitle(sprintf("%s: log_2 (IPD) vs. log_2 (IPD predicted using the PacBio software) per base", name), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    p <- p + xlab(x.label) + ylab(y.label)
    print(p)
    #print(p + scale_fill_continuous(trans = "log10"))
    subdata <- subdata[coverage >= coverage_thres]
    #subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, r_square = Inf, rmse = Inf, sample = name, coverage_thres = coverage_thres), by = base]
    #for (b in bases) {
    #    base_data <- subdata[base == b]
    #    if (base_data[,.N] >= 3L) {
    #        cor.ret <- cor.test(base_data[,log2(tMean)], base_data[,log2(modelPrediction)], method = "pearson")
    #        base_r_square <- base_data[, 1 - sum((log2tMean - log2modelPrediction)^2) / sum((log2tMean - mean(log2tMean))^2)]
    #        base_rmse <- base_data[, sqrt(sum((log2tMean - log2modelPrediction)^2) / .N)]
    #        subdata_stat[base == b, c("cor.value", "p.value", "r_square", "rmse") := .(cor.ret$estimate, cor.ret$p.value, base_r_square, base_rmse)]
    #    }
    #}
    #fwrite(subdata_stat, file = csv_name, append = TRUE)
    p <- ggplot(subdata, aes(tMean, tMean - modelPrediction)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both)
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    #p <- p + ggtitle(sprintf("%s (coverage >= %g): log_2 (IPD) vs. log_2 (IPD predicted using the PacBio software) per base", name, coverage_thres), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    x.label <- bquote("IPD observed in" ~ .(name) ~ "with its depth" >= ~ .(coverage_thres))
    p <- p + xlab(x.label) + ylab(y.label)
    print(p)
    #print(p + scale_fill_continuous(trans = "log10"))
    cat(sprintf("Printed plots for %s\n", name))
}

source("load_jun2018_data.saved.fst.basic_wga.R")

setkey(ab_data, refName, base)
setkey(cd_data, refName, base)
#setkey(k_normBy_ab_data, refName, base)
#setkey(l_normBy_cd_data, refName, base)

pdf(sprintf("validate_pacbio_ipd.residual.linear.coverage%g.pdf", coverage_thres), width = 12, height = 12)
plot_ipd_vs_pacbio_estimate(ab_data, "VC2010+OP50(WGA)")
plot_ipd_vs_pacbio_estimate(ab_data[refName != ecoli_chr], "VC2010+OP50(WGA) in C. elegans")
plot_ipd_vs_pacbio_estimate(ab_data[refName == ecoli_chr], "VC2010+OP50(WGA) in E. coli")
plot_ipd_vs_pacbio_estimate(cd_data, "VC2010(WGA)")
plot_ipd_vs_pacbio_estimate(cd_data[refName != ecoli_chr], "VC2010(WGA) in C. elegans")
plot_ipd_vs_pacbio_estimate(cd_data[refName == ecoli_chr], "VC2010(WGA) in E. coli")
#plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data, "VC2010+OP50")
#plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data[refName != ecoli_chr], "VC2010+OP50 in C. elegans")
#plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data[refName == ecoli_chr], "VC2010+OP50 in E. coli")
#plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data, "VC2010")
#plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data[refName != ecoli_chr], "VC2010 in C. elegans")
#plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data[refName == ecoli_chr], "VC2010 in E. coli")
#plot_ipd_vs_pacbio_estimate(ab_data, "AB_WGA_VC2010_OP50")
#plot_ipd_vs_pacbio_estimate(ab_data[refName != ecoli_chr], "AB_WGA_VC2010_OP50 in C. elegans")
#plot_ipd_vs_pacbio_estimate(ab_data[refName == ecoli_chr], "AB_WGA_VC2010_OP50 in E. coli")
#plot_ipd_vs_pacbio_estimate(cd_data, "CD_WGA_VC2010")
#plot_ipd_vs_pacbio_estimate(cd_data[refName != ecoli_chr], "CD_WGA_VC2010 in C. elegans")
#plot_ipd_vs_pacbio_estimate(cd_data[refName == ecoli_chr], "CD_WGA_VC2010 in E. coli")
#plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data, "k_normBy_ab (VC2010+OP50)")
#plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data[refName != ecoli_chr], "k_normBy_ab (VC2010+OP50) in C. elegans")
#plot_ipdratio_with_or_without_pacbio_estimate(k_normBy_ab_data[refName == ecoli_chr], "k_normBy_ab (VC2010+OP50) in E. coli")
#plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data, "l_normBy_cd (VC2010)")
#plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data[refName != ecoli_chr], "l_normBy_cd (VC2010) in C. elegans")
#plot_ipdratio_with_or_without_pacbio_estimate(l_normBy_cd_data[refName == ecoli_chr], "l_normBy_cd (VC2010) in E. coli")
invisible(dev.off())
