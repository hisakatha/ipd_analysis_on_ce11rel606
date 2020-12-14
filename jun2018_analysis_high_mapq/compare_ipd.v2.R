library(data.table)

score_thres <- 20
ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
coverage_thres <- 25
bases <- c("A", "C", "G", "T")

csv_name <- sprintf("compare_ipd.v2.coverage%g.csv", coverage_thres)

library(ggplot2)

plot_control_ipd_vs_native_ipd <- function(data, name1, name2) {
    # name1: name of the control sample
    # name2: name of the native sample
    x.label <- bquote(log[2] ~ "(IPD observed in" ~ .(name1) ~ ")")
    y.label <- bquote(log[2] ~ "(IPD observed in" ~ .(name2) ~ ")")
    subdata <- data[control_tMean > 0 & native_tMean > 0]
    subdata[, c("log2control_tMean", "log2native_tMean") := .(log2(control_tMean), log2(native_tMean))]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, r_square = Inf, rmse = Inf, sample = paste0(name1, ":", name2), coverage_thres = 0), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2control_tMean], base_data[,log2native_tMean], method = "pearson")
            base_r_square <- base_data[, 1 - sum((log2control_tMean - log2native_tMean)^2) / sum((log2control_tMean - mean(log2control_tMean))^2)]
            base_rmse <- base_data[, sqrt(sum((log2control_tMean - log2native_tMean)^2) / .N)]
            subdata_stat[base == b, c("cor.value", "p.value", "r_square", "rmse") := .(cor.ret$estimate, cor.ret$p.value, base_r_square, base_rmse)]
        }
    }
    fwrite(subdata_stat, file = csv_name, append = TRUE)
    p <- ggplot(subdata, aes(log2control_tMean, log2native_tMean)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nRMSE = %.3g\nR2 = %.3g\ncor = %.3g (p = %.3g)", N, rmse, r_square, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + xlab(x.label) + ylab(y.label)
    #print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    subdata <- subdata[control_coverage >= coverage_thres & native_coverage >= coverage_thres]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, r_square = Inf, rmse = Inf, sample = paste0(name1, ":", name2), coverage_thres = coverage_thres), by = base]
    for (b in bases) {
        base_data <- subdata[base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2control_tMean], base_data[,log2native_tMean], method = "pearson")
            base_r_square <- base_data[, 1 - sum((log2control_tMean - log2native_tMean)^2) / sum((log2control_tMean - mean(log2control_tMean))^2)]
            base_rmse <- base_data[, sqrt(sum((log2control_tMean - log2native_tMean)^2) / .N)]
            subdata_stat[base == b, c("cor.value", "p.value", "r_square", "rmse") := .(cor.ret$estimate, cor.ret$p.value, base_r_square, base_rmse)]
        }
    }
    fwrite(subdata_stat, file = csv_name, append = TRUE)
    p <- ggplot(subdata, aes(log2control_tMean, log2native_tMean)) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nRMSE = %.3g\nR2 = %.3g\ncor = %.3g (p = %.3g)", N, rmse, r_square, cor.value, p.value)), hjust = "inward", vjust = "inward") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    #p <- p + ggtitle(sprintf("%s (coverage >= %g): log_2 (IPD) vs. log_2 (IPD predicted using the PacBio software) per base", name, coverage_thres), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    x.label <- bquote(log[2] ~ "(IPD observed in" ~ .(name1) ~ "with its depth" >= ~ .(coverage_thres) ~ ")")
    y.label <- bquote(log[2] ~ "(IPD observed in" ~ .(name2) ~ "with its depth" >= ~ .(coverage_thres) ~ ")")
    p <- p + xlab(x.label) + ylab(y.label)
    #print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    cat(sprintf("Printed plots for %s and %s\n", name1, name2))
}

source("load_jun2018_data.saved.fst.basic_normalized.R")

setkey(k_normBy_ab_data, refName, base)
setkey(l_normBy_cd_data, refName, base)

pdf(sprintf("compare_ipd.v2.coverage%g.pdf", coverage_thres), width = 12, height = 12)
plot_control_ipd_vs_native_ipd(k_normBy_ab_data, "replicate 1 (WGA)", "replicate 1 (native)")
plot_control_ipd_vs_native_ipd(k_normBy_ab_data[refName != ecoli_chr], "replicate 1 (WGA) in C. elegans", "replicate 1 (native) in C. elegans")
plot_control_ipd_vs_native_ipd(k_normBy_ab_data[refName == ecoli_chr], "replicate 1 (WGA) in E. coli", "replicate 1 (native) in E. coli")
plot_control_ipd_vs_native_ipd(l_normBy_cd_data, "replicate 2 (WGA)", "replicate 2 (native)")
plot_control_ipd_vs_native_ipd(l_normBy_cd_data[refName != ecoli_chr], "replicate 2 (WGA) in C. elegans", "replicate 2 (native) in C. elegans")
plot_control_ipd_vs_native_ipd(l_normBy_cd_data[refName == ecoli_chr], "replicate 2 (WGA) in E. coli", "replicate 2 (native) in E. coli")
invisible(dev.off())
