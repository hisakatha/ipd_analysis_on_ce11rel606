library(data.table)
library(ggplot2)
library(cowplot)

qq_plot_per_base <- function(data, title_text){
    # data: a data.table with columns: base, value, orig_base; sorted by value
    data_dup <- copy(data)
    data_dup[, "base" := .("all")]
    data_dup <- rbind(data_dup, data)
    data_dup[, "theoretical" := .(qnorm(ppoints(.N))), by = base]
    data_dup$base <- factor(data_dup$base, levels = c("A", "C", "G", "T", "all"))
    data_dup_quantiles <- data_dup[, .(theoretical_q = quantile(theoretical, probs = c(0.25, 0.75), names = FALSE), value_q = quantile(value, probs = c(0.25, 0.75), names = FALSE)), by = base][, .(x1 = theoretical_q[1], x2 = theoretical_q[2], y1 = value_q[1], y2 = value_q[2]), by = base]
    ggplot(data_dup, aes(x = theoretical, y = value)) + geom_point(aes(color = orig_base)) +
        geom_abline(data = data_dup_quantiles, mapping = aes(slope = (y2 - y1) / (x2 - x1), intercept = y1 - ((y2 - y1) / (x2 - x1)) * x1)) +
        facet_wrap(~ base, labeller = label_both, ncol = 2) +
        ggtitle(title_text, subtitle = paste0("Sampled ", data[, .N], " points")) +
        labs(x = "Theoretical quantile", y = "Observed quantile", color = "base")
}

source("load_jun2018_data.saved.fst.basic.R")

coverage_thres <- 25

sample_size_max <- 20000
set.seed(1234)
ab_data <- ab_data[coverage >= coverage_thres, .(base, orig_base = base, value = tMean)][sample.int(min(.N, sample_size_max), size = sample_size_max)]
setorder(ab_data, value)
cd_data <- cd_data[coverage >= coverage_thres, .(base, orig_base = base, value = tMean)][sample.int(min(.N, sample_size_max), size = sample_size_max)]
setorder(cd_data, value)
k_normBy_ab_data <- k_normBy_ab_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres, .(base, orig_base = base, value = ipdRatio)][sample.int(min(.N, sample_size_max), size = sample_size_max)]
setorder(k_normBy_ab_data, value)
l_normBy_cd_data <- l_normBy_cd_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres, .(base, orig_base = base, value = ipdRatio)][sample.int(min(.N, sample_size_max), size = sample_size_max)]
setorder(l_normBy_cd_data, value)

pdf("kinetics_qq_plot.pdf", width = 9, height = 12)
qq_plot_per_base(ab_data, "Normal Q-Q plot of IPDs of AB (WGA, VC2010+OP50)")
qq_plot_per_base(ab_data[, "value" := .(log2(value))], "Normal Q-Q plot of log2 IPDs of AB (WGA, VC2010+OP50)")
qq_plot_per_base(cd_data, "Normal Q-Q plot of IPDs of CD (WGA, VC2010)")
qq_plot_per_base(cd_data[, "value" := .(log2(value))], "Normal Q-Q plot of log2 IPDs of CD (WGA, VC2010)")
qq_plot_per_base(k_normBy_ab_data, "Normal Q-Q plot of IPD ratios of K normalized by AB (VC2010+OP50)")
qq_plot_per_base(k_normBy_ab_data[, "value" := .(log2(value))], "Normal Q-Q plot of log2 IPD ratios of K normalized by AB (VC2010+OP50)")
qq_plot_per_base(l_normBy_cd_data, "Normal Q-Q plot of IPD ratios of L normalized by CD (VC2010)")
qq_plot_per_base(l_normBy_cd_data[, "value" := .(log2(value))], "Normal Q-Q plot of log2 IPD ratios of L normalized by CD (VC2010)")
invisible(dev.off())
