library(data.table)

#score_thres <- 20
#ipdRatio_thres <- 2.0
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
library(ggExtra)
library(cowplot)

plot_ipd_for_2sets <- function(data, name, a.name, b.name) {
    ret_list <- NULL
    binwd <- 0.4
    x.label <- sprintf("%s log2 (IPD)", a.name)
    y.label <- sprintf("%s log2 (IPD)", b.name)
    # plot all the bases
    subdata <- data[a.tMean > 0 & b.tMean > 0]
    subdata[, c("a.log2_tMean", "b.log2_tMean") := .(log2(a.tMean), log2(b.tMean))]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, mean_a.log2_tMean = mean(a.log2_tMean), mean_b.log2_tMean = mean(b.log2_tMean)), by = a.base]
    p_linear_list <- NULL
    p_log_list <- NULL
    for (b in bases) {
        base_data <- subdata[a.base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,a.log2_tMean], base_data[,b.log2_tMean], method = "pearson")
            subdata_stat[a.base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
        p_linear <- ggplot(subdata, aes(a.log2_tMean, b.log2_tMean)) + geom_hex(binwidth = binwd) +
            geom_text(data = subdata_stat[a.base == b], mapping = aes(x = -Inf, y = Inf,
                label = sprintf("N = %d, cor = %.3g (p = %.3g)\nmean_x = %.3g, mean_y = %.3g", N, cor.value, p.value, mean_a.log2_tMean, mean_b.log2_tMean)), hjust = "inward", vjust = "inward") +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") + ggtitle(sprintf("Base: %s", b)) + xlab(x.label) + ylab(y.label)
        if(subdata[,.N] == 0){ p_linear <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
        p_log <- p_linear + scale_fill_continuous(trans = "log10")
        # TODO: replace ggMarginal
        p_linear <- ggMarginal(p_linear, type = "histogram", binwidth = binwd)
        p_log <- ggMarginal(p_log, type = "histogram", binwidth = binwd)
        p_linear_list <- c(p_linear_list, list(p_linear))
        p_log_list <- c(p_log_list, list(p_log))
    }
    p_linear_main <- plot_grid(plotlist = p_linear_list, align = "hv", axis = "lbtr", ncol = 2)
    p_log_main <- plot_grid(plotlist = p_log_list, align = "hv", axis = "lbtr", ncol = 2)
    p_title <- ggdraw() + draw_label(sprintf("%s: log2 (IPD) per base with Pearson's correlation tests", name), fontface = "bold")
    ret_list <- c(ret_list, list(plot_grid(p_title, p_linear_main, ncol = 1, rel_heights = c(0.1, 1))))
    ret_list <- c(ret_list, list(plot_grid(p_title, p_log_main, ncol = 1, rel_heights = c(0.1, 1))))

    # plot bases with a coverage threshold
    subdata <- subdata[a.coverage >= coverage_thres & b.coverage >= coverage_thres]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf, mean_a.log2_tMean = mean(a.log2_tMean), mean_b.log2_tMean = mean(b.log2_tMean)), by = a.base]
    p_linear_list <- NULL
    p_log_list <- NULL
    for (b in bases) {
        base_data <- subdata[a.base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,a.log2_tMean], base_data[,b.log2_tMean], method = "pearson")
            subdata_stat[a.base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
        p_linear <- ggplot(subdata, aes(a.log2_tMean, b.log2_tMean)) + geom_hex(binwidth = binwd) +
            geom_text(data = subdata_stat[a.base == b], mapping = aes(x = -Inf, y = Inf,
                label = sprintf("N = %d, cor = %.3g (p = %.3g)\nmean_x = %.3g, mean_y = %.3g", N, cor.value, p.value, mean_a.log2_tMean, mean_b.log2_tMean)), hjust = "inward", vjust = "inward") +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") + ggtitle(sprintf("Base: %s", b)) + xlab(x.label) + ylab(y.label)
        if(subdata[,.N] == 0){ p_linear <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
        p_log <- p_linear + scale_fill_continuous(trans = "log10")
        # TODO: replace ggMarginal
        p_linear <- ggMarginal(p_linear, type = "histogram", binwidth = binwd)
        p_log <- ggMarginal(p_log, type = "histogram", binwidth = binwd)
        p_linear_list <- c(p_linear_list, list(p_linear))
        p_log_list <- c(p_log_list, list(p_log))
    }
    p_linear_main <- plot_grid(plotlist = p_linear_list, align = "hv", axis = "lbtr", ncol = 2)
    p_log_main <- plot_grid(plotlist = p_log_list, align = "hv", axis = "lbtr", ncol = 2)
    p_title <- ggdraw() + draw_label(sprintf("%s: log2 (IPD) per base with coverage >= %g with Pearson's correlation tests", name, coverage_thres), fontface = "bold")
    ret_list <- c(ret_list, list(plot_grid(p_title, p_linear_main, ncol = 1, rel_heights = c(0.1, 1))))
    ret_list <- c(ret_list, list(plot_grid(p_title, p_log_main, ncol = 1, rel_heights = c(0.1, 1))))
    cat(sprintf("Generated plots for %s\n", name))
    return(ret_list)
}

#source("load_jun2018_data.saved.fst.wga.R")
source("load_jun2018_data.saved.fst.R")

#union_data <- cbind(a = ab_data, b = cd_data)
union_data1 <- data.table(a.refName = ab_data[,refName], a.base = ab_data[,base], a.tMean = ab_data[,tMean], a.coverage = ab_data[,coverage], b.tMean = cd_data[,tMean], b.coverage = cd_data[,coverage])
union_data2 <- data.table(a.refName = ab_data[,refName], a.base = ab_data[,base], a.tMean = ab_data[,tMean], a.coverage = ab_data[,coverage], b.tMean = k_data[,tMean], b.coverage = k_data[,coverage])
union_data3 <- data.table(a.refName = cd_data[,refName], a.base = cd_data[,base], a.tMean = cd_data[,tMean], a.coverage = cd_data[,coverage], b.tMean = l_data[,tMean], b.coverage = l_data[,coverage])
setkey(union_data1, a.refName, a.base)
setkey(union_data2, a.refName, a.base)
setkey(union_data3, a.refName, a.base)

pdf_width <- 14
pdf_height <- 14

pdf(file = NULL, width = pdf_width, height = pdf_height)
p_all <- NULL
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data1, "Intersection of AB_WGA_VC2010_OP50 and CD_WGA_VC2010", "AB_WGA_VC2010_OP50", "CD_WGA_VC2010")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data1[a.refName != ecoli_chr], "Intersection of AB_WGA_VC2010_OP50 and CD_WGA_VC2010 in C. elegans", "AB_WGA_VC2010_OP50", "CD_WGA_VC2010")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data1[a.refName == ecoli_chr], "Intersection of AB_WGA_VC2010_OP50 and CD_WGA_VC2010 in E. coli", "AB_WGA_VC2010_OP50", "CD_WGA_VC2010")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data2, "Intersection of AB_WGA_VC2010_OP50 and K_native_VC2010_OP50", "AB_WGA_VC2010_OP50", "K_native_VC2010_OP50")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data2[a.refName != ecoli_chr], "Intersection of AB_WGA_VC2010_OP50 and K_native_VC2010_OP50 in C. elegans", "AB_WGA_VC2010_OP50", "K_native_VC2010_OP50")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data2[a.refName == ecoli_chr], "Intersection of AB_WGA_VC2010_OP50 and K_native_VC2010_OP50 in E. coli", "AB_WGA_VC2010_OP50", "K_native_VC2010_OP50")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data3, "Intersection of CD_WGA_VC2010 and L_native_VC2010", "CD_WGA_VC2010", "L_native_VC2010")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data3[a.refName != ecoli_chr], "Intersection of CD_WGA_VC2010 and L_native_VC2010 in C. elegans", "CD_WGA_VC2010", "L_native_VC2010")))
p_all <- c(p_all, list(plot_ipd_for_2sets(union_data3[a.refName == ecoli_chr], "Intersection of CD_WGA_VC2010 and L_native_VC2010 in E. coli", "CD_WGA_VC2010", "L_native_VC2010")))
invisible(dev.off())
pdf(sprintf("compare_ipd.coverage%g.pdf", coverage_thres), width = pdf_width, height = pdf_height)
for (p in p_all){
    for (psub in p) { print(psub) }
}
invisible(dev.off())

