library(data.table)
library(ggplot2)
library(cowplot)

plot_rmse_cor_vs_coverage <- function(data1, name1, data2, name2) {
    # 'coverage' means valid IPD counts
    merged_data <- data.table(COV1 = data1[,coverage], COV2 = data2[,coverage], LOG2IPD1 = log2(data1[,tMean]), LOG2IPD2 = log2(data2[,tMean]))
    merged_data <- merged_data[is.finite(LOG2IPD1) & is.finite(LOG2IPD2), .(lower_cov = pmin(COV1, COV2), log2ipdr = LOG2IPD2 - LOG2IPD1, LOG2IPD1, LOG2IPD2)]
    max_lower_cov <- merged_data[,max(lower_cov)]
    stats_vs_coverage <- data.table(coverage_thres = 1:max_lower_cov, ipd_cor = 0, ipd_rmse = 0)
    for (i in 1:max_lower_cov) {
        current_data <- merged_data[lower_cov >= i]
        current_ipd_cor <- current_data[,cor(LOG2IPD1, LOG2IPD2)]
        current_ipd_rmse <- current_data[,sqrt(sum((LOG2IPD1 - LOG2IPD2)^2) / .N)]
        stats_vs_coverage[coverage_thres == i, c("ipd_cor", "ipd_rmse") := .(current_ipd_cor, current_ipd_rmse)]
    }
    if (stats_vs_coverage[,.N] > 0) {
        p <- ggplot(stats_vs_coverage, aes(coverage_thres, ipd_cor)) + geom_line() +
            ylab(bquote("Pearson's correlation coefficient of" ~ log[2] ~ "of IPDs")) +
            xlab(paste0("Lower threshold of the lower count of valid IPDs per strand\nbetween ", name1, " and ", name2)) +
            coord_cartesian(ylim = c(0, 1)) + theme_classic()
        print(p)
        p <- ggplot(stats_vs_coverage, aes(coverage_thres, ipd_rmse)) + geom_line() +
            ylab(bquote("Root mean squared error of" ~ log[2] ~ "of IPDs")) +
            xlab(paste0("Lower threshold of the lower count of valid IPDs per strand\nbetween ", name1, " and ", name2)) +
            coord_cartesian(ylim = c(0, stats_vs_coverage[, max(ipd_rmse)])) + theme_classic()
        print(p)
    }
}

source("load_jun2018_data.saved.fst.basic_wga.R")
ab_name <- "VC2010+OP50(WGA)"
cd_name <- "VC2010(WGA)"

pdf_width <- 6;
pdf_height <- 4;

pdf("ipdRatio.rmse_cor_vs_coverage.ab_cd.pdf", width = pdf_width, height = pdf_height, onefile = TRUE)
plot_rmse_cor_vs_coverage(ab_data, ab_name, cd_data, cd_name)
invisible(dev.off())
