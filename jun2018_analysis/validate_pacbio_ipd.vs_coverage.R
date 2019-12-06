library(data.table)
library(ggplot2)
library(cowplot)

plot_pacbio_rmse_cor_vs_coverage <- function(data1, name1) {
    # 'coverage' means valid IPD counts
    max_cov <- data1[,max(coverage)]
    log2_data <- data1[tMean > 0 & is.finite(tMean), .(coverage, base, log2tMean = log2(tMean), log2modelPrediction = log2(modelPrediction))]
    stats_vs_coverage <- data.table(coverage_thres = 1:max_cov, ipd_cor = 0, ipd_rmse = 0)
    for (i in 1:max_cov) {
        current_data <- log2_data[coverage >= i]
        current_ipd_cor <- current_data[,cor(log2tMean, log2modelPrediction)]
        current_ipd_rmse <- current_data[,sqrt(sum((log2tMean - log2modelPrediction)^2) / .N)]
        stats_vs_coverage[coverage_thres == i, c("ipd_cor", "ipd_rmse") := .(current_ipd_cor, current_ipd_rmse)]
    }
    if (stats_vs_coverage[,.N] > 0) {
        p <- ggplot(stats_vs_coverage, aes(coverage_thres, ipd_cor)) + geom_line() +
            ylab(bquote("Pearson's correlation coefficient of" ~ log[2] ~ "of IPDs")) +
            xlab(paste0("Lower threshold of the valid IPD count per strand in ", name1)) +
            coord_cartesian(ylim = c(0, 1)) + theme_classic()
        print(p)
        p <- ggplot(stats_vs_coverage, aes(coverage_thres, ipd_rmse)) + geom_line() +
            ylab(bquote("Root mean squared error of" ~ log[2] ~ "of IPDs")) +
            xlab(paste0("Lower threshold of the valid IPD count per strand in ", name1)) +
            coord_cartesian(ylim = c(0, stats_vs_coverage[, max(ipd_rmse)])) + theme_classic()
        print(p)
    }
}

source("load_jun2018_data.saved.fst.basic_wga.R")
ab_name <- "VC2010+OP50(WGA)"
cd_name <- "VC2010(WGA)"

pdf_width <- 6;
pdf_height <- 4;

pdf("validate_pacbio_ipd.vs_coverage.pdf", width = pdf_width, height = pdf_height, onefile = TRUE)
plot_pacbio_rmse_cor_vs_coverage(ab_data, ab_name)
plot_pacbio_rmse_cor_vs_coverage(cd_data, cd_name)
invisible(dev.off())
