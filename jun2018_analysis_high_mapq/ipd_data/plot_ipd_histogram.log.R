library(data.table)
library(ggplot2)

plot_ipd_histogram <- function(data, sample, xthres) {
    data[, "log2tMean" := .(log2(tMean))]
    xthres <- log2(xthres)
    binwidth1 <- 0.2
    #data_stat <- data[,.(data_size = .N, data_mean = mean(log2tMean)), by = base]
    #data_stat <- rbind(data_stat, data[,.(data_size = .N, data_mean = mean(log2tMean), base = "all")])
    #p_all <- ggplot(data) + geom_histogram(aes(log2tMean), binwidth = binwidth1) + xlab(bquote(log[2]~IPD)) +
    #    geom_text(data = data_stat[base=="all"], mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nMean = %.3g", data_size, data_mean)), hjust = "inward", vjust = "inward")
    #p_base <- p_all + facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
    #    geom_text(data = data_stati[base!="all"], mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nMean = %.3g", data_size, data_mean)), hjust = "inward", vjust = "inward")
    p_all <- ggplot(data) + geom_histogram(aes(log2tMean), binwidth = binwidth1) + xlab(bquote(log[2]~IPD)) + coord_cartesian(xlim = c(-4, 6)) +
        theme(panel.grid.minor = element_blank()) + geom_vline(xintercept = xthres, linetype = "dashed")
    p_base <- p_all + facet_wrap(vars(base = base), ncol = 2, labeller = label_both)
    pdf(sprintf("ipd_histogram.log.all.%s.pdf", sample), height = 2, width = 3)
    if(data[,.N] > 0){
        print(p_all)
        print(p_all + scale_y_continuous(trans = "log10"))
    }else{
        ggplot() + annotate("text", x = 1, y = 1, label = "No data")
    }
    invisible(dev.off())
    pdf(sprintf("ipd_histogram.log.base.%s.pdf", sample), height = 3, width = 3)
    if(data[,.N] > 0){
        print(p_base)
        print(p_base + scale_y_continuous(trans = "log10"))
    }else{
        ggplot() + annotate("text", x = 1, y = 1, label = "No data")
    }
    invisible(dev.off())
}

plot_modelPrediction_histogram <- function(data, sample) {
    data[, "log2modelPrediction" := .(log2(modelPrediction))]
    binwidth1 <- 0.2
    #data_stat <- data[,.(data_size = .N, data_mean = mean(log2modelPrediction)), by = base]
    #data_stat <- rbind(data_stat, data[,.(data_size = .N, data_mean = mean(log2modelPrediction), base = "all")])
    #p_all <- ggplot(data) + geom_histogram(aes(log2modelPrediction), binwidth = binwidth1) + xlab(bquote(log[2]~IPD~estimate)) +
    #    geom_text(data = data_stat[base=="all"], mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nMean = %.3g", data_size, data_mean)), hjust = "inward", vjust = "inward")
    #p_base <- p_all + facet_wrap(vars(base = base), ncol = 2, labeller = label_both) +
    #    geom_text(data = data_stat[base!="all"], mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d\nMean = %.3g", data_size, data_mean)), hjust = "inward", vjust = "inward")
    p_all <- ggplot(data) + geom_histogram(aes(log2modelPrediction), binwidth = binwidth1) + xlab(bquote(log[2]~IPD~estimate)) + coord_cartesian(xlim = c(-4, 6)) + theme(panel.grid.minor = element_blank())
    p_base <- p_all + facet_wrap(vars(base = base), ncol = 2, labeller = label_both)
    pdf(sprintf("model_prediction_histogram.all.%s.pdf", sample), height = 2, width = 3)
    if(data[,.N] > 0){
        print(p_all)
    }else{
        ggplot() + annotate("text", x = 1, y = 1, label = "No data")
    }
    invisible(dev.off())
    pdf(sprintf("model_prediction_histogram.base.%s.pdf", sample), height = 3, width = 3)
    if(data[,.N] > 0){
        print(p_base)
    }else{
        ggplot() + annotate("text", x = 1, y = 1, label = "No data")
    }
    invisible(dev.off())
}

source("../load_jun2018_data.saved.fst.basic_wga.R", chdir = TRUE)
source("../load_and_save_PD2182_data.saved.fst.PD2182.R", chdir = TRUE)
source("../load_and_save_PD2182sequel_data.saved.fst.PD2182sequel.R", chdir = TRUE)

# threshold of valid IPD counts (coverage)
thres <- 25
# raw IPD values at the top 1% and the bottom 1%, including reads on E. coli
ipd_ab_top1 <- 3.08223
ipd_cd_top1 <- 3.05611
ipd_PD2182_top1 <- 3.71169
ipd_PD2182sequel_top1 <- 4.17400
ipd_ab_bottom1 <- 0.30158
ipd_cd_bottom1 <- 0.29231
ipd_PD2182_bottom1 <- 0.38372
ipd_PD2182sequel_bottom1 <- 0.35110
# max of IPD with coverage >= 25, which was in sample cd
ipd_max <- 61.35193

setkey(ab_data, refName, base)
setkey(cd_data, refName, base)
setkey(PD2182_data, refName, base)
setkey(PD2182sequel_data, refName, base)

ab_data1 <- ab_data[refName != ecoli_chr & coverage >= thres]
cd_data1 <- cd_data[refName != ecoli_chr & coverage >= thres]
PD2182_data1 <- PD2182_data[refName != ecoli_chr & coverage >= thres]
PD2182sequel_data1 <- PD2182sequel_data[refName != ecoli_chr & coverage >= thres]
plot_ipd_histogram(ab_data1, "ab_deep_celegans", c(ipd_ab_bottom1, ipd_ab_top1))
plot_ipd_histogram(cd_data1, "cd_deep_celegans", c(ipd_cd_bottom1, ipd_cd_top1))
plot_ipd_histogram(PD2182_data1, "PD2182_deep_celegans", c(ipd_PD2182_bottom1, ipd_PD2182_top1))
plot_ipd_histogram(PD2182sequel_data1, "PD2182sequel_deep_celegans", c(ipd_PD2182sequel_bottom1, ipd_PD2182sequel_top1))
#plot_modelPrediction_histogram(ab_data1, "ab_deep_celegans")
#plot_modelPrediction_histogram(cd_data1, "cd_deep_celegans")
#plot_modelPrediction_histogram(PD2182_data1, "PD2182_deep_celegans")
#plot_modelPrediction_histogram(PD2182sequel_data1, "PD2182sequel_deep_celegans")
