library(data.table)
library(ggplot2)
library(cowplot)

ipdRatio_vs_coverage <- function(data1, name1, data2, name2) {
    # 'coverage' means valid IPD counts
    merged_data <- data.table(COV1 = data1[,coverage], COV2 = data2[,coverage], LOG2IPD1 = log2(data1[,tMean]), LOG2IPD2 = log2(data2[,tMean]))
    plot_data <- merged_data[is.finite(LOG2IPD1) & is.finite(LOG2IPD2), .(lower_cov = pmin(COV1, COV2), log2ipdr = LOG2IPD2 - LOG2IPD1)]
    if (merged_data[,.N] > 0) {
        p <- ggplot(plot_data, aes(lower_cov, log2ipdr)) + geom_hex(binwidth=c(5,0.2)) +
            ylab(bquote(log[2] ~ "(IPD in" ~ .(name2) ~ ") / (IPD in" ~ .(name1) ~ ")")) +
            xlab(paste0("The lower count of valid IPDs per strand between ", name1, " and ", name2)) +
            geom_hline(yintercept = 0, lty = 2) + theme(legend.position = "left") +
            scale_fill_continuous(trans = "log10")
        marginal_size <- 0.2
        p_xhist <- axis_canvas(p, axis = "x") + geom_histogram(data = plot_data, aes(lower_cov), size = marginal_size)
        p_yhist <- axis_canvas(p, axis = "y", coord_flip = TRUE) + geom_histogram(data = plot_data, aes(log2ipdr), size = marginal_size) + coord_flip()
        p_tmp <- insert_xaxis_grob(p, p_xhist, grid::unit(marginal_size, "null"), position = "top")
        p2 <- insert_yaxis_grob(p_tmp, p_yhist, grid::unit(marginal_size, "null"), position = "right")
    } else {
        p2 <- NULL;
    }
    return(p2)
}

source("load_jun2018_data.saved.fst.basic_wga.R")
ab_name <- "VC2010+OP50(WGA)"
cd_name <- "VC2010(WGA)"

pdf_width <- 7;
pdf_height <- 7;

stopifnot(ab_data[refName == "", .N] == 0, cd_data[refName == "", .N] == 0, all(ab_data[, refName] == cd_data[, refName]))
ab_data <- ab_data[refName != ecoli_chr]
cd_data <- cd_data[refName != ecoli_chr]
pdf("ipdRatio.vs_coverage.ab_cd.c_elegans.pdf", width = pdf_width, height = pdf_height, onefile = TRUE);
result1 <- ipdRatio_vs_coverage(ab_data, ab_name, cd_data, cd_name);
ggdraw(result1);
invisible(dev.off());
