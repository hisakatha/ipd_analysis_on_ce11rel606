library(data.table)
library(ggplot2)

plot_means <- function(data, ylab_text, title_text, ylim) {
    ggplot(data, aes(base, mean)) + geom_point() +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
        ggtitle(title_text) +
        ylab(ylab_text) +
        theme(panel.grid = element_blank()) +
        scale_y_continuous(limits = ylim)
}

title_ab <- "VC2010+OP50/WGA"
title_cd <- "VC2010/WGA"
title_k <- "VC2010+OP50/native"
title_l <- "VC2010/native"
title_PD2182 <- "PD2182 (PacBio RS II)"
title_PD2182sequel <- "PD2182 (PacBio Sequel)"
ylab_text <- "IPD (mean +/- standard deviation)"
type_text <- "IPD"
pdf_height <- 3
pdf_width <- 3

data_celegans <- fread("kinetics_stats.c_elegans.csv")
ylim <- c(0, 2.5)
pdf("kinetics_stats.plot.c_elegans.pdf", height = pdf_height, width = pdf_width)
plot_means(data_celegans[sample == "ab_deep" & base != "all" & type == type_text], ylab_text, title_ab, ylim)
plot_means(data_celegans[sample == "cd_deep" & base != "all" & type == type_text], ylab_text, title_cd, ylim)
plot_means(data_celegans[sample == "k_deep" & base != "all" & type == type_text], ylab_text, title_k, ylim)
plot_means(data_celegans[sample == "l_deep" & base != "all" & type == type_text], ylab_text, title_l, ylim)
plot_means(data_celegans[sample == "PD2182_deep" & base != "all" & type == type_text], ylab_text, title_PD2182, ylim)
plot_means(data_celegans[sample == "PD2182sequel_deep" & base != "all" & type == type_text], ylab_text, title_PD2182sequel, ylim)
invisible(dev.off())

data_ecoli <- fread("kinetics_stats.e_coli.csv")
ylim <- c(-0.6, 4)
pdf("kinetics_stats.plot.e_coli.pdf", height = pdf_height, width = pdf_width)
plot_means(data_ecoli[sample == "ab_deep" & base != "all" & type == type_text], ylab_text, title_ab, ylim)
plot_means(data_ecoli[sample == "cd_deep" & base != "all" & type == type_text], ylab_text, title_cd, ylim)
plot_means(data_ecoli[sample == "k_deep" & base != "all" & type == type_text], ylab_text, title_k, ylim)
plot_means(data_ecoli[sample == "l_deep" & base != "all" & type == type_text], ylab_text, title_l, ylim)
plot_means(data_ecoli[sample == "PD2182_deep" & base != "all" & type == type_text], ylab_text, title_PD2182, ylim)
plot_means(data_ecoli[sample == "PD2182sequel_deep" & base != "all" & type == type_text], ylab_text, title_PD2182sequel, ylim)
invisible(dev.off())
