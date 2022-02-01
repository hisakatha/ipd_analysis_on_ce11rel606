library(tidyverse)
library(data.table)

d_celegans <- fread("motif_ipd.jun2021_pd1074_ccs.c_elegans.csv")
d_ecoli <- fread("motif_ipd.jun2021_pd1074_ccs.e_coli.csv")

target_celegans <- d_celegans[position == 22 & strand == "+" & value > 0]
target_ecoli <- d_ecoli[position == 22 & strand == "+" & value > 0]

fwrite(target_celegans[order(-value)], file = "motif_ipd.jun2021_pd1074_ccs.c_elegans.sorted.csv")
fwrite(target_ecoli[order(-value)], file = "motif_ipd.jun2021_pd1074_ccs.e_coli.sorted.csv")

# Plots for IPD
title1 <- "IPD at G(A)TC in C. elegans"
title2 <- "IPD at G(A)TC in E. coli"
labelx <- "IPD"
log_breaks <- c(0, 10, 100, 1000, 10000, 100000)
xlim1 <- c(0, 23)

p1 <- ggplot(target_celegans) + geom_histogram(aes(value), binwidth = 1, boundary = 0) + xlab(labelx) + coord_cartesian(xlim = xlim1) + ggtitle(title1)
p1.2 <- p1 + scale_y_continuous(trans = "log1p", breaks = log_breaks)
p1.3 <- ggplot(target_celegans) + geom_histogram(aes(value, fill = factor(ref_strand)), position = "identity", alpha = 0.7, binwidth = 1, boundary = 0) + xlab(labelx) + coord_cartesian(xlim = xlim1) + ggtitle(title1) + theme_bw() + scale_y_continuous(trans = "log1p", breaks = log_breaks)

p2 <- ggplot(target_ecoli) + geom_histogram(aes(value), binwidth = 1, boundary = 0) + xlab(labelx) + coord_cartesian(xlim = xlim1) + ggtitle(title2)
p2.2 <- p2 + scale_y_continuous(trans = "log1p", breaks = log_breaks)
p2.3 <- ggplot(target_ecoli) + geom_histogram(aes(value, fill = factor(ref_strand)), position = "identity", alpha = 0.7, binwidth = 1, boundary = 0) + xlab(labelx) + coord_cartesian(xlim = xlim1) + ggtitle(title2) + theme_bw() + scale_y_continuous(trans = "log1p", breaks = log_breaks)

# Plots for coverage
title1 <- "CCS coverage at G(A)TC in C. elegans"
title2 <- "CCS coverage at G(A)TC in E. coli"
labelx <- "CCS coverage"
log_breaks <- c(0, 10, 100, 1000, 10000, 100000)
xlim1 <- c(0, NA)

p.cov.c <- ggplot(target_celegans) + geom_histogram(aes(coverage), binwidth = 5, center = 0) + xlab(labelx) + coord_cartesian(xlim = xlim1) + ggtitle(title1)
p.cov.c.2 <- p.cov.c + scale_y_continuous(trans = "log1p", breaks = log_breaks)

p.cov.e <- ggplot(target_ecoli) + geom_histogram(aes(coverage), binwidth = 1, center = 0) + xlab(labelx) + coord_cartesian(xlim = xlim1) + ggtitle(title2)
p.cov.e.2 <- p.cov.e + scale_y_continuous(trans = "log1p", breaks = log_breaks)

pdf_path <- "plot_target_base.jun2021_pd1074_ccs.pdf"
cairo_pdf(pdf_path, height = 4, width = 4, onefile = TRUE)
print(p1)
print(p1.2)
print(p1.3)
print(p2)
print(p2.2)
print(p2.3)
print(p.cov.c)
print(p.cov.c.2)
print(p.cov.e)
print(p.cov.e.2)
invisible(dev.off())
