library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
d <- fread(data_path)

title1 <- args[2]
xlab1 <- bquote("IPD in" ~ log[2] ~ "scale")
xbreaks <- -3:5
xlabels <- sprintf("%.3g", 2 ^ xbreaks)
d[base == "N", "base" := .("All")]
d$base <- factor(d$base, levels = c("A", "C", "G", "T", "All"))

pdf(paste0(data_path, ".pdf"), height = 4.5, width = 4.5)
ggplot(d, aes((lower_exclusive + upper_inclusive) / 2, y = count, color = base, shape = base)) +
    geom_point() + geom_line() +
    ylab("Frequency") + xlab(xlab1) +
    ggtitle(title1) +
    theme_classic() +
    scale_x_continuous(breaks = xbreaks, labels = xlabels, limits = c(min(xbreaks), max(xbreaks))) +
    theme(legend.position = c(0.98, 0.98), legend.title = element_blank(), legend.justification = c("right", "top")) +
    scale_color_manual(values = c("red3", "darkblue", "green3", "yellow3", "gray60")) +
    theme(plot.title = element_text(size = 10), legend.background = element_blank())
ggplot(d, aes((lower_exclusive + upper_inclusive) / 2, y = count, color = base, shape = base)) +
    geom_point() + geom_line() +
    ylab(bquote("Frequency in" ~ log[10] ~ "scale")) + xlab(xlab1) +
    ggtitle(title1) +
    theme_classic() +
    scale_x_continuous(breaks = xbreaks, labels = xlabels, limits = c(min(xbreaks), max(xbreaks))) +
    theme(legend.position = c(0.98, 0.98), legend.title = element_blank(), legend.justification = c("right", "top")) +
    scale_color_manual(values = c("red3", "darkblue", "green3", "yellow3", "gray60")) +
    scale_y_continuous(trans = "log10") +
    theme(plot.title = element_text(size = 10), legend.background = element_blank())
invisible(dev.off())
