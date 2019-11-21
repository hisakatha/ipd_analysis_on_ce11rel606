library(data.table)
#library(ggplot2)
library(corrplot)

mean_log2value <- function(kinetics) {
    if (kinetics[, .N] == 0) {
        data.table(position = numeric(0), strand = character(0), m = numeric(0))
    } else {
        kinetics[value > 0 & is.finite(value)][, .(position, strand, log2value = log2(value))][, .(m = mean(log2value)), by = .(position, strand)]
    }
}

ab_ipd <- fread("motif_ipd.ab.csv")
cd_ipd <- fread("motif_ipd.cd.csv")
PD2182_ipd <- fread("motif_ipd.PD2182.csv")
PD2182sequel_ipd <- fread("motif_ipd.PD2182sequel.csv")
k_normBy_ab_ipdratio <- fread("motif_ipdratio.k_normBy_ab.csv")
l_normBy_cd_ipdratio <- fread("motif_ipdratio.l_normBy_cd.csv")

ab_mean <- mean_log2value(ab_ipd)[, .(position, strand, ab = m)]
cd_mean <- mean_log2value(cd_ipd)[, .(position, strand, cd = m)]
PD2182_mean <- mean_log2value(PD2182_ipd)[, .(position, strand, PD2182 = m)]
PD2182sequel_mean <- mean_log2value(PD2182sequel_ipd)[, .(position, strand, PD2182sequel = m)]
k_normBy_ab_mean <- mean_log2value(k_normBy_ab_ipdratio)[, .(position, strand, k_normBy_ab = m)]
l_normBy_cd_mean <- mean_log2value(l_normBy_cd_ipdratio)[, .(position, strand, l_normBy_cd = m)]

# Rename columns
setnames(ab_mean, "ab", "VC2010+OP50\n(WGA)")
setnames(cd_mean, "cd", "VC2010\n(WGA)")
setnames(PD2182_mean, "PD2182", "PD2182\n(PacBio RS II)")
setnames(PD2182sequel_mean, "PD2182sequel", "PD2182\n(PacBio Sequel)")
setnames(k_normBy_ab_mean, "k_normBy_ab", "VC2010+OP50")
setnames(l_normBy_cd_mean, "l_normBy_cd", "VC2010")

ipd_table <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(ab_mean, cd_mean, PD2182_mean, PD2182sequel_mean))
ipd_table[, c("position", "strand") := NULL]
ipd_table2 <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(ab_mean, cd_mean))
ipd_table2[, c("position", "strand") := NULL]
ipdratio_table <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(k_normBy_ab_mean, l_normBy_cd_mean))
ipdratio_table[, c("position", "strand") := NULL]

corrplot.na <- function (input_table) {
    if (nrow(input_table) == 0) {
        plot.new()
        title(sprintf("No data in %s", paste0(colnames(input_table), collapse = ", ")))
        return()
    }
    cors <- cor(input_table)
    if (any(is.na(input_table))) {
        corrplot(cors, method = "color", type = "upper")
        corrplot(cors, method = "number", type = "upper")
    } else {
        pvalues <- cor.mtest(input_table)$p
        corrplot(cors, p.mat = pvalues, method = "color", type = "upper",
                 sig.level = c(.001, .01, .05), pch.cex = .9,
                 insig = "label_sig", pch.col = "white", tl.col = "black")
        pvalues[lower.tri(pvalues)] <- NA
        corrplot.mixed(cors, p.mat = pvalues, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", lower = "number", upper = "color",
                       tl.col = "black", tl.cex = 0.8, pch.cex = 0.9, pch.col = "white")
    }
}

pdf("plot_kinetics_correlation.pdf", onefile = TRUE, width = 6, height = 6)
corrplot.na(ipd_table)
corrplot.na(ipd_table2)
corrplot.na(ipdratio_table)
invisible(dev.off())
