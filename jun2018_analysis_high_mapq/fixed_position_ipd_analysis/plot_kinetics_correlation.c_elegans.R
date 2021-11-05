library(data.table)

mean_log2value <- function(kinetics) {
    if (kinetics[, .N] == 0) {
        data.table(position = numeric(0), strand = character(0), m = numeric(0))
    } else {
        # Use only the motif region
        kinetics[value > 0 & is.finite(value) & substr(label, 1, 1) == "m"][, .(position, strand, log2value = log2(value))][, .(m = mean(log2value)), by = .(position, strand)]
    }
}

ab_ipd <- fread("motif_ipd.ab.c_elegans.csv")
cd_ipd <- fread("motif_ipd.cd.c_elegans.csv")
PD2182_ipd <- fread("motif_ipd.PD2182.c_elegans.csv")
PD2182sequel_ipd <- fread("motif_ipd.PD2182sequel.c_elegans.csv")
k_normBy_ab_ipdratio <- fread("motif_ipdratio.k_normBy_ab.c_elegans.csv")
l_normBy_cd_ipdratio <- fread("motif_ipdratio.l_normBy_cd.c_elegans.csv")
k_ipd <- fread("motif_ipd.k.c_elegans.csv")
l_ipd <- fread("motif_ipd.l.c_elegans.csv")

ab_mean <- mean_log2value(ab_ipd)[, .(position, strand, ab = m)]
cd_mean <- mean_log2value(cd_ipd)[, .(position, strand, cd = m)]
PD2182_mean <- mean_log2value(PD2182_ipd)[, .(position, strand, PD2182 = m)]
PD2182sequel_mean <- mean_log2value(PD2182sequel_ipd)[, .(position, strand, PD2182sequel = m)]
k_normBy_ab_mean <- mean_log2value(k_normBy_ab_ipdratio)[, .(position, strand, k_normBy_ab = m)]
l_normBy_cd_mean <- mean_log2value(l_normBy_cd_ipdratio)[, .(position, strand, l_normBy_cd = m)]
k_mean <- mean_log2value(k_ipd)[, .(position, strand, k = m)]
l_mean <- mean_log2value(l_ipd)[, .(position, strand, l = m)]

# Rename columns
setnames(ab_mean, "ab", "Replicate 1\n/WGA")
setnames(cd_mean, "cd", "Replicate 2\n/WGA")
setnames(PD2182_mean, "PD2182", "PD2182\n(RS II)\n/native")
setnames(PD2182sequel_mean, "PD2182sequel", "PD2182\n/native")
setnames(k_normBy_ab_mean, "k_normBy_ab", "Replicate 1")
setnames(l_normBy_cd_mean, "l_normBy_cd", "Replicate 2")
setnames(k_mean, "k", "Replicate 1\n/native")
setnames(l_mean, "l", "Replicate 2\n/native")

ipd_table <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(ab_mean, cd_mean, k_mean, l_mean, PD2182_mean, PD2182sequel_mean))
ipd_table[, c("position", "strand") := NULL]
ipd_table2 <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(ab_mean, cd_mean))
ipd_table2[, c("position", "strand") := NULL]
ipdratio_table <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(k_normBy_ab_mean, l_normBy_cd_mean))
ipdratio_table[, c("position", "strand") := NULL]

# function: corrplot.na
source("../plot_kinetics_correlation.corrplot.na.R")

pdf("plot_kinetics_correlation.c_elegans.pdf", onefile = TRUE, width = 8, height = 8)
corrplot.na(ipd_table)
corrplot.na(ipd_table2)
corrplot.na(ipdratio_table)
invisible(dev.off())
