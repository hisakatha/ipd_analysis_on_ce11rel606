library(data.table)
#library(ggplot2)
library(corrplot)

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
PD2182sequel_ipd <- fread("motif_ipd.PD2182sequel.c_elegans.csv")
k_normBy_ab_ipdratio <- fread("motif_ipdratio.k_normBy_ab.c_elegans.csv")
l_normBy_cd_ipdratio <- fread("motif_ipdratio.l_normBy_cd.c_elegans.csv")
k_ipd <- fread("motif_ipd.k.c_elegans.csv")
l_ipd <- fread("motif_ipd.l.c_elegans.csv")
ab_ipd_ecoli <- fread("motif_ipd.ab.e_coli.csv")
cd_ipd_ecoli <- fread("motif_ipd.cd.e_coli.csv")
k_normBy_ab_ipdratio_ecoli <- fread("motif_ipdratio.k_normBy_ab.e_coli.csv")
l_normBy_cd_ipdratio_ecoli <- fread("motif_ipdratio.l_normBy_cd.e_coli.csv")
k_ipd_ecoli <- fread("motif_ipd.k.e_coli.csv")
l_ipd_ecoli <- fread("motif_ipd.l.e_coli.csv")

ab_mean <- mean_log2value(ab_ipd)[, .(position, strand, ab = m)]
cd_mean <- mean_log2value(cd_ipd)[, .(position, strand, cd = m)]
PD2182sequel_mean <- mean_log2value(PD2182sequel_ipd)[, .(position, strand, PD2182sequel = m)]
k_normBy_ab_mean <- mean_log2value(k_normBy_ab_ipdratio)[, .(position, strand, k_normBy_ab = m)]
l_normBy_cd_mean <- mean_log2value(l_normBy_cd_ipdratio)[, .(position, strand, l_normBy_cd = m)]
k_mean <- mean_log2value(k_ipd)[, .(position, strand, k = m)]
l_mean <- mean_log2value(l_ipd)[, .(position, strand, l = m)]
ab_mean_ecoli <- mean_log2value(ab_ipd_ecoli)[, .(position, strand, ab = m)]
cd_mean_ecoli <- mean_log2value(cd_ipd_ecoli)[, .(position, strand, cd = m)]
k_normBy_ab_mean_ecoli <- mean_log2value(k_normBy_ab_ipdratio_ecoli)[, .(position, strand, k_normBy_ab = m)]
l_normBy_cd_mean_ecoli <- mean_log2value(l_normBy_cd_ipdratio_ecoli)[, .(position, strand, l_normBy_cd = m)]
k_mean_ecoli <- mean_log2value(k_ipd_ecoli)[, .(position, strand, k = m)]
l_mean_ecoli <- mean_log2value(l_ipd_ecoli)[, .(position, strand, l = m)]

# Rename columns
setnames(ab_mean, "ab", "VC2010+OP50\n/WGA\n/C. elegans")
setnames(cd_mean, "cd", "VC2010\n/WGA\n/C. elegans")
setnames(PD2182sequel_mean, "PD2182sequel", "PD2182\n(PacBio Sequel)\n/C. elegans")
setnames(k_normBy_ab_mean, "k_normBy_ab", "VC2010+OP50\n/C. elegans")
setnames(l_normBy_cd_mean, "l_normBy_cd", "VC2010\n/C. elegans")
setnames(k_mean, "k", "VC2010+OP50\n/native\n/C. elegans")
setnames(l_mean, "l", "VC2010\n/native\n/C. elegans")
setnames(ab_mean_ecoli, "ab", "VC2010+OP50\n/WGA\n/E. coli")
setnames(cd_mean_ecoli, "cd", "VC2010\n/WGA\n/E. coli")
setnames(k_normBy_ab_mean_ecoli, "k_normBy_ab", "VC2010+OP50\n/E. coli")
setnames(l_normBy_cd_mean_ecoli, "l_normBy_cd", "VC2010\n/E. coli")
setnames(k_mean_ecoli, "k", "VC2010+OP50\n/native\n/E. coli")
setnames(l_mean_ecoli, "l", "VC2010\n/native\n/E. coli")

ipd_table <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(ab_mean, cd_mean, k_mean, l_mean, PD2182sequel_mean, ab_mean_ecoli, cd_mean_ecoli, k_mean_ecoli, l_mean_ecoli))
ipd_table[, c("position", "strand") := NULL]
ipd_table2 <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(ab_mean, cd_mean, ab_mean_ecoli, cd_mean_ecoli))
ipd_table2[, c("position", "strand") := NULL]
ipdratio_table <- Reduce(function(x, y) merge(x, y, all = TRUE, by = c("position", "strand")), list(k_normBy_ab_mean, l_normBy_cd_mean, k_normBy_ab_mean_ecoli, l_normBy_cd_mean_ecoli))
ipdratio_table[, c("position", "strand") := NULL]

corrplot.na <- function (input_table) {
    if (nrow(input_table) == 0) {
        plot.new()
        title(sprintf("No data\nin %s", paste0(colnames(input_table), collapse = ", ")), line = -12)
        return()
    }
    cors <- cor(input_table)
    if (all(is.na(cors[upper.tri(cors)]))) {
        plot.new()
        title(sprintf("Not enough data for correlation matrix\nin %s", paste0(colnames(input_table), collapse = ", ")), line = -12)
        return()
    }
    if (any(is.na(input_table))) {
        corrplot(cors, method = "color", type = "upper", tl.col = "black", na.label = "NA")
        corrplot.mixed(cors, lower = "number", upper = "color",
            tl.col = "black", tl.cex = 0.6, lower.col = "black", na.label = "NA")
    } else {
        pvalues <- cor.mtest(input_table)$p
        corrplot(cors, p.mat = pvalues, method = "color", type = "upper",
            sig.level = c(.001, .01, .05), pch.cex = .9,
            insig = "label_sig", pch.col = "white", tl.col = "black", na.label = "NA")
        pvalues[lower.tri(pvalues)] <- NA
        corrplot.mixed(cors, p.mat = pvalues, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", lower = "number", upper = "color",
            tl.col = "black", tl.cex = 0.6, pch.cex = 0.9, pch.col = "white", lower.col = "black", na.label = "NA")
    }
}

pdf("plot_kinetics_correlation.celegans_ecoli_subset1.pdf", onefile = TRUE, width = 11, height = 11)
corrplot.na(ipd_table)
corrplot.na(ipd_table2)
corrplot.na(ipdratio_table)
invisible(dev.off())
