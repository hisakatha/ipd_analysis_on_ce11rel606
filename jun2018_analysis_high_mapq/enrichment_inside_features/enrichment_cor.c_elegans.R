library(data.table)
library(tidyverse)
library(corrplot)

d_high <- fread("high_ipd_enrichment.c_elegans.csv")
d_low <- fread("low_ipd_enrichment.c_elegans.csv")
d_mod <- fread("modification_enrichment.c_elegans.csv")

# Tables of enrichment by sample
d_high_table <- spread(d_high[region != "ALL" & base != "ALL", .(sample, region, base, enrichment)], sample, enrichment)
d_low_table <- spread(d_low[region != "ALL" & base != "ALL", .(sample, region, base, enrichment)], sample, enrichment)
d_mod_table <- spread(d_mod[region != "ALL" & base != "ALL", .(sample, region, base, enrichment)], sample, enrichment)

# Rename columns
old_cols <- c("ab", "cd", "k", "l", "PD2182sequel")
new_cols <- c("replicate 1\n/WGA\n/C. elegans", "replicate 2\n/WGA\n/C. elegans", "replicate 1\n/native\n/C. elegans", "replicate 2\n/native\n/C. elegans", "PD2182\n/native\n/C. elegans")
d_high_table <- d_high_table[, ..old_cols]
d_low_table <- d_low_table[, ..old_cols]
setnames(d_high_table, old = old_cols, new = new_cols)
setnames(d_low_table, old = old_cols, new = new_cols)

old_cols <- c("k_normBy_ab", "l_normBy_cd")
new_cols <- c("replicate 1\n/C. elegans", "replicate 2\n/C. elegans")
d_mod_table <- d_mod_table[, ..old_cols]
setnames(d_mod_table, old = old_cols, new = new_cols)

corrplot.na <- function (input_table) {
    if (nrow(input_table) == 0) {
        plot.new()
        title(sprintf("No data in %s", paste0(colnames(input_table), collapse = ", ")), line = -12)
        return()
    }
    cex1 <- 0.25 + 4 / ncol(input_table)
    cors <- cor(input_table)
    if (any(is.na(input_table))) {
        corrplot(cors, method = "color", type = "upper", tl.col = "black", tl.cex = cex1, na.label = "NA")
        corrplot.mixed(cors, lower = "number", upper = "color",
            tl.col = "black", tl.cex = cex1, lower.col = "black", na.label = "NA")
    } else {
        pvalues <- cor.mtest(input_table)$p
        corrplot(cors, p.mat = pvalues, method = "color", type = "upper",
                 sig.level = c(.001, .01, .05), pch.cex = cex1 * 1.5,
                 insig = "label_sig", pch.col = "white", tl.col = "black", tl.cex = cex1, na.label = "NA")
        pvalues[lower.tri(pvalues)] <- NA
        corrplot.mixed(cors, p.mat = pvalues, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", lower = "number", upper = "color",
            tl.col = "black", tl.cex = cex1, pch.cex = cex1 * 1.5, number.cex = cex1 * 1.5, pch.col = "white", lower.col = "black", na.label = "NA")
    }
}

width1 <- 8
height1 <- 8
pdf("enrichment_cor.c_elegans.high.pdf", onefile = TRUE, width = width1, height = height1)
corrplot.na(d_high_table)
invisible(dev.off())
pdf("enrichment_cor.c_elegans.low.pdf", onefile = TRUE, width = width1, height = height1)
corrplot.na(d_low_table)
invisible(dev.off())
pdf("enrichment_cor.c_elegans.modification.pdf", onefile = TRUE, width = width1, height = height1)
corrplot.na(d_mod_table)
invisible(dev.off())
