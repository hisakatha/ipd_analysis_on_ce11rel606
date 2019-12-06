library(data.table)
library(tidyverse)
library(corrplot)

d_high <- fread("high_ipd_enrichment.csv")
d_low <- fread("low_ipd_enrichment.csv")
d_mod <- fread("modification_enrichment.csv")

# Tables of enrichment by sample
d_high_table <- spread(d_high[region != "ALL" & base != "ALL", .(sample, region, base, enrichment)], sample, enrichment)
d_low_table <- spread(d_low[region != "ALL" & base != "ALL", .(sample, region, base, enrichment)], sample, enrichment)
d_mod_table <- spread(d_mod[region != "ALL" & base != "ALL", .(sample, region, base, enrichment)], sample, enrichment)

# Rename columns
old_cols <- c("ab", "cd", "PD2182", "PD2182sequel")
new_cols <- c("VC2010+OP50\n(WGA)", "VC2010\n(WGA)", "PD2182\n(PacBio RS II)", "PD2182\n(PacBio Sequel)")
d_high_table <- d_high_table[, ..old_cols]
d_low_table <- d_low_table[, ..old_cols]
setnames(d_high_table, old = old_cols, new = new_cols)
setnames(d_low_table, old = old_cols, new = new_cols)

old_cols <- c("k_normBy_ab", "l_normBy_cd")
new_cols <- c("VC2010+OP50", "VC2010")
d_mod_table <- d_mod_table[, ..old_cols]
setnames(d_mod_table, old = old_cols, new = new_cols)

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

pdf("enrichment_cor.high.pdf", onefile = TRUE, width = 6, height = 6)
corrplot.na(d_high_table)
invisible(dev.off())
pdf("enrichment_cor.low.pdf", onefile = TRUE, width = 6, height = 6)
corrplot.na(d_low_table)
invisible(dev.off())
pdf("enrichment_cor.modification.pdf", onefile = TRUE, width = 6, height = 6)
corrplot.na(d_mod_table)
invisible(dev.off())
