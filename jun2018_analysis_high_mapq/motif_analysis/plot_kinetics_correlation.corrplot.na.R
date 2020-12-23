library(data.table)
library(corrplot)

corrplot.na <- function (input_table) {
    if (nrow(input_table) == 0) {
        plot.new()
        title(sprintf("No data\nin %s", paste0(colnames(input_table), collapse = ", ")), line = -12)
        return()
    }
    cex1 <- 0.25 + 4 / ncol(input_table)
    cors <- cor(input_table)
    if (all(is.na(cors[upper.tri(cors)]))) {
        plot.new()
        title(sprintf("Not enough data for correlation matrix\nin %s", paste0(colnames(input_table), collapse = ", ")), line = -12)
        return()
    }
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
