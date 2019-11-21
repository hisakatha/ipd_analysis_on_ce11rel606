source("load_jun2018_data.saved.R")

library(data.table)
export_as_gff3 <- function(data, file) {
    cat("##gff-version 3\n", file = file, append = FALSE)
    formatted_data <- data.table(
        seqid = data$refName,
        source = "call_modification.R",
        type = "modified_DNA_base",
        start = data$tpl,
        end = data$tpl,
        score = data$score,
        strand = ifelse(data$strand == 0, "+", "-"),
        phase = ".",
        attributes = sprintf("IPDRatio=%.5f;native_coverage=%d;control_coverage=%d", data$ipdRatio, data$native_coverage, data$control_coverage))
    fwrite(formatted_data, file, append = TRUE, sep = "\t", col.names = FALSE)
}

export_intersection_as_gff3 <- function(data_list, name_list, file) {
    data1 <- data_list[[1]]
    cat("##gff-version 3\n", file = file, append = FALSE)
    formatted_data <- data.table(
        seqid = data1$refName,
        source = "call_modification.R",
        type = "modified_DNA_base",
        start = data1$tpl,
        end = data1$tpl,
        score = rowMeans(sapply(data_list, function(data){ data$score })),
        strand = ifelse(data1$strand == 0, "+", "-"),
        phase = ".",
        attributes = apply(mapply(function(data, name){sprintf("%s_IPDRatio=%.5f;%s_native_coverage=%d;%s_control_coverage=%d", name, data$ipdRatio, name, data$native_coverage, name, data$control_coverage)}, data_list, name_list), 1, paste0, collapse = ";"))
    fwrite(formatted_data, file, append = TRUE, sep = "\t", col.names = FALSE)
}

score_thres <- 20
ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
coverage_thres <- 25
bases <- c("A", "C", "G", "T")

k_normBy_ab_data[, deviated_flags := list(score >= score_thres)]
export_as_gff3(k_normBy_ab_data[deviated_flags == TRUE], sprintf("modifications.k_normBy_ab.score%g.gff", score_thres))
k_normBy_ab_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]
export_as_gff3(k_normBy_ab_data[modified_flags == TRUE], sprintf("modifications.k_normBy_ab.score%g_ipdratio%g_coverage%g.gff", score_thres, ipdRatio_thres, coverage_thres))
for (b in bases) { export_as_gff3(k_normBy_ab_data[modified_flags == TRUE & base == b], sprintf("modifications.k_normBy_ab.score%g_ipdratio%g_coverage%g.%s.gff", score_thres, ipdRatio_thres, coverage_thres, b)) }

l_normBy_cd_data[, deviated_flags := list(score >= score_thres)]
export_as_gff3(l_normBy_cd_data[deviated_flags == TRUE], sprintf("modifications.l_normBy_cd.score%g.gff", score_thres))
l_normBy_cd_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]
export_as_gff3(l_normBy_cd_data[modified_flags == TRUE], sprintf("modifications.l_normBy_cd.score%g_ipdratio%g_coverage%g.gff", score_thres, ipdRatio_thres, coverage_thres))
for (b in bases) { export_as_gff3(l_normBy_cd_data[modified_flags == TRUE & base == b], sprintf("modifications.l_normBy_cd.score%g_ipdratio%g_coverage%g.%s.gff", score_thres, ipdRatio_thres, coverage_thres, b)) }

kl_normBy_abcd_data[, deviated_flags := list(score >= score_thres)]
export_as_gff3(kl_normBy_abcd_data[deviated_flags == TRUE], sprintf("modifications.kl_normBy_abcd.score%g.gff", score_thres))
kl_normBy_abcd_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]
export_as_gff3(kl_normBy_abcd_data[modified_flags == TRUE], sprintf("modifications.kl_normBy_abcd.score%g_ipdratio%g_coverage%g.gff", score_thres, ipdRatio_thres, coverage_thres))
for (b in bases) { export_as_gff3(kl_normBy_abcd_data[modified_flags == TRUE & base == b], sprintf("modifications.kl_normBy_abcd.score%g_ipdratio%g_coverage%g.%s.gff", score_thres, ipdRatio_thres, coverage_thres, b)) }


union_data <- cbind(a = k_normBy_ab_data, b = l_normBy_cd_data)

deviated_flags_intersection <- k_normBy_ab_data$deviated_flags & l_normBy_cd_data$deviated_flags
#deviated_intersection_data <- union_data[deviated_flags_intersection]
export_intersection_as_gff3(list(k_normBy_ab_data[deviated_flags_intersection], l_normBy_cd_data[deviated_flags_intersection]), list("k_normBy_ab", "l_normBy_cd"), sprintf("modifications.intersection.score%g.gff", score_thres))

modified_flags_intersection <- k_normBy_ab_data$modified_flags & l_normBy_cd_data$modified_flags
#modified_intersection_data <- union_data[modified_flags_intersection]
export_intersection_as_gff3(list(k_normBy_ab_data[modified_flags_intersection], l_normBy_cd_data[modified_flags_intersection]), list("k_normBy_ab", "l_normBy_cd"), sprintf("modifications.intersection.score%g_ipdratio%g_coverage%g.gff", score_thres, ipdRatio_thres, coverage_thres))
for (b in bases) { export_intersection_as_gff3(list(k_normBy_ab_data[modified_flags_intersection & base == b], l_normBy_cd_data[modified_flags_intersection & base == b]), list("k_normBy_ab", "l_normBy_cd"), sprintf("modifications.intersection.score%g_ipdratio%g_coverage%g.%s.gff", score_thres, ipdRatio_thres, coverage_thres, b)) }

library(ggplot2)

plot_score_distribution <- function(data, name) {
    p <- ggplot(data[score > 0], aes(score, color = ipdRatio > 1)) + geom_freqpoly(binwidth = 1) + facet_wrap(~ base, ncol=2) + ggtitle(sprintf("%s: modification score per base", name))
    print(p)
    print(p + scale_y_log10(limits = c(1,NA)))
    p <- ggplot(data[modified_flags == TRUE], aes(score, color = ipdRatio > 1)) + geom_freqpoly(binwidth = 1) + facet_wrap(~ base, ncol=2) + ggtitle(sprintf("%s (score >= %g, ipdRatio >= %g, coverage >= %g): modification score per base", name, score_thres, ipdRatio_thres, coverage_thres))
    print(p)
    print(p + scale_y_log10(limits = c(1,NA)))
    cat(sprintf("Printed plots for %s\n", name))
}

plot_score_distribution_for_2sets <- function(data, name, a.name, b.name) {
    a.label <- sprintf("%s IPD ratio > 1", a.name)
    x.label <- sprintf("%s score", a.name)
    b.label <- sprintf("%s IPD ratio > 1", b.name)
    y.label <- sprintf("%s score", b.name)
    subdata <- data[a.score > 0 & b.score > 0][, c("a.isLargeIpdRatio", "b.isLargeIpdRatio") := .(a.ipdRatio > 1, b.ipdRatio > 1)]
    subdata_count <- subdata[,.N, by = .(a.isLargeIpdRatio, b.isLargeIpdRatio, a.base)]
    p <- ggplot(subdata, aes(a.score, b.score)) + geom_hex(binwidth = 2) +
        facet_grid(rows = vars(!!a.label := a.isLargeIpdRatio, !!b.label := b.isLargeIpdRatio), cols = vars(base = a.base), labeller = label_both) +
        geom_text(data = subdata_count, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d", N)), hjust = "inward", vjust = "inward")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s: modification score per base", name)) + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    subdata <- data[a.modified_flags == TRUE & b.modified_flags == TRUE][, c("a.isLargeIpdRatio", "b.isLargeIpdRatio") := .(a.ipdRatio > 1, b.ipdRatio > 1)]
    subdata_count <- subdata[,.N, by = .(a.isLargeIpdRatio, b.isLargeIpdRatio, a.base)]
    p <- ggplot(subdata, aes(a.score, b.score)) + geom_hex(binwidth = 2) +
        facet_grid(rows = vars(!!a.label := a.isLargeIpdRatio, !!b.label := b.isLargeIpdRatio), cols = vars(base = a.base), labeller = label_both) +
        geom_text(data = subdata_count, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d", N)), hjust = "inward", vjust = "inward")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s (score >= %g, ipdRatio >= %g, coverage >= %g): modification score per base", name, score_thres, ipdRatio_thres, coverage_thres)) + xlab(x.label) + ylab(y.label)
    print(p)
    cat(sprintf("Printed plots for %s\n", name))
}

plot_ipdratio_distribution_for_2sets <- function(data, name, a.name, b.name) {
    x.label <- sprintf("%s log2 (IPD ratio)", a.name)
    y.label <- sprintf("%s log2 (IPD ratio)", b.name)
    subdata <- data[a.score > 0 & b.score > 0 & a.ipdRatio > 0 & b.ipdRatio > 0]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf), by = a.base]
    for (b in bases) {
        base_data <- subdata[a.base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2(a.ipdRatio)], base_data[,log2(b.ipdRatio)], method = "pearson")
            subdata_stat[a.base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
    }
    p <- ggplot(subdata, aes(log2(a.ipdRatio), log2(b.ipdRatio))) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = a.base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d, cor = %.3g (p = %.3g)", N, cor.value, p.value)), hjust = "inward", vjust = "inward")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s: log_2 (IPD ratio) per base", name), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    subdata <- data[a.modified_flags == TRUE & b.modified_flags == TRUE]
    subdata_stat <- subdata[, .(.N, cor.value = Inf, p.value = Inf), by = a.base]
    for (b in bases) {
        base_data <- subdata[a.base == b]
        if (base_data[,.N] >= 3L) {
            cor.ret <- cor.test(base_data[,log2(a.ipdRatio)], base_data[,log2(b.ipdRatio)], method = "pearson")
            subdata_stat[a.base == b, c("cor.value", "p.value") := .(cor.ret$estimate, cor.ret$p.value)]
        }
    }
    p <- ggplot(subdata, aes(log2(a.ipdRatio), log2(b.ipdRatio))) + geom_hex(binwidth = 0.2) +
        facet_wrap(vars(base = a.base), ncol = 2, labeller = label_both) +
        geom_text(data = subdata_stat, mapping = aes(x = -Inf, y = Inf, label = sprintf("N = %d, cor = %.3g (p = %.3g)", N, cor.value, p.value)), hjust = "inward", vjust = "inward")
    if(subdata[,.N] == 0){ p <- ggplot() + annotate("text", x = 1, y = 1, label = "No data") }
    p <- p + ggtitle(sprintf("%s (score >= %g, ipdRatio >= %g, coverage >= %g): log_2 (IPD ratio) per base", name, score_thres, ipdRatio_thres, coverage_thres), subtitle = "with Pearson's correlation tests") + xlab(x.label) + ylab(y.label)
    print(p)
    print(p + scale_fill_continuous(trans = "log10"))
    cat(sprintf("Printed plots for %s\n", name))
}

pdf(sprintf("call_modification.score%g_ipdratio%g_coverage%g.pdf", score_thres, ipdRatio_thres, coverage_thres), width = 16, height = 16)
plot_score_distribution(k_normBy_ab_data, "k_normBy_ab")
plot_score_distribution(k_normBy_ab_data[refName != ecoli_chr], "k_normBy_ab C. elegans")
plot_score_distribution(k_normBy_ab_data[refName == ecoli_chr], "k_normBy_ab E. coli")
plot_score_distribution(l_normBy_cd_data, "l_normBy_cd")
plot_score_distribution(l_normBy_cd_data[refName != ecoli_chr], "l_normBy_cd C. elegans")
plot_score_distribution(l_normBy_cd_data[refName == ecoli_chr], "l_normBy_cd E. coli")
plot_score_distribution(kl_normBy_abcd_data, "kl_normBy_abcd")
plot_score_distribution(kl_normBy_abcd_data[refName != ecoli_chr], "kl_normBy_abcd C. elegans")
plot_score_distribution(kl_normBy_abcd_data[refName == ecoli_chr], "kl_normBy_abcd E. coli")

plot_score_distribution_for_2sets(union_data, "intersection of k_normBy_ab and l_normBy_cd", "k_normBy_ab", "l_normBy_cd")
plot_score_distribution_for_2sets(union_data[a.refName != ecoli_chr], "intersection of k_normBy_ab and l_normBy_cd in C. elegans", "k_normBy_ab", "l_normBy_cd")
plot_score_distribution_for_2sets(union_data[a.refName == ecoli_chr], "intersection of k_normBy_ab and l_normBy_cd in E. coli", "k_normBy_ab", "l_normBy_cd")
plot_ipdratio_distribution_for_2sets(union_data, "intersection of k_normBy_ab and l_normBy_cd", "k_normBy_ab", "l_normBy_cd")
plot_ipdratio_distribution_for_2sets(union_data[a.refName != ecoli_chr], "intersection of k_normBy_ab and l_normBy_cd in C. elegans", "k_normBy_ab", "l_normBy_cd")
plot_ipdratio_distribution_for_2sets(union_data[a.refName == ecoli_chr], "intersection of k_normBy_ab and l_normBy_cd in E. coli", "k_normBy_ab", "l_normBy_cd")
invisible(dev.off())

