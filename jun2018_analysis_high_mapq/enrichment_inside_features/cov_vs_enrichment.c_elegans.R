library(data.table)
library(ggplot2)
d25 <- fread("cov100_enrichment.c_elegans.csv")[base == "ALL", .(region, num = callableNum, fraction = callableNum / callableNum[which(region == "ALL")]), by = sample]
d100 <- fread("cov100_enrichment.c_elegans.simple.csv")
d200 <- fread("cov200_enrichment.c_elegans.simple.csv")
d300 <- fread("cov300_enrichment.c_elegans.simple.csv")
d400 <- fread("cov400_enrichment.c_elegans.simple.csv")
d500 <- fread("cov500_enrichment.c_elegans.simple.csv")
d600 <- fread("cov600_enrichment.c_elegans.simple.csv")

d25[, "cov_thres" := .(25)]
d100[, "cov_thres" := .(100)]
d200[, "cov_thres" := .(200)]
d300[, "cov_thres" := .(300)]
d400[, "cov_thres" := .(400)]
d500[, "cov_thres" := .(500)]
d600[, "cov_thres" := .(600)]

d <- rbindlist(list(d25, d100, d200, d300, d400, d500, d600))
d <- d[sample != "PD2182"]
d[, sample := ifelse(sample == "ab", "VC2010+OP50/WGA", ifelse(sample == "cd", "VC2010/WGA", ifelse(sample == "k", "VC2010+OP50/native",
            ifelse(sample == "l", "VC2010/native", ifelse(sample == "PD2182sequel", "PD2182/native", ifelse(sample == "PD2182", "PD2182(RS II)/native", "UNKNOWN"))))))]
d$sample <- factor(d$sample, levels = c("VC2010+OP50/WGA", "VC2010/WGA", "VC2010+OP50/native", "VC2010/native", "PD2182/native", "PD2182(RS II)/native", "UNKNOWN"))

num_region <- d[, length(unique(region))]
g_all <- ggplot(d[region != "ALL"], aes(cov_thres, fraction)) +
    geom_point(aes(color = region, shape = region)) + geom_line(aes(color = region)) +
    facet_grid(sample ~ .) +
    xlab("Minimum threshold of coverage depth") +
    ylab("Fraction of bases in each region") +
    scale_shape_manual(values = 0:(num_region - 1))

g_sub <- ggplot(d[region == "exon" | region == "tandem repeat"], aes(cov_thres, fraction)) +
    geom_point(aes(color = region, shape = region)) + geom_line(aes(color = region)) +
    facet_grid(sample ~ .) +
    xlab("Minimum threshold of coverage depth") +
    ylab("Fraction of bases in each region")

pdf("cov_vs_enrichment.c_elegans.pdf")
print(g_all)
print(g_sub)
invisible(dev.off())
