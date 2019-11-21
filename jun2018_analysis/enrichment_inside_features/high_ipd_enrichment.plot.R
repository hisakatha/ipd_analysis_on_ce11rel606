library(data.table)
library(ggplot2)
source("plot.functions.R")

d <- fread("high_ipd_enrichment.csv")
d <- d[region != "ALL"]
#d <- d[base != "ALL"]
#d$region <- factor(d$region, levels = c("ALL", d[region != "ALL", unique(region)]))

rename_region(d)

pdf("high_ipd_enrichment.plot.pdf", height = 5, width = 12, onefile = TRUE)

ylab_text <- "Fold change of fraction of high IPD bases"
plot.enrichment(d, "ab", ylab_text)
#plot.intron_details(d, "ab")
plot.enrichment(d, "cd", ylab_text)
#plot.intron_details(d, "cd")
plot.enrichment(d, "abcd", ylab_text)
#plot.intron_details(d, "abcd")
#plot.enrichment(d, "intersection of callable k_normBy_ab and l_normBy_cd")
plot.enrichment(d, "PD2182", ylab_text)
plot.enrichment(d, "PD2182sequel", ylab_text)

plot.enrichment.dual_facet(d, "ab", ylab_text)
plot.enrichment.dual_facet(d, "cd", ylab_text)
plot.enrichment.dual_facet(d, "abcd", ylab_text)
plot.enrichment.dual_facet(d, "PD2182", ylab_text)
plot.enrichment.dual_facet(d, "PD2182sequel", ylab_text)

invisible(dev.off())
pdf("high_ipd_enrichment.plot.all_samples.pdf", height = 8, width = 12, onefile = TRUE)

d <- d[sample != "abcd"]
plot.enrichment.dual_facet.all_sample(d, ylab_text)

invisible(dev.off())

