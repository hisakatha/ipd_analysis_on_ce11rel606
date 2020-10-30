library(data.table)
library(ggplot2)
source("plot.functions.R")

d <- fread("cov500_enrichment.c_elegans.csv")
d <- d[region != "ALL"]
#d <- d[base != "ALL"]
#d$region <- factor(d$region, levels = c("ALL", d[region != "ALL", unique(region)]))
set_region_label(d)

#pdf("high_ipd_enrichment.c_elegans.plot.pdf", height = 5, width = 12, onefile = TRUE)
ylab_text <- "Fold change of fraction of bases with coverage 500 or more in each region"
#plot.enrichment(d, "ab", ylab_text)
##plot.intron_details(d, "ab")
#plot.enrichment(d, "cd", ylab_text)
##plot.intron_details(d, "cd")
#plot.enrichment(d, "abcd", ylab_text)
##plot.intron_details(d, "abcd")
##plot.enrichment(d, "intersection of callable k_normBy_ab and l_normBy_cd")
#plot.enrichment(d, "PD2182", ylab_text)
#plot.enrichment(d, "PD2182sequel", ylab_text)
#
#plot.enrichment.dual_facet(d, "ab", ylab_text)
#plot.enrichment.dual_facet(d, "cd", ylab_text)
#plot.enrichment.dual_facet(d, "abcd", ylab_text)
#plot.enrichment.dual_facet(d, "PD2182", ylab_text)
#plot.enrichment.dual_facet(d, "PD2182sequel", ylab_text)
#invisible(dev.off())

d <- d[sample != "abcd"]
pdf("cov500_enrichment.c_elegans.plot.all_samples.pdf", height = 8, width = 15, onefile = TRUE)
plot.enrichment.dual_facet.all_sample.v2(d, ylab_text)
invisible(dev.off())

pdf("cov500_enrichment.c_elegans.plot.all_samples.subset1.pdf", height = 6, width = 6, onefile = TRUE)
plot.enrichment.dual_facet.all_sample.subset1.v2(d, ylab_text)
plot.enrichment.dual_facet.all_sample.subset2.v2(d, ylab_text)
invisible(dev.off())
