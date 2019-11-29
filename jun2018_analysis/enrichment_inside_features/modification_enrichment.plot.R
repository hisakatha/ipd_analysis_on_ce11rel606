library(data.table)
library(ggplot2)
source("plot.functions.R")

d <- fread("modification_enrichment.csv")
d <- d[region != "ALL"]
d <- d[base != "ALL"]
#d$region <- factor(d$region, levels = c("ALL", d[region != "ALL", unique(region)]))
set_region_label(d)

pdf("modification_enrichment.plot.pdf", height = 5, width = 12, onefile = TRUE)
ylab_text <- "Fold change of modification fraction in each region"
plot.enrichment(d, "k_normBy_ab", ylab_text)
#plot.intron_details(d, "k_normBy_ab")
plot.enrichment(d, "l_normBy_cd", ylab_text)
#plot.intron_details(d, "l_normBy_cd")
plot.enrichment(d, "kl_normBy_abcd", ylab_text)
#plot.intron_details(d, "kl_normBy_abcd")
#plot.enrichment(d, "intersection of callable k_normBy_ab and l_normBy_cd")

plot.enrichment.dual_facet(d, "k_normBy_ab", ylab_text)
plot.enrichment.dual_facet(d, "l_normBy_cd", ylab_text)
plot.enrichment.dual_facet(d, "kl_normBy_abcd", ylab_text)
invisible(dev.off())

d <- d[sample != "kl_normBy_abcd"]
pdf("modification_enrichment.plot.all_samples.pdf", height = 8, width = 12, onefile = TRUE)
plot.enrichment.dual_facet.all_sample(d, ylab_text)
invisible(dev.off())

pdf("modification_enrichment.plot.all_samples.subset1.pdf", height = 6, width = 5, onefile = TRUE)
plot.enrichment.dual_facet.all_sample.subset3(d, ylab_text)
invisible(dev.off())
