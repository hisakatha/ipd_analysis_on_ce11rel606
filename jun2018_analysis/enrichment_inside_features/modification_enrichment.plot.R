library(data.table)
library(ggplot2)
source("plot.functions.R")

d <- fread("modification_enrichment.csv")
d <- d[region != "ALL"]
#d$region <- factor(d$region, levels = c("ALL", d[region != "ALL", unique(region)]))

rename_region(d)

pdf("modification_enrichment.plot.pdf", height = 5, width = 12, onefile = TRUE)

plot.enrichment(d, "k_normBy_ab", "Fold change of modification fraction")
#plot.intron_details(d, "k_normBy_ab")
plot.enrichment(d, "l_normBy_cd", "Fold change of modification fraction")
#plot.intron_details(d, "l_normBy_cd")
plot.enrichment(d, "kl_normBy_abcd", "Fold change of modification fraction")
#plot.intron_details(d, "kl_normBy_abcd")
#plot.enrichment(d, "intersection of callable k_normBy_ab and l_normBy_cd")

invisible(dev.off())

