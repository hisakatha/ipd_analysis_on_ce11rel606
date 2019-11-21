library(data.table)
library(Biostrings)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

cat("Start loading load_jun2018_data.RData\n")
load("/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_analysis/load_jun2018_data.RData")
cat("Loaded ab_data, cd_data, k_data, l_data, abcd_data, kl_data, k_normBy_ab_data, l_normBy_cd_data, and kl_normBy_abcd_data from load_jun2018_data.RData\n")

