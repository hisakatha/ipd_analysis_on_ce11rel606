library(data.table)
library(Biostrings)
library(fst)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

#prev_wd <- getwd()
#setwd("/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_analysis")
cat("Start loading load_jun2018_data.*.fst\n")
ab_data <- read_fst("load_jun2018_data.ab_data.fst", as.data.table = TRUE)
cd_data <- read_fst("load_jun2018_data.cd_data.fst", as.data.table = TRUE)
#abcd_data <- read_fst("load_jun2018_data.abcd_data.fst", as.data.table = TRUE)
k_normBy_ab_data <- read_fst("load_jun2018_data.k_normBy_ab_data.fst", as.data.table = TRUE)
l_normBy_cd_data <- read_fst("load_jun2018_data.l_normBy_cd_data.fst", as.data.table = TRUE)
#kl_normBy_abcd_data <- read_fst("load_jun2018_data.kl_normBy_abcd_data.fst", as.data.table = TRUE)
cat("Loaded ab_data, cd_data, k_normBy_ab_data, and l_normBy_cd_data from load_jun2018_data.*.fst\n")
#setwd(prev_wd)
