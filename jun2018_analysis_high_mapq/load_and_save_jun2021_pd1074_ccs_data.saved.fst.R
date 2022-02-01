library(data.table)
library(Biostrings)
library(fst)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

cat("Start loading fst data\n")
jun2021_pd1074_ccs_data <- read_fst("load_and_save_jun2021_pd1074_ccs_data.jun2021_pd1074_ccs.fst", as.data.table = TRUE)
cat("Loaded jun2021_pd1074_ccs_data from load_and_save_jun2021_pd1074_ccs_data.jun2021_pd1074_ccs.fst\n")
