library(data.table)
library(Biostrings)
library(fst)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

cat("Start loading load_jun2018_data.k_data.fst\n")
k_data <- read_fst("load_jun2018_data.k_data.fst", as.data.table = TRUE)
cat("Loaded k_data from load_jun2018_data.*.fst\n")

