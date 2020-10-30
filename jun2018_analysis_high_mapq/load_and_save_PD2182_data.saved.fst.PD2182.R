library(data.table)
library(Biostrings)
library(fst)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

cat("Start loading fst data\n")
PD2182_data <- read_fst("load_and_save_PD2182_data.PD2182_data.fst", as.data.table = TRUE)
cat("Loaded PD2182_data from load_and_save_PD2182_data.PD2182_data.fst\n")
