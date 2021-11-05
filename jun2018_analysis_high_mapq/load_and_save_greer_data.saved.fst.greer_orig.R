library(data.table)
library(Biostrings)
library(fst)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

cat("Start loading fst data\n")
greer_orig_data <- read_fst("load_and_save_greer_data.greer_orig.fst", as.data.table = TRUE)
cat("Loaded greer_orig_data from load_and_save_greer_data.greer_orig.fst\n")
