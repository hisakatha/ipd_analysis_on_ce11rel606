library(data.table)
library(Biostrings)
library(hdf5r)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

parse_kinetics_csv <- function(file) {
    # TODO: This may not work for a large genome
    data <- fread(file, header = TRUE)
    return(data)
}

data_path <- "../jun2021_pd1074_ccs/ipd_summary.csv"
jun2021_pd1074_ccs_data <- parse_kinetics_csv(data_path)
cat("Loaded data jun2021_pd1074_ccs\n")

#save(jun2021_pd1074_ccs_data, file = "load_and_save_jun2021_pd1074_ccs_data.Rdata")

library(fst)

write_fst(jun2021_pd1074_ccs_data, "load_and_save_jun2021_pd1074_ccs_data.jun2021_pd1074_ccs.fst", compress = 100)

