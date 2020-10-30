library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
input <- args[1]
output <- args[2]
d_in <- fread(input)
d_out <- d_in[base == "ALL", .(region, num = modifiedNum, fraction = modifiedNum / modifiedNum[which(region == "ALL")]), by = sample]
fwrite(d_out, file = output)
