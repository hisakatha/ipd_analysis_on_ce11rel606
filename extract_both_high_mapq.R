library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
input <- args[1]
output <- args[2]
d <- fread(input)
# header: read_name,top_mapq_celegans,count_celegans,top_mapq_ecoli,count_ecoli,abs_diff_mapq
d_out <- d[top_mapq_celegans == 254 & top_mapq_ecoli == 254, .(read_name)]
fwrite(d_out, file = output, col.names = FALSE)
