library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
file_path <- args[1]
d <- fread(file_path)
d[, "mean_ratio" := .(mean / genome_mean)]
fwrite(d)
