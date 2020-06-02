library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
data.path <- args[1]
occ.path <- args[2]
out.path <- args[3]
chr.ecoli <- "E._coli_REL606"

data <- fread(data.path)
if (ncol(data) == 0) {
    data <- data.table(position=integer(0), strand=character(0), value=numeric(0), label=character(0), src=integer(0))
}

occ <- fread(occ.path, sep = " ", header = FALSE, col.names = c("chr", "start", "strand"))
if (ncol(occ) == 0) {
    occ <- data.table(chr=character(0), start=integer(0), strand=character(0))
}

src.celegans <- occ[chr != chr.ecoli, .I]
data.celegans <- data[src %in% src.celegans]
fwrite(data.celegans, file = out.path)
