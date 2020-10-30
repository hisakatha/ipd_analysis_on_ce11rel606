library(data.table)
source("../load_jun2018_data.saved.fst.ab_cd_k_l.R", chdir = TRUE)
source("../load_and_save_PD2182_data.saved.fst.PD2182.R", chdir = TRUE)
source("../load_and_save_PD2182sequel_data.saved.fst.PD2182sequel.R", chdir = TRUE)

# threshold of valid IPD counts (coverage)
thres <- 25

setkey(ab_data, refName, base)
setkey(cd_data, refName, base)
setkey(k_data, refName, base)
setkey(l_data, refName, base)
setkey(PD2182_data, refName, base)
setkey(PD2182sequel_data, refName, base)

ab_data1 <- ab_data[coverage >= thres, .(refName, log2tMean = log2(tMean))]
cd_data1 <- cd_data[coverage >= thres, .(refName, log2tMean = log2(tMean))]
k_data1 <- k_data[coverage >= thres, .(refName, log2tMean = log2(tMean))]
l_data1 <- l_data[coverage >= thres, .(refName, log2tMean = log2(tMean))]
PD2182_data1 <- PD2182_data[coverage >= thres, .(refName, log2tMean = log2(tMean))]
PD2182sequel_data1 <- PD2182sequel_data[coverage >= thres, .(refName, log2tMean = log2(tMean))]

fwrite(ab_data1[refName != ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.ab_celegans.csv", thres), col.names = FALSE)
fwrite(cd_data1[refName != ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.cd_celegans.csv", thres), col.names = FALSE)
fwrite(k_data1[refName != ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.k_celegans.csv", thres), col.names = FALSE)
fwrite(l_data1[refName != ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.l_celegans.csv", thres), col.names = FALSE)
fwrite(PD2182_data1[refName != ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.PD2182_celegans.csv", thres), col.names = FALSE)
fwrite(PD2182sequel_data1[refName != ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.PD2182sequel_celegans.csv", thres), col.names = FALSE)
fwrite(ab_data1[refName == ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.ab_ecoli.csv", thres), col.names = FALSE)
fwrite(cd_data1[refName == ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.cd_ecoli.csv", thres), col.names = FALSE)
fwrite(k_data1[refName == ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.k_ecoli.csv", thres), col.names = FALSE)
fwrite(l_data1[refName == ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.l_ecoli.csv", thres), col.names = FALSE)
fwrite(PD2182_data1[refName == ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.PD2182_ecoli.csv", thres), col.names = FALSE)
fwrite(PD2182sequel_data1[refName == ecoli_chr, .(log2tMean)], file = sprintf("log2ipd_coverage%d.PD2182sequel_ecoli.csv", thres), col.names = FALSE)
