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

#TODO: add counts per base
ab_data1 <- ab_data[coverage >= thres, .(refName, log2tMean = log2(tMean), base)]
cd_data1 <- cd_data[coverage >= thres, .(refName, log2tMean = log2(tMean), base)]
k_data1 <- k_data[coverage >= thres, .(refName, log2tMean = log2(tMean), base)]
l_data1 <- l_data[coverage >= thres, .(refName, log2tMean = log2(tMean), base)]
PD2182_data1 <- PD2182_data[coverage >= thres, .(refName, log2tMean = log2(tMean), base)]
PD2182sequel_data1 <- PD2182sequel_data[coverage >= thres, .(refName, log2tMean = log2(tMean), base)]

cut_step <- 0.2
cut_breaks <- seq(-4.1, 6.1, by = cut_step)
ab_data1[, bin := cut(log2tMean, cut_breaks, labels = FALSE)]
cd_data1[, bin := cut(log2tMean, cut_breaks, labels = FALSE)]
k_data1[, bin := cut(log2tMean, cut_breaks, labels = FALSE)]
l_data1[, bin := cut(log2tMean, cut_breaks, labels = FALSE)]
PD2182_data1[, bin := cut(log2tMean, cut_breaks, labels = FALSE)]
PD2182sequel_data1[, bin := cut(log2tMean, cut_breaks, labels = FALSE)]

bases <- c("A", "C", "G", "T")
base_table <- data.table(lower_exclusive = rep(cut_breaks[-length(cut_breaks)], each = length(bases)), upper_inclusive = rep(cut_breaks[-1], each = length(bases)), base = rep(bases, length(cut_breaks) - 1))

ab_count_celegans <- merge(ab_data1[refName != ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
ab_count_celegans[is.na(count), count := 0]
ab_count_celegans <- rbindlist(list(ab_count_celegans, ab_count_celegans[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
cd_count_celegans <- merge(cd_data1[refName != ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
cd_count_celegans[is.na(count), count := 0]
cd_count_celegans <- rbindlist(list(cd_count_celegans, cd_count_celegans[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
k_count_celegans <- merge(k_data1[refName != ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
k_count_celegans[is.na(count), count := 0]
k_count_celegans <- rbindlist(list(k_count_celegans, k_count_celegans[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
l_count_celegans <- merge(l_data1[refName != ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
l_count_celegans[is.na(count), count := 0]
l_count_celegans <- rbindlist(list(l_count_celegans, l_count_celegans[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
PD2182_count_celegans <- merge(PD2182_data1[refName != ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
PD2182_count_celegans[is.na(count), count := 0]
PD2182_count_celegans <- rbindlist(list(PD2182_count_celegans, PD2182_count_celegans[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
PD2182sequel_count_celegans <- merge(PD2182sequel_data1[refName != ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
PD2182sequel_count_celegans[is.na(count), count := 0]
PD2182sequel_count_celegans <- rbindlist(list(PD2182sequel_count_celegans, PD2182sequel_count_celegans[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))

ab_count_ecoli <- merge(ab_data1[refName == ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
ab_count_ecoli[is.na(count), count := 0]
ab_count_ecoli <- rbindlist(list(ab_count_ecoli, ab_count_ecoli[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
cd_count_ecoli <- merge(cd_data1[refName == ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
cd_count_ecoli[is.na(count), count := 0]
cd_count_ecoli <- rbindlist(list(cd_count_ecoli, cd_count_ecoli[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
k_count_ecoli <- merge(k_data1[refName == ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
k_count_ecoli[is.na(count), count := 0]
k_count_ecoli <- rbindlist(list(k_count_ecoli, k_count_ecoli[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
l_count_ecoli <- merge(l_data1[refName == ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
l_count_ecoli[is.na(count), count := 0]
l_count_ecoli <- rbindlist(list(l_count_ecoli, l_count_ecoli[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
PD2182_count_ecoli <- merge(PD2182_data1[refName == ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
PD2182_count_ecoli[is.na(count), count := 0]
PD2182_count_ecoli <- rbindlist(list(PD2182_count_ecoli, PD2182_count_ecoli[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))
PD2182sequel_count_ecoli <- merge(PD2182sequel_data1[refName == ecoli_chr, .(count = .N), keyby = .(bin, base)][, .(lower_exclusive = cut_breaks[bin], upper_inclusive = cut_breaks[bin + 1], base, count)], base_table, all = TRUE)
PD2182sequel_count_ecoli[is.na(count), count := 0]
PD2182sequel_count_ecoli <- rbindlist(list(PD2182sequel_count_ecoli, PD2182sequel_count_ecoli[, .(base = "N", count = sum(count)), keyby = .(lower_exclusive, upper_inclusive)]))

fwrite(ab_count_celegans, file = sprintf("log2ipd_count_coverage%d.ab_celegans.csv", thres))
fwrite(cd_count_celegans, file = sprintf("log2ipd_count_coverage%d.cd_celegans.csv", thres))
fwrite(k_count_celegans, file = sprintf("log2ipd_count_coverage%d.k_celegans.csv", thres))
fwrite(l_count_celegans, file = sprintf("log2ipd_count_coverage%d.l_celegans.csv", thres))
fwrite(PD2182_count_celegans, file = sprintf("log2ipd_count_coverage%d.PD2182_celegans.csv", thres))
fwrite(PD2182sequel_count_celegans, file = sprintf("log2ipd_count_coverage%d.PD2182sequel_celegans.csv", thres))

fwrite(ab_count_ecoli, file = sprintf("log2ipd_count_coverage%d.ab_ecoli.csv", thres))
fwrite(cd_count_ecoli, file = sprintf("log2ipd_count_coverage%d.cd_ecoli.csv", thres))
fwrite(k_count_ecoli, file = sprintf("log2ipd_count_coverage%d.k_ecoli.csv", thres))
fwrite(l_count_ecoli, file = sprintf("log2ipd_count_coverage%d.l_ecoli.csv", thres))
fwrite(PD2182_count_ecoli, file = sprintf("log2ipd_count_coverage%d.PD2182_ecoli.csv", thres))
fwrite(PD2182sequel_count_ecoli, file = sprintf("log2ipd_count_coverage%d.PD2182sequel_ecoli.csv", thres))
