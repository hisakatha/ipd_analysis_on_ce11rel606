library(data.table)

source("load_jun2018_data.saved.fst.R")
source("load_and_save_PD2182_data.saved.fst.PD2182.R")
source("load_and_save_PD2182sequel_data.saved.fst.PD2182sequel.R")

get_stats_per_base <- function(data, sample_name, value_type_name){
    # data: a data.table with columns: value, base
    ret1 <- data[, .(sample = sample_name, type = value_type_name, mean = mean(value), sd = sd(value), .N), keyby = base]
    ret2 <- data[, .(sample = sample_name, type = value_type_name, mean = mean(value), sd = sd(value), .N, base = "all")]
    rbind(ret1, ret2)
}

coverage_thres <- 25
ab_data_deep <- ab_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
cd_data_deep <- cd_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
k_data_deep <- k_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
l_data_deep <- l_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
abcd_data_deep <- abcd_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
kl_data_deep <- kl_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
PD2182_data_deep <- PD2182_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
PD2182sequel_data_deep <- PD2182sequel_data[coverage >= coverage_thres][, "value" := .(log2(tMean))]
k_normBy_ab_data_deep <- k_normBy_ab_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres][, "value" := .(log2(ipdRatio))]
l_normBy_cd_data_deep <- l_normBy_cd_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres][, "value" := .(log2(ipdRatio))]
kl_normBy_abcd_data_deep <- kl_normBy_abcd_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres][, "value" := .(log2(ipdRatio))]

ab_stats <- get_stats_per_base(ab_data_deep, "ab_deep", "log2_IPD")
cd_stats <- get_stats_per_base(cd_data_deep, "cd_deep", "log2_IPD")
k_stats <- get_stats_per_base(k_data_deep, "k_deep", "log2_IPD")
l_stats <- get_stats_per_base(l_data_deep, "l_deep", "log2_IPD")
abcd_stats <- get_stats_per_base(abcd_data_deep, "abcd_deep", "log2_IPD")
kl_stats <- get_stats_per_base(kl_data_deep, "kl_deep", "log2_IPD")
PD2182_stats <- get_stats_per_base(PD2182_data_deep, "PD2182_deep", "log2_IPD")
PD2182sequel_stats <- get_stats_per_base(PD2182sequel_data_deep, "PD2182sequel_deep", "log2_IPD")
k_normBy_ab_stats <- get_stats_per_base(k_normBy_ab_data_deep, "k_normBy_ab_deep", "log2_IPD_ratio")
l_normBy_cd_stats <- get_stats_per_base(l_normBy_cd_data_deep, "l_normBy_cd_deep", "log2_IPD_ratio")
kl_normBy_abcd_stats <- get_stats_per_base(kl_normBy_abcd_data_deep, "kl_normBy_abcd_deep", "log2_IPD_ratio")
ab_stats2 <- get_stats_per_base(ab_data_deep[, "value" := .(log2(ipdRatio))], "ab_deep", "log2_IPD_ratio")
cd_stats2 <- get_stats_per_base(cd_data_deep[, "value" := .(log2(ipdRatio))], "cd_deep", "log2_IPD_ratio")
PD2182_stats2 <- get_stats_per_base(PD2182_data_deep[, "value" := .(log2(ipdRatio))], "PD2182_deep", "log2_IPD_ratio")
PD2182sequel_stats2 <- get_stats_per_base(PD2182sequel_data_deep[, "value" := .(log2(ipdRatio))], "PD2182sequel_deep", "log2_IPD_ratio")

all_stats <- rbindlist(list(ab_stats, cd_stats, k_stats, l_stats, abcd_stats, kl_stats, PD2182_stats, PD2182sequel_stats, k_normBy_ab_stats, l_normBy_cd_stats, kl_normBy_abcd_stats, ab_stats2, cd_stats2, PD2182_stats2, PD2182sequel_stats2))
fwrite(all_stats, file = "kinetics_stats.csv")
