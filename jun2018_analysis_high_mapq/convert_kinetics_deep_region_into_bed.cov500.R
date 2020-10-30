library(data.table)
kinetics_deep_region_to_bed <- function (data, coverage_thres, out) {
    fwrite(data[coverage >= coverage_thres, .(refName, tpl - 1, tpl, "deep_kinetics", coverage, ifelse(strand == 0, "+", "-"))] , sep = "\t", file = out, col.names = FALSE)
}

coverage_thres <- 500

source("load_jun2018_data.saved.fst.ab_cd_k_l.R")
kinetics_deep_region_to_bed(ab_data, coverage_thres, "deep_kinetics_region_cov500.ab.bed")
kinetics_deep_region_to_bed(cd_data, coverage_thres, "deep_kinetics_region_cov500.cd.bed")
kinetics_deep_region_to_bed(k_data, coverage_thres, "deep_kinetics_region_cov500.k.bed")
kinetics_deep_region_to_bed(l_data, coverage_thres, "deep_kinetics_region_cov500.l.bed")
source("load_jun2018_data.saved.fst.abcd.R")
kinetics_deep_region_to_bed(abcd_data, coverage_thres, "deep_kinetics_region_cov500.abcd.bed")
source("load_jun2018_data.saved.fst.kl.R")
kinetics_deep_region_to_bed(kl_data, coverage_thres, "deep_kinetics_region_cov500.kl.bed")

source("load_and_save_PD2182_data.saved.fst.PD2182.R")
kinetics_deep_region_to_bed(PD2182_data, coverage_thres, "deep_kinetics_region_cov500.PD2182.bed")

source("load_and_save_PD2182sequel_data.saved.fst.PD2182sequel.R")
kinetics_deep_region_to_bed(PD2182sequel_data, coverage_thres, "deep_kinetics_region_cov500.PD2182sequel.bed")
