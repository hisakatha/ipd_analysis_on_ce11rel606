library(data.table)

#score_thres <- 20
#ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
coverage_thres <- 25
collect_frac <- 0.01
bases <- c("A", "C", "G", "T")

#k_normBy_ab_data[, deviated_flags := list(score >= score_thres)]
#k_normBy_ab_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]

#l_normBy_cd_data[, deviated_flags := list(score >= score_thres)]
#l_normBy_cd_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]

#kl_normBy_abcd_data[, deviated_flags := list(score >= score_thres)]
#kl_normBy_abcd_data[, modified_flags := list(deviated_flags == TRUE & native_coverage >= coverage_thres & control_coverage >= coverage_thres & ipdRatio >= ipdRatio_thres)]


#union_data <- cbind(a = k_normBy_ab_data, b = l_normBy_cd_data)

#deviated_flags_intersection <- k_normBy_ab_data$deviated_flags & l_normBy_cd_data$deviated_flags

#modified_flags_intersection <- k_normBy_ab_data$modified_flags & l_normBy_cd_data$modified_flags

library(ggplot2)
library(ggExtra)
library(cowplot)

export_as_gff3 <- function(data, file) {
    cat("##gff-version 3\n", file = file, append = FALSE)
    formatted_data <- data.table(
        seqid = data$refName,
        source = "call_extreme_ipd.R",
        type = "base",
        start = data$tpl,
        end = data$tpl,
        score = data$score,
        strand = ifelse(data$strand == 0, "+", "-"),
        phase = ".",
        attributes = sprintf("tMean=%.5f;tErr=%.5f;modelPrediction=%.5f;IPDRatio=%.5f;coverage=%d", data$tMean, data$tErr, data$modelPrediction, data$ipdRatio, data$coverage))
    fwrite(formatted_data, file, append = TRUE, sep = "\t", col.names = FALSE)
}

export_intersection_as_gff3 <- function(data_list, name_list, file) {
    data1 <- data_list[[1]]
    cat("##gff-version 3\n", file = file, append = FALSE)
    formatted_data <- data.table(
        seqid = data1$refName,
        source = "call_extreme_ipd.R",
        type = "base",
        start = data1$tpl,
        end = data1$tpl,
        score = rowMeans(sapply(data_list, function(data){ data$score })),
        strand = ifelse(data1$strand == 0, "+", "-"),
        phase = ".",
        attributes = apply(mapply(function(data, name){sprintf("%s_tMean=%.5f;%s_tErr=%.5f;%s_modelPrediction=%.5f;%s_IPDRatio=%.5f;%s_coverage=%d", name, data$tMean, name, data$tErr, name, data$modelPrediction, name, data$ipdRatio, name, data$coverage)}, data_list, name_list), 1, paste0, collapse = ";"))
    fwrite(formatted_data, file, append = TRUE, sep = "\t", col.names = FALSE)
}

save_extreme_ipd <- function(data, name){
    data_high <- data[min(ceiling((1 - collect_frac) * .N) + 1, .N):.N]
    export_as_gff3(data_high, sprintf("extreme_ipd.high.%s.gff", name))
    data_low <- data[1:max(floor(collect_frac * .N), 1)]
    export_as_gff3(data_low, sprintf("extreme_ipd.low.%s.gff", name))
    data <- data[coverage >= coverage_thres]
    data_high <- data[min(ceiling((1 - collect_frac) * .N) + 1, .N):.N]
    export_as_gff3(data_high, sprintf("extreme_ipd.high.%s.coverage%g.gff", name, coverage_thres))
    data_low <- data[1:max(floor(collect_frac * .N), 1)]
    export_as_gff3(data_low, sprintf("extreme_ipd.low.%s.coverage%g.gff", name, coverage_thres))
}

save_intersection <- function(data1, data_orig1, name1, thres1, data2, data_orig2, name2, thres2){
    both_high <- (data_orig1[,tMean] >= thres1[1]) & (data_orig2[,tMean] >= thres2[1])
    both_low <- (data_orig1[,tMean] > 0) & (data_orig2[,tMean] > 0) & (data_orig1[,tMean] <= thres1[2]) & (data_orig2[,tMean] <= thres2[2])
    export_intersection_as_gff3(list(data_orig1[both_high], data_orig2[both_high]), list(name1, name2), sprintf("extreme_ipd.high.intersection_%s_%s.gff", name1, name2))
    export_intersection_as_gff3(list(data_orig1[both_low], data_orig2[both_low]), list(name1, name2), sprintf("extreme_ipd.low.intersection_%s_%s.gff", name1, name2))

    both_high_cov <- (data_orig1[,coverage] >= coverage_thres) & (data_orig2[,coverage] >= coverage_thres)
    data_orig1 <- data_orig1[both_high_cov]
    data_orig2 <- data_orig2[both_high_cov]
    both_high_with_cov <- (data_orig1[,tMean] >= thres1[3]) & (data_orig2[,tMean] >= thres2[3])
    both_low_with_cov <- (data_orig1[,tMean] <= thres1[4]) & (data_orig2[,tMean] <= thres2[4])
    export_intersection_as_gff3(list(data_orig1[both_high_with_cov], data_orig2[both_high_with_cov]), list(name1, name2), sprintf("extreme_ipd.high.intersection_%s_%s.coverage%g.gff", name1, name2, coverage_thres))
    export_intersection_as_gff3(list(data_orig1[both_low_with_cov], data_orig2[both_low_with_cov]), list(name1, name2), sprintf("extreme_ipd.low.intersection_%s_%s.coverage%g.gff", name1, name2, coverage_thres))
}

#source("load_jun2018_data.saved.fst.wga.R")
source("load_jun2018_data.saved.fst.ab_cd_k_l.R")

ab_data_orig <- copy(ab_data)
cd_data_orig <- copy(cd_data)
k_data_orig <- copy(k_data)
l_data_orig <- copy(l_data)

ab_data <- ab_data[tMean > 0]
cd_data <- cd_data[tMean > 0]
k_data <- k_data[tMean >0]
l_data <- l_data[tMean >0]

setkey(ab_data, tMean)
setkey(cd_data, tMean)
setkey(k_data, tMean)
setkey(l_data, tMean)

save_extreme_ipd(ab_data, "ab")
save_extreme_ipd(cd_data, "cd")
save_extreme_ipd(k_data, "k")
save_extreme_ipd(l_data, "l")

get_thres <- function(data){
    high_thres <- data[min(ceiling((1 - collect_frac) * .N) + 1, .N), tMean]
    low_thres <- data[max(floor(collect_frac * .N), 1), tMean]
    data <- data[coverage >= coverage_thres]
    high_thres_with_cov <- data[min(ceiling((1 - collect_frac) * .N) + 1, .N), tMean]
    low_thres_with_cov <- data[max(floor(collect_frac * .N), 1), tMean]
    return(c(high_thres, low_thres, high_thres_with_cov, low_thres_with_cov))
}

ab_thres <- get_thres(ab_data)
cd_thres <- get_thres(cd_data)
k_thres <- get_thres(k_data)
l_thres <- get_thres(l_data)

save_intersection(ab_data, ab_data_orig, "ab", ab_thres, cd_data, cd_data_orig, "cd", cd_thres)
save_intersection(ab_data, ab_data_orig, "ab", ab_thres, k_data, k_data_orig, "k", k_thres)
save_intersection(cd_data, cd_data_orig, "cd", cd_thres, l_data, l_data_orig, "l", l_thres)

