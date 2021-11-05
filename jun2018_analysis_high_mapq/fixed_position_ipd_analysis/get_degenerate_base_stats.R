library(data.table)

# WIP
# TODO: Deal with the case N == 0
get_degenerate_base_stats <- function(kinetics, sample1, type1){
    target <- kinetics[sample == sample1 & type == type1]
    print("")
    print(sample1)
    print(type1)
    print(target)
    stopifnot(target[base == "A", N] > 0)
    stopifnot(target[base == "C", N] > 0)
    stopifnot(target[base == "G", N] > 0)
    stopifnot(target[base == "T", N] > 0)
    kinetics_new <- rbindlist(list(
        data.table(base = "A", sample = sample1, type = type1, mean = target[base == "A", mean], var = target[base == "A", sd ^ 2], N = target[base == "A", N]),
        data.table(base = "C", sample = sample1, type = type1, mean = target[base == "C", mean], var = target[base == "C", sd ^ 2], N = target[base == "C", N]),
        data.table(base = "G", sample = sample1, type = type1, mean = target[base == "G", mean], var = target[base == "G", sd ^ 2], N = target[base == "G", N]),
        data.table(base = "T", sample = sample1, type = type1, mean = target[base == "T", mean], var = target[base == "T", sd ^ 2], N = target[base == "T", N]),
        data.table(base = "W", sample = sample1, type = type1,
            mean = (target[base == "A", mean * N] + target[base == "T", mean * N]) / (target[base == "A", N] + target[base == "T", N]),
            var = (target[base == "A", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "T", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "A", N] + target[base == "T", N]) *
                   ((target[base == "A", N * mean] + target[base == "T", N * mean]) / (target[base == "A", N] + target[base == "T", N])) ^ 2) /
                   (target[base == "A", N] + target[base == "T", N] - 1),
            N = target[base == "A", N] + target[base == "T", N]),
        data.table(base = "S", sample = sample1, type = type1,
            mean = (target[base == "C", mean * N] + target[base == "G", mean * N]) / (target[base == "C", N] + target[base == "G", N]),
            var = (target[base == "C", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "G", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "C", N] + target[base == "G", N]) *
                   ((target[base == "C", N * mean] + target[base == "G", N * mean]) / (target[base == "C", N] + target[base == "G", N])) ^ 2) /
                   (target[base == "C", N] + target[base == "G", N] - 1),
            N = target[base == "C", N] + target[base == "G", N]),
        data.table(base = "M", sample = sample1, type = type1,
            mean = (target[base == "A", mean * N] + target[base == "C", mean * N]) / (target[base == "A", N] + target[base == "C", N]),
            var = (target[base == "A", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "C", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "A", N] + target[base == "C", N]) *
                   ((target[base == "A", N * mean] + target[base == "C", N * mean]) / (target[base == "A", N] + target[base == "C", N])) ^ 2) /
                   (target[base == "A", N] + target[base == "C", N] - 1),
            N = target[base == "A", N] + target[base == "C", N]),
        data.table(base = "K", sample = sample1, type = type1,
            mean = (target[base == "G", mean * N] + target[base == "T", mean * N]) / (target[base == "G", N] + target[base == "T", N]),
            var = (target[base == "G", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "T", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "G", N] + target[base == "T", N]) *
                   ((target[base == "G", N * mean] + target[base == "T", N * mean]) / (target[base == "G", N] + target[base == "T", N])) ^ 2) /
                   (target[base == "G", N] + target[base == "T", N] - 1),
            N = target[base == "G", N] + target[base == "T", N]),
        data.table(base = "R", sample = sample1, type = type1,
            mean = (target[base == "A", mean * N] + target[base == "G", mean * N]) / (target[base == "A", N] + target[base == "G", N]),
            var = (target[base == "A", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "G", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "A", N] + target[base == "G", N]) *
                   ((target[base == "A", N * mean] + target[base == "G", N * mean]) / (target[base == "A", N] + target[base == "G", N])) ^ 2) /
                   (target[base == "A", N] + target[base == "G", N] - 1),
            N = target[base == "A", N] + target[base == "G", N]),
        data.table(base = "Y", sample = sample1, type = type1,
            mean = (target[base == "C", mean * N] + target[base == "T", mean * N]) / (target[base == "C", N] + target[base == "T", N]),
            var = (target[base == "C", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "T", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "C", N] + target[base == "T", N]) *
                   ((target[base == "C", N * mean] + target[base == "T", N * mean]) / (target[base == "C", N] + target[base == "T", N])) ^ 2) /
                   (target[base == "C", N] + target[base == "T", N] - 1),
            N = target[base == "C", N] + target[base == "T", N]),
        data.table(base = "B", sample = sample1, type = type1,
            mean = (target[base == "C", mean * N] + target[base == "G", mean * N] + target[base == "T", mean * N]) / (target[base == "C", N] + target[base == "G", N] + target[base == "T", N]),
            var = (target[base == "C", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "G", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "T", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "C", N] + target[base == "G", N] + target[base == "T", N]) *
                   ((target[base == "C", N * mean] + target[base == "G", N * mean] + target[base == "T", N * mean]) / (target[base == "C", N] + target[base == "G", N] + target[base == "T", N])) ^ 2) /
                   (target[base == "C", N] + target[base == "G", N] + target[base == "T", N] - 1),
            N = target[base == "C", N] + target[base == "G", N] + target[base == "T", N]),
        data.table(base = "D", sample = sample1, type = type1,
            mean = (target[base == "A", mean * N] + target[base == "G", mean * N] + target[base == "T", mean * N]) / (target[base == "A", N] + target[base == "G", N] + target[base == "T", N]),
            var = (target[base == "A", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "G", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "T", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "A", N] + target[base == "G", N] + target[base == "T", N]) *
                   ((target[base == "A", N * mean] + target[base == "G", N * mean] + target[base == "T", N * mean]) / (target[base == "A", N] + target[base == "G", N] + target[base == "T", N])) ^ 2) /
                   (target[base == "A", N] + target[base == "G", N] + target[base == "T", N] - 1),
            N = target[base == "A", N] + target[base == "G", N] + target[base == "T", N]),
        data.table(base = "H", sample = sample1, type = type1,
            mean = (target[base == "A", mean * N] + target[base == "C", mean * N] + target[base == "T", mean * N]) / (target[base == "A", N] + target[base == "C", N] + target[base == "T", N]),
            var = (target[base == "A", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "C", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "T", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "A", N] + target[base == "C", N] + target[base == "T", N]) *
                   ((target[base == "A", N * mean] + target[base == "C", N * mean] + target[base == "T", N * mean]) / (target[base == "A", N] + target[base == "C", N] + target[base == "T", N])) ^ 2) /
                   (target[base == "A", N] + target[base == "C", N] + target[base == "T", N] - 1),
            N = target[base == "A", N] + target[base == "C", N] + target[base == "T", N]),
        data.table(base = "V", sample = sample1, type = type1,
            mean = (target[base == "A", mean * N] + target[base == "C", mean * N] + target[base == "G", mean * N]) / (target[base == "A", N] + target[base == "C", N] + target[base == "G", N]),
            var = (target[base == "A", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "C", (N - 1) * (sd ^ 2) + N * mean ^ 2] +
                   target[base == "G", (N - 1) * (sd ^ 2) + N * mean ^ 2] -
                   (target[base == "A", N] + target[base == "C", N] + target[base == "G", N]) *
                   ((target[base == "A", N * mean] + target[base == "C", N * mean] + target[base == "G", N * mean]) / (target[base == "A", N] + target[base == "C", N] + target[base == "G", N])) ^ 2) /
                   (target[base == "A", N] + target[base == "C", N] + target[base == "G", N] - 1),
            N = target[base == "A", N] + target[base == "C", N] + target[base == "G", N]),
        data.table(base = "N", sample = sample1, type = type1, mean = target[base == "all", mean], var = target[base == "all", sd ^ 2], N = target[base == "all", N])
        ))
    return(kinetics_new)
}

kinetics_stats_celegans <- fread("../kinetics_stats.c_elegans.csv")
kinetics_stats_celegans_degenerate <- rbindlist(list(get_degenerate_base_stats(kinetics_stats_celegans, "ab_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "cd_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "k_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "l_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "abcd_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "kl_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182sequel_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "k_normBy_ab_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "l_normBy_cd_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "kl_normBy_abcd_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "ab_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "cd_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182sequel_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "ab_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "cd_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182sequel_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "ab_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "cd_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "k_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "l_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "abcd_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "kl_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_celegans, "PD2182sequel_deep", "IPD")
                                                     ))
fwrite(kinetics_stats_celegans_degenerate, file = "kinetics_stats.c_elegans.degenerate.csv")

kinetics_stats_ecoli <- fread("../kinetics_stats.e_coli.csv")
kinetics_stats_ecoli_degenerate <- rbindlist(list(get_degenerate_base_stats(kinetics_stats_ecoli, "ab_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "cd_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "k_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "l_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "abcd_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "kl_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182sequel_deep", "log2_IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "k_normBy_ab_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "l_normBy_cd_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "kl_normBy_abcd_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "ab_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "cd_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182sequel_deep", "log2_IPD_ratio"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "ab_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "cd_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182sequel_deep", "log2_modelPrediction"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "ab_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "cd_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "k_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "l_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "abcd_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "kl_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182_deep", "IPD"),
                                                     get_degenerate_base_stats(kinetics_stats_ecoli, "PD2182sequel_deep", "IPD")
                                                     ))
fwrite(kinetics_stats_ecoli_degenerate, file = "kinetics_stats.e_coli.degenerate.csv")


