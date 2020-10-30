source("load_jun2018_data.saved.R")

library(data.table)
#score_thres <- 20
#ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
#coverage_thres <- 25
bases <- c("A", "C", "G", "T")

k_normBy_ab_data[, isCallable := list(is.finite(ipdRatio) & ipdRatio > 0.0)]
l_normBy_cd_data[, isCallable := list(is.finite(ipdRatio) & ipdRatio > 0.0)]
kl_normBy_abcd_data[, isCallable := list(is.finite(ipdRatio) & ipdRatio > 0.0)]
callableIntersection <- k_normBy_ab_data$isCallable & l_normBy_cd_data$isCallable
bothCallableNum <- sum(callableIntersection)

calc_jaccard <- function(score_thres, ipdRatio_thres, coverage_thres){
    k_normBy_ab_data[, isDeep := list(native_coverage >= coverage_thres & control_coverage >= coverage_thres)]
    k_normBy_ab_data[, isDeviated := list(score >= score_thres)]
    k_normBy_ab_data[, isModified := list(isDeviated & isDeep & ipdRatio >= ipdRatio_thres)]
    if(k_normBy_ab_data[isDeviated == TRUE, !all(isCallable)]){ warning("k_normBy_ab_data: There is (isDeviated AND not isCallable)") }

    l_normBy_cd_data[, isDeep := list(native_coverage >= coverage_thres & control_coverage >= coverage_thres)]
    l_normBy_cd_data[, isDeviated := list(score >= score_thres)]
    l_normBy_cd_data[, isModified := list(isDeviated & isDeep & ipdRatio >= ipdRatio_thres)]
    if(l_normBy_cd_data[isDeviated == TRUE, !all(isCallable)]){ warning("l_normBy_cd_data: There is (isDeviated AND not isCallable)") }
    
    kl_normBy_abcd_data[, isDeep := list(native_coverage >= coverage_thres & control_coverage >= coverage_thres)]
    kl_normBy_abcd_data[, isDeviated := list(score >= score_thres)]
    kl_normBy_abcd_data[, isModified := list(isDeviated & isDeep & ipdRatio >= ipdRatio_thres)]
    if(kl_normBy_abcd_data[isDeviated == TRUE, !all(isCallable)]){ warning("kl_normBy_abcd_data: There is (isDeviated AND not isCallable)") }

    deepIntersection <- k_normBy_ab_data$isDeep & l_normBy_cd_data$isDeep
    bothDeepNum <- sum(deepIntersection)
    deviatedIntersection <- k_normBy_ab_data$isDeviated & l_normBy_cd_data$isDeviated
    modifiedIntersection <- k_normBy_ab_data$isModified & l_normBy_cd_data$isModified
    union <- cbind(a = k_normBy_ab_data, b = l_normBy_cd_data)

    satisfactoryIntersection <- callableIntersection & deepIntersection
    bothSatisfactoryNum <- sum(satisfactoryIntersection)
    mod_stats <- data.table("Condition" = "Callable (With valid IPD ratios)",
                            "k_normBy_ab" = union[a.isCallable == TRUE, .N],
                            "l_normBy_cd" = union[b.isCallable == TRUE, .N],
                            "k_normBy_ab in both the satisfactory regions" = bothCallableNum,
                            "l_normBy_cd in both the satisfactory regions" = bothCallableNum,
                            "Intersection" = bothCallableNum,
                            "Jaccard index" = 1.0,
                            "kl_normBy_abcd" = kl_normBy_abcd_data[isCallable == TRUE, .N])

    mod_stats <- rbind(mod_stats, list(sprintf("Deep (coverage >= %g)", coverage_thres), union[a.isDeep == TRUE, .N], union[b.isDeep == TRUE, .N],
                                       bothDeepNum, bothDeepNum, bothDeepNum, 1.0, kl_normBy_abcd_data[isDeep == TRUE, .N]))

    mod_stats <- rbind(mod_stats, list("Satisfactory (Callable & Deep)", union[a.isCallable & a.isDeep, .N], union[b.isCallable & b.isDeep, .N],
                                       bothSatisfactoryNum, bothSatisfactoryNum, bothSatisfactoryNum, 1.0, kl_normBy_abcd_data[isCallable & isDeep, .N]))

    a <- union[satisfactoryIntersection & a.isDeviated, .N]
    b <- union[satisfactoryIntersection & b.isDeviated, .N]
    i <- sum(satisfactoryIntersection & deviatedIntersection)
    mod_stats <- rbind(mod_stats, list(sprintf("score >= %g", score_thres), union[a.isDeviated == TRUE, .N], union[b.isDeviated == TRUE, .N], a, b, i, i / (a + b - i),
                                       kl_normBy_abcd_data[isDeviated == TRUE, .N]))
    
    a <- union[satisfactoryIntersection & a.isModified, .N]
    b <- union[satisfactoryIntersection & b.isModified, .N]
    i <- sum(modifiedIntersection)
    mod_stats <- rbind(mod_stats, list(sprintf("Satisfactory & score >= %g & ipdRatio >= %g", score_thres, ipdRatio_thres),
        union[a.isModified == TRUE, .N], union[b.isModified == TRUE, .N], a, b, i, i / (a + b - i), kl_normBy_abcd_data[isModified == TRUE, .N]))

    for (BASE in bases) {
        a <- union[satisfactoryIntersection & a.base == BASE, .N]
        b <- a #union[satisfactoryIntersection & b.base == BASE, .N]
        i <- a #union[satisfactoryIntersection & a.base == BASE & b.base == BASE, .N]
        mod_stats <- rbind(mod_stats, list(sprintf("Satisfactory & base = %s", BASE),
            union[a.isCallable & a.isDeep & a.base == BASE, .N], union[b.isCallable & b.isDeep & b.base == BASE, .N], a, b, i, 1.0, kl_normBy_abcd_data[isCallable & isDeep & base == BASE, .N]))
        a <- union[satisfactoryIntersection & a.isModified & a.base == BASE, .N]
        b <- union[satisfactoryIntersection & b.isModified & b.base == BASE, .N]
        i <- union[modifiedIntersection & a.base == BASE, .N]
        mod_stats <- rbind(mod_stats, list(sprintf("Satisfactory & score >= %g & ipdRatio >= %g & base = %s", score_thres, ipdRatio_thres, BASE),
            union[a.isModified & a.base == BASE, .N], union[b.isModified & b.base == BASE, .N], a, b, i, i / (a + b - i), kl_normBy_abcd_data[isModified & base == BASE, .N]))
    }
    return(mod_stats)
}


d1 <- calc_jaccard(score_thres = 20, ipdRatio_thres = 2.0, coverage_thres = 25)
d1
d2 <- calc_jaccard(score_thres = 20, ipdRatio_thres = 4.0, coverage_thres = 25)
d2

d <- rbind(d1, d2[5:.N])
fwrite(d, file = "call_modification.calc_jaccard.csv")

