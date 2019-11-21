source("load_jun2018_data.saved.fst.wga.R")

library(data.table)
#score_thres <- 20
#ipdRatio_thres <- 2.0
#ipdRatio_thres <- 4.0
coverage_thres <- 25
collect_frac <- 0.01
bases <- c("A", "C", "G", "T")

ab_data[, isCallable := list(is.finite(tMean) & tMean > 0.0)]
cd_data[, isCallable := list(is.finite(tMean) & tMean > 0.0)]
abcd_data[, isCallable := list(is.finite(tMean) & tMean > 0.0)]
callableIntersection <- ab_data$isCallable & cd_data$isCallable
bothCallableNum <- sum(callableIntersection)

get_thres <- function(data_orig){
    data <- copy(data_orig)
    data <- data[tMean > 0]
    setkey(data, tMean)
    high_thres <- data[min(ceiling((1 - collect_frac) * .N) + 1, .N), tMean]
    low_thres <- data[max(floor(collect_frac * .N), 1), tMean]
    data <- data[coverage >= coverage_thres]
    high_thres_with_cov <- data[min(ceiling((1 - collect_frac) * .N) + 1, .N), tMean]
    low_thres_with_cov <- data[max(floor(collect_frac * .N), 1), tMean]
    return(c(high_thres, low_thres, high_thres_with_cov, low_thres_with_cov))
}

calc_jaccard <- function(){
    ab_thres <- get_thres(ab_data)
    ab_data[, isDeep := list(coverage >= coverage_thres)]
    ab_data[, isHigh := list(isCallable & tMean >= ab_thres[1])]
    ab_data[, isLow := list(isCallable & tMean <= ab_thres[2])]
    ab_data[, isHighWithCov := list(isDeep & tMean >= ab_thres[3])]
    ab_data[, isLowWithCov := list(isDeep & tMean <= ab_thres[4])]
    if(ab_data[isDeep == TRUE, !all(isCallable)]){ warning("ab_data: There is (isDeep AND !isCallable)") }

    cd_thres <- get_thres(cd_data)
    cd_data[, isDeep := list(coverage >= coverage_thres)]
    cd_data[, isHigh := list(isCallable & tMean >= cd_thres[1])]
    cd_data[, isLow := list(isCallable & tMean <= cd_thres[2])]
    cd_data[, isHighWithCov := list(isDeep & tMean >= cd_thres[3])]
    cd_data[, isLowWithCov := list(isDeep & tMean <= cd_thres[4])]
    if(cd_data[isDeep == TRUE, !all(isCallable)]){ warning("cd_data: There is (isDeep AND !isCallable)") }
    
    abcd_thres <- get_thres(abcd_data)
    abcd_data[, isDeep := list(coverage >= coverage_thres)]
    abcd_data[, isHigh := list(isCallable & tMean >= abcd_thres[1])]
    abcd_data[, isLow := list(isCallable & tMean <= abcd_thres[2])]
    abcd_data[, isHighWithCov := list(isDeep & tMean >= abcd_thres[3])]
    abcd_data[, isLowWithCov := list(isDeep & tMean <= abcd_thres[4])]
    if(abcd_data[isDeep == TRUE, !all(isCallable)]){ warning("abcd_data: There is (isDeep AND !isCallable)") }

    deepIntersection <- ab_data$isDeep & cd_data$isDeep
    bothDeepNum <- sum(deepIntersection)
    highIntersection <- ab_data$isHigh & cd_data$isHigh
    lowIntersection <- ab_data$isLow & cd_data$isLow
    highWithCovIntersection <- ab_data$isHighWithCov & cd_data$isHighWithCov
    lowWithCovIntersection <- ab_data$isLowWithCov & cd_data$isLowWithCov
    union <- cbind(a = ab_data, b = cd_data)

    satisfactoryIntersection <- callableIntersection & deepIntersection
    bothSatisfactoryNum <- sum(satisfactoryIntersection)
    mod_stats <- data.table("Condition" = "Callable (With valid IPD)",
                            "ab" = union[a.isCallable == TRUE, .N],
                            "cd" = union[b.isCallable == TRUE, .N],
                            "ab in both the satisfactory regions" = bothCallableNum,
                            "cd in both the satisfactory regions" = bothCallableNum,
                            "Intersection" = bothCallableNum,
                            "Jaccard index" = 1.0,
                            "abcd" = abcd_data[isCallable == TRUE, .N])

    mod_stats <- rbind(mod_stats, list(sprintf("Deep (coverage >= %g)", coverage_thres), union[a.isDeep == TRUE, .N], union[b.isDeep == TRUE, .N],
                                       bothDeepNum, bothDeepNum, bothDeepNum, 1.0, abcd_data[isDeep == TRUE, .N]))

    mod_stats <- rbind(mod_stats, list("Satisfactory (Callable & Deep)", union[a.isCallable & a.isDeep, .N], union[b.isCallable & b.isDeep, .N],
                                       bothSatisfactoryNum, bothSatisfactoryNum, bothSatisfactoryNum, 1.0, abcd_data[isCallable & isDeep, .N]))

    a <- union[satisfactoryIntersection & a.isHigh, .N]
    b <- union[satisfactoryIntersection & b.isHigh, .N]
    i <- sum(satisfactoryIntersection & highIntersection)
    mod_stats <- rbind(mod_stats, list(sprintf("high IPD (thres:ab=%.3g,cd=%.3g,abcd=%.3g)", ab_thres[1], cd_thres[1], abcd_thres[1]), union[a.isHigh == TRUE, .N], union[b.isHigh == TRUE, .N], a, b, i, i / (a + b - i),
                                       abcd_data[isHigh == TRUE, .N]))
    
    a <- union[satisfactoryIntersection & a.isLow, .N]
    b <- union[satisfactoryIntersection & b.isLow, .N]
    i <- sum(satisfactoryIntersection & lowIntersection)
    mod_stats <- rbind(mod_stats, list(sprintf("low IPD (thres:ab=%.3g,cd=%.3g,abcd=%.3g)", ab_thres[2], cd_thres[2], abcd_thres[2]), union[a.isLow == TRUE, .N], union[b.isLow == TRUE, .N], a, b, i, i / (a + b - i),
                                       abcd_data[isLow == TRUE, .N]))

    a <- union[satisfactoryIntersection & a.isHighWithCov, .N]
    b <- union[satisfactoryIntersection & b.isHighWithCov, .N]
    i <- sum(satisfactoryIntersection & highWithCovIntersection)
    mod_stats <- rbind(mod_stats, list(sprintf("high IPD with coverage >= %g (thres:ab=%.3g,cd=%.3g,abcd=%.3g)", coverage_thres, ab_thres[3], cd_thres[3], abcd_thres[3]), union[a.isHighWithCov == TRUE, .N], union[b.isHighWithCov == TRUE, .N], a, b, i, i / (a + b - i),
                                       abcd_data[isHighWithCov == TRUE, .N]))
    
    a <- union[satisfactoryIntersection & a.isLowWithCov, .N]
    b <- union[satisfactoryIntersection & b.isLowWithCov, .N]
    i <- sum(satisfactoryIntersection & lowWithCovIntersection)
    mod_stats <- rbind(mod_stats, list(sprintf("low IPD with coverage >= %g (thres:ab=%.3g,cd=%.3g,abcd=%.3g)", coverage_thres, ab_thres[4], cd_thres[4], abcd_thres[4]), union[a.isLowWithCov == TRUE, .N], union[b.isLowWithCov == TRUE, .N], a, b, i, i / (a + b - i),
                                       abcd_data[isLowWithCov == TRUE, .N]))
    
#    for (BASE in bases) {
#        a <- union[satisfactoryIntersection & a.base == BASE, .N]
#        b <- a #union[satisfactoryIntersection & b.base == BASE, .N]
#        i <- a #union[satisfactoryIntersection & a.base == BASE & b.base == BASE, .N]
#        mod_stats <- rbind(mod_stats, list(sprintf("Satisfactory & base = %s", BASE),
#            union[a.isCallable & a.isDeep & a.base == BASE, .N], union[b.isCallable & b.isDeep & b.base == BASE, .N], a, b, i, 1.0, kl_normBy_abcd_data[isCallable & isDeep & base == BASE, .N]))
#        a <- union[satisfactoryIntersection & a.isModified & a.base == BASE, .N]
#        b <- union[satisfactoryIntersection & b.isModified & b.base == BASE, .N]
#        i <- union[modifiedIntersection & a.base == BASE, .N]
#        mod_stats <- rbind(mod_stats, list(sprintf("Satisfactory & score >= %g & ipdRatio >= %g & base = %s", score_thres, ipdRatio_thres, BASE),
#            union[a.isModified & a.base == BASE, .N], union[b.isModified & b.base == BASE, .N], a, b, i, i / (a + b - i), kl_normBy_abcd_data[isModified & base == BASE, .N]))
#    }
    return(mod_stats)
}


d1 <- calc_jaccard()
d1

d <- rbind(d1)
fwrite(d, file = "call_extreme_ipd.calc_jaccard.csv")

