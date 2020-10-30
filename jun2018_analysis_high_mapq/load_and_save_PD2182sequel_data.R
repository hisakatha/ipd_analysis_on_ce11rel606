library(data.table)
library(Biostrings)
library(hdf5r)

reference <- readDNAStringSet("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa")

chrs <- names(reference)
ecoli_chr <- "E._coli_REL606"

parse_kinetics_h5 <- function(file, chromosomes) {
    data_hdf5 <- H5File$new(file, mode = "r")
    data <- NULL
    for (chr in chromosomes) {
        data_chr <- data.table(
            refName = chr,
            tpl = data_hdf5[[paste0(chr, "/tpl")]][],
            strand = data_hdf5[[paste0(chr, "/strand")]][],
            base = data_hdf5[[paste0(chr, "/base")]][],
            score = data_hdf5[[paste0(chr, "/score")]][],
            tMean = data_hdf5[[paste0(chr, "/tMean")]][],
            tErr = data_hdf5[[paste0(chr, "/tErr")]][],
            modelPrediction = data_hdf5[[paste0(chr, "/modelPrediction")]][],
            ipdRatio = data_hdf5[[paste0(chr, "/ipdRatio")]][],
            coverage = data_hdf5[[paste0(chr, "/coverage")]][],
            frac = data_hdf5[[paste0(chr, "/frac")]][],
            fracLow = data_hdf5[[paste0(chr, "/fracLow")]][],
            fracUp = data_hdf5[[paste0(chr, "/fracUp")]][])
        data <- rbind(data, data_chr)
    }
    data_hdf5$close_all()
    return(data)
}

# AB, CD
source("load_jun2018_data.saved.fst.basic_wga.R")
cat("Loaded data AB\n")
cat("Loaded data CD\n")

# ABCD
#cat("Loaded data ABCD\n")

# PD2182-sequel
data_path <- "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/PD2182-sequel_no_chunk/ipd_summary_high_mapq.h5"
PD2182sequel_data <- parse_kinetics_h5(data_path, chrs)
cat("Loaded data PD2182sequel\n")

calculate_score <- function(native_tMean, native_tErr, native_coverage, control_tMean, control_tErr, control_coverage) {
    # tErr calculated by kineticsTools is based on population variance (denominator N) as of Aug 2nd, 2018.
    # tErr = sqrt( population variance ) / sqrt(N)
    native_tErr_sq <- ifelse(native_tMean > 0 & native_tErr > 0 & native_coverage > 1, (native_tErr ^ 2) * native_coverage / (native_coverage - 1), NA)
    control_tErr_sq <- ifelse(control_tMean > 0 & control_tErr > 0 & control_coverage > 1, (control_tErr ^ 2) * control_coverage / (control_coverage - 1), NA)
    t_stat <- (native_tMean - control_tMean) / sqrt(native_tErr_sq + control_tErr_sq)
    dof <- (native_tErr_sq + control_tErr_sq)^2 / ((native_tErr_sq ^ 2) / (native_coverage - 1) + (control_tErr_sq ^ 2) / (control_coverage - 1))
    # Two-sided p-value
    pvalue <- pt(-abs(t_stat), dof) * 2
    return(ifelse(is.na(pvalue), 0, -10 * log10(pvalue)))
}

normalize_native_by_control <- function(native, control) {
    mod_score = calculate_score(native$tMean, native$tErr, native$coverage, control$tMean, control$tErr, control$coverage)
    data.table(
        refName = native$refName,
        tpl = native$tpl,
        strand = native$strand,
        base = native$base,
        score = mod_score,
        ipdRatio = native$tMean / control$tMean,
        native_score = native$score,
        control_score = control$score,
        native_tMean = native$tMean,
        control_tMean = control$tMean,
        native_tErr = native$tErr,
        control_tErr = control$tErr,
        native_coverage = native$coverage,
        control_coverage = control$coverage,
        modelPrediction = native$modelPrediction)
}

PD2182sequel_normBy_ab_data <- normalize_native_by_control(PD2182sequel_data, ab_data)
cat("Generated PD2182sequel_normBy_ab_data\n")
PD2182sequel_normBy_cd_data <- normalize_native_by_control(PD2182sequel_data, cd_data)
cat("Generated PD2182sequel_normBy_cd_data\n")

save(PD2182sequel_data, PD2182sequel_normBy_ab_data, PD2182sequel_normBy_cd_data, file = "load_and_save_PD2182sequel_data.Rdata")

library(fst)

write_fst(PD2182sequel_data, "load_and_save_PD2182sequel_data.PD2182sequel_data.fst", compress = 100)
write_fst(PD2182sequel_normBy_ab_data, "load_and_save_PD2182sequel_data.PD2182sequel_normBy_ab_data.fst", compress = 100)
write_fst(PD2182sequel_normBy_cd_data, "load_and_save_PD2182sequel_data.PD2182sequel_normBy_cd_data.fst", compress = 100)

