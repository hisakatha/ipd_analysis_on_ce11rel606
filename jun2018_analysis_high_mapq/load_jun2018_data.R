library(data.table)
library(Biostrings)
library(hdf5r)
library(fst)

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

# AB
data_path <- "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_ab_pcr_vc2010_op50_no_chunk/ipd_summary_high_mapq.h5"
ab_data <- parse_kinetics_h5(data_path, chrs)
cat("Loaded data AB\n")

# CD
data_path <- "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_cd_pcr_vc2010_no_chunk/ipd_summary_high_mapq.h5"
cd_data <- parse_kinetics_h5(data_path, chrs)
cat("Loaded data CD\n")

# K
data_path <- "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_k_nopcr_vc2010_op50_no_chunk/ipd_summary_high_mapq.h5"
k_data <- parse_kinetics_h5(data_path, chrs)
cat("Loaded data K\n")

# L
data_path <- "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_l_nopcr_vc2010_no_chunk/ipd_summary_high_mapq.h5"
l_data <- parse_kinetics_h5(data_path, chrs)
cat("Loaded data L\n")

# ABCD
#data_path <- "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_abcd_pcr/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5"
#abcd_data <- parse_kinetics_h5(data_path, chrs)
#cat("Loaded data ABCD\n")

# KL
#data_path <- "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_kl_nopcr/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5"
#kl_data <- parse_kinetics_h5(data_path, chrs)
#cat("Loaded data KL\n")

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

k_normBy_ab_data <- normalize_native_by_control(k_data, ab_data)
cat("Generated k_normBy_ab_data\n")
l_normBy_cd_data <- normalize_native_by_control(l_data, cd_data)
cat("Generated l_normBy_cd_data\n")
#kl_normBy_abcd_data <- normalize_native_by_control(kl_data, abcd_data)
#cat("Generated kl_normBy_abcd_data\n")

#save(ab_data, cd_data, k_data, l_data, abcd_data, kl_data, k_normBy_ab_data, l_normBy_cd_data, kl_normBy_abcd_data, file = "load_jun2018_data.RData")
save(ab_data, cd_data, k_data, l_data, k_normBy_ab_data, l_normBy_cd_data, file = "load_jun2018_data.RData")

write_fst(ab_data, "load_jun2018_data.ab_data.fst", compress = 100)
write_fst(cd_data, "load_jun2018_data.cd_data.fst", compress = 100)
write_fst(k_data, "load_jun2018_data.k_data.fst", compress = 100)
write_fst(l_data, "load_jun2018_data.l_data.fst", compress = 100)
#write_fst(abcd_data, "load_jun2018_data.abcd_data.fst", compress = 100)
#write_fst(kl_data, "load_jun2018_data.kl_data.fst", compress = 100)
write_fst(k_normBy_ab_data, "load_jun2018_data.k_normBy_ab_data.fst", compress = 100)
write_fst(l_normBy_cd_data, "load_jun2018_data.l_normBy_cd_data.fst", compress = 100)
#write_fst(kl_normBy_abcd_data, "load_jun2018_data.kl_normBy_abcd_data.fst", compress = 100)
