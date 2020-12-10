library(data.table)
library(ggplot2)
library(cowplot)

convert_sum_data <- function(data) {
    # data: a data.table with columns:
    # kmer_string,kmer_number,position,chromosome,file_index,ipd_sum,ipd_sq_sum,log2_ipd_sum,log2_ipd_sq_sum,prediction_sum,prediction_sq_sum,log2_prediction_sum,log2_prediction_sq_sum,count
    # merge by chromosome and file_index
    data[, lapply(.SD, sum), keyby = .(kmer_string, kmer_number, position), .SDcols = !c("chromosome", "file_index")][,
        .(kmer_string, kmer_number, position, strand = "+",
          ipd_mean = ipd_sum / count, ipd_var = ipd_sq_sum / (count - 1) - ipd_sum ^ 2 / count / (count - 1),
          log2_ipd_mean = log2_ipd_sum / count, log2_ipd_var = log2_ipd_sq_sum / (count - 1) - log2_ipd_sum ^ 2 / count / (count - 1),
          prediction_mean = prediction_sum / count, prediction_var = prediction_sq_sum / (count - 1) - prediction_sum ^ 2 / count / (count - 1),
          log2_prediction_mean = log2_prediction_sum / count, log2_prediction_var = log2_prediction_sq_sum / count - log2_prediction_sum ^ 2 / count / (count - 1), count)]
}

plot_value_per_kmer <- function(summary_data, title_text, ylab_text, y_expr) {
    xlab_text <- "4-mer"
    y_expr <- enquo(y_expr)
    plot_data <- summary_data[substr(kmer_string, 2, 2) == "A" & position == 2]
    plot_data[, "isGATC" := .(ifelse(kmer_string == "GATC", "GATC", "not GATC"))]
    p1 <- ggplot(plot_data, aes(kmer_string, !!y_expr, color = isGATC)) + geom_point() +
        xlab(xlab_text) + ylab(ylab_text) +
        ggtitle(title_text) +
        labs(color = "Motif") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
    print(p1)
}

ylab_text_ipd <- "Mean IPD"
ylab_text_pred <- "Mean predicted IPD"
ylab_text_mean_ratio <- "(Mean observed IPD) / (Mean predicted IPD)"
plot_multiple_yvalues <- function(data, data_text) {
    plot_value_per_kmer(data, data_text, ylab_text_ipd, ipd_mean)
    plot_value_per_kmer(data, data_text, ylab_text_pred, prediction_mean)
    plot_value_per_kmer(data, data_text, ylab_text_mean_ratio, ipd_mean / prediction_mean)
}

stats_ab <- fread("tetramer_ipd_stats.ab.csv")
stats_cd <- fread("tetramer_ipd_stats.cd.csv")
stats_k <- fread("tetramer_ipd_stats.k.csv")
stats_l <- fread("tetramer_ipd_stats.l.csv")
stats_PD2182sequel <- fread("tetramer_ipd_stats.PD2182sequel.csv")
stats_abcd <- fread("tetramer_ipd_stats.abcd.csv")
stats_kl <- fread("tetramer_ipd_stats.kl.csv")

ecoli_chr <- "E._coli_REL606"
mito_chr <- "chrM"
summary_ab_celegans <- convert_sum_data(stats_ab[chromosome != ecoli_chr & chromosome != mito_chr])
summary_cd_celegans <- convert_sum_data(stats_cd[chromosome != ecoli_chr & chromosome != mito_chr])
summary_k_celegans <- convert_sum_data(stats_k[chromosome != ecoli_chr & chromosome != mito_chr])
summary_l_celegans <- convert_sum_data(stats_l[chromosome != ecoli_chr & chromosome != mito_chr])
summary_PD2182sequel_celegans <- convert_sum_data(stats_PD2182sequel[chromosome != ecoli_chr & chromosome != mito_chr])
summary_abcd_celegans <- convert_sum_data(stats_abcd[chromosome != ecoli_chr & chromosome != mito_chr])
summary_kl_celegans <- convert_sum_data(stats_kl[chromosome != ecoli_chr & chromosome != mito_chr])
summary_ab_ecoli <- convert_sum_data(stats_ab[chromosome == ecoli_chr])
summary_cd_ecoli <- convert_sum_data(stats_cd[chromosome == ecoli_chr])
summary_k_ecoli <- convert_sum_data(stats_k[chromosome == ecoli_chr])
summary_l_ecoli <- convert_sum_data(stats_l[chromosome == ecoli_chr])
summary_PD2182sequel_ecoli <- convert_sum_data(stats_PD2182sequel[chromosome == ecoli_chr])
summary_abcd_ecoli <- convert_sum_data(stats_abcd[chromosome == ecoli_chr])
summary_kl_ecoli <- convert_sum_data(stats_kl[chromosome == ecoli_chr])
summary_ab_mito <- convert_sum_data(stats_ab[chromosome == mito_chr])
summary_cd_mito <- convert_sum_data(stats_cd[chromosome == mito_chr])
summary_k_mito <- convert_sum_data(stats_k[chromosome == mito_chr])
summary_l_mito <- convert_sum_data(stats_l[chromosome == mito_chr])
summary_PD2182sequel_mito <- convert_sum_data(stats_PD2182sequel[chromosome == mito_chr])
summary_abcd_mito <- convert_sum_data(stats_abcd[chromosome == mito_chr])
summary_kl_mito <- convert_sum_data(stats_kl[chromosome == mito_chr])

normalize_summarized_data <- function(data_native, data_wga) {
    merge(data_native, data_wga, by = c("kmer_string", "kmer_number", "position", "strand"), all = TRUE, suffixes = c(".native", ".wga"))[
        , "normalized_ipd_mean" := .(ipd_mean.native / ipd_mean.wga)]
}

summary_k_celegans_normBy_ab_celegans <- normalize_summarized_data(summary_k_celegans, summary_ab_celegans)
summary_l_celegans_normBy_cd_celegans <- normalize_summarized_data(summary_l_celegans, summary_cd_celegans)
summary_kl_celegans_normBy_abcd_celegans <- normalize_summarized_data(summary_kl_celegans, summary_abcd_celegans)
summary_k_mito_normBy_ab_mito <- normalize_summarized_data(summary_k_mito, summary_ab_mito)
summary_l_mito_normBy_cd_mito <- normalize_summarized_data(summary_l_mito, summary_cd_mito)
summary_kl_mito_normBy_abcd_mito <- normalize_summarized_data(summary_kl_mito, summary_abcd_mito)
summary_k_ecoli_normBy_ab_ecoli <- normalize_summarized_data(summary_k_ecoli, summary_ab_ecoli)
summary_l_ecoli_normBy_cd_ecoli <- normalize_summarized_data(summary_l_ecoli, summary_cd_ecoli)
summary_kl_ecoli_normBy_abcd_ecoli <- normalize_summarized_data(summary_kl_ecoli, summary_abcd_ecoli)

ylab_text_ipd_normalized <- "(Mean IPD in native) / (Mean IPD in WGA)"
pdf("ipd_at_second_adenine.pdf", height = 4, width = 10)
plot_value_per_kmer(summary_k_celegans_normBy_ab_celegans, "VC2010+OP50/C. elegans nuclear", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_l_celegans_normBy_cd_celegans, "VC2010/C. elegans nuclear", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_kl_celegans_normBy_abcd_celegans, "Merged/C. elegans nuclear", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_k_mito_normBy_ab_mito, "VC2010+OP50/C. elegans mitochondria", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_l_mito_normBy_cd_mito, "VC2010/C. elegans mitochondria", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_kl_mito_normBy_abcd_mito, "Merged/C. elegans mitochondria", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_k_ecoli_normBy_ab_ecoli, "VC2010+OP50/E. coli", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_l_ecoli_normBy_cd_ecoli, "VC2010/E. coli", ylab_text_ipd_normalized, normalized_ipd_mean)
plot_value_per_kmer(summary_kl_ecoli_normBy_abcd_ecoli, "Merged/E. coli", ylab_text_ipd_normalized, normalized_ipd_mean)

plot_multiple_yvalues(summary_ab_celegans, "VC2010+OP50/WGA/C. elegans nuclear")
plot_multiple_yvalues(summary_cd_celegans, "VC2010/WGA/C. elegans nuclear")
plot_multiple_yvalues(summary_k_celegans, "VC2010+OP50/native/C. elegans nuclear")
plot_multiple_yvalues(summary_l_celegans, "VC2010/native/C. elegans nuclear")
plot_multiple_yvalues(summary_PD2182sequel_celegans, "PD2182/native/C. elegans nuclear")
plot_multiple_yvalues(summary_abcd_celegans, "Merged/WGA/C. elegans nuclear")
plot_multiple_yvalues(summary_kl_celegans, "Merged/native/C. elegans nuclear")
plot_multiple_yvalues(summary_ab_mito, "VC2010+OP50/WGA/C. elegans mitochondria")
plot_multiple_yvalues(summary_cd_mito, "VC2010/WGA/C. elegans mitochondria")
plot_multiple_yvalues(summary_k_mito, "VC2010+OP50/native/C. elegans mitochondria")
plot_multiple_yvalues(summary_l_mito, "VC2010/native/C. elegans mitochondria")
plot_multiple_yvalues(summary_PD2182sequel_mito, "PD2182/native/C. elegans mitochondria")
plot_multiple_yvalues(summary_abcd_mito, "Merged/WGA/C. elegans mitochondria")
plot_multiple_yvalues(summary_kl_mito, "Merged/native/C. elegans mitochondria")
plot_multiple_yvalues(summary_ab_ecoli, "VC2010+OP50/WGA/E. coli")
plot_multiple_yvalues(summary_cd_ecoli, "VC2010/WGA/E. coli")
plot_multiple_yvalues(summary_k_ecoli, "VC2010+OP50/native/E. coli")
plot_multiple_yvalues(summary_l_ecoli, "VC2010/native/E. coli")
plot_multiple_yvalues(summary_PD2182sequel_ecoli, "PD2182/native/E. coli")
plot_multiple_yvalues(summary_abcd_ecoli, "Merged/WGA/E. coli")
plot_multiple_yvalues(summary_kl_ecoli, "Merged/native/E. coli")
invisible(dev.off())
