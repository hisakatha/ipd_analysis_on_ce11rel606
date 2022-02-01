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

plot_value_per_kmer.wga_vs_native <- function(summary_data_wga, summary_data_native, title_text, ylab_text, y_expr) {
    xlab_text <- "4-mer"
    y_expr <- enquo(y_expr)
    type1 <- "WGA"
    type2 <- "native"
    summary_data_wga2 <- copy(summary_data_wga)
    summary_data_native2 <- copy(summary_data_native)
    summary_data_wga2[, "type" := .(type1)]
    summary_data_native2[, "type" := .(type2)]
    summary_data <- rbindlist(list(summary_data_wga2, summary_data_native2))
    summary_data$type <- factor(summary_data$type, levels = c(type1, type2))
    plot_data <- summary_data[substr(kmer_string, 2, 2) == "A" & position == 2]
    p1 <- ggplot(plot_data, aes(kmer_string, !!y_expr, color = type)) + geom_point(alpha = 0.8) +
        xlab(xlab_text) + ylab(ylab_text) +
        ggtitle(title_text) +
        labs(color = "Sample") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, family = "mono", face = "bold"))
    print(p1)
}

ylab_text_ipd <- "Mean IPD"

stats_ab <- fread("tetramer_ipd_stats.ab.csv")
stats_cd <- fread("tetramer_ipd_stats.cd.csv")
stats_k <- fread("tetramer_ipd_stats.k.csv")
stats_l <- fread("tetramer_ipd_stats.l.csv")
stats_abcd <- fread("tetramer_ipd_stats.abcd.csv")
stats_kl <- fread("tetramer_ipd_stats.kl.csv")

ecoli_chr <- "E._coli_REL606"
mito_chr <- "chrM"
summary_ab_celegans <- convert_sum_data(stats_ab[chromosome != ecoli_chr])
summary_cd_celegans <- convert_sum_data(stats_cd[chromosome != ecoli_chr])
summary_k_celegans <- convert_sum_data(stats_k[chromosome != ecoli_chr])
summary_l_celegans <- convert_sum_data(stats_l[chromosome != ecoli_chr])
summary_abcd_celegans <- convert_sum_data(stats_abcd[chromosome != ecoli_chr])
summary_kl_celegans <- convert_sum_data(stats_kl[chromosome != ecoli_chr])
summary_ab_celegans_nuclear <- convert_sum_data(stats_ab[chromosome != ecoli_chr & chromosome != mito_chr])
summary_cd_celegans_nuclear <- convert_sum_data(stats_cd[chromosome != ecoli_chr & chromosome != mito_chr])
summary_k_celegans_nuclear <- convert_sum_data(stats_k[chromosome != ecoli_chr & chromosome != mito_chr])
summary_l_celegans_nuclear <- convert_sum_data(stats_l[chromosome != ecoli_chr & chromosome != mito_chr])
summary_abcd_celegans_nuclear <- convert_sum_data(stats_abcd[chromosome != ecoli_chr & chromosome != mito_chr])
summary_kl_celegans_nuclear <- convert_sum_data(stats_kl[chromosome != ecoli_chr & chromosome != mito_chr])
summary_ab_ecoli <- convert_sum_data(stats_ab[chromosome == ecoli_chr])
summary_cd_ecoli <- convert_sum_data(stats_cd[chromosome == ecoli_chr])
summary_k_ecoli <- convert_sum_data(stats_k[chromosome == ecoli_chr])
summary_l_ecoli <- convert_sum_data(stats_l[chromosome == ecoli_chr])
summary_abcd_ecoli <- convert_sum_data(stats_abcd[chromosome == ecoli_chr])
summary_kl_ecoli <- convert_sum_data(stats_kl[chromosome == ecoli_chr])
summary_ab_mito <- convert_sum_data(stats_ab[chromosome == mito_chr])
summary_cd_mito <- convert_sum_data(stats_cd[chromosome == mito_chr])
summary_k_mito <- convert_sum_data(stats_k[chromosome == mito_chr])
summary_l_mito <- convert_sum_data(stats_l[chromosome == mito_chr])
summary_abcd_mito <- convert_sum_data(stats_abcd[chromosome == mito_chr])
summary_kl_mito <- convert_sum_data(stats_kl[chromosome == mito_chr])

cairo_pdf("ipd_at_second_adenine.wga_vs_native.pdf", height = 2.5, width = 9, onefile = TRUE)
plot_value_per_kmer.wga_vs_native(summary_ab_ecoli, summary_k_ecoli, "Replicate 1 / E. coli", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_ab_celegans, summary_k_celegans, "Replicate 1 / C. elegans", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_ab_celegans_nuclear, summary_k_celegans_nuclear, "Replicate 1 / C. elegans nuclear", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_ab_mito, summary_k_mito, "Replicate 1 / C. elegans mitochondria", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_cd_ecoli, summary_l_ecoli, "Replicate 2 / E. coli", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_cd_celegans, summary_l_celegans, "Replicate 2 / C. elegans", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_cd_celegans_nuclear, summary_l_celegans_nuclear, "Replicate 2 / C. elegans nuclear", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_cd_mito, summary_l_mito, "Replicate 2 / C. elegans mitochondria", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_abcd_ecoli, summary_kl_ecoli, "Merged / E. coli", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_abcd_celegans, summary_kl_celegans, "Merged / C. elegans", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_abcd_celegans_nuclear, summary_kl_celegans_nuclear, "Merged / C. elegans nuclear", ylab_text_ipd, ipd_mean)
plot_value_per_kmer.wga_vs_native(summary_abcd_mito, summary_kl_mito, "Merged / C. elegans mitochondria", ylab_text_ipd, ipd_mean)
invisible(dev.off())
