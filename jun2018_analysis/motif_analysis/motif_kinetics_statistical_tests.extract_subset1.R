library(data.table)
library(doParallel)
d_celegans <- fread("motif_kinetics_statistical_tests.c_elegans.csv")
d_ecoli <- fread("motif_kinetics_statistical_tests.e_coli.csv")

subset1_motifs <- c("GAGG", "AGAA", "GATC", "ACGCRTG", "ATCAGCTG", "(GGN)4",
    "AGCTATAT", "CAGYTG", "CRACGAS", "DCGAGACC", "GAAGGATC", "GATATRGY",
    "GCGACCTA", "GCGCGCGC", "GGHGGY", "GTAGATCA", "GTATCGTA", "TGACGTCA",
    "TGGTGSA", "CGGYTTGA", "GCGCGTCA")
subset1_positions <- c(2, 3, 2, 4, 5, 3,
    3, 3, 6, 2, 2, 2,
    4, 2, 3, 1, 1, 3,
    3, 4, 4)

stopifnot(length(subset1_motifs) == length(subset1_positions))

extract_data <- function(data1, data2, motif1, position1) {
    rbindlist(list(
        data1[motif_name == motif1 & position == position1 & strand == "+" & genome_stats_type == "log2_IPD", .(motif_name, position, sample_name, pvalue, ks_normal_pvalue)],
        data2[motif_name == motif1 & position == position1 & strand == "+" & genome_stats_type == "log2_IPD", .(motif_name, position, sample_name, pvalue, ks_normal_pvalue)]))
}

d_celegans <- d_celegans[sample_id != "PD2182_deep"]
d_ecoli <- d_ecoli[sample_id != "PD2182_deep"]

extracted_tables <- foreach(i = 1:length(subset1_motifs)) %do% { extract_data(d_celegans, d_ecoli, subset1_motifs[i], subset1_positions[i]) }
extracted_data <- rbindlist(extracted_tables)

extracted_data <- extracted_data[is.finite(pvalue)]
#n_test <- extracted_data[, .N]
#cat(sprintf("n_test: %d\n", n_test))
#extracted_data[, "pvalue_bonferroni_adjusted" := .(pvalue * n_test)]

fwrite(extracted_data, file = "motif_kinetics_statistical_tests.subset1.csv")
