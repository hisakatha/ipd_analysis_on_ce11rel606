library(data.table)
# unique motifs
motifs_high_ipd <- c("AAAABT", "AGAGAGTA", "AGCTATAT", "AGGCAGGC", "ATGGGAYA", "ATTGTTAC",
"CAAACTAC", "CAGYTG", "CCAATCAG", "CRACGAS",
"DCGAGACC", "DGCTTC",
"GAAGGATC", "GATATRGY", "GCGACCTA", "GCGCGCGC", "GGHGGY", "GTAGATCA", "GTATCGTA",
"TGGTGSA", "YGGAR")
motifs_low_ipd <- c("AATMAATA", "AGACGCAG",
"CGGYTTGA", "CTKCAA",
"GCGCGTCA", "GTGTGTGY",
"TACCCCKA", "TACCTTGA")
motifs_high_ipd <- c(c("ACGCRTG", "ATCAGCTG", "GGNGGNGGNGGN", "RGTA"), motifs_high_ipd)
motifs_low_ipd <- c(c("ACATMTGG", "CTGDAR", "TACTGTAG"), motifs_low_ipd)
# common motifs
motifs_high_and_low_ipd <- "TGACGTCA"

motifs_select1 <- c("GAGG", "AGAA", "GATC")

motif_types <- c(rep("high_ipd_motif", length(motifs_high_ipd)),
                 rep("high_ipd_motif|low_ipd_motif", length(motifs_high_and_low_ipd)),
                 rep("low_ipd_motif", length(motifs_low_ipd)),
                 rep("6mA_motif", length(motifs_select1)))
names(motif_types) <- c(motifs_high_ipd, motifs_high_and_low_ipd, motifs_low_ipd, motifs_select1)

occ_path_prefix <- "count_occ.celegans_and_ecoli.binary_fisher.set1"
occ_data <- fread(sprintf("%s.csv", occ_path_prefix))
occ_data[, "motif_type" := .(motif_types[motif])]
# Shift motif_type to the fourth column
setcolorder(occ_data, c("sample_id", "chr", "motif", "motif_type"))
fwrite(occ_data[, lapply(.SD, as.character)], file = sprintf("%s.v2.csv", occ_path_prefix))
