library(data.table)
#library(ggplot2)
#library(cowplot)

welch.t.pvalue <- function(a_mean, a_variance, a_size, b_mean, b_variance, b_size) {
    a_sq <- (a_variance ^ 2) / a_size
    b_sq <- (b_variance ^ 2) / b_size
    t_stat <- (a_mean - b_mean) / sqrt(a_sq + b_sq)
    dof <- ((a_sq + b_sq) ^ 2) / ((a_sq ^ 2) / (a_size - 1) + (b_sq ^ 2) / (b_size - 1))
    # Two-sided p-value
    pt(-abs(t_stat), dof) * 2
}

complement_bases <- list("A" = "T", "C" = "G", "G" = "C", "T" = "A", "U" = "A",
    "W" = "W", "S" = "S", "M" = "K", "K" = "M", "R" = "Y", "Y" = "R",
    "B" = "V", "D" = "H", "H" = "D", "V" = "B", "N" = "N")

test_motif_kinetics <- function(kinetics, motif_name, motif_string, outside_length, sample_id, sample_name, type1, genome_stats, output_path, print_header){
    kinetics <- kinetics[value > 0 & is.finite(value)]
    if (kinetics[, .N] == 0) {
        return(data.table())
    }
    if (type1 == "IPD") {
        kinetics_summary <- kinetics[, .(mean = mean(value), var = var(value), .N, ks_normal_pvalue = ks.test(value, "pnorm", mean = mean(value), sd = sd(value))$p.value), by = .(position, strand, label)]
    } else if (type1 == "log2_IPD") {
        kinetics_summary <- kinetics[, "log2_value" := .(log2(value))][, .(mean = mean(log2_value), var = var(log2_value), .N, ks_normal_pvalue = ks.test(log2_value, "pnorm", mean = mean(log2_value), sd = sd(log2_value))$p.value), by = .(position, strand, label)]
    } else {
        stop("Unexpected genome_stats type")
    }
    kinetics_summary[, "position" := .(position - outside_length)]
    kinetics_summary[, c("motif_name", "motif_string", "sample_id", "sample_name", "genome_stats_type") := .(motif_name, motif_string, sample_id, sample_name, type1)]
    genome_stats_target <- genome_stats[sample == sample_id & type == type1]
    for (i in 1:nrow(kinetics_summary)) {
        region_label <- substr(kinetics_summary[i, label], 1, 1)
        if (region_label == "s" || region_label == "e") {
            current_char <- "N"
        } else if (region_label == "m") {
            current_position <- kinetics_summary[i, position]
            current_char <- substr(motif_string, current_position, current_position)
            if (kinetics_summary[i, strand] == "+") {
                # do nothing
            } else if (kinetics_summary[i, strand] == "-") {
                current_char <- complement_bases[[current_char]]
            } else {
                stop("Unexpected strand")
            }
        } else {
            stop("Unexpected region")
        }
        if (genome_stats_target[base == current_char, !(is.finite(N) & N > 0)]) {
            genome_stats_target[base == current_char, c("mean", "sd", "N") := .(NA, NA, 0)]
        }
        kinetics_summary[i, c("char", "genome_mean", "genome_var", "genome_N", "pvalue") := .(current_char,
            genome_stats_target[base == current_char, mean],
            genome_stats_target[base == current_char, sd ^ 2],
            genome_stats_target[base == current_char, N],
            welch.t.pvalue(mean, var, N, genome_stats_target[base == current_char, mean], genome_stats_target[base == current_char, sd ^ 2], genome_stats_target[base == current_char, N]))]
    }
    fwrite(kinetics_summary, file = output_path, append = TRUE, col.names = print_header)
    return(kinetics_summary)
}

motifs_high_ipd <- c("AAAABT", "AGAGAGTA", "AGCTATAT", "AGGCAGGC", "ATGGGAYA", "ATTGTTAC",
"CAAACTAC", "CAGYTG", "CCAATCAG", "CRACGAS",
"DCGAGACC", "DGCTTC",
"GAAGGATC", "GATATRGY", "GCGACCTA", "GCGCGCGC", "GGHGGY", "GTAGATCA", "GTATCGTA",
"TGACGTCA", "TGGTGSA", "YGGAR")

motifs_low_ipd <- c("AATMAATA", "AGACGCAG",
"CGGYTTGA", "CTKCAA",
"GCGCGTCA", "GTGTGTGY",
"TACCCCKA", "TACCTTGA", "TGACGTCA")

motifs_high_ipd <- c(c("ACGCRTG", "ATCAGCTG", "GGN_4", "RGTA"), motifs_high_ipd)
motifs_low_ipd <- c(c("ACATMTGG", "CTGDAR", "TACTGTAG"), motifs_low_ipd)
#motifs_select1 <- c("GAGG", "BGAGG", "AGAA", "ASAAD")
motifs_select1 <- c("GAGG", "AGAA")
motifs_select2 <- "GATC"

motifs_extreme1 <- c("ACGCRTG", "ATCAGCTG", "GGN_4", "AGCTATAT", "CAGYTG", "CRACGAS", "DCGAGACC",
"GAAGGATC", "GATATRGY", "GCGACCTA", "GCGCGCGC", "GGHGGY", "GTAGATCA", "GTATCGTA", "TGACGTCA", "TGGTGSA",
"CGGYTTGA", "GCGCGTCA")

#motif_dirs <- c(c("ACGCRTG", "ATCAGCTG", "GGN_4", "RGTA"), c("ACATMTGG", "CTGDAR", "TACTGTAG"))
motif_dirs <- unique(c(motifs_high_ipd, motifs_low_ipd, motifs_select1, motifs_select2))
#motif_titles <- list("ACATMTGG", "ACGCRTG", "AGGCTT_4" = "(AGGCTT)4", "AGGY_7" = "(AGGY)7", "ATCAGCTG", "ATWTTB", "CGRAACCCG", "CTGDAR", "CTGRAA", "GAGD", "GCGC", "GGN_10" = "(GGN)10", "GGN_4" = "(GGN)4", "GGTCTCGC", "GKCGGY", "GTAS", "HKGAGCA", "HRGTA", "RGTAB", "RTAB", "TACTGTAG", "tandem_repeat_s4_01down", "tandem_repeat_s4_01up", "tandem_repeat_s4_02down", "tandem_repeat_s4_02up", "tandem_repeat_s4_03down", "tandem_repeat_s4_03up", "tandem_repeat_s4_04down", "tandem_repeat_s4_04up", "tandem_repeat_s4_05down", "tandem_repeat_s4_05up", "TATAA", "TTTT", "YGGWR")
motif_titles <- list("AGGCTT_4" = "(AGGCTT)4", "AGGY_7" = "(AGGY)7", "GGN_10" = "(GGN)10", "GGN_4" = "(GGN)4")
motif_strings <- list("AGGCTT_4" = "AGGCTTAGGCTTAGGCTTAGGCTT", "AGGY_7" = "AGGYAGGYAGGYAGGYAGGYAGGYAGGY", "GGN_10" = "GGNGGNGGNGGNGGNGGNGGNGGNGGNGGN", "GGN_4" = "GGNGGNGGNGGN")

motifs_select1 <- c(motifs_select1, motifs_select2)

lookup_motif_title <- function(motif_dir){
    if(is.null(motif_titles[[motif_dir]])){
        return(motif_dir)
    }else{
        return(motif_titles[[motif_dir]])
    }
}

lookup_motif_string <- function(motif_dir){
    if(is.null(motif_strings[[motif_dir]])){
        return(motif_dir)
    }else{
        return(motif_strings[[motif_dir]])
    }
}

kinetics_stats_celegans_degenerate <- fread("../kinetics_stats.c_elegans.degenerate.csv")
kinetics_stats_ecoli_degenerate <- fread("../kinetics_stats.e_coli.degenerate.csv")

output_celegans <- "motif_kinetics_statistical_tests.c_elegans.csv"
invisible(file.remove(output_celegans))
output_ecoli <- "motif_kinetics_statistical_tests.e_coli.csv"
invisible(file.remove(output_ecoli))

outside_length <- 20
print_header_celegans <- TRUE
print_header_ecoli <- TRUE

for (motif_dir in motif_dirs) {
    cat(sprintf("Start processing motif %s\n", motif_dir))
    ab_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.c_elegans.csv"))
    test_motif_kinetics(ab_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep", "VC2010+OP50/WGA/C. elegans", "IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    
    print_header_celegans <- FALSE
    test_motif_kinetics(ab_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep", "VC2010+OP50/WGA/C. elegans", "log2_IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)

    cd_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.c_elegans.csv"))
    test_motif_kinetics(cd_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep", "VC2010/WGA/C. elegans", "IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    test_motif_kinetics(cd_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep", "VC2010/WGA/C. elegans", "log2_IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    
    PD2182_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.c_elegans.csv"))
    test_motif_kinetics(PD2182_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182_deep", "PD2182 (PacBio RS II)/native/C. elegans", "IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    test_motif_kinetics(PD2182_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182_deep", "PD2182 (PacBio RS II)/native/C. elegans", "log2_IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    
    PD2182sequel_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.c_elegans.csv"))
    test_motif_kinetics(PD2182sequel_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182sequel_deep", "PD2182 (PacBio Sequel)/native/C. elegans", "IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    test_motif_kinetics(PD2182sequel_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182sequel_deep", "PD2182 (PacBio Sequel)/native/C. elegans", "log2_IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    
    k_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.c_elegans.csv"))
    test_motif_kinetics(k_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "k_deep", "VC2010+OP50/native/C. elegans", "IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    test_motif_kinetics(k_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "k_deep", "VC2010+OP50/native/C. elegans", "log2_IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    
    l_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.c_elegans.csv"))
    test_motif_kinetics(l_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "l_deep", "VC2010/native/C. elegans", "IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    test_motif_kinetics(l_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "l_deep", "VC2010/native/C. elegans", "log2_IPD", kinetics_stats_celegans_degenerate, output_celegans, print_header_celegans)
    
    ab_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.e_coli.csv"))
    test_motif_kinetics(ab_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep", "VC2010+OP50/WGA/E. coli", "IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)

    print_header_ecoli <- FALSE
    test_motif_kinetics(ab_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep", "VC2010+OP50/WGA/E. coli", "log2_IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)

    cd_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.e_coli.csv"))
    test_motif_kinetics(cd_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep", "VC2010/WGA/E. coli", "IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)
    test_motif_kinetics(cd_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep", "VC2010/WGA/E. coli", "log2_IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)

    k_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.e_coli.csv"))
    test_motif_kinetics(k_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "k_deep", "VC2010+OP50/native/E. coli", "IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)
    test_motif_kinetics(k_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "k_deep", "VC2010+OP50/native/E. coli", "log2_IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)

    l_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.e_coli.csv"))
    test_motif_kinetics(l_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "l_deep", "VC2010/native/E. coli", "IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)
    test_motif_kinetics(l_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "l_deep", "VC2010/native/E. coli", "log2_IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)

    PD2182_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.e_coli.csv"))
    test_motif_kinetics(PD2182_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182_deep", "PD2182 (PacBio RS II)/native/E. coli", "IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)
    test_motif_kinetics(PD2182_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182_deep", "PD2182 (PacBio RS II)/native/E. coli", "log2_IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)

    PD2182sequel_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.e_coli.csv"))
    test_motif_kinetics(PD2182sequel_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182sequel_deep", "PD2182 (PacBio Sequel)/native/E. coli", "IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)
    test_motif_kinetics(PD2182sequel_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182sequel_deep", "PD2182 (PacBio Sequel)/native/E. coli", "log2_IPD", kinetics_stats_ecoli_degenerate, output_ecoli, print_header_ecoli)
}
