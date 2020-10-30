library(data.table)
library(doParallel)

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

test_motif_kinetics <- function(kinetics, motif_name, motif_string, outside_length, sample_id, sample_name, type1, genome_kinetics_list, output_path, print_header){
    kinetics <- kinetics[value > 0 & is.finite(value)]
    if (kinetics[, .N] == 0) {
        return(data.table())
    }
    if (type1 == "IPD") {
        kinetics_summary <- kinetics[, .(mean = mean(value), var = var(value), .N), by = .(position, strand, label)]
    } else if (type1 == "log2_IPD") {
        kinetics_summary <- kinetics[, "log2_value" := .(log2(value))][, .(mean = mean(log2_value), var = var(log2_value), .N), by = .(position, strand, label)]
    } else {
        stop("Unexpected genome_stats type")
    }
    kinetics_summary[, "position" := .(position - outside_length)]
    kinetics_summary[, c("motif_name", "motif_string", "sample_id", "sample_name", "genome_stats_type") := .(motif_name, motif_string, sample_id, sample_name, type1)]
    stopifnot(kinetics_summary[, .N] > 0)
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
        current_genome_kinetics <- genome_kinetics_list[[current_char]]
        current_kinetics_values <- kinetics[position == current_position + outside_length & strand == kinetics_summary[i, strand], value]
        stopifnot(length(current_kinetics_values) > 0)
        kinetics_summary[i, c("char", "genome_mean", "genome_var", "genome_N", "pvalue") := .(current_char,
            mean(current_genome_kinetics),
            var(current_genome_kinetics),
            length(current_genome_kinetics),
            ifelse(length(current_genome_kinetics) >= 1, wilcox.test(current_kinetics_values, current_genome_kinetics)$p.value, NA))]
    }
    fwrite(kinetics_summary, file = output_path, append = TRUE, col.names = print_header)
    return(kinetics_summary)
}

#motifs_high_ipd <- c("AAAABT", "AGAGAGTA", "AGCTATAT", "AGGCAGGC", "ATGGGAYA", "ATTGTTAC",
#"CAAACTAC", "CAGYTG", "CCAATCAG", "CRACGAS",
#"DCGAGACC", "DGCTTC",
#"GAAGGATC", "GATATRGY", "GCGACCTA", "GCGCGCGC", "GGHGGY", "GTAGATCA", "GTATCGTA",
#"TGACGTCA", "TGGTGSA", "YGGAR")

#motifs_low_ipd <- c("AATMAATA", "AGACGCAG",
#"CGGYTTGA", "CTKCAA",
#"GCGCGTCA", "GTGTGTGY",
#"TACCCCKA", "TACCTTGA", "TGACGTCA")

motif_titles <- list("AGGCTT_4" = "(AGGCTT)4", "AGGY_7" = "(AGGY)7", "GGN_10" = "(GGN)10", "GGN_4" = "(GGN)4")
motif_strings <- list("AGGCTT_4" = "AGGCTTAGGCTTAGGCTTAGGCTT", "AGGY_7" = "AGGYAGGYAGGYAGGYAGGYAGGYAGGY", "GGN_10" = "GGNGGNGGNGGNGGNGGNGGNGGNGGNGGN", "GGN_4" = "GGNGGNGGNGGN")


subset1_motifs <- c("ANNGATC", "CNNGATC", "GNNGATC", "TNNGATC", "GATCGATC",
                    "ANNNATCAGCTG", "CNNNATCAGCTG", "GNNNATCAGCTG", "TNNNATCAGCTG")
subset1_positions <- c(1, 1, 1, 1, 2,
                       1, 1, 1, 1)
subset1_strand <- "+"

stopifnot(length(subset1_motifs) == length(subset1_positions))

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

get_values_per_base <- function(data){
    # data: a data.table with columns: value, base
    setkey(data, base)
    list(
        "A" = data[base == "A", value],
        "C" = data[base == "C", value],
        "G" = data[base == "G", value],
        "T" = data[base == "T", value],
        "W" = data[base == "A" | base == "T", value],
        "S" = data[base == "C" | base == "G", value],
        "M" = data[base == "A" | base == "C", value],
        "K" = data[base == "G" | base == "T", value],
        "R" = data[base == "A" | base == "G", value],
        "Y" = data[base == "C" | base == "T", value],
        "B" = data[base == "C" | base == "G" | base == "T", value],
        "D" = data[base == "A" | base == "G" | base == "T", value],
        "H" = data[base == "A" | base == "C" | base == "T", value],
        "V" = data[base == "A" | base == "C" | base == "G", value],
        "N" = data[, value])
}

source("../load_jun2018_data.saved.fst.ab_cd_k_l.R", chdir = TRUE)
#source("../load_and_save_PD2182_data.saved.fst.PD2182.R", chdir = TRUE)
source("../load_and_save_PD2182sequel_data.saved.fst.PD2182sequel.R", chdir = TRUE)

coverage_thres <- 25
ecoli_chr <- "E._coli_REL606"
ab_data_deep_celegans <- ab_data[coverage >= coverage_thres & refName != ecoli_chr]
cd_data_deep_celegans <- cd_data[coverage >= coverage_thres & refName != ecoli_chr]
k_data_deep_celegans <- k_data[coverage >= coverage_thres & refName != ecoli_chr]
l_data_deep_celegans <- l_data[coverage >= coverage_thres & refName != ecoli_chr]
#abcd_data_deep_celegans <- abcd_data[coverage >= coverage_thres & refName != ecoli_chr]
#kl_data_deep_celegans <- kl_data[coverage >= coverage_thres & refName != ecoli_chr]
#PD2182_data_deep_celegans <- PD2182_data[coverage >= coverage_thres & refName != ecoli_chr]
PD2182sequel_data_deep_celegans <- PD2182sequel_data[coverage >= coverage_thres & refName != ecoli_chr]

ab_data_deep_ecoli <- ab_data[coverage >= coverage_thres & refName == ecoli_chr]
cd_data_deep_ecoli <- cd_data[coverage >= coverage_thres & refName == ecoli_chr]
k_data_deep_ecoli <- k_data[coverage >= coverage_thres & refName == ecoli_chr]
l_data_deep_ecoli <- l_data[coverage >= coverage_thres & refName == ecoli_chr]
#abcd_data_deep_ecoli <- abcd_data[coverage >= coverage_thres & refName == ecoli_chr]
#kl_data_deep_ecoli <- kl_data[coverage >= coverage_thres & refName == ecoli_chr]
#PD2182_data_deep_ecoli <- PD2182_data[coverage >= coverage_thres & refName == ecoli_chr]
PD2182sequel_data_deep_ecoli <- PD2182sequel_data[coverage >= coverage_thres & refName == ecoli_chr]

cat("Start processing get_values_per_base for C. elegans\n")
ab_celegans_ipd_list <- get_values_per_base(ab_data_deep_celegans[, "value" := .(tMean)])
cd_celegans_ipd_list <- get_values_per_base(cd_data_deep_celegans[, "value" := .(tMean)])
k_celegans_ipd_list <- get_values_per_base(k_data_deep_celegans[, "value" := .(tMean)])
l_celegans_ipd_list <- get_values_per_base(l_data_deep_celegans[, "value" := .(tMean)])
#abcd_celegans_ipd_list <- get_values_per_base(abcd_data_deep_celegans[, "value" := .(tMean)])
#kl_celegans_ipd_list <- get_values_per_base(kl_data_deep_celegans[, "value" := .(tMean)])
#PD2182_celegans_ipd_list <- get_values_per_base(PD2182_data_deep_celegans[, "value" := .(tMean)])
PD2182sequel_celegans_ipd_list <- get_values_per_base(PD2182sequel_data_deep_celegans[, "value" := .(tMean)])

cat("Start processing get_values_per_base for E. coli\n")
ab_ecoli_ipd_list <- get_values_per_base(ab_data_deep_ecoli[, "value" := .(tMean)])
cd_ecoli_ipd_list <- get_values_per_base(cd_data_deep_ecoli[, "value" := .(tMean)])
k_ecoli_ipd_list <- get_values_per_base(k_data_deep_ecoli[, "value" := .(tMean)])
l_ecoli_ipd_list <- get_values_per_base(l_data_deep_ecoli[, "value" := .(tMean)])
#abcd_ecoli_ipd_list <- get_values_per_base(abcd_data_deep_ecoli[, "value" := .(tMean)])
#kl_ecoli_ipd_list <- get_values_per_base(kl_data_deep_ecoli[, "value" := .(tMean)])
#PD2182_ecoli_ipd_list <- get_values_per_base(PD2182_data_deep_ecoli[, "value" := .(tMean)])
PD2182sequel_ecoli_ipd_list <- get_values_per_base(PD2182sequel_data_deep_ecoli[, "value" := .(tMean)])


#output_celegans <- "motif_kinetics_statistical_tests.rank_sum.subset1.c_elegans.csv"
#invisible(file.remove(output_celegans))
#output_ecoli <- "motif_kinetics_statistical_tests.rank_sum.subset1.e_coli.csv"
#invisible(file.remove(output_ecoli))

outside_length <- 20
#print_header_celegans <- TRUE
#print_header_ecoli <- TRUE

#registerDoParallel()

foreach_ret <- foreach(i = 1:length(subset1_motifs)) %do% {
    motif_dir <- subset1_motifs[i]
    output <- sprintf("motif_kinetics_statistical_tests.rank_sum.subset2.%s.csv", motif_dir)
    invisible(file.remove(output))
    print_header <- TRUE
    target_position <- subset1_positions[i] + outside_length
    target_strand <- subset1_strand
    cat(sprintf("Start processing motif %s\n", motif_dir))
    ab_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(ab_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep_celegans", "VC2010+OP50/WGA/C. elegans", "IPD", ab_celegans_ipd_list, output, print_header)
    
    print_header <- FALSE
    #print_header_celegans <- FALSE

    cd_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(cd_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep_celegans", "VC2010/WGA/C. elegans", "IPD", cd_celegans_ipd_list, output, print_header)
    
    #PD2182_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.c_elegans.csv"))[position == target_position & strand == target_strand]
    #test_motif_kinetics(PD2182_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182_deep_celegans", "PD2182 (PacBio RS II)/native/C. elegans", "IPD", PD2182_celegans_ipd_list, output, print_header)
    
    PD2182sequel_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(PD2182sequel_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182sequel_deep_celegans", "PD2182 (PacBio Sequel)/native/C. elegans", "IPD", PD2182sequel_celegans_ipd_list, output, print_header)
    
    k_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(k_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "k_deep_celegans", "VC2010+OP50/native/C. elegans", "IPD", k_celegans_ipd_list, output, print_header)
    
    l_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(l_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "l_deep_celegans", "VC2010/native/C. elegans", "IPD", l_celegans_ipd_list, output, print_header)
    
    ab_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(ab_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep_ecoli", "VC2010+OP50/WGA/E. coli", "IPD", ab_ecoli_ipd_list, output, print_header)

    #print_header_ecoli <- FALSE

    cd_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(cd_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep_ecoli", "VC2010/WGA/E. coli", "IPD", cd_ecoli_ipd_list, output, print_header)

    k_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(k_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "k_deep_ecoli", "VC2010+OP50/native/E. coli", "IPD", k_ecoli_ipd_list, output, print_header)

    l_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(l_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "l_deep_ecoli", "VC2010/native/E. coli", "IPD", l_ecoli_ipd_list, output, print_header)

    #PD2182_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.e_coli.csv"))[position == target_position & strand == target_strand]
    #test_motif_kinetics(PD2182_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182_deep_ecoli", "PD2182 (PacBio RS II)/native/E. coli", "IPD", PD2182_ecoli_ipd_list, output, print_header)

    PD2182sequel_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(PD2182sequel_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "PD2182sequel_deep_ecoli", "PD2182 (PacBio Sequel)/native/E. coli", "IPD", PD2182sequel_ecoli_ipd_list, output, print_header)

    fread(output)
}

stopifnot(length(foreach_ret) == length(subset1_motifs))
fwrite(rbindlist(foreach_ret), file = "motif_kinetics_statistical_tests.rank_sum.subset2.csv")
