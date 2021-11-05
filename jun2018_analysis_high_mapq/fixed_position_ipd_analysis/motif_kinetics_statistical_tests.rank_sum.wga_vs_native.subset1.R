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

test_motif_kinetics <- function(kinetics1, kinetics2, motif_name, motif_string, outside_length, sample1_id, sample1_name, sample2_id, sample2_name, type1, output_path, print_header){
    kinetics1 <- kinetics1[value > 0 & is.finite(value)]
    kinetics2 <- kinetics2[value > 0 & is.finite(value)]
    if (kinetics1[, .N] == 0 || kinetics2[, .N] == 0) {
        return(data.table())
    }
    if (type1 == "IPD") {
        kinetics1_summary <- kinetics1[, .(mean = mean(value), var = var(value), .N), by = .(position, strand, label)]
        kinetics2_summary <- kinetics2[, .(mean = mean(value), var = var(value), .N), by = .(position, strand, label)]
    } else if (type1 == "log2_IPD") {
        kinetics1_summary <- kinetics1[, "log2_value" := .(log2(value))][, .(mean = mean(log2_value), var = var(log2_value), .N), by = .(position, strand, label)]
        kinetics2_summary <- kinetics2[, "log2_value" := .(log2(value))][, .(mean = mean(log2_value), var = var(log2_value), .N), by = .(position, strand, label)]
    } else {
        stop("Unexpected genome_stats type")
    }
    kinetics1_summary[, "position" := .(position - outside_length)]
    kinetics2_summary[, "position" := .(position - outside_length)]
    #stopifnot(kinetics1_summary[, .N] > 0)
    #stopifnot(kinetics2_summary[, .N] > 0)
    test_summary <- merge(kinetics1_summary, kinetics2_summary, by = c("position", "strand", "label"), all = TRUE, suffixes = c("_sample1", "_sample2"))
    test_summary[, c("motif_name", "motif_string", "id_sample1", "name_sample1", "id_sample2", "name_sample2", "stats_type") :=
                 .(motif_name, motif_string, sample1_id, sample1_name, sample2_id, sample2_name, type1)]
    for (i in 1:nrow(test_summary)) {
        region_label <- substr(test_summary[i, label], 1, 1)
        current_position <- test_summary[i, position]
        current_strand <- test_summary[i, strand]
        if (region_label == "s" || region_label == "e") {
            current_char <- "N"
        } else if (region_label == "m") {
            current_char <- substr(motif_string, current_position, current_position)
            if (current_strand == "+") {
                # do nothing
            } else if (current_strand == "-") {
                current_char <- complement_bases[[current_char]]
            } else {
                stop("Unexpected strand")
            }
        } else {
            stop("Unexpected region")
        }
        #current_genome_kinetics <- genome_kinetics_list[[current_char]]
        current_kinetics1_values <- kinetics1[position == current_position + outside_length & strand == current_strand, value]
        current_kinetics2_values <- kinetics2[position == current_position + outside_length & strand == current_strand, value]
        #stopifnot(length(current_kinetics1_values) > 0)
        test_summary[i, c("char", "pvalue") :=
            .(char = current_char,
            pvalue = ifelse(length(current_kinetics1_values) >= 1 && length(current_kinetics2_values) >= 1, wilcox.test(current_kinetics1_values, current_kinetics2_values)$p.value, NA))]
    }
    fwrite(test_summary, file = output_path, append = TRUE, col.names = print_header)
    return(test_summary)
}

motif_titles <- list("greer_pacbio_m6A_cov50" = "Greer et al. P6-C4 data 6mA position", "greer_pacbio_orig_m6A_cov50" = "Greer et al. mixed data 6mA position")
motif_strings <- list("greer_pacbio_m6A_cov50" = "A", "greer_pacbio_orig_m6A_cov50" = "A")

subset1_motifs <- c("greer_pacbio_m6A_cov50", "greer_pacbio_orig_m6A_cov50")
subset1_positions <- c(1, 1)
subset1_strand <- c("+", "+")

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

outside_length <- 20

#registerDoParallel()

foreach_ret <- foreach(i = 1:length(subset1_motifs)) %do% {
    motif_dir <- subset1_motifs[i]
    target_position <- subset1_positions[i] + outside_length
    target_strand <- subset1_strand[i]
    output <- sprintf("motif_kinetics_statistical_tests.rank_sum.wga_vs_native.subset1.%s_%d_%d.csv", motif_dir, target_position, ifelse(target_strand == "+", 0, 1))
    invisible(file.remove(output))
    print_header <- TRUE
    cat(sprintf("Start processing motif %s\n", motif_dir))
    ab_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.c_elegans.csv"))[position == target_position & strand == target_strand]
    k_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(ab_ipd, k_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep_celegans", "VC2010+OP50/WGA/C. elegans", "k_deep_celegans", "VC2010+OP50/native/C. elegans", "IPD", output, print_header)

    print_header <- FALSE

    ab_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.e_coli.csv"))[position == target_position & strand == target_strand]
    k_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(ab_ecoli_ipd, k_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "ab_deep_ecoli", "VC2010+OP50/WGA/E. coli", "k_deep_ecoli", "VC2010+OP50/native/E. coli", "IPD", output, print_header)

    cd_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.c_elegans.csv"))[position == target_position & strand == target_strand]
    l_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(cd_ipd, l_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep_celegans", "VC2010/WGA/C. elegans", "l_deep_celegans", "VC2010/native/C. elegans", "IPD", output, print_header)

    cd_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.e_coli.csv"))[position == target_position & strand == target_strand]
    l_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(cd_ecoli_ipd, l_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "cd_deep_ecoli", "VC2010/WGA/E. coli", "l_deep_ecoli", "VC2010/native/E. coli", "IPD", output, print_header)

    abcd_ipd <- fread(paste0(motif_dir, "/motif_ipd.abcd.c_elegans.csv"))[position == target_position & strand == target_strand]
    kl_ipd <- fread(paste0(motif_dir, "/motif_ipd.kl.c_elegans.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(abcd_ipd, kl_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "abcd_deep_celegans", "Merged/WGA/C. elegans", "kl_deep_celegans", "Merged/native/C. elegans", "IPD", output, print_header)

    abcd_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.abcd.e_coli.csv"))[position == target_position & strand == target_strand]
    kl_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.kl.e_coli.csv"))[position == target_position & strand == target_strand]
    test_motif_kinetics(abcd_ecoli_ipd, kl_ecoli_ipd, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), outside_length, "abcd_deep_ecoli", "Merged/WGA/E. coli", "kl_deep_ecoli", "Merged/native/E. coli", "IPD", output, print_header)
    fread(output)
}

stopifnot(length(foreach_ret) == length(subset1_motifs))
fwrite(rbindlist(foreach_ret), file = "motif_kinetics_statistical_tests.rank_sum.wga_vs_native.subset1.csv")
