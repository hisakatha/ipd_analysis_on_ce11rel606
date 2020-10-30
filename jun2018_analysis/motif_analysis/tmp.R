library(data.table)
library(ggplot2)
library(cowplot)

welch.t.pvalue <- function(a_mean, a_variance, a_size, b_mean, b_variance, b_size) {
    a_sq <- (a_variance ^ 2) / a_size
    b_sq <- (b_variance ^ 2) / b_size
    t_stat <- (a_mean - b_mean) / sqrt(a_sq + b_sq)
    dof <- ((a_sq + b_sq) ^ 2) / ((a_sq ^ 2) / (a_size - 1) + (b_sq ^ 2) / (b_size - 1))
    # Two-sided p-value
    pt(-abs(t_stat), dof) * 2
}

get_integer_breaks <- function(values) {
    max_values <- max(values)
    min_values <- min(values)
    values_range <- max_values - min_values + 1
    fixed_point <- 1
    if (values_range <= 10) {
        return(min_values:max_values)
    } else if (values_range <= 25) {
        step <- 3
    } else {
        step <- 5
    }
    return(seq(min_values + (fixed_point - min_values) %% step, max_values, by = step))
}

plot_motif_kinetics <- function(kinetics, estimate, title_text, ylab_text, global_mean, global_var, global_size){
    if (kinetics[, .N] == 0 | estimate[, .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "No Data"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("ipd" = g1, "ipd_and_estimate" = g1, "ipd_outside20" = g1, "ipd_outside10" = g1,
                    "ipd_and_estimate_outside20" = g1, "ipd_and_estimate_outside10" = g1))
    }
    kinetics_summary <- kinetics[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = "observed"), by = .(position, strand, label)]
    #kinetics_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics_summary$strand <- factor(kinetics_summary$strand, levels = c("+", "-"))
    kinetics_summary$region <- factor(kinetics_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics_summary[,"position" := .(position - 20)]
    kinetics_occ <- ifelse(kinetics[, .N] == 0, 0, kinetics[, max(src)])
    #cat(sprintf("kinetics: max = %.3g\tmin = %.3g\t(%s)\n", kinetics_summary[, max(mean)], kinetics_summary[, min(mean)], title_text))
    estimate_summary <- estimate[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = "estimated"), by = .(position, strand, label)]
    #estimate_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    estimate_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    estimate_summary$strand <- factor(estimate_summary$strand, levels = c("+", "-"))
    estimate_summary$region <- factor(estimate_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    estimate_summary[,"position" := .(position - 20)]
    estimate_occ <- ifelse(estimate[, .N] == 0, 0, estimate[, max(src)])
    if (kinetics_summary[region == "Motif", .N] == 0 | estimate_summary[region == "Motif", .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "No Data"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("ipd" = g1, "ipd_and_estimate" = g1, "ipd_outside20" = g1, "ipd_outside10" = g1,
                    "ipd_and_estimate_outside20" = g1, "ipd_and_estimate_outside10" = g1))
    }
    #cat(sprintf("estimate: max = %.3g\tmin = %.3g\t(%s)\n", estimate_summary[, max(mean)], estimate_summary[, min(mean)], title_text))
    merged_summary <- rbind(kinetics_summary, estimate_summary)
    merged_summary$type <- factor(merged_summary$type, levels = c("observed", "estimated"))
    if (kinetics_summary[region == "Motif", .N] > 3 & estimate_summary[region == "Motif", .N] > 3) {
        merged_summary_table <- merge(kinetics_summary, estimate_summary, by = c("position", "strand", "label", "region"), all = TRUE, suffixes = c(".kinetics", ".estimate"))
        mean_cor <- cor(merged_summary_table[region == "Motif", mean.kinetics], merged_summary_table[region == "Motif", mean.estimate], use = "pair")
    }else{
        mean_cor <- NA
    }
    if(is.na(mean_cor)){
        mean_cor_p <- NA
    }else{
        mean_cor_p <- cor.test(merged_summary_table[region == "Motif", mean.kinetics], merged_summary_table[region == "Motif", mean.estimate])$p.value
    }
    mean_min <- min(kinetics_summary[region == "Motif", mean], estimate_summary[region == "Motif", mean])
    mean_max <- max(kinetics_summary[region == "Motif", mean], estimate_summary[region == "Motif", mean])
    stopifnot(length(mean_max) == 1)
    ylim_lower <- ifelse(mean_min > (2 ** -1.9), 2 ** -2, 0)
    ylim_upper <- ifelse(mean_max < (2 ** 2.2), 2 ** 2.3, ifelse(mean_max < (2 ** 3.3), 2 ** 3.4, mean_max + (2 ** 0.2)))
    plot_data <- kinetics_summary[region == "Motif"]
    g1 <- ggplot(plot_data, aes(position, mean)) + geom_point(alpha = 0.5) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N))) +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank(), plot.subtitle = element_text(size = 6)) +
        ggtitle(title_text, subtitle = sprintf("#occ. = %d", kinetics_occ)) + ylab(ylab_text) +
        geom_hline(yintercept = global_mean, linetype = "dashed") +
        scale_x_continuous(breaks = get_integer_breaks(plot_data[, position])) +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper))
        #geom_text(aes(y = Inf, label = format(pvalue, digits = 3, scientific = T)), show.legend = F, angle = 90, hjust = 1.2, size = 2, color = "black")
        #scale_x_continuous(breaks = kinetics_summary$position, labels = kinetics_summary$label)
        #scale_x_continuous(breaks = 1:kinetics_summary[,max(position)], labels = kinetics_summary$label)
    end_pos <- kinetics[, max(position)] - 20
    g1_outside20 <- g1 %+% kinetics_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed") + scale_x_continuous(breaks = get_integer_breaks(kinetics_summary[, position]))
    plot_data <- kinetics_summary[-9 <= position & position <= end_pos - 10]
    g1_outside10 <- g1_outside20 %+% plot_data + scale_x_continuous(breaks = get_integer_breaks(plot_data[, position]))
    plot_data <- merged_summary[region == "Motif"]
    g2 <- ggplot(plot_data, aes(position, mean, color = type)) + geom_point(alpha = 0.75) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), alpha = 0.75) +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank(), plot.subtitle = element_text(size = 6)) +
        ggtitle(title_text, subtitle = sprintf("r = %.3g (p = %.3g), #occ (obs) = %d, #occ (est) = %d", mean_cor, mean_cor_p, kinetics_occ, estimate_occ)) +
        ylab(ylab_text) +
        geom_hline(yintercept = global_mean, linetype = "dashed") +
        scale_x_continuous(breaks = get_integer_breaks(plot_data[, position])) + labs(color = "Data source") +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper))
    g2_outside20 <- g2 %+% merged_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed") + scale_x_continuous(breaks = get_integer_breaks(merged_summary[, position]))
    plot_data <- merged_summary[-9 <= position & position <= end_pos - 10]
    g2_outside10 <- g2_outside20 %+% plot_data + scale_x_continuous(breaks = get_integer_breaks(plot_data[, position]))
    return(list("ipd" = g1, "ipd_and_estimate" = g2, "ipd_outside20" = g1_outside20, "ipd_outside10" = g1_outside10,
                "ipd_and_estimate_outside20" = g2_outside20, "ipd_and_estimate_outside10" = g2_outside10))
}

plot_motif_kinetics_comparison <- function(kinetics1, kinetics2, type1, type2, title_text, ylab_text, global_mean, global_var, global_size){
    if (kinetics1[, .N] == 0 & kinetics2[, .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "No Data"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("comparison" = g1, "comparison_outside20" = g1, "comparison_outside10" = g1))
    }
    kinetics1_summary <- kinetics1[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = type1), by = .(position, strand, label)]
    #kinetics1_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics1_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics1_summary$strand <- factor(kinetics1_summary$strand, levels = c("+", "-"))
    kinetics1_summary$region <- factor(kinetics1_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics1_summary[,"position" := .(position - 20)]
    kinetics1_occ <- ifelse(kinetics1[, .N] == 0, 0, kinetics1[, max(src)])
    #kinetics1_summary[,"position" := .(position - 20 - 0.15)]
    kinetics2_summary <- kinetics2[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = type2), by = .(position, strand, label)]
    #kinetics2_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics2_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics2_summary$strand <- factor(kinetics2_summary$strand, levels = c("+", "-"))
    kinetics2_summary$region <- factor(kinetics2_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics2_summary[,"position" := .(position - 20)]
    kinetics2_occ <- ifelse(kinetics2[, .N] == 0, 0, kinetics2[, max(src)])
    #kinetics2_summary[,"position" := .(position - 20 + 0.15)]
    if (kinetics1_summary[region == "Motif", .N] == 0 & kinetics2_summary[region == "Motif", .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "No Data"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("comparison" = g1, "comparison_outside20" = g1, "comparison_outside10" = g1))
    }
    merged_summary <- rbind(kinetics1_summary, kinetics2_summary)
    merged_summary$type <- factor(merged_summary$type, levels = c(type1, type2))
    if (kinetics1_summary[region == "Motif", .N] > 3 & kinetics2_summary[region == "Motif", .N] > 3) {
        merged_summary_table <- merge(kinetics1_summary, kinetics2_summary, by = c("position", "strand", "label", "region"), all = TRUE, suffixes = c(".kinetics1", ".kinetics2"))
        mean_cor <- cor(merged_summary_table[region == "Motif", mean.kinetics1], merged_summary_table[region == "Motif", mean.kinetics2], use = "pair")
    }else{
        mean_cor <- NA
    }
    if(is.na(mean_cor)){
        mean_cor_p <- NA
    }else{
        mean_cor_p <- cor.test(merged_summary_table[region == "Motif", mean.kinetics1], merged_summary_table[region == "Motif", mean.kinetics2])$p.value
    }
    mean_min <- min(kinetics1_summary[region == "Motif", mean], kinetics2_summary[region == "Motif", mean])
    mean_max <- max(kinetics1_summary[region == "Motif", mean], kinetics2_summary[region == "Motif", mean])
    stopifnot(length(mean_max) == 1)
    ylim_lower <- ifelse(mean_min > (2 ** -1.9), 2 ** -2, 0)
    ylim_upper <- ifelse(mean_max < (2 ** 2.2), 2 ** 2.3, ifelse(mean_max < (2 ** 3.3), 2 ** 3.4, mean_max + (2 ** 0.2)))
    plot_data <- merged_summary[region == "Motif"]
    #g2 <- ggplot(plot_data, aes(position, mean, color = type)) +
    g2 <- ggplot(plot_data, aes(ifelse(type == type1, position - 0.15, position + 0.15), mean, color = type)) +
        geom_hline(yintercept = global_mean, linetype = "dashed", size = 0.1) +
        geom_point(alpha = 0.5, size = 0.1) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), alpha = 0.5, size = 0.2) +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "gray"), panel.background = element_rect(fill = NA), plot.subtitle = element_text(size = 6)) +
        ggtitle(title_text, subtitle = sprintf("r = %.3g (p = %.3g), #occ (%s) = %d, #occ (%s) = %d", mean_cor, mean_cor_p, type1, kinetics1_occ, type2, kinetics2_occ)) +
        ylab(ylab_text) + xlab("position") +
        scale_x_continuous(breaks = get_integer_breaks(plot_data[, position])) +
        labs(color = "Data source") +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper))
    end_pos <- max(kinetics1[, position], kinetics2[, position]) - 20
    g2_outside20 <- g2 %+% merged_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) + scale_x_continuous(breaks = get_integer_breaks(merged_summary[, position]))
    plot_data <- merged_summary[-9 <= position & position <= end_pos - 10]
    g2_outside10 <- g2_outside20 %+% plot_data + scale_x_continuous(breaks = get_integer_breaks(plot_data[, position]))
    return(list("comparison" = g2, "comparison_outside20" = g2_outside20, "comparison_outside10" = g2_outside10))
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

#motifs_high_ipd <- c(c("ACGCRTG", "ATCAGCTG", "GGN_4", "RGTA"), motifs_high_ipd)
#motifs_low_ipd <- c(c("ACATMTGG", "CTGDAR", "TACTGTAG"), motifs_low_ipd)
motifs_high_ipd <- c(c("ACGCRTG", "ATCAGCTG", "GGN_4"))
motifs_low_ipd <- c(c("ACATMTGG", "CTGDAR"))

#motif_dirs <- c(c("ACGCRTG", "ATCAGCTG", "GGN_4", "RGTA"), c("ACATMTGG", "CTGDAR", "TACTGTAG"))
motif_dirs <- unique(c(motifs_high_ipd, motifs_low_ipd))
#motif_titles <- list("ACATMTGG", "ACGCRTG", "AGGCTT_4" = "(AGGCTT)4", "AGGY_7" = "(AGGY)7", "ATCAGCTG", "ATWTTB", "CGRAACCCG", "CTGDAR", "CTGRAA", "GAGD", "GCGC", "GGN_10" = "(GGN)10", "GGN_4" = "(GGN)4", "GGTCTCGC", "GKCGGY", "GTAS", "HKGAGCA", "HRGTA", "RGTAB", "RTAB", "TACTGTAG", "tandem_repeat_s4_01down", "tandem_repeat_s4_01up", "tandem_repeat_s4_02down", "tandem_repeat_s4_02up", "tandem_repeat_s4_03down", "tandem_repeat_s4_03up", "tandem_repeat_s4_04down", "tandem_repeat_s4_04up", "tandem_repeat_s4_05down", "tandem_repeat_s4_05up", "TATAA", "TTTT", "YGGWR")
motif_titles <- list("AGGCTT_4" = "(AGGCTT)4", "AGGY_7" = "(AGGY)7", "GGN_10" = "(GGN)10", "GGN_4" = "(GGN)4")

lookup_motif_title <- function(motif_dir){
    if(is.null(motif_titles[[motif_dir]])){
        return(motif_dir)
    }else{
        return(motif_titles[[motif_dir]])
    }
}

# These constants were derived from deep kinetics regions (valid IPD count >= 25) in C. elegans (kinetics_stats.c_elegans.csv).
#ab_log2ipd_mean <- -0.188928879061703
#ab_log2ipd_var <- 0.497446
#ab_log2ipd_size <- 26968723
#cd_log2ipd_mean <- -0.19021812891614
#cd_log2ipd_var <- 0.533584
#cd_log2ipd_size <- 45958624
#PD2182_log2ipd_mean <- -0.144740199751773
#PD2182_log2ipd_var <- 0.416585
#PD2182_log2ipd_size <- 1356981
#PD2182sequel_log2ipd_mean <- -0.192455782396496
#PD2182sequel_log2ipd_var <- 0.510954
#PD2182sequel_log2ipd_size <- 166271875
#k_normBy_ab_log2ipdratio_mean <- -0.0083460611408774
#k_normBy_ab_log2ipdratio_var <- 0.313728
#k_normBy_ab_log2ipdratio_size <- 19537938
#l_normBy_cd_log2ipdratio_mean <- -0.00172669273254656
#l_normBy_cd_log2ipdratio_var <- 0.292948
#l_normBy_cd_log2ipdratio_size <- 35513307
#k_log2ipd_mean <- -0.204852396784493
#k_log2ipd_var <- 0.6066397
#k_log2ipd_size <- 134024867
#l_log2ipd_mean <- -0.201510486629003
#l_log2ipd_var <- 0.59264
#l_log2ipd_size <- 140652495
all_ab_deep_IPD_mean <- 1.00122651085507
all_ab_deep_IPD_var <- 0.525581
all_ab_deep_IPD_size <- 26968723
all_cd_deep_IPD_mean <- 1.00646974138881
all_cd_deep_IPD_var <- 0.49649
all_cd_deep_IPD_size <- 45958624
all_k_deep_IPD_mean <- 1.0097679025642
all_k_deep_IPD_var <- 0.457066
all_k_deep_IPD_size <- 134024867
all_l_deep_IPD_mean <- 1.00886151291594
all_l_deep_IPD_var <- 0.453217
all_l_deep_IPD_size <- 140652495
all_abcd_deep_IPD_mean <- 1.00205937097585
all_abcd_deep_IPD_var <- 0.476666
all_abcd_deep_IPD_size <- 64798351
all_kl_deep_IPD_mean <- 1.00402923278361
all_kl_deep_IPD_var <- 0.403352
all_kl_deep_IPD_size <- 188829575
all_PD2182_deep_IPD_mean <- 1.02115179403736
all_PD2182_deep_IPD_var <- 0.542253
all_PD2182_deep_IPD_size <- 1356981
all_PD2182sequel_deep_IPD_mean <- 1.01455192554966
all_PD2182sequel_deep_IPD_var <- 0.621856
all_PD2182sequel_deep_IPD_size <- 166271875

ab_plots <- list()
cd_plots <- list()
PD2182_plots <- list()
PD2182sequel_plots <- list()
k_plots <- list()
l_plots <- list()
ab_plots_with_estimate <- list()
cd_plots_with_estimate <- list()
PD2182_plots_with_estimate <- list()
PD2182sequel_plots_with_estimate <- list()
k_plots_with_estimate <- list()
l_plots_with_estimate <- list()
ab_plots_outside20 <- list()
cd_plots_outside20 <- list()
PD2182_plots_outside20 <- list()
PD2182sequel_plots_outside20 <- list()
k_plots_outside20 <- list()
l_plots_outside20 <- list()
ab_plots_outside10 <- list()
cd_plots_outside10 <- list()
PD2182_plots_outside10 <- list()
PD2182sequel_plots_outside10 <- list()
k_plots_outside10 <- list()
l_plots_outside10 <- list()
ab_plots_with_estimate_outside20 <- list()
cd_plots_with_estimate_outside20 <- list()
PD2182_plots_with_estimate_outside20 <- list()
PD2182sequel_plots_with_estimate_outside20 <- list()
k_plots_with_estimate_outside20 <- list()
l_plots_with_estimate_outside20 <- list()
ab_plots_with_estimate_outside10 <- list()
cd_plots_with_estimate_outside10 <- list()
PD2182_plots_with_estimate_outside10 <- list()
PD2182sequel_plots_with_estimate_outside10 <- list()
k_plots_with_estimate_outside10 <- list()
l_plots_with_estimate_outside10 <- list()
ab_ecoli_plots <- list()
ab_ecoli_plots_with_estimate <- list()
ab_ecoli_plots_outside20 <- list()
ab_ecoli_plots_outside10 <- list()
ab_ecoli_plots_with_estimate_outside20 <- list()
ab_ecoli_plots_with_estimate_outside10 <- list()
cd_ecoli_plots <- list()
cd_ecoli_plots_with_estimate <- list()
cd_ecoli_plots_outside20 <- list()
cd_ecoli_plots_outside10 <- list()
cd_ecoli_plots_with_estimate_outside20 <- list()
cd_ecoli_plots_with_estimate_outside10 <- list()
k_ecoli_plots <- list()
k_ecoli_plots_with_estimate <- list()
k_ecoli_plots_outside20 <- list()
k_ecoli_plots_outside10 <- list()
k_ecoli_plots_with_estimate_outside20 <- list()
k_ecoli_plots_with_estimate_outside10 <- list()
l_ecoli_plots <- list()
l_ecoli_plots_with_estimate <- list()
l_ecoli_plots_outside20 <- list()
l_ecoli_plots_outside10 <- list()
l_ecoli_plots_with_estimate_outside20 <- list()
l_ecoli_plots_with_estimate_outside10 <- list()
PD2182_ecoli_plots <- list()
PD2182_ecoli_plots_with_estimate <- list()
PD2182_ecoli_plots_outside20 <- list()
PD2182_ecoli_plots_outside10 <- list()
PD2182_ecoli_plots_with_estimate_outside20 <- list()
PD2182_ecoli_plots_with_estimate_outside10 <- list()
PD2182sequel_ecoli_plots <- list()
PD2182sequel_ecoli_plots_with_estimate <- list()
PD2182sequel_ecoli_plots_outside20 <- list()
PD2182sequel_ecoli_plots_outside10 <- list()
PD2182sequel_ecoli_plots_with_estimate_outside20 <- list()
PD2182sequel_ecoli_plots_with_estimate_outside10 <- list()
# Comparison of C. elegans IPDs and E. coli IPDs
ab_comparison_plots <- list()
ab_comparison_plots_outside20 <- list()
ab_comparison_plots_outside10 <- list()
cd_comparison_plots <- list()
cd_comparison_plots_outside20 <- list()
cd_comparison_plots_outside10 <- list()
k_comparison_plots <- list()
k_comparison_plots_outside20 <- list()
k_comparison_plots_outside10 <- list()
l_comparison_plots <- list()
l_comparison_plots_outside20 <- list()
l_comparison_plots_outside10 <- list()
PD2182_comparison_plots <- list()
PD2182_comparison_plots_outside20 <- list()
PD2182_comparison_plots_outside10 <- list()
PD2182sequel_comparison_plots <- list()
PD2182sequel_comparison_plots_outside20 <- list()
PD2182sequel_comparison_plots_outside10 <- list()

legend_plot <- NULL
legend_plot_with_estimate <- NULL
legend_plot_comparison <- NULL
ylab_expr2 <- bquote(log[2]~IPD)
ylab_expr <- "IPD"

for (motif_dir in motif_dirs) {
    ab_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.c_elegans.csv"))
    ab_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.ab.c_elegans.csv"))
    ab_plots_tmp <- plot_motif_kinetics(ab_ipd, ab_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_ab_deep_IPD_mean, all_ab_deep_IPD_var, all_ab_deep_IPD_size)
    ab_plots[[motif_dir]] <- ab_plots_tmp[["ipd"]] + theme(legend.position = "none")
    ab_plots_with_estimate[[motif_dir]] <- ab_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    ab_plots_outside20[[motif_dir]] <- ab_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    ab_plots_outside10[[motif_dir]] <- ab_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    ab_plots_with_estimate_outside20[[motif_dir]] <- ab_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    ab_plots_with_estimate_outside10[[motif_dir]] <- ab_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")

    if (is.null(legend_plot_with_estimate)) {
        legend_plot <- get_legend(ab_plots_tmp[["ipd"]])
        legend_plot_with_estimate <- get_legend(ab_plots_tmp[["ipd_and_estimate"]])
    }

    cd_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.c_elegans.csv"))
    cd_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.cd.c_elegans.csv"))
    cd_plots_tmp <- plot_motif_kinetics(cd_ipd, cd_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_cd_deep_IPD_mean, all_cd_deep_IPD_var, all_cd_deep_IPD_size)
    cd_plots[[motif_dir]] <- cd_plots_tmp[["ipd"]] + theme(legend.position = "none")
    cd_plots_with_estimate[[motif_dir]] <- cd_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    cd_plots_outside20[[motif_dir]] <- cd_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    cd_plots_outside10[[motif_dir]] <- cd_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    cd_plots_with_estimate_outside20[[motif_dir]] <- cd_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    cd_plots_with_estimate_outside10[[motif_dir]] <- cd_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")
    
    PD2182_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.c_elegans.csv"))
    PD2182_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182.c_elegans.csv"))
    PD2182_plots_tmp <- plot_motif_kinetics(PD2182_ipd, PD2182_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_PD2182_deep_IPD_mean, all_PD2182_deep_IPD_var, all_PD2182_deep_IPD_size)
    PD2182_plots[[motif_dir]] <- PD2182_plots_tmp[["ipd"]] + theme(legend.position = "none")
    PD2182_plots_with_estimate[[motif_dir]] <- PD2182_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    PD2182_plots_outside20[[motif_dir]] <- PD2182_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    PD2182_plots_outside10[[motif_dir]] <- PD2182_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    PD2182_plots_with_estimate_outside20[[motif_dir]] <- PD2182_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    PD2182_plots_with_estimate_outside10[[motif_dir]] <- PD2182_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")
    
    PD2182sequel_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.c_elegans.csv"))
    PD2182sequel_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182sequel.c_elegans.csv"))
    PD2182sequel_plots_tmp <- plot_motif_kinetics(PD2182sequel_ipd, PD2182sequel_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_PD2182sequel_deep_IPD_mean, all_PD2182sequel_deep_IPD_var, all_PD2182sequel_deep_IPD_size)
    PD2182sequel_plots[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd"]] + theme(legend.position = "none")
    PD2182sequel_plots_with_estimate[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    PD2182sequel_plots_outside20[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    PD2182sequel_plots_outside10[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    PD2182sequel_plots_with_estimate_outside20[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    PD2182sequel_plots_with_estimate_outside10[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")
    
    k_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.c_elegans.csv"))
    k_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.k.c_elegans.csv"))
    k_plots_tmp <- plot_motif_kinetics(k_ipd, k_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_k_deep_IPD_mean, all_k_deep_IPD_var, all_k_deep_IPD_size)
    k_plots[[motif_dir]] <- k_plots_tmp[["ipd"]] + theme(legend.position = "none")
    k_plots_with_estimate[[motif_dir]] <- k_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    k_plots_outside20[[motif_dir]] <- k_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    k_plots_outside10[[motif_dir]] <- k_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    k_plots_with_estimate_outside20[[motif_dir]] <- k_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    k_plots_with_estimate_outside10[[motif_dir]] <- k_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")
    
    l_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.c_elegans.csv"))
    l_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.l.c_elegans.csv"))
    l_plots_tmp <- plot_motif_kinetics(l_ipd, l_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_l_deep_IPD_mean, all_l_deep_IPD_var, all_l_deep_IPD_size)
    l_plots[[motif_dir]] <- l_plots_tmp[["ipd"]] + theme(legend.position = "none")
    l_plots_with_estimate[[motif_dir]] <- l_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    l_plots_outside20[[motif_dir]] <- l_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    l_plots_outside10[[motif_dir]] <- l_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    l_plots_with_estimate_outside20[[motif_dir]] <- l_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    l_plots_with_estimate_outside10[[motif_dir]] <- l_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")
    
    ab_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.e_coli.csv"))
    ab_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.ab.e_coli.csv"))
    ab_ecoli_plots_tmp <- plot_motif_kinetics(ab_ecoli_ipd, ab_ecoli_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_ab_deep_IPD_mean, all_ab_deep_IPD_var, all_ab_deep_IPD_size)
    ab_ecoli_plots[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd"]] + theme(legend.position = "none")
    ab_ecoli_plots_with_estimate[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    ab_ecoli_plots_outside20[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    ab_ecoli_plots_outside10[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    ab_ecoli_plots_with_estimate_outside20[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    ab_ecoli_plots_with_estimate_outside10[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")

    cd_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.e_coli.csv"))
    cd_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.cd.e_coli.csv"))
    cd_ecoli_plots_tmp <- plot_motif_kinetics(cd_ecoli_ipd, cd_ecoli_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_cd_deep_IPD_mean, all_cd_deep_IPD_var, all_cd_deep_IPD_size)
    cd_ecoli_plots[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd"]] + theme(legend.position = "none")
    cd_ecoli_plots_with_estimate[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    cd_ecoli_plots_outside20[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    cd_ecoli_plots_outside10[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    cd_ecoli_plots_with_estimate_outside20[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    cd_ecoli_plots_with_estimate_outside10[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")

    k_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.e_coli.csv"))
    k_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.k.e_coli.csv"))
    k_ecoli_plots_tmp <- plot_motif_kinetics(k_ecoli_ipd, k_ecoli_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_k_deep_IPD_mean, all_k_deep_IPD_var, all_k_deep_IPD_size)
    k_ecoli_plots[[motif_dir]] <- k_ecoli_plots_tmp[["ipd"]] + theme(legend.position = "none")
    k_ecoli_plots_with_estimate[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    k_ecoli_plots_outside20[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    k_ecoli_plots_outside10[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    k_ecoli_plots_with_estimate_outside20[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    k_ecoli_plots_with_estimate_outside10[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")

    l_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.e_coli.csv"))
    l_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.l.e_coli.csv"))
    l_ecoli_plots_tmp <- plot_motif_kinetics(l_ecoli_ipd, l_ecoli_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_l_deep_IPD_mean, all_l_deep_IPD_var, all_l_deep_IPD_size)
    l_ecoli_plots[[motif_dir]] <- l_ecoli_plots_tmp[["ipd"]] + theme(legend.position = "none")
    l_ecoli_plots_with_estimate[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    l_ecoli_plots_outside20[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    l_ecoli_plots_outside10[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    l_ecoli_plots_with_estimate_outside20[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    l_ecoli_plots_with_estimate_outside10[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")

    PD2182_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.e_coli.csv"))
    PD2182_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182.e_coli.csv"))
    PD2182_ecoli_plots_tmp <- plot_motif_kinetics(PD2182_ecoli_ipd, PD2182_ecoli_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_PD2182_deep_IPD_mean, all_PD2182_deep_IPD_var, all_PD2182_deep_IPD_size)
    PD2182_ecoli_plots[[motif_dir]] <- PD2182_ecoli_plots_tmp[["ipd"]] + theme(legend.position = "none")
    PD2182_ecoli_plots_with_estimate[[motif_dir]] <- PD2182_ecoli_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    PD2182_ecoli_plots_outside20[[motif_dir]] <- PD2182_ecoli_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    PD2182_ecoli_plots_outside10[[motif_dir]] <- PD2182_ecoli_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    PD2182_ecoli_plots_with_estimate_outside20[[motif_dir]] <- PD2182_ecoli_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    PD2182_ecoli_plots_with_estimate_outside10[[motif_dir]] <- PD2182_ecoli_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")

    PD2182sequel_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.e_coli.csv"))
    PD2182sequel_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182sequel.e_coli.csv"))
    PD2182sequel_ecoli_plots_tmp <- plot_motif_kinetics(PD2182sequel_ecoli_ipd, PD2182sequel_ecoli_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_PD2182sequel_deep_IPD_mean, all_PD2182sequel_deep_IPD_var, all_PD2182sequel_deep_IPD_size)
    PD2182sequel_ecoli_plots[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd"]] + theme(legend.position = "none")
    PD2182sequel_ecoli_plots_with_estimate[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none")
    PD2182sequel_ecoli_plots_outside20[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_outside20"]] + theme(legend.position = "none")
    PD2182sequel_ecoli_plots_outside10[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_outside10"]] + theme(legend.position = "none")
    PD2182sequel_ecoli_plots_with_estimate_outside20[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_and_estimate_outside20"]] + theme(legend.position = "none")
    PD2182sequel_ecoli_plots_with_estimate_outside10[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_and_estimate_outside10"]] + theme(legend.position = "none")

    type1_text <- "C. elegans"
    type2_text <- "E. coli"
    ab_comparison_plots_tmp <- plot_motif_kinetics_comparison(ab_ipd, ab_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), ylab_expr, all_ab_deep_IPD_mean, all_ab_deep_IPD_var, all_ab_deep_IPD_size)
    ab_comparison_plots[[motif_dir]] <- ab_comparison_plots_tmp[["comparison"]] + theme(legend.position = "none")
    ab_comparison_plots_outside20[[motif_dir]] <- ab_comparison_plots_tmp[["comparison_outside20"]] + theme(legend.position = "none")
    ab_comparison_plots_outside10[[motif_dir]] <- ab_comparison_plots_tmp[["comparison_outside10"]] + theme(legend.position = "none")
    cd_comparison_plots_tmp <- plot_motif_kinetics_comparison(cd_ipd, cd_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), ylab_expr, all_cd_deep_IPD_mean, all_cd_deep_IPD_var, all_cd_deep_IPD_size)
    cd_comparison_plots[[motif_dir]] <- cd_comparison_plots_tmp[["comparison"]] + theme(legend.position = "none")
    cd_comparison_plots_outside20[[motif_dir]] <- cd_comparison_plots_tmp[["comparison_outside20"]] + theme(legend.position = "none")
    cd_comparison_plots_outside10[[motif_dir]] <- cd_comparison_plots_tmp[["comparison_outside10"]] + theme(legend.position = "none")
    k_comparison_plots_tmp <- plot_motif_kinetics_comparison(k_ipd, k_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), ylab_expr, all_k_deep_IPD_mean, all_k_deep_IPD_var, all_k_deep_IPD_size)
    k_comparison_plots[[motif_dir]] <- k_comparison_plots_tmp[["comparison"]] + theme(legend.position = "none")
    k_comparison_plots_outside20[[motif_dir]] <- k_comparison_plots_tmp[["comparison_outside20"]] + theme(legend.position = "none")
    k_comparison_plots_outside10[[motif_dir]] <- k_comparison_plots_tmp[["comparison_outside10"]] + theme(legend.position = "none")
    l_comparison_plots_tmp <- plot_motif_kinetics_comparison(l_ipd, l_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), ylab_expr, all_l_deep_IPD_mean, all_l_deep_IPD_var, all_l_deep_IPD_size)
    l_comparison_plots[[motif_dir]] <- l_comparison_plots_tmp[["comparison"]] + theme(legend.position = "none")
    l_comparison_plots_outside20[[motif_dir]] <- l_comparison_plots_tmp[["comparison_outside20"]] + theme(legend.position = "none")
    l_comparison_plots_outside10[[motif_dir]] <- l_comparison_plots_tmp[["comparison_outside10"]] + theme(legend.position = "none")
    PD2182_comparison_plots_tmp <- plot_motif_kinetics_comparison(PD2182_ipd, PD2182_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), ylab_expr, all_PD2182_deep_IPD_mean, all_PD2182_deep_IPD_var, all_PD2182_deep_IPD_size)
    PD2182_comparison_plots[[motif_dir]] <- PD2182_comparison_plots_tmp[["comparison"]] + theme(legend.position = "none")
    PD2182_comparison_plots_outside20[[motif_dir]] <- PD2182_comparison_plots_tmp[["comparison_outside20"]] + theme(legend.position = "none")
    PD2182_comparison_plots_outside10[[motif_dir]] <- PD2182_comparison_plots_tmp[["comparison_outside10"]] + theme(legend.position = "none")
    PD2182sequel_comparison_plots_tmp <- plot_motif_kinetics_comparison(PD2182sequel_ipd, PD2182sequel_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), ylab_expr, all_PD2182sequel_deep_IPD_mean, all_PD2182sequel_deep_IPD_var, all_PD2182sequel_deep_IPD_size)
    PD2182sequel_comparison_plots[[motif_dir]] <- PD2182sequel_comparison_plots_tmp[["comparison"]] + theme(legend.position = "none")
    PD2182sequel_comparison_plots_outside20[[motif_dir]] <- PD2182sequel_comparison_plots_tmp[["comparison_outside20"]] + theme(legend.position = "none")
    PD2182sequel_comparison_plots_outside10[[motif_dir]] <- PD2182sequel_comparison_plots_tmp[["comparison_outside10"]] + theme(legend.position = "none")
    if (is.null(legend_plot_comparison)) {
        legend_plot_comparison <- get_legend(ab_comparison_plots_tmp[["comparison"]])
    }
}

ab_plots[["legend"]] <- legend_plot
cd_plots[["legend"]] <- legend_plot
PD2182_plots[["legend"]] <- legend_plot
PD2182sequel_plots[["legend"]] <- legend_plot
k_plots[["legend"]] <- legend_plot
l_plots[["legend"]] <- legend_plot
ab_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
cd_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
PD2182_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
PD2182sequel_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
k_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
l_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
ab_plots_outside20[["legend"]] <- legend_plot
cd_plots_outside20[["legend"]] <- legend_plot
PD2182_plots_outside20[["legend"]] <- legend_plot
PD2182sequel_plots_outside20[["legend"]] <- legend_plot
k_plots_outside20[["legend"]] <- legend_plot
l_plots_outside20[["legend"]] <- legend_plot
ab_plots_outside10[["legend"]] <- legend_plot
cd_plots_outside10[["legend"]] <- legend_plot
PD2182_plots_outside10[["legend"]] <- legend_plot
PD2182sequel_plots_outside10[["legend"]] <- legend_plot
k_plots_outside10[["legend"]] <- legend_plot
l_plots_outside10[["legend"]] <- legend_plot
ab_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
cd_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
PD2182_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
PD2182sequel_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
k_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
l_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
ab_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
cd_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
PD2182_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
PD2182sequel_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
k_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
l_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
ab_ecoli_plots[["legend"]] <- legend_plot
ab_ecoli_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
ab_ecoli_plots_outside20[["legend"]] <- legend_plot
ab_ecoli_plots_outside10[["legend"]] <- legend_plot
ab_ecoli_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
ab_ecoli_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
cd_ecoli_plots[["legend"]] <- legend_plot
cd_ecoli_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
cd_ecoli_plots_outside20[["legend"]] <- legend_plot
cd_ecoli_plots_outside10[["legend"]] <- legend_plot
cd_ecoli_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
cd_ecoli_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
k_ecoli_plots[["legend"]] <- legend_plot
k_ecoli_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
k_ecoli_plots_outside20[["legend"]] <- legend_plot
k_ecoli_plots_outside10[["legend"]] <- legend_plot
k_ecoli_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
k_ecoli_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
l_ecoli_plots[["legend"]] <- legend_plot
l_ecoli_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
l_ecoli_plots_outside20[["legend"]] <- legend_plot
l_ecoli_plots_outside10[["legend"]] <- legend_plot
l_ecoli_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
l_ecoli_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
PD2182_ecoli_plots[["legend"]] <- legend_plot
PD2182_ecoli_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
PD2182_ecoli_plots_outside20[["legend"]] <- legend_plot
PD2182_ecoli_plots_outside10[["legend"]] <- legend_plot
PD2182_ecoli_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
PD2182_ecoli_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
PD2182sequel_ecoli_plots[["legend"]] <- legend_plot
PD2182sequel_ecoli_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
PD2182sequel_ecoli_plots_outside20[["legend"]] <- legend_plot
PD2182sequel_ecoli_plots_outside10[["legend"]] <- legend_plot
PD2182sequel_ecoli_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
PD2182sequel_ecoli_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
ab_comparison_plots[["legend"]] <- legend_plot_comparison
ab_comparison_plots_outside20[["legend"]] <- legend_plot_comparison
ab_comparison_plots_outside10[["legend"]] <- legend_plot_comparison
cd_comparison_plots[["legend"]] <- legend_plot_comparison
cd_comparison_plots_outside20[["legend"]] <- legend_plot_comparison
cd_comparison_plots_outside10[["legend"]] <- legend_plot_comparison
k_comparison_plots[["legend"]] <- legend_plot_comparison
k_comparison_plots_outside20[["legend"]] <- legend_plot_comparison
k_comparison_plots_outside10[["legend"]] <- legend_plot_comparison
l_comparison_plots[["legend"]] <- legend_plot_comparison
l_comparison_plots_outside20[["legend"]] <- legend_plot_comparison
l_comparison_plots_outside10[["legend"]] <- legend_plot_comparison
PD2182_comparison_plots[["legend"]] <- legend_plot_comparison
PD2182_comparison_plots_outside20[["legend"]] <- legend_plot_comparison
PD2182_comparison_plots_outside10[["legend"]] <- legend_plot_comparison
PD2182sequel_comparison_plots[["legend"]] <- legend_plot_comparison
PD2182sequel_comparison_plots_outside20[["legend"]] <- legend_plot_comparison
PD2182sequel_comparison_plots_outside10[["legend"]] <- legend_plot_comparison

height_per_motif <- 2.4
n_motifs_high_ipd <- length(motifs_high_ipd)
n_motifs_low_ipd <- length(motifs_low_ipd)

n_sample <- 6
pdf_height_high2 <- n_sample * height_per_motif
pdf_height_low2 <- n_sample * height_per_motif

p_text_ab <- ggdraw() + draw_label("VC2010+OP50/WGA", angle = 90)
p_text_cd <- ggdraw() + draw_label("VC2010/WGA", angle = 90)
p_text_PD2182 <- ggdraw() + draw_label("PD2182 (PacBio RS II)", angle = 90)
p_text_PD2182sequel <- ggdraw() + draw_label("PD2182 (PacBio Sequel)", angle = 90)
p_text_k <- ggdraw() + draw_label("VC2010+OP50/native", angle = 90)
p_text_l <- ggdraw() + draw_label("VC2010/native", angle = 90)
#rel_widths2 <- c(0.1, rep(1, ncol2))

# All sample arrangement (comparison of C. elegans and E. coli)
pdf_width <- 13
ncol2 <- 5
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
pdf("tmp.v5_linear.all.high_ecoli_comparison.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots[p_motifs], list(NULL), list(p_text_cd), cd_comparison_plots[p_motifs], list(NULL),
               list(p_text_k), k_comparison_plots[p_motifs], list(NULL), list(p_text_l), l_comparison_plots[p_motifs], list(NULL),
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots[p_motifs], list(NULL), list(p_text_PD2182), PD2182_comparison_plots[p_motifs], ab_comparison_plots["legend"])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 2, rel_widths = c(0.1, rep(1, ncol2), 0.5)))
}
invisible(dev.off())
pdf("tmp.v5_linear.all.low_ecoli_comparison.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots[p_motifs], list(NULL), list(p_text_cd), cd_comparison_plots[p_motifs], list(NULL),
               list(p_text_k), k_comparison_plots[p_motifs], list(NULL), list(p_text_l), l_comparison_plots[p_motifs], list(NULL),
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots[p_motifs], list(NULL), list(p_text_PD2182), PD2182_comparison_plots[p_motifs], ab_comparison_plots["legend"])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 2, rel_widths = c(0.1, rep(1, ncol2), 0.5)))
}
invisible(dev.off())

ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
pdf("tmp.v5_linear.all.high_ecoli_comparison_outside20.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside20[p_motifs], list(NULL), list(p_text_cd), cd_comparison_plots_outside20[p_motifs], list(NULL),
               list(p_text_k), k_comparison_plots_outside20[p_motifs], list(NULL), list(p_text_l), l_comparison_plots_outside20[p_motifs], list(NULL),
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside20[p_motifs], list(NULL), list(p_text_PD2182), PD2182_comparison_plots_outside20[p_motifs],
               ab_comparison_plots_outside20["legend"])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 2, rel_widths = c(0.1, rep(1, ncol2), 0.5)))
}
invisible(dev.off())
pdf("tmp.v5_linear.all.low_ecoli_comparison_outside20.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside20[p_motifs], list(NULL), list(p_text_cd), cd_comparison_plots_outside20[p_motifs], list(NULL),
               list(p_text_k), k_comparison_plots_outside20[p_motifs], list(NULL), list(p_text_l), l_comparison_plots_outside20[p_motifs], list(NULL),
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside20[p_motifs], list(NULL), list(p_text_PD2182), PD2182_comparison_plots_outside20[p_motifs],
               ab_comparison_plots_outside20["legend"])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 2, rel_widths = c(0.1, rep(1, ncol2), 0.5)))
}
invisible(dev.off())
pdf("tmp.v5_linear.all.high_ecoli_comparison_outside10.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside10[p_motifs], list(NULL), list(p_text_cd), cd_comparison_plots_outside10[p_motifs], list(NULL),
               list(p_text_k), k_comparison_plots_outside10[p_motifs], list(NULL), list(p_text_l), l_comparison_plots_outside10[p_motifs], list(NULL),
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside10[p_motifs], list(NULL), list(p_text_PD2182), PD2182_comparison_plots_outside10[p_motifs],
               ab_comparison_plots_outside10["legend"])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 2, rel_widths = c(0.1, rep(1, ncol2), 0.5)))
}
invisible(dev.off())
pdf("tmp.v5_linear.all.low_ecoli_comparison_outside10.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside10[p_motifs], list(NULL), list(p_text_cd), cd_comparison_plots_outside10[p_motifs], list(NULL),
               list(p_text_k), k_comparison_plots_outside10[p_motifs], list(NULL), list(p_text_l), l_comparison_plots_outside10[p_motifs], list(NULL),
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside10[p_motifs], list(NULL), list(p_text_PD2182), PD2182_comparison_plots_outside10[p_motifs],
               ab_comparison_plots_outside10["legend"])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 2, rel_widths = c(0.1, rep(1, ncol2), 0.5)))
}
invisible(dev.off())
