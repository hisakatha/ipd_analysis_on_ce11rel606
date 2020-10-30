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

plot_motif_kinetics <- function(kinetics, estimate, title_text, motif_string, sample_name, ylab_text, global_mean, global_var, global_size){
    kinetics_occ <- ifelse(kinetics[, .N] == 0, 0, kinetics[strand == "+", max(src)])
    estimate_occ <- ifelse(estimate[, .N] == 0, 0, estimate[strand == "+", max(src)])
    occ_threshold <- 3
    if (kinetics_occ < occ_threshold | estimate_occ < occ_threshold) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "NA"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("ipd" = g1, "ipd_and_estimate" = g1, "ipd_outside20" = g1, "ipd_outside10" = g1,
                    "ipd_and_estimate_outside20" = g1, "ipd_and_estimate_outside10" = g1))
    }
    kinetics <- kinetics[strand == "+"]
    estimate <- estimate[strand == "+"]
    kinetics_summary <- kinetics[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = "observed"), by = .(position, strand, label)]
    #kinetics_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics_summary$strand <- factor(kinetics_summary$strand, levels = c("+", "-"))
    kinetics_summary$region <- factor(kinetics_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics_summary[,"position" := .(position - 20)]
    #cat(sprintf("kinetics: max = %.3g\tmin = %.3g\t(%s)\n", kinetics_summary[, max(mean)], kinetics_summary[, min(mean)], title_text))
    estimate_summary <- estimate[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = "estimated"), by = .(position, strand, label)]
    #estimate_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    estimate_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    estimate_summary$strand <- factor(estimate_summary$strand, levels = c("+", "-"))
    estimate_summary$region <- factor(estimate_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    estimate_summary[,"position" := .(position - 20)]
    if (kinetics_summary[region == "Motif", .N] == 0 | estimate_summary[region == "Motif", .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "NA"), size = 12) +
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
    #ylim_lower <- ifelse(mean_min > (2 ** -1.9), 2 ** -2, 0)
    #ylim_upper <- ifelse(mean_max < (2 ** 2.2), 2 ** 2.3, ifelse(mean_max < (2 ** 3.3), 2 ** 3.4, mean_max + (2 ** 0.2)))
    ylim_lower <- 0
    ylim_upper <- ifelse(mean_max < 2.5, 3, mean_max + 1)
    plot_data <- kinetics_summary[region == "Motif"]
    g1 <- ggplot(plot_data, aes(position, mean)) +
        geom_point(size = 1) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), size = 1) +
        theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "gray"), panel.background = element_rect(fill = NA), plot.subtitle = element_text(size = 6)) +
        theme(axis.title.x = element_blank()) +
        ylab(ylab_text) +
        geom_hline(yintercept = global_mean, linetype = "dashed", size = 0.1) +
        scale_x_continuous(breaks = 1:nchar(motif_string), labels = unlist(strsplit(motif_string, ""))) +
        theme(axis.ticks.x = element_blank()) +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper))
        #geom_text(aes(y = Inf, label = format(pvalue, digits = 3, scientific = T)), show.legend = F, angle = 90, hjust = 1.2, size = 2, color = "black")
        #scale_x_continuous(breaks = kinetics_summary$position, labels = kinetics_summary$label)
        #scale_x_continuous(breaks = 1:kinetics_summary[,max(position)], labels = kinetics_summary$label)
    annot_x_rate <- 0.7
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g1_ret <- g1 + annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\n#occ = %d", title_text, sample_name, kinetics_occ), size = 2)
    # g1_outside20
    end_pos <- kinetics[, max(position)] - 20
    annot_x <- quantile(c(kinetics_summary[, min(position)], kinetics_summary[, max(position)]), annot_x_rate, names = FALSE)
    g1_outside20 <- g1 %+% kinetics_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\n#occ = %d", title_text, sample_name, kinetics_occ), size = 2) +
        theme(axis.text.x = element_text(size = 5))
    # g1_outside10
    plot_data <- kinetics_summary[-9 <= position & position <= end_pos - 10]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g1_outside10 <- g1 %+% plot_data + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\n#occ = %d", title_text, sample_name, kinetics_occ), size = 2)
    # g2
    plot_data <- merged_summary[region == "Motif"]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g2 <- ggplot(plot_data, aes(ifelse(type == "observed", position - 0.15, position + 0.15), mean, color = type)) +
        geom_hline(yintercept = global_mean, linetype = "dashed", size = 0.1) +
        geom_point(alpha = 0.8, size = 0.1) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), alpha = 0.8, size = 0.2) +
        theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "gray"), panel.background = element_rect(fill = NA), plot.subtitle = element_text(size = 6)) +
        theme(legend.position = c(0,1), legend.justification = c("left","top"), legend.background = element_rect(fill = NA, color = "gray")) +
        ylab(ylab_text) +
        scale_x_continuous(breaks = 1:nchar(motif_string), labels = unlist(strsplit(motif_string, ""))) +
        theme(axis.ticks.x = element_blank(), legend.text = element_text(size = 4), legend.title = element_text(size = 5), legend.key.size = unit(0.5, "lines")) +
        theme(axis.title.x = element_blank()) +
        labs(color = "Data source") +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper)) +
        scale_colour_brewer(palette = "Set1", limits = c("observed", "estimated"), labels = c("observed", "estimated"), drop = FALSE)
    g2_ret <- g2 + annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ = %d", title_text, sample_name, mean_cor, mean_cor_p, kinetics_occ), size = 2)
    # g2_outside20
    annot_x <- quantile(c(merged_summary[, min(position)], merged_summary[, max(position)]), annot_x_rate, names = FALSE)
    g2_outside20 <- g2 %+% merged_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ = %d", title_text, sample_name, mean_cor, mean_cor_p, kinetics_occ), size = 2) +
        theme(axis.text.x = element_text(size = 5))
    # g2_outside10
    plot_data <- merged_summary[-9 <= position & position <= end_pos - 10]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g2_outside10 <- g2 %+% plot_data + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ = %d", title_text, sample_name, mean_cor, mean_cor_p, kinetics_occ), size = 2)
    return(list("ipd" = g1_ret, "ipd_and_estimate" = g2_ret, "ipd_outside20" = g1_outside20, "ipd_outside10" = g1_outside10,
                "ipd_and_estimate_outside20" = g2_outside20, "ipd_and_estimate_outside10" = g2_outside10))
}

plot_motif_kinetics_comparison <- function(kinetics1, kinetics2, type1, type2, title_text, motif_string, sample_name, ylab_text, global_mean, global_var, global_size){
    kinetics1_occ <- ifelse(kinetics1[, .N] == 0, 0, kinetics1[strand == "+", max(src)])
    kinetics2_occ <- ifelse(kinetics2[, .N] == 0, 0, kinetics2[strand == "+", max(src)])
    occ_threshold <- 3
    if (kinetics1_occ < occ_threshold & kinetics2_occ < occ_threshold) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "NA"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("comparison" = g1, "comparison_outside20" = g1, "comparison_outside10" = g1))
    }
    if (kinetics1_occ < occ_threshold) {
        kinetics1 <- kinetics1[FALSE]
        kinetics1_occ <- NA
    } else {
        kinetics1 <- kinetics1[strand == "+"]
    }
    if (kinetics2_occ < occ_threshold) {
        kinetics2 <- kinetics2[FALSE]
        kinetics2_occ <- NA
    } else {
        kinetics2 <- kinetics2[strand == "+"]
    }
    kinetics1_summary <- kinetics1[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = type1), by = .(position, strand, label)]
    #kinetics1_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics1_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics1_summary$strand <- factor(kinetics1_summary$strand, levels = c("+", "-"))
    kinetics1_summary$region <- factor(kinetics1_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics1_summary[,"position" := .(position - 20)]
    kinetics2_summary <- kinetics2[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, type = type2), by = .(position, strand, label)]
    #kinetics2_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics2_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics2_summary$strand <- factor(kinetics2_summary$strand, levels = c("+", "-"))
    kinetics2_summary$region <- factor(kinetics2_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics2_summary[,"position" := .(position - 20)]
    if (kinetics1_summary[region == "Motif", .N] == 0 & kinetics2_summary[region == "Motif", .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "NA"), size = 12) +
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
    #ylim_lower <- ifelse(mean_min > (2 ** -1.9), 2 ** -2, 0)
    #ylim_upper <- ifelse(mean_max < (2 ** 2.2), 2 ** 2.3, ifelse(mean_max < (2 ** 3.3), 2 ** 3.4, mean_max + (2 ** 0.2)))
    ylim_lower <- 0
    ylim_upper <- ifelse(mean_max < 2.5, 3, mean_max + 1)
    plot_data <- merged_summary[region == "Motif"]
    annot_x_rate <- 0.7
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g2 <- ggplot(plot_data, aes(ifelse(type == type1, position - 0.15, position + 0.15), mean, color = type)) +
        geom_hline(yintercept = global_mean, linetype = "dashed", size = 0.1) +
        geom_point(alpha = 0.8, size = 0.1) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), alpha = 0.8, size = 0.2) +
        theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "gray"), panel.background = element_rect(fill = NA), plot.subtitle = element_text(size = 6)) +
        theme(legend.position = c(0,1), legend.justification = c("left","top"), legend.background = element_rect(fill = NA, color = "gray"), legend.key.size = unit(0.5, "lines")) +
        ylab(ylab_text) +
        scale_x_continuous(breaks = 1:nchar(motif_string), labels = unlist(strsplit(motif_string, ""))) +
        theme(axis.ticks.x = element_blank(), legend.text = element_text(size = 4), legend.title = element_text(size = 5), legend.key.size = unit(0.5, "lines")) +
        theme(axis.title.x = element_blank()) +
        labs(color = "Data source") +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper)) +
        scale_colour_brewer(palette = "Set1", limits = c(type1, type2), labels = c(type1, type2), drop = FALSE)
        #scale_x_continuous(breaks = get_integer_breaks(plot_data[, position])) +
    g2_ret <- g2 + annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ (%s) = %d,\n#occ (%s) = %d", title_text, sample_name, mean_cor, mean_cor_p, type1, kinetics1_occ, type2, kinetics2_occ), size = 2)
    # g2_outside20
    annot_x <- quantile(c(merged_summary[, min(position)], merged_summary[, max(position)]), annot_x_rate, names = FALSE)
    end_pos <- max(kinetics1[, position], kinetics2[, position]) - 20
    g2_outside20 <- g2 %+% merged_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ (%s) = %d,\n#occ (%s) = %d", title_text, sample_name, mean_cor, mean_cor_p, type1, kinetics1_occ, type2, kinetics2_occ), size = 2) +
        theme(axis.text.x = element_text(size = 5))
    # g2_outside10
    plot_data <- merged_summary[-9 <= position & position <= end_pos - 10]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g2_outside10 <- g2 %+% plot_data + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ (%s) = %d,\n#occ (%s) = %d", title_text, sample_name, mean_cor, mean_cor_p, type1, kinetics1_occ, type2, kinetics2_occ), size = 2)
    return(list("comparison" = g2_ret, "comparison_outside20" = g2_outside20, "comparison_outside10" = g2_outside10))
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
motifs_select3 <- c("ANNNATCAGCTG", "CNNNATCAGCTG", "GNNNATCAGCTG", "TNNNATCAGCTG", "ATCAGCTGATCAGCTG",
"ANNGATC", "CNNGATC", "GNNGATC", "TNNGATC", "GATCGATC")
motifs_select4 <- c("ATGCAT", "TGANNNNNNNNTGCT", "AATT", "ATCGAT", "ATGCAT", "ATTAAT", "CATG", "CTCGAG", "GANTC", "RAATTY", "TTAA", "AGCT", "CTAG")

motifs_extreme1 <- c("ACGCRTG", "ATCAGCTG", "GGN_4", "AGCTATAT", "CAGYTG", "CRACGAS", "DCGAGACC",
"GAAGGATC", "GATATRGY", "GCGACCTA", "GCGCGCGC", "GGHGGY", "GTAGATCA", "GTATCGTA", "TGACGTCA", "TGGTGSA",
"CGGYTTGA", "GCGCGTCA")

#motif_dirs <- c(c("ACGCRTG", "ATCAGCTG", "GGN_4", "RGTA"), c("ACATMTGG", "CTGDAR", "TACTGTAG"))
motif_dirs <- unique(c(motifs_high_ipd, motifs_low_ipd, motifs_select1, motifs_select2, motifs_select3, motifs_select4))

motif_titles <- list("AGGCTT_4" = "(AGGCTT)4", "AGGY_7" = "(AGGY)7", "GGN_10" = "(GGN)10", "GGN_4" = "(GGN)4")
motif_strings <- list("AGGCTT_4" = "AGGCTTAGGCTTAGGCTTAGGCTT", "AGGY_7" = "AGGYAGGYAGGYAGGYAGGYAGGYAGGY", "GGN_10" = "GGNGGNGGNGGNGGNGGNGGNGGNGGNGGN", "GGN_4" = "GGNGGNGGNGGN")

# Assuming ncol is 3, preparing arrangement of select1 and select2
motifs_select1 <- c(motifs_select1, "null", motifs_select2, "null", "null", motifs_select1, motifs_select2, motifs_select3, motifs_select4)

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

# These constants were derived from deep kinetics regions (valid IPD count >= 25) in C. elegans (kinetics_stats.c_elegans.degenerate.csv).
N_ab_deep_IPD_mean <- 1.0016125417229
N_ab_deep_IPD_var <- 0.530361
N_ab_deep_IPD_size <- 24896768
N_cd_deep_IPD_mean <- 1.00697069898825
N_cd_deep_IPD_var <- 0.499876
N_cd_deep_IPD_size <- 43864792
N_k_deep_IPD_mean <- 1.00987841064884
N_k_deep_IPD_var <- 0.457901
N_k_deep_IPD_size <- 132281513
N_l_deep_IPD_mean <- 1.00898633147954
N_l_deep_IPD_var <- 0.453772
N_l_deep_IPD_size <- 138935492
N_PD2182_deep_IPD_mean <- 1.02215998241509
N_PD2182_deep_IPD_var <- 0.551352
N_PD2182_deep_IPD_size <- 1223970
N_PD2182sequel_deep_IPD_mean <- 1.01456759843303
N_PD2182sequel_deep_IPD_var <- 0.622626
N_PD2182sequel_deep_IPD_size <- 165136774

ab_plots <- list()
cd_plots <- list()
PD2182sequel_plots <- list()
k_plots <- list()
l_plots <- list()
ab_plots_with_estimate <- list()
cd_plots_with_estimate <- list()
PD2182sequel_plots_with_estimate <- list()
k_plots_with_estimate <- list()
l_plots_with_estimate <- list()
ab_plots_outside20 <- list()
cd_plots_outside20 <- list()
PD2182sequel_plots_outside20 <- list()
k_plots_outside20 <- list()
l_plots_outside20 <- list()
ab_plots_outside10 <- list()
cd_plots_outside10 <- list()
PD2182sequel_plots_outside10 <- list()
k_plots_outside10 <- list()
l_plots_outside10 <- list()
ab_plots_with_estimate_outside20 <- list()
cd_plots_with_estimate_outside20 <- list()
PD2182sequel_plots_with_estimate_outside20 <- list()
k_plots_with_estimate_outside20 <- list()
l_plots_with_estimate_outside20 <- list()
ab_plots_with_estimate_outside10 <- list()
cd_plots_with_estimate_outside10 <- list()
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
PD2182sequel_comparison_plots <- list()
PD2182sequel_comparison_plots_outside20 <- list()
PD2182sequel_comparison_plots_outside10 <- list()
# Comparison of C. elegans WGA and native
ab_k_comparison_plots <- list()
ab_k_comparison_plots_outside20 <- list()
ab_k_comparison_plots_outside10 <- list()
cd_l_comparison_plots <- list()
cd_l_comparison_plots_outside20 <- list()
cd_l_comparison_plots_outside10 <- list()
null_PD2182sequel_comparison_plots <- list()
null_PD2182sequel_comparison_plots_outside20 <- list()
null_PD2182sequel_comparison_plots_outside10 <- list()

legend_plot <- NULL
legend_plot_with_estimate <- NULL
legend_plot_comparison <- NULL
legend_plot_wga_native_comparison <- NULL
ylab_expr2 <- bquote(log[2]~IPD)
ylab_expr <- "IPD"

#position,strand,value,label,src
null_ipd <- data.table(position = integer(0), strand = character(0), value = numeric(0), label = character(0), src = integer(0))

for (motif_dir in motif_dirs) {
    cat(sprintf("Start processing motif %s\n", motif_dir))
    ab_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.c_elegans.csv"))
    ab_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.ab.c_elegans.csv"))
    ab_plots_tmp <- plot_motif_kinetics(ab_ipd, ab_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010+OP50/WGA\n/C. elegans", ylab_expr, N_ab_deep_IPD_mean, N_ab_deep_IPD_var, N_ab_deep_IPD_size)
    ab_plots[[motif_dir]] <- ab_plots_tmp[["ipd"]]
    ab_plots_with_estimate[[motif_dir]] <- ab_plots_tmp[["ipd_and_estimate"]]
    ab_plots_outside20[[motif_dir]] <- ab_plots_tmp[["ipd_outside20"]]
    ab_plots_outside10[[motif_dir]] <- ab_plots_tmp[["ipd_outside10"]]
    ab_plots_with_estimate_outside20[[motif_dir]] <- ab_plots_tmp[["ipd_and_estimate_outside20"]]
    ab_plots_with_estimate_outside10[[motif_dir]] <- ab_plots_tmp[["ipd_and_estimate_outside10"]]

    if (is.null(legend_plot_with_estimate)) {
        legend_plot <- get_legend(ab_plots_tmp[["ipd"]])
        legend_plot_with_estimate <- get_legend(ab_plots_tmp[["ipd_and_estimate"]])
    }

    cd_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.c_elegans.csv"))
    cd_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.cd.c_elegans.csv"))
    cd_plots_tmp <- plot_motif_kinetics(cd_ipd, cd_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010/WGA\n/C. elegans", ylab_expr, N_cd_deep_IPD_mean, N_cd_deep_IPD_var, N_cd_deep_IPD_size)
    cd_plots[[motif_dir]] <- cd_plots_tmp[["ipd"]]
    cd_plots_with_estimate[[motif_dir]] <- cd_plots_tmp[["ipd_and_estimate"]]
    cd_plots_outside20[[motif_dir]] <- cd_plots_tmp[["ipd_outside20"]]
    cd_plots_outside10[[motif_dir]] <- cd_plots_tmp[["ipd_outside10"]]
    cd_plots_with_estimate_outside20[[motif_dir]] <- cd_plots_tmp[["ipd_and_estimate_outside20"]]
    cd_plots_with_estimate_outside10[[motif_dir]] <- cd_plots_tmp[["ipd_and_estimate_outside10"]]
    
    PD2182sequel_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.c_elegans.csv"))
    PD2182sequel_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182sequel.c_elegans.csv"))
    PD2182sequel_plots_tmp <- plot_motif_kinetics(PD2182sequel_ipd, PD2182sequel_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "PD2182/native\n/C. elegans", ylab_expr, N_PD2182sequel_deep_IPD_mean, N_PD2182sequel_deep_IPD_var, N_PD2182sequel_deep_IPD_size)
    PD2182sequel_plots[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd"]]
    PD2182sequel_plots_with_estimate[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_and_estimate"]]
    PD2182sequel_plots_outside20[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_outside20"]]
    PD2182sequel_plots_outside10[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_outside10"]]
    PD2182sequel_plots_with_estimate_outside20[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_and_estimate_outside20"]]
    PD2182sequel_plots_with_estimate_outside10[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_and_estimate_outside10"]]
    
    k_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.c_elegans.csv"))
    k_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.k.c_elegans.csv"))
    k_plots_tmp <- plot_motif_kinetics(k_ipd, k_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010+OP50/native\n/C. elegans", ylab_expr, N_k_deep_IPD_mean, N_k_deep_IPD_var, N_k_deep_IPD_size)
    k_plots[[motif_dir]] <- k_plots_tmp[["ipd"]]
    k_plots_with_estimate[[motif_dir]] <- k_plots_tmp[["ipd_and_estimate"]]
    k_plots_outside20[[motif_dir]] <- k_plots_tmp[["ipd_outside20"]]
    k_plots_outside10[[motif_dir]] <- k_plots_tmp[["ipd_outside10"]]
    k_plots_with_estimate_outside20[[motif_dir]] <- k_plots_tmp[["ipd_and_estimate_outside20"]]
    k_plots_with_estimate_outside10[[motif_dir]] <- k_plots_tmp[["ipd_and_estimate_outside10"]]
    
    l_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.c_elegans.csv"))
    l_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.l.c_elegans.csv"))
    l_plots_tmp <- plot_motif_kinetics(l_ipd, l_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010/native\n/C. elegans", ylab_expr, N_l_deep_IPD_mean, N_l_deep_IPD_var, N_l_deep_IPD_size)
    l_plots[[motif_dir]] <- l_plots_tmp[["ipd"]]
    l_plots_with_estimate[[motif_dir]] <- l_plots_tmp[["ipd_and_estimate"]]
    l_plots_outside20[[motif_dir]] <- l_plots_tmp[["ipd_outside20"]]
    l_plots_outside10[[motif_dir]] <- l_plots_tmp[["ipd_outside10"]]
    l_plots_with_estimate_outside20[[motif_dir]] <- l_plots_tmp[["ipd_and_estimate_outside20"]]
    l_plots_with_estimate_outside10[[motif_dir]] <- l_plots_tmp[["ipd_and_estimate_outside10"]]
    
    ab_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.e_coli.csv"))
    ab_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.ab.e_coli.csv"))
    ab_ecoli_plots_tmp <- plot_motif_kinetics(ab_ecoli_ipd, ab_ecoli_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010+OP50/WGA\n/E. coli", ylab_expr, N_ab_deep_IPD_mean, N_ab_deep_IPD_var, N_ab_deep_IPD_size)
    ab_ecoli_plots[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd"]]
    ab_ecoli_plots_with_estimate[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_and_estimate"]]
    ab_ecoli_plots_outside20[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_outside20"]]
    ab_ecoli_plots_outside10[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_outside10"]]
    ab_ecoli_plots_with_estimate_outside20[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_and_estimate_outside20"]]
    ab_ecoli_plots_with_estimate_outside10[[motif_dir]] <- ab_ecoli_plots_tmp[["ipd_and_estimate_outside10"]]

    cd_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.e_coli.csv"))
    cd_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.cd.e_coli.csv"))
    cd_ecoli_plots_tmp <- plot_motif_kinetics(cd_ecoli_ipd, cd_ecoli_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010/WGA\n/E. coli", ylab_expr, N_cd_deep_IPD_mean, N_cd_deep_IPD_var, N_cd_deep_IPD_size)
    cd_ecoli_plots[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd"]]
    cd_ecoli_plots_with_estimate[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_and_estimate"]]
    cd_ecoli_plots_outside20[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_outside20"]]
    cd_ecoli_plots_outside10[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_outside10"]]
    cd_ecoli_plots_with_estimate_outside20[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_and_estimate_outside20"]]
    cd_ecoli_plots_with_estimate_outside10[[motif_dir]] <- cd_ecoli_plots_tmp[["ipd_and_estimate_outside10"]]

    k_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.k.e_coli.csv"))
    k_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.k.e_coli.csv"))
    k_ecoli_plots_tmp <- plot_motif_kinetics(k_ecoli_ipd, k_ecoli_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010+OP50/native\n/E. coli", ylab_expr, N_k_deep_IPD_mean, N_k_deep_IPD_var, N_k_deep_IPD_size)
    k_ecoli_plots[[motif_dir]] <- k_ecoli_plots_tmp[["ipd"]]
    k_ecoli_plots_with_estimate[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_and_estimate"]]
    k_ecoli_plots_outside20[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_outside20"]]
    k_ecoli_plots_outside10[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_outside10"]]
    k_ecoli_plots_with_estimate_outside20[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_and_estimate_outside20"]]
    k_ecoli_plots_with_estimate_outside10[[motif_dir]] <- k_ecoli_plots_tmp[["ipd_and_estimate_outside10"]]

    l_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.l.e_coli.csv"))
    l_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.l.e_coli.csv"))
    l_ecoli_plots_tmp <- plot_motif_kinetics(l_ecoli_ipd, l_ecoli_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010/native\n/E. coli", ylab_expr, N_l_deep_IPD_mean, N_l_deep_IPD_var, N_l_deep_IPD_size)
    l_ecoli_plots[[motif_dir]] <- l_ecoli_plots_tmp[["ipd"]]
    l_ecoli_plots_with_estimate[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_and_estimate"]]
    l_ecoli_plots_outside20[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_outside20"]]
    l_ecoli_plots_outside10[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_outside10"]]
    l_ecoli_plots_with_estimate_outside20[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_and_estimate_outside20"]]
    l_ecoli_plots_with_estimate_outside10[[motif_dir]] <- l_ecoli_plots_tmp[["ipd_and_estimate_outside10"]]

    PD2182sequel_ecoli_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.e_coli.csv"))
    PD2182sequel_ecoli_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182sequel.e_coli.csv"))
    PD2182sequel_ecoli_plots_tmp <- plot_motif_kinetics(PD2182sequel_ecoli_ipd, PD2182sequel_ecoli_modelPrediction, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "PD2182/native\n/E. coli", ylab_expr, N_PD2182sequel_deep_IPD_mean, N_PD2182sequel_deep_IPD_var, N_PD2182sequel_deep_IPD_size)
    PD2182sequel_ecoli_plots[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd"]]
    PD2182sequel_ecoli_plots_with_estimate[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_and_estimate"]]
    PD2182sequel_ecoli_plots_outside20[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_outside20"]]
    PD2182sequel_ecoli_plots_outside10[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_outside10"]]
    PD2182sequel_ecoli_plots_with_estimate_outside20[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_and_estimate_outside20"]]
    PD2182sequel_ecoli_plots_with_estimate_outside10[[motif_dir]] <- PD2182sequel_ecoli_plots_tmp[["ipd_and_estimate_outside10"]]

    type1_text <- "C. elegans"
    type2_text <- "E. coli"
    ab_comparison_plots_tmp <- plot_motif_kinetics_comparison(ab_ipd, ab_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010+OP50/WGA", ylab_expr, N_ab_deep_IPD_mean, N_ab_deep_IPD_var, N_ab_deep_IPD_size)
    ab_comparison_plots[[motif_dir]] <- ab_comparison_plots_tmp[["comparison"]]
    ab_comparison_plots_outside20[[motif_dir]] <- ab_comparison_plots_tmp[["comparison_outside20"]]
    ab_comparison_plots_outside10[[motif_dir]] <- ab_comparison_plots_tmp[["comparison_outside10"]]
    cd_comparison_plots_tmp <- plot_motif_kinetics_comparison(cd_ipd, cd_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010/WGA", ylab_expr, N_cd_deep_IPD_mean, N_cd_deep_IPD_var, N_cd_deep_IPD_size)
    cd_comparison_plots[[motif_dir]] <- cd_comparison_plots_tmp[["comparison"]]
    cd_comparison_plots_outside20[[motif_dir]] <- cd_comparison_plots_tmp[["comparison_outside20"]]
    cd_comparison_plots_outside10[[motif_dir]] <- cd_comparison_plots_tmp[["comparison_outside10"]]
    k_comparison_plots_tmp <- plot_motif_kinetics_comparison(k_ipd, k_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010+OP50/native", ylab_expr, N_k_deep_IPD_mean, N_k_deep_IPD_var, N_k_deep_IPD_size)
    k_comparison_plots[[motif_dir]] <- k_comparison_plots_tmp[["comparison"]]
    k_comparison_plots_outside20[[motif_dir]] <- k_comparison_plots_tmp[["comparison_outside20"]]
    k_comparison_plots_outside10[[motif_dir]] <- k_comparison_plots_tmp[["comparison_outside10"]]
    l_comparison_plots_tmp <- plot_motif_kinetics_comparison(l_ipd, l_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010/native", ylab_expr, N_l_deep_IPD_mean, N_l_deep_IPD_var, N_l_deep_IPD_size)
    l_comparison_plots[[motif_dir]] <- l_comparison_plots_tmp[["comparison"]]
    l_comparison_plots_outside20[[motif_dir]] <- l_comparison_plots_tmp[["comparison_outside20"]]
    l_comparison_plots_outside10[[motif_dir]] <- l_comparison_plots_tmp[["comparison_outside10"]]
    PD2182sequel_comparison_plots_tmp <- plot_motif_kinetics_comparison(PD2182sequel_ipd, PD2182sequel_ecoli_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "PD2182/native", ylab_expr, N_PD2182sequel_deep_IPD_mean, N_PD2182sequel_deep_IPD_var, N_PD2182sequel_deep_IPD_size)
    PD2182sequel_comparison_plots[[motif_dir]] <- PD2182sequel_comparison_plots_tmp[["comparison"]]
    PD2182sequel_comparison_plots_outside20[[motif_dir]] <- PD2182sequel_comparison_plots_tmp[["comparison_outside20"]]
    PD2182sequel_comparison_plots_outside10[[motif_dir]] <- PD2182sequel_comparison_plots_tmp[["comparison_outside10"]]
    if (is.null(legend_plot_comparison)) {
        legend_plot_comparison <- get_legend(ab_comparison_plots_tmp[["comparison"]])
    }

    type1_text <- "WGA"
    type2_text <- "native"
    ab_k_comparison_plots_tmp <- plot_motif_kinetics_comparison(ab_ipd, k_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010+OP50/C. elegans", ylab_expr, N_ab_deep_IPD_mean, N_ab_deep_IPD_var, N_ab_deep_IPD_size)
    ab_k_comparison_plots[[motif_dir]] <- ab_k_comparison_plots_tmp[["comparison"]]
    ab_k_comparison_plots_outside20[[motif_dir]] <- ab_k_comparison_plots_tmp[["comparison_outside20"]]
    ab_k_comparison_plots_outside10[[motif_dir]] <- ab_k_comparison_plots_tmp[["comparison_outside10"]]
    cd_l_comparison_plots_tmp <- plot_motif_kinetics_comparison(cd_ipd, l_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "VC2010/C. elegans", ylab_expr, N_cd_deep_IPD_mean, N_cd_deep_IPD_var, N_cd_deep_IPD_size)
    cd_l_comparison_plots[[motif_dir]] <- cd_l_comparison_plots_tmp[["comparison"]]
    cd_l_comparison_plots_outside20[[motif_dir]] <- cd_l_comparison_plots_tmp[["comparison_outside20"]]
    cd_l_comparison_plots_outside10[[motif_dir]] <- cd_l_comparison_plots_tmp[["comparison_outside10"]]
    null_PD2182sequel_comparison_plots_tmp <- plot_motif_kinetics_comparison(null_ipd, PD2182sequel_ipd, type1_text, type2_text, lookup_motif_title(motif_dir), lookup_motif_string(motif_dir), "PD2182/C. elegans", ylab_expr, N_PD2182sequel_deep_IPD_mean, N_PD2182sequel_deep_IPD_var, N_PD2182sequel_deep_IPD_size)
    null_PD2182sequel_comparison_plots[[motif_dir]] <- null_PD2182sequel_comparison_plots_tmp[["comparison"]]
    null_PD2182sequel_comparison_plots_outside20[[motif_dir]] <- null_PD2182sequel_comparison_plots_tmp[["comparison_outside20"]]
    null_PD2182sequel_comparison_plots_outside10[[motif_dir]] <- null_PD2182sequel_comparison_plots_tmp[["comparison_outside10"]]
    if (is.null(legend_plot_wga_native_comparison)) {
        legend_plot_wga_native_comparison <- get_legend(ab_k_comparison_plots_tmp[["comparison"]])
    }
}

ab_plots[["legend"]] <- legend_plot
cd_plots[["legend"]] <- legend_plot
PD2182sequel_plots[["legend"]] <- legend_plot
k_plots[["legend"]] <- legend_plot
l_plots[["legend"]] <- legend_plot
ab_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
cd_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
PD2182sequel_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
k_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
l_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
ab_plots_outside20[["legend"]] <- legend_plot
cd_plots_outside20[["legend"]] <- legend_plot
PD2182sequel_plots_outside20[["legend"]] <- legend_plot
k_plots_outside20[["legend"]] <- legend_plot
l_plots_outside20[["legend"]] <- legend_plot
ab_plots_outside10[["legend"]] <- legend_plot
cd_plots_outside10[["legend"]] <- legend_plot
PD2182sequel_plots_outside10[["legend"]] <- legend_plot
k_plots_outside10[["legend"]] <- legend_plot
l_plots_outside10[["legend"]] <- legend_plot
ab_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
cd_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
PD2182sequel_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
k_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
l_plots_with_estimate_outside20[["legend"]] <- legend_plot_with_estimate
ab_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
cd_plots_with_estimate_outside10[["legend"]] <- legend_plot_with_estimate
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
PD2182sequel_comparison_plots[["legend"]] <- legend_plot_comparison
PD2182sequel_comparison_plots_outside20[["legend"]] <- legend_plot_comparison
PD2182sequel_comparison_plots_outside10[["legend"]] <- legend_plot_comparison
ab_k_comparison_plots[["legend"]] <- legend_plot_wga_native_comparison
ab_k_comparison_plots_outside20[["legend"]] <- legend_plot_wga_native_comparison
ab_k_comparison_plots_outside10[["legend"]] <- legend_plot_wga_native_comparison
cd_l_comparison_plots[["legend"]] <- legend_plot_wga_native_comparison
cd_l_comparison_plots_outside20[["legend"]] <- legend_plot_wga_native_comparison
cd_l_comparison_plots_outside10[["legend"]] <- legend_plot_wga_native_comparison
null_PD2182sequel_comparison_plots[["legend"]] <- legend_plot_wga_native_comparison
null_PD2182sequel_comparison_plots_outside20[["legend"]] <- legend_plot_wga_native_comparison
null_PD2182sequel_comparison_plots_outside10[["legend"]] <- legend_plot_wga_native_comparison

# All sample arrangement (C. elegans)
pdf_width <- 12
n_motifs_high_ipd <- length(motifs_high_ipd)
n_motifs_low_ipd <- length(motifs_low_ipd)
n_motifs_select1 <- length(motifs_select1)
n_motifs_extreme1 <- length(motifs_extreme1)
n_sample <- 5
ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
height_per_motif <- 2.4
pdf_height_high2 <- n_sample * height_per_motif
pdf_height_low2 <- n_sample * height_per_motif
pdf_height_select1 <- n_sample * height_per_motif
pdf_height_extreme1 <- n_sample * height_per_motif

p_text_ab <- ggdraw() + draw_label("VC2010+OP50/WGA", angle = 90)
p_text_cd <- ggdraw() + draw_label("VC2010/WGA", angle = 90)
p_text_PD2182sequel <- ggdraw() + draw_label("PD2182 (PacBio Sequel)", angle = 90)
p_text_k <- ggdraw() + draw_label("VC2010+OP50/native", angle = 90)
p_text_l <- ggdraw() + draw_label("VC2010/native", angle = 90)
rel_widths2 <- c(0.1, rep(1, ncol2))

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    # ab_plots["null"] returns NULL, because there is no such element
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots[p_motifs], list(p_text_cd), cd_plots[p_motifs],
               list(p_text_k), k_plots[p_motifs], list(p_text_l), l_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots[p_motifs], list(p_text_cd), cd_plots[p_motifs],
               list(p_text_k), k_plots[p_motifs], list(p_text_l), l_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_outside20.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside20[p_motifs], list(p_text_cd), cd_plots_outside20[p_motifs],
               list(p_text_k), k_plots_outside20[p_motifs], list(p_text_l), l_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_outside20.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside20[p_motifs], list(p_text_cd), cd_plots_outside20[p_motifs],
               list(p_text_k), k_plots_outside20[p_motifs], list(p_text_l), l_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_outside10.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside10[p_motifs], list(p_text_cd), cd_plots_outside10[p_motifs],
               list(p_text_k), k_plots_outside10[p_motifs], list(p_text_l), l_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_outside10.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside10[p_motifs], list(p_text_cd), cd_plots_outside10[p_motifs],
               list(p_text_k), k_plots_outside10[p_motifs], list(p_text_l), l_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

pdf_width <- 13
ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_with_estimate.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate[p_motifs], list(p_text_cd), cd_plots_with_estimate[p_motifs],
               list(p_text_k), k_plots_with_estimate[p_motifs], list(p_text_l), l_plots_with_estimate[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_with_estimate.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate[p_motifs], list(p_text_cd), cd_plots_with_estimate[p_motifs],
               list(p_text_k), k_plots_with_estimate[p_motifs], list(p_text_l), l_plots_with_estimate[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_with_estimate_outside20.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside20[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside20[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside20[p_motifs], list(p_text_l), l_plots_with_estimate_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_with_estimate_outside20.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside20[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside20[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside20[p_motifs], list(p_text_l), l_plots_with_estimate_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_with_estimate_outside10.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside10[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside10[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside10[p_motifs], list(p_text_l), l_plots_with_estimate_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_with_estimate_outside10.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside10[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside10[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside10[p_motifs], list(p_text_l), l_plots_with_estimate_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

# All sample arrangement (comparison of C. elegans and E. coli)
pdf_width <- 13
ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_ecoli_comparison.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots[p_motifs], list(p_text_cd), cd_comparison_plots[p_motifs],
               list(p_text_k), k_comparison_plots[p_motifs], list(p_text_l), l_comparison_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_ecoli_comparison.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots[p_motifs], list(p_text_cd), cd_comparison_plots[p_motifs],
               list(p_text_k), k_comparison_plots[p_motifs], list(p_text_l), l_comparison_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.select1_ecoli_comparison.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots[p_motifs], list(p_text_cd), cd_comparison_plots[p_motifs],
               list(p_text_k), k_comparison_plots[p_motifs], list(p_text_l), l_comparison_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.extreme1_ecoli_comparison.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots[p_motifs], list(p_text_cd), cd_comparison_plots[p_motifs],
               list(p_text_k), k_comparison_plots[p_motifs], list(p_text_l), l_comparison_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())


ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_ecoli_comparison_outside20.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside20[p_motifs], list(p_text_cd), cd_comparison_plots_outside20[p_motifs],
               list(p_text_k), k_comparison_plots_outside20[p_motifs], list(p_text_l), l_comparison_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_ecoli_comparison_outside20.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside20[p_motifs], list(p_text_cd), cd_comparison_plots_outside20[p_motifs],
               list(p_text_k), k_comparison_plots_outside20[p_motifs], list(p_text_l), l_comparison_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.select1_ecoli_comparison_outside20.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside20[p_motifs], list(p_text_cd), cd_comparison_plots_outside20[p_motifs],
               list(p_text_k), k_comparison_plots_outside20[p_motifs], list(p_text_l), l_comparison_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.extreme1_ecoli_comparison_outside20.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside20[p_motifs], list(p_text_cd), cd_comparison_plots_outside20[p_motifs],
               list(p_text_k), k_comparison_plots_outside20[p_motifs], list(p_text_l), l_comparison_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_ecoli_comparison_outside10.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside10[p_motifs], list(p_text_cd), cd_comparison_plots_outside10[p_motifs],
               list(p_text_k), k_comparison_plots_outside10[p_motifs], list(p_text_l), l_comparison_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_ecoli_comparison_outside10.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside10[p_motifs], list(p_text_cd), cd_comparison_plots_outside10[p_motifs],
               list(p_text_k), k_comparison_plots_outside10[p_motifs], list(p_text_l), l_comparison_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.select1_ecoli_comparison_outside10.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside10[p_motifs], list(p_text_cd), cd_comparison_plots_outside10[p_motifs],
               list(p_text_k), k_comparison_plots_outside10[p_motifs], list(p_text_l), l_comparison_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.extreme1_ecoli_comparison_outside10.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_comparison_plots_outside10[p_motifs], list(p_text_cd), cd_comparison_plots_outside10[p_motifs],
               list(p_text_k), k_comparison_plots_outside10[p_motifs], list(p_text_l), l_comparison_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

# All sample arrangement (C. elegans and E. coli)
pdf_width <- 12
n_motifs_high_ipd <- length(motifs_high_ipd)
n_motifs_low_ipd <- length(motifs_low_ipd)
n_motifs_select1 <- length(motifs_select1)
n_motifs_extreme1 <- length(motifs_extreme1)
n_sample <- 9
ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
height_per_motif <- 2.3
pdf_height_high2 <- n_sample * height_per_motif
pdf_height_low2 <- n_sample * height_per_motif
pdf_height_select1 <- n_sample * height_per_motif
pdf_height_extreme1 <- n_sample * height_per_motif

p_text_ab <- ggdraw() + draw_label("VC2010+OP50/WGA\nC. elegans", angle = 90)
p_text_cd <- ggdraw() + draw_label("VC2010/WGA\nC. elegans", angle = 90)
p_text_PD2182 <- ggdraw() + draw_label("PD2182 (PacBio RS II)\nC. elegans", angle = 90)
p_text_PD2182sequel <- ggdraw() + draw_label("PD2182 (PacBio Sequel)\nC. elegans", angle = 90)
p_text_k <- ggdraw() + draw_label("VC2010+OP50/native\nC. elegans", angle = 90)
p_text_l <- ggdraw() + draw_label("VC2010/native\nC. elegans", angle = 90)
p_text_ab_ecoli <- ggdraw() + draw_label("VC2010+OP50/WGA\nE. coli", angle = 90)
p_text_cd_ecoli <- ggdraw() + draw_label("VC2010/WGA\nE. coli", angle = 90)
p_text_PD2182sequel_ecoli <- ggdraw() + draw_label("PD2182 (PacBio Sequel)\nE. coli", angle = 90)
p_text_k_ecoli <- ggdraw() + draw_label("VC2010+OP50/native\nE. coli", angle = 90)
p_text_l_ecoli <- ggdraw() + draw_label("VC2010/native\nE. coli", angle = 90)
rel_widths2 <- c(0.15, rep(1, ncol2))

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.high.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    # ab_plots["null"] returns NULL, because there is no such element
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots[p_motifs], list(p_text_cd), cd_plots[p_motifs],
               list(p_text_k), k_plots[p_motifs], list(p_text_l), l_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots[p_motifs], list(p_text_l_ecoli), l_ecoli_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.low.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots[p_motifs], list(p_text_cd), cd_plots[p_motifs],
               list(p_text_k), k_plots[p_motifs], list(p_text_l), l_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots[p_motifs], list(p_text_l_ecoli), l_ecoli_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.select1.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots[p_motifs], list(p_text_cd), cd_plots[p_motifs],
               list(p_text_k), k_plots[p_motifs], list(p_text_l), l_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots[p_motifs], list(p_text_l_ecoli), l_ecoli_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.extreme1.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots[p_motifs], list(p_text_cd), cd_plots[p_motifs],
               list(p_text_k), k_plots[p_motifs], list(p_text_l), l_plots[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots[p_motifs], list(p_text_l_ecoli), l_ecoli_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())

ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.high_outside20.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside20[p_motifs], list(p_text_cd), cd_plots_outside20[p_motifs],
               list(p_text_k), k_plots_outside20[p_motifs], list(p_text_l), l_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.low_outside20.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside20[p_motifs], list(p_text_cd), cd_plots_outside20[p_motifs],
               list(p_text_k), k_plots_outside20[p_motifs], list(p_text_l), l_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.select1_outside20.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside20[p_motifs], list(p_text_cd), cd_plots_outside20[p_motifs],
               list(p_text_k), k_plots_outside20[p_motifs], list(p_text_l), l_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.extreme1_outside20.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside20[p_motifs], list(p_text_cd), cd_plots_outside20[p_motifs],
               list(p_text_k), k_plots_outside20[p_motifs], list(p_text_l), l_plots_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.high_outside10.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside10[p_motifs], list(p_text_cd), cd_plots_outside10[p_motifs],
               list(p_text_k), k_plots_outside10[p_motifs], list(p_text_l), l_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.low_outside10.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside10[p_motifs], list(p_text_cd), cd_plots_outside10[p_motifs],
               list(p_text_k), k_plots_outside10[p_motifs], list(p_text_l), l_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.select1_outside10.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside10[p_motifs], list(p_text_cd), cd_plots_outside10[p_motifs],
               list(p_text_k), k_plots_outside10[p_motifs], list(p_text_l), l_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.extreme1_outside10.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_outside10[p_motifs], list(p_text_cd), cd_plots_outside10[p_motifs],
               list(p_text_k), k_plots_outside10[p_motifs], list(p_text_l), l_plots_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())

pdf_width <- 13
ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
rel_widths2 <- c(0.15, rep(1, ncol2))
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.high_with_estimate.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate[p_motifs], list(p_text_cd), cd_plots_with_estimate[p_motifs],
               list(p_text_k), k_plots_with_estimate[p_motifs], list(p_text_l), l_plots_with_estimate[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.low_with_estimate.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate[p_motifs], list(p_text_cd), cd_plots_with_estimate[p_motifs],
               list(p_text_k), k_plots_with_estimate[p_motifs], list(p_text_l), l_plots_with_estimate[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.select1_with_estimate.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate[p_motifs], list(p_text_cd), cd_plots_with_estimate[p_motifs],
               list(p_text_k), k_plots_with_estimate[p_motifs], list(p_text_l), l_plots_with_estimate[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.extreme1_with_estimate.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate[p_motifs], list(p_text_cd), cd_plots_with_estimate[p_motifs],
               list(p_text_k), k_plots_with_estimate[p_motifs], list(p_text_l), l_plots_with_estimate[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())

ncol2 <- 3
rel_widths2 <- c(0.15, rep(1, ncol2))
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.high_with_estimate_outside20.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside20[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside20[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside20[p_motifs], list(p_text_l), l_plots_with_estimate_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.low_with_estimate_outside20.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside20[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside20[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside20[p_motifs], list(p_text_l), l_plots_with_estimate_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.select1_with_estimate_outside20.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside20[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside20[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside20[p_motifs], list(p_text_l), l_plots_with_estimate_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.extreme1_with_estimate_outside20.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside20[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside20[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside20[p_motifs], list(p_text_l), l_plots_with_estimate_outside20[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside20[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside20[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside20[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.high_with_estimate_outside10.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside10[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside10[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside10[p_motifs], list(p_text_l), l_plots_with_estimate_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.low_with_estimate_outside10.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside10[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside10[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside10[p_motifs], list(p_text_l), l_plots_with_estimate_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.select1_with_estimate_outside10.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside10[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside10[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside10[p_motifs], list(p_text_l), l_plots_with_estimate_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all_with_ecoli.extreme1_with_estimate_outside10.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab), ab_plots_with_estimate_outside10[p_motifs], list(p_text_cd), cd_plots_with_estimate_outside10[p_motifs],
               list(p_text_k), k_plots_with_estimate_outside10[p_motifs], list(p_text_l), l_plots_with_estimate_outside10[p_motifs],
               list(p_text_PD2182sequel), PD2182sequel_plots_with_estimate_outside10[p_motifs],
               list(p_text_ab_ecoli), ab_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_cd_ecoli), cd_ecoli_plots_with_estimate_outside10[p_motifs],
               list(p_text_k_ecoli), k_ecoli_plots_with_estimate_outside10[p_motifs], list(p_text_l_ecoli), l_ecoli_plots_with_estimate_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = rel_widths2))
}
invisible(dev.off())

# All sample arrangement (comparison of C. elegans WGA and native)
p_text_ab_k <- ggdraw() + draw_label("VC2010+OP50\nC. elegans", angle = 90)
p_text_cd_l <- ggdraw() + draw_label("VC2010\nC. elegans", angle = 90)
p_text_null_PD2182sequel <- ggdraw() + draw_label("PD2182\nC. elegans", angle = 90)
n_motifs_high_ipd <- length(motifs_high_ipd)
n_motifs_low_ipd <- length(motifs_low_ipd)
n_motifs_select1 <- length(motifs_select1)
n_motifs_extreme1 <- length(motifs_extreme1)
n_sample <- 3
ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
height_per_motif <- 2.3
pdf_height_high2 <- n_sample * height_per_motif
pdf_height_low2 <- n_sample * height_per_motif
pdf_height_select1 <- n_sample * height_per_motif
pdf_height_extreme1 <- n_sample * height_per_motif
pdf_width <- 13
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_wga_native_comparison.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots[p_motifs], list(p_text_cd_l), cd_l_comparison_plots[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_wga_native_comparison.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots[p_motifs], list(p_text_cd_l), cd_l_comparison_plots[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.select1_wga_native_comparison.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots[p_motifs], list(p_text_cd_l), cd_l_comparison_plots[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.extreme1_wga_native_comparison.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots[p_motifs], list(p_text_cd_l), cd_l_comparison_plots[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

ncol2 <- 3
n_page_high_ipd <- ((n_motifs_high_ipd - 1) %/% ncol2) + 1
n_page_low_ipd <- ((n_motifs_low_ipd - 1) %/% ncol2) + 1
n_page_select1 <- ((n_motifs_select1 - 1) %/% ncol2) + 1
n_page_extreme1 <- ((n_motifs_extreme1 - 1) %/% ncol2) + 1
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_wga_native_comparison_outside20.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside20[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside20[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_wga_native_comparison_outside20.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside20[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside20[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.select1_wga_native_comparison_outside20.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside20[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside20[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.extreme1_wga_native_comparison_outside20.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside20[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside20[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside20[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.high_wga_native_comparison_outside10.pdf", width = pdf_width, height = pdf_height_high2)
for(i in 1:n_page_high_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_high_ipd, i * ncol2))
    p_motifs <- c(motifs_high_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside10[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside10[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.low_wga_native_comparison_outside10.pdf", width = pdf_width, height = pdf_height_low2)
for(i in 1:n_page_low_ipd){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_low_ipd, i * ncol2))
    p_motifs <- c(motifs_low_ipd[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside10[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside10[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.select1_wga_native_comparison_outside10.pdf", width = pdf_width, height = pdf_height_select1)
for(i in 1:n_page_select1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_select1, i * ncol2))
    p_motifs <- c(motifs_select1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside10[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside10[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v8_linear_positive_strand.low_occ_filter.all.extreme1_wga_native_comparison_outside10.pdf", width = pdf_width, height = pdf_height_extreme1)
for(i in 1:n_page_extreme1){
    p_range <- (1 + (i - 1) * ncol2):(min(n_motifs_extreme1, i * ncol2))
    p_motifs <- c(motifs_extreme1[p_range], rep("null", ncol2 - length(p_range)))
    plots <- c(list(p_text_ab_k), ab_k_comparison_plots_outside10[p_motifs], list(p_text_cd_l), cd_l_comparison_plots_outside10[p_motifs],
        list(p_text_null_PD2182sequel), null_PD2182sequel_comparison_plots_outside10[p_motifs])
    print(plot_grid(plotlist = plots, align = "none", ncol = ncol2 + 1, rel_widths = c(0.1, rep(1, ncol2))))
}
invisible(dev.off())
