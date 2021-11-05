plot_motif_kinetics_comparison <- function(kinetics1_summary, kinetics2_summary, type1, type2, title_text, motif_string, sample_name, ylab_text, global_mean, global_var, global_size){
    kinetics1_occ <- ifelse(kinetics1_summary[, .N] == 0, 0, kinetics1_summary[, max(motif_occ)])
    kinetics2_occ <- ifelse(kinetics2_summary[, .N] == 0, 0, kinetics2_summary[, max(motif_occ)])
    occ_threshold <- 3
    if (kinetics1_occ < occ_threshold & kinetics2_occ < occ_threshold) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "NA"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("comparison" = g1, "comparison_outside20" = g1, "comparison_outside10" = g1))
    }
    if (kinetics1_occ < occ_threshold) {
        kinetics1_summary <- kinetics1_summary[FALSE]
        kinetics1_occ <- NA
    } else {
        kinetics1_summary <- kinetics1_summary[strand == "+"]
    }
    if (kinetics2_occ < occ_threshold) {
        kinetics2_summary <- kinetics2_summary[FALSE]
        kinetics2_occ <- NA
    } else {
        kinetics2_summary <- kinetics2_summary[strand == "+"]
    }
    kinetics1_summary[, "type" := .(type1)]
    #kinetics1_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics1_summary$strand <- factor(kinetics1_summary$strand, levels = c("+", "-"))
    kinetics1_summary$region <- factor(kinetics1_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics2_summary[, "type" := .(type2)]
    #kinetics2_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics2_summary$strand <- factor(kinetics2_summary$strand, levels = c("+", "-"))
    kinetics2_summary$region <- factor(kinetics2_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
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
    legend_face <- ifelse(type1 == "C. elegans", "italic", "plain")
    g2 <- ggplot(plot_data, aes(ifelse(type == type1, position - 0.15, position + 0.15), mean, color = type)) +
        geom_hline(yintercept = global_mean, linetype = "dashed", size = 0.1) +
        geom_point(alpha = 0.8, size = 0.1) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), alpha = 0.8, size = 0.2) +
        theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "gray"), panel.background = element_rect(fill = NA), plot.subtitle = element_text(size = 6)) +
        theme(legend.position = c(0,1), legend.justification = c(0,0.85), legend.background = element_blank()) +
        ylab(ylab_text) +
        scale_x_continuous(breaks = 1:nchar(motif_string), labels = unlist(strsplit(motif_string, ""))) +
        theme(axis.ticks.x = element_blank(), legend.text = element_text(size = 9, face = legend_face), legend.key.size = unit(0.9, "lines")) +
        theme(axis.title.x = element_blank()) +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper)) +
        scale_colour_brewer(palette = "Set1", limits = c(type1, type2), labels = c(type1, type2), drop = FALSE, name = NULL)
        #scale_x_continuous(breaks = get_integer_breaks(plot_data[, position])) +
    #g2_ret <- g2 + annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ (%s) = %d,\n#occ (%s) = %d", title_text, sample_name, mean_cor, mean_cor_p, type1, kinetics1_occ, type2, kinetics2_occ), size = 2)
    fwrite(data.table(motif_string = title_text, sample_name = sample_name, mean_cor = mean_cor, mean_cor_p = mean_cor_p, type1 = type1, kinetics1_occ = kinetics1_occ, type2 = type2, kinetics2_occ = kinetics2_occ), file = correlation_output2, append = TRUE)
    g2_ret <- g2
    # g2_outside20
    annot_x <- quantile(c(merged_summary[, min(position)], merged_summary[, max(position)]), annot_x_rate, names = FALSE)
    end_pos <- max(kinetics1_summary[, position], kinetics2_summary[, position])
    g2_outside20 <- g2 %+% merged_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        theme(axis.text.x = element_text(size = 5))
    # g2_outside10
    plot_data <- merged_summary[-9 <= position & position <= end_pos - 10]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g2_outside10 <- g2 %+% plot_data + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1)
    return(list("comparison" = g2_ret, "comparison_outside20" = g2_outside20, "comparison_outside10" = g2_outside10))
}
