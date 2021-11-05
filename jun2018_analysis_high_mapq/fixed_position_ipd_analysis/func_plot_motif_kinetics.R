plot_motif_kinetics <- function(kinetics_summary, estimate_summary, title_text, motif_string, sample_name, ylab_text, global_mean, global_var, global_size){
    kinetics_occ <- ifelse(kinetics_summary[, .N] == 0, 0, kinetics_summary[, max(motif_occ)])
    estimate_occ <- ifelse(estimate_summary[, .N] == 0, 0, estimate_summary[, max(motif_occ)])
    occ_threshold <- 3
    if (kinetics_occ < occ_threshold | estimate_occ < occ_threshold) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "NA"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("ipd" = g1, "ipd_and_estimate" = g1, "ipd_outside20" = g1, "ipd_outside10" = g1,
                    "ipd_and_estimate_outside20" = g1, "ipd_and_estimate_outside10" = g1))
    }
    kinetics_summary <- kinetics_summary[strand == "+"]
    estimate_summary <- estimate_summary[strand == "+"]
    kinetics_summary[, "type" := .("observed")]
    #kinetics_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics_summary$strand <- factor(kinetics_summary$strand, levels = c("+", "-"))
    kinetics_summary$region <- factor(kinetics_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    #cat(sprintf("kinetics: max = %.3g\tmin = %.3g\t(%s)\n", kinetics_summary[, max(mean)], kinetics_summary[, min(mean)], title_text))
    estimate_summary[, "type" := .("estimated")]
    #estimate_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    estimate_summary$strand <- factor(estimate_summary$strand, levels = c("+", "-"))
    estimate_summary$region <- factor(estimate_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
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
    #g1_ret <- g1 + annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\n#occ = %d", title_text, sample_name, kinetics_occ), size = 2)
    g1_ret <- g1
    # g1_outside20
    end_pos <- max(kinetics_summary[, position], estimate_summary[, position])
    annot_x <- quantile(c(kinetics_summary[, min(position)], kinetics_summary[, max(position)]), annot_x_rate, names = FALSE)
    g1_outside20 <- g1 %+% kinetics_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        theme(axis.text.x = element_text(size = 5))
    # g1_outside10
    plot_data <- kinetics_summary[-9 <= position & position <= end_pos - 10]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g1_outside10 <- g1 %+% plot_data + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1)
    # g2
    plot_data <- merged_summary[region == "Motif"]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g2 <- ggplot(plot_data, aes(ifelse(type == "observed", position - 0.15, position + 0.15), mean, color = type)) +
        geom_hline(yintercept = global_mean, linetype = "dashed", size = 0.1) +
        geom_point(alpha = 0.8, size = 0.1) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), alpha = 0.8, size = 0.2) +
        theme(panel.grid = element_blank(), panel.border = element_rect(fill = NA, color = "gray"), panel.background = element_rect(fill = NA), plot.subtitle = element_text(size = 6)) +
        theme(legend.position = c(0,1), legend.justification = c(0,0.85), legend.background = element_blank()) +
        ylab(ylab_text) +
        scale_x_continuous(breaks = 1:nchar(motif_string), labels = unlist(strsplit(motif_string, ""))) +
        theme(axis.ticks.x = element_blank(), legend.text = element_text(size = 9), legend.key.size = unit(0.9, "lines")) +
        theme(axis.title.x = element_blank()) +
        coord_cartesian(ylim = c(ylim_lower, ylim_upper)) +
        scale_colour_brewer(palette = "Set1", limits = c("observed", "estimated"), labels = c("observed", "estimated"), drop = FALSE, name = NULL)
    #g2_ret <- g2 + annotate("text", x = annot_x, y = Inf, hjust = "left", vjust = 1.1, label = sprintf("%s,\n%s,\nr = %.3g (p = %.3g),\n#occ = %d", title_text, sample_name, mean_cor, mean_cor_p, kinetics_occ), size = 2)
    fwrite(data.table(motif_string = title_text, sample_name = sample_name, mean_cor = mean_cor, mean_cor_p = mean_cor_p, kinetics_occ = kinetics_occ), file = correlation_output1, append = TRUE)
    g2_ret <- g2
    # g2_outside20
    annot_x <- quantile(c(merged_summary[, min(position)], merged_summary[, max(position)]), annot_x_rate, names = FALSE)
    g2_outside20 <- g2 %+% merged_summary + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1) +
        theme(axis.text.x = element_text(size = 5))
    # g2_outside10
    plot_data <- merged_summary[-9 <= position & position <= end_pos - 10]
    annot_x <- quantile(c(plot_data[, min(position)], plot_data[, max(position)]), annot_x_rate, names = FALSE)
    g2_outside10 <- g2 %+% plot_data + geom_vline(xintercept = c(0.5, end_pos - 19.5), linetype = "dashed", size = 0.1)
    return(list("ipd" = g1_ret, "ipd_and_estimate" = g2_ret, "ipd_outside20" = g1_outside20, "ipd_outside10" = g1_outside10,
                "ipd_and_estimate_outside20" = g2_outside20, "ipd_and_estimate_outside10" = g2_outside10))
}
