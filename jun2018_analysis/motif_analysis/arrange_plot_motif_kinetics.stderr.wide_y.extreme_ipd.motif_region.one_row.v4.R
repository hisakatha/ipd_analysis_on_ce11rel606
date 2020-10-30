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

get_integer_breaks <- function(limits) {
    max_limits <- max(limits)
    if (max_limits <= 10) {
        return(1:max_limits)
    } else if (max_limits <= 30) {
        return(seq(1, max_limits - 1, by = 3))
    } else {
        return(seq(1, max_limits - 1, by = 4))
    }
}

plot_motif_kinetics <- function(kinetics, estimate, title_text, ylab_text, global_mean, global_var, global_size){
    if (kinetics[, .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "No Data"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(list("ipd" = g1, "ipd_and_estimate" = g1))
    }
    kinetics_summary <- kinetics[value > 0 & is.finite(value)][, .(position, strand, label, log2value = log2(value))][, .(mean = mean(log2value), var = var(log2value), .N, type = "observed"), by = .(position, strand, label)]
    #kinetics_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics_summary$strand <- factor(kinetics_summary$strand, levels = c("+", "-"))
    kinetics_summary$region <- factor(kinetics_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics_summary[,"position" := .(position - 20)]
    cat(sprintf("kinetics: max = %.3g\tmin = %.3g\t(%s)\n", kinetics_summary[, max(mean)], kinetics_summary[, min(mean)], title_text))
    estimate_summary <- estimate[value > 0 & is.finite(value)][, .(position, strand, label, log2value = log2(value))][, .(mean = mean(log2value), var = var(log2value), .N, type = "estimated"), by = .(position, strand, label)]
    #estimate_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    estimate_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    estimate_summary$strand <- factor(estimate_summary$strand, levels = c("+", "-"))
    estimate_summary$region <- factor(estimate_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    estimate_summary[,"position" := .(position - 20)]
    cat(sprintf("estimate: max = %.3g\tmin = %.3g\t(%s)\n", estimate_summary[, max(mean)], estimate_summary[, min(mean)], title_text))
    merged_summary <- rbind(kinetics_summary, estimate_summary)
    merged_summary$type <- factor(merged_summary$type, levels = c("observed", "estimated"))
    merged_summary_table <- merge(kinetics_summary, estimate_summary, by = c("position", "strand", "label", "region"), all = TRUE, suffixes = c(".kinetics", ".estimate"))
    mean_cor <- cor(merged_summary_table[region == "Motif", mean.kinetics], merged_summary_table[region == "Motif", mean.estimate], use = "pair")
    if(is.na(mean_cor)){
        mean_cor_p <- NA
    }else{
        mean_cor_p <- cor.test(merged_summary_table[region == "Motif", mean.kinetics], merged_summary_table[region == "Motif", mean.estimate])$p.value
    }
    g1 <- ggplot(kinetics_summary[region == "Motif"], aes(position, mean)) + geom_point(alpha = 0.5) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N))) +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank()) +
        ggtitle(title_text) + ylab(ylab_text) +
        geom_hline(yintercept = global_mean, linetype = "dashed") +
        scale_x_continuous(breaks = get_integer_breaks)
        #geom_text(aes(y = Inf, label = format(pvalue, digits = 3, scientific = T)), show.legend = F, angle = 90, hjust = 1.2, size = 2, color = "black")
        #scale_x_continuous(breaks = kinetics_summary$position, labels = kinetics_summary$label)
        #scale_x_continuous(breaks = 1:kinetics_summary[,max(position)], labels = kinetics_summary$label)
    g2 <- ggplot(merged_summary[region == "Motif"], aes(position, mean, color = type)) + geom_point(alpha = 0.75) +
        geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N)), alpha = 0.75) +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank()) +
        ggtitle(title_text, subtitle = sprintf("r = %.3g (p = %.3g)", mean_cor, mean_cor_p)) + ylab(ylab_text) +
        geom_hline(yintercept = global_mean, linetype = "dashed") +
        scale_x_continuous(breaks = get_integer_breaks) + labs(color = "Data source")
    return(list("ipd" = g1, "ipd_and_estimate" = g2))
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
ab_log2ipd_mean <- -0.188928879061703
ab_log2ipd_var <- 0.497446
ab_log2ipd_size <- 26968723
cd_log2ipd_mean <- -0.19021812891614
cd_log2ipd_var <- 0.533584
cd_log2ipd_size <- 45958624
PD2182_log2ipd_mean <- -0.144740199751773
PD2182_log2ipd_var <- 0.416585
PD2182_log2ipd_size <- 1356981
PD2182sequel_log2ipd_mean <- -0.192455782396496
PD2182sequel_log2ipd_var <- 0.510954
PD2182sequel_log2ipd_size <- 166271875
k_normBy_ab_log2ipdratio_mean <- -0.0083460611408774
k_normBy_ab_log2ipdratio_var <- 0.313728
k_normBy_ab_log2ipdratio_size <- 19537938
l_normBy_cd_log2ipdratio_mean <- -0.00172669273254656
l_normBy_cd_log2ipdratio_var <- 0.292948
l_normBy_cd_log2ipdratio_size <- 35513307

ab_plots <- list()
cd_plots <- list()
PD2182_plots <- list()
PD2182sequel_plots <- list()
ab_plots_with_estimate <- list()
cd_plots_with_estimate <- list()
PD2182_plots_with_estimate <- list()
PD2182sequel_plots_with_estimate <- list()
#k_normBy_ab_plots <- list()
#l_normBy_cd_plots <- list()
legend_plot <- NULL
legend_plot_with_estimate <- NULL
ylim_ipd <- c(-2, 2.3)
#ylim_ipdratio <- c(-2.5, 3)
ylab_expr <- bquote(log[2]~IPD)

for (motif_dir in motif_dirs) {
    ab_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.c_elegans.csv"))
    ab_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.ab.c_elegans.csv"))
    ab_plots_tmp <- plot_motif_kinetics(ab_ipd, ab_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, ab_log2ipd_mean, ab_log2ipd_var, ab_log2ipd_size)
    ab_plots[[motif_dir]] <- ab_plots_tmp[["ipd"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    ab_plots_with_estimate[[motif_dir]] <- ab_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)

    cd_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.c_elegans.csv"))
    cd_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.cd.c_elegans.csv"))
    cd_plots_tmp <- plot_motif_kinetics(cd_ipd, cd_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, cd_log2ipd_mean, cd_log2ipd_var, cd_log2ipd_size)
    cd_plots[[motif_dir]] <- cd_plots_tmp[["ipd"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    cd_plots_with_estimate[[motif_dir]] <- cd_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    
    PD2182_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.c_elegans.csv"))
    PD2182_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182.c_elegans.csv"))
    PD2182_plots_tmp <- plot_motif_kinetics(PD2182_ipd, PD2182_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, PD2182_log2ipd_mean, PD2182_log2ipd_var, PD2182_log2ipd_size)
    PD2182_plots[[motif_dir]] <- PD2182_plots_tmp[["ipd"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    PD2182_plots_with_estimate[[motif_dir]] <- PD2182_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    
    PD2182sequel_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.c_elegans.csv"))
    PD2182sequel_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182sequel.c_elegans.csv"))
    PD2182sequel_plots_tmp <- plot_motif_kinetics(PD2182sequel_ipd, PD2182sequel_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, PD2182sequel_log2ipd_mean, PD2182sequel_log2ipd_var, PD2182sequel_log2ipd_size)
    PD2182sequel_plots[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    PD2182sequel_plots_with_estimate[[motif_dir]] <- PD2182sequel_plots_tmp[["ipd_and_estimate"]] + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    #k_normBy_ab_ipdratio <- fread(paste0(motif_dir, "/motif_ipdratio.k_normBy_ab.csv"))
    #k_normBy_ab_plots[[motif_dir]] <- plot_motif_kinetics(k_normBy_ab_ipdratio, motif_dir, "log2 IPD ratio", k_normBy_ab_log2ipdratio_mean, k_normBy_ab_log2ipdratio_var, k_normBy_ab_log2ipdratio_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipdratio)
    #l_normBy_cd_ipdratio <- fread(paste0(motif_dir, "/motif_ipdratio.l_normBy_cd.csv"))
    #l_normBy_cd_plots[[motif_dir]] <- plot_motif_kinetics(l_normBy_cd_ipdratio, motif_dir, "log2 IPD ratio", l_normBy_cd_log2ipdratio_mean, l_normBy_cd_log2ipdratio_var, l_normBy_cd_log2ipdratio_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipdratio)
    if (is.null(legend_plot_with_estimate)) {
        legend_plot_tmp <- plot_motif_kinetics(ab_ipd, ab_modelPrediction, motif_dir, ylab_expr, ab_log2ipd_mean, ab_log2ipd_var, ab_log2ipd_size)
        legend_plot <- get_legend(legend_plot_tmp[["ipd"]])
        legend_plot_with_estimate <- get_legend(legend_plot_tmp[["ipd_and_estimate"]])
    }
}

ab_plots[["legend"]] <- legend_plot
cd_plots[["legend"]] <- legend_plot
PD2182_plots[["legend"]] <- legend_plot
PD2182sequel_plots[["legend"]] <- legend_plot
#k_normBy_ab_plots[["legend"]] <- legend_plot
#l_normBy_cd_plots[["legend"]] <- legend_plot
ab_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
cd_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
PD2182_plots_with_estimate[["legend"]] <- legend_plot_with_estimate
PD2182sequel_plots_with_estimate[["legend"]] <- legend_plot_with_estimate

pdf_width <- 12
# num = 26
pdf_height_high <- 6 * 2.4
# num = 12
pdf_height_low <- 3 * 2.4
ncol1 <- 5
ylim2 <- c(-2, 3.4)
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.ab.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = ab_plots[motifs_high_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(ab_plots[motifs_high_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.cd.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = cd_plots[motifs_high_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(cd_plots[motifs_high_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = PD2182_plots[motifs_high_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182_plots[motifs_high_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182sequel.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = PD2182sequel_plots[motifs_high_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182sequel_plots[motifs_high_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.ab.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = ab_plots[motifs_low_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(ab_plots[motifs_low_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.cd.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = cd_plots[motifs_low_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(cd_plots[motifs_low_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = PD2182_plots[motifs_low_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182_plots[motifs_low_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182sequel.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = PD2182sequel_plots[motifs_low_ipd], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182sequel_plots[motifs_low_ipd], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.ab.high_with_estimate.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = ab_plots_with_estimate[c(motifs_high_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(ab_plots_with_estimate[c(motifs_high_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.cd.high_with_estimate.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = cd_plots_with_estimate[c(motifs_high_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(cd_plots_with_estimate[c(motifs_high_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182.high_with_estimate.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = PD2182_plots_with_estimate[c(motifs_high_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182_plots_with_estimate[c(motifs_high_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182sequel.high_with_estimate.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = PD2182sequel_plots_with_estimate[c(motifs_high_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182sequel_plots_with_estimate[c(motifs_high_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())

pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.ab.low_with_estimate.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = ab_plots_with_estimate[c(motifs_low_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(ab_plots_with_estimate[c(motifs_low_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.cd.low_with_estimate.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = cd_plots_with_estimate[c(motifs_low_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(cd_plots_with_estimate[c(motifs_low_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182.low_with_estimate.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = PD2182_plots_with_estimate[c(motifs_low_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182_plots_with_estimate[c(motifs_low_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.motif_region.one_row.v4.PD2182sequel.low_with_estimate.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = PD2182sequel_plots_with_estimate[c(motifs_low_ipd, "legend")], align = "none", ncol = ncol1)
plot_grid(plotlist = lapply(PD2182sequel_plots_with_estimate[c(motifs_low_ipd, "legend")], function(p){ p + coord_cartesian(ylim = ylim2) }), align = "none", ncol = ncol1)
invisible(dev.off())
