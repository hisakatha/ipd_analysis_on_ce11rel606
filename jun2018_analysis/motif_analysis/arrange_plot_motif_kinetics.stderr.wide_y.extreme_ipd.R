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

plot_motif_kinetics <- function(kinetics, title_text, ylab_text, global_mean, global_var, global_size){
    if (kinetics[, .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "No Data"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank())
        return(g1)
    }
    kinetics_summary <- kinetics[value > 0 & is.finite(value)][, .(position, strand, label, log2value = log2(value))][, .(mean = mean(log2value), var = var(log2value), .N), by = .(position, strand, label)]
    #kinetics_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics_summary$strand <- factor(kinetics_summary$strand, levels = c("+", "-"))
    kinetics_summary$region <- factor(kinetics_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    g1 <- ggplot(kinetics_summary, aes(position, mean, color = region)) + geom_point(alpha = 0.5) + geom_errorbar(aes(ymin = mean - 2 * sqrt(var / N), ymax = mean + 2 * sqrt(var / N))) +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank()) +
        ggtitle(title_text) + ylab(ylab_text) +
        geom_hline(yintercept = global_mean, linetype = "dashed") +
        scale_color_discrete(breaks = c("Upstream", "Motif", "Downstream"), drop = FALSE)
        #geom_text(aes(y = Inf, label = format(pvalue, digits = 3, scientific = T)), show.legend = F, angle = 90, hjust = 1.2, size = 2, color = "black")
        #scale_x_continuous(breaks = kinetics_summary$position, labels = kinetics_summary$label)
        #scale_x_continuous(breaks = 1:kinetics_summary[,max(position)], labels = kinetics_summary$label)
    return(g1)
}

motif_dirs <- c("ACATMTGG", "ACGCRTG", "AGGCTT_4", "AGGY_7", "ATCAGCTG", "ATWTTB", "CGRAACCCG", "CTGDAR", "CTGRAA", "GAGD", "GCGC", "GGN_10", "GGN_4", "GGTCTCGC", "GKCGGY", "GTAS", "HKGAGCA", "HRGTA", "RGTAB", "RTAB", "TACTGTAG", "tandem_repeat_s4_01down", "tandem_repeat_s4_01up", "tandem_repeat_s4_02down", "tandem_repeat_s4_02up", "tandem_repeat_s4_03down", "tandem_repeat_s4_03up", "tandem_repeat_s4_04down", "tandem_repeat_s4_04up", "tandem_repeat_s4_05down", "tandem_repeat_s4_05up", "TATAA", "TTTT", "YGGWR")

#ab_ipd <- fread("motif_ipd.ab.csv")
#cd_ipd <- fread("motif_ipd.cd.csv")
#PD2182_ipd <- fread("motif_ipd.PD2182.csv")
#PD2182sequel_ipd <- fread("motif_ipd.PD2182sequel.csv")
#k_normBy_ab_ipdratio <- fread("motif_ipdratio.k_normBy_ab.csv")
#l_normBy_cd_ipdratio <- fread("motif_ipdratio.l_normBy_cd.csv")

# These constants were derived from deep regions (coverage >= 25).
ab_log2ipd_mean <- -0.1889734
ab_log2ipd_var <- 0.4976316
ab_log2ipd_size <- 27095150
cd_log2ipd_mean <- -0.1902569
cd_log2ipd_var <- 0.5337555
cd_log2ipd_size <- 46121818
PD2182_log2ipd_mean <- -0.144740199751773
PD2182_log2ipd_var <- 0.4165849
PD2182_log2ipd_size <- 1356981
PD2182sequel_log2ipd_mean <- -0.192455782396496
PD2182sequel_log2ipd_var <- 0.8454646
PD2182sequel_log2ipd_size <- 166271875
k_normBy_ab_log2ipdratio_mean <- -0.008352887
k_normBy_ab_log2ipdratio_var <- 0.3137371
k_normBy_ab_log2ipdratio_size <- 19539347
l_normBy_cd_log2ipdratio_mean <- -0.001731761
l_normBy_cd_log2ipdratio_var <- 0.2929514
l_normBy_cd_log2ipdratio_size <- 35514409

ab_plots <- list()
cd_plots <- list()
PD2182_plots <- list()
PD2182sequel_plots <- list()
k_normBy_ab_plots <- list()
l_normBy_cd_plots <- list()
legend_plot <- NULL
ylim_ipd <- c(-2, 3.4)
#ylim_ipdratio <- c(-2.5, 3)

for (motif_dir in motif_dirs) {
    ab_ipd <- fread(paste0(motif_dir, "/motif_ipd.ab.csv"))
    ab_plots[[motif_dir]] <- plot_motif_kinetics(ab_ipd, motif_dir, "log2 IPD", ab_log2ipd_mean, ab_log2ipd_var, ab_log2ipd_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    cd_ipd <- fread(paste0(motif_dir, "/motif_ipd.cd.csv"))
    cd_plots[[motif_dir]] <- plot_motif_kinetics(cd_ipd, motif_dir, "log2 IPD", cd_log2ipd_mean, cd_log2ipd_var, cd_log2ipd_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    PD2182_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182.csv"))
    PD2182_plots[[motif_dir]] <- plot_motif_kinetics(PD2182_ipd, motif_dir, "log2 IPD", PD2182_log2ipd_mean, PD2182_log2ipd_var, PD2182_log2ipd_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    PD2182sequel_ipd <- fread(paste0(motif_dir, "/motif_ipd.PD2182sequel.csv"))
    PD2182sequel_plots[[motif_dir]] <- plot_motif_kinetics(PD2182sequel_ipd, motif_dir, "log2 IPD", PD2182sequel_log2ipd_mean, PD2182sequel_log2ipd_var, PD2182sequel_log2ipd_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    #k_normBy_ab_ipdratio <- fread(paste0(motif_dir, "/motif_ipdratio.k_normBy_ab.csv"))
    #k_normBy_ab_plots[[motif_dir]] <- plot_motif_kinetics(k_normBy_ab_ipdratio, motif_dir, "log2 IPD ratio", k_normBy_ab_log2ipdratio_mean, k_normBy_ab_log2ipdratio_var, k_normBy_ab_log2ipdratio_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipdratio)
    #l_normBy_cd_ipdratio <- fread(paste0(motif_dir, "/motif_ipdratio.l_normBy_cd.csv"))
    #l_normBy_cd_plots[[motif_dir]] <- plot_motif_kinetics(l_normBy_cd_ipdratio, motif_dir, "log2 IPD ratio", l_normBy_cd_log2ipdratio_mean, l_normBy_cd_log2ipdratio_var, l_normBy_cd_log2ipdratio_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipdratio)
    if (is.null(legend_plot)) {
        legend_plot <- get_legend(plot_motif_kinetics(ab_ipd, motif_dir, "log2 IPD", ab_log2ipd_mean, ab_log2ipd_var, ab_log2ipd_size))
    }
}

ab_plots[["legend"]] <- legend_plot
cd_plots[["legend"]] <- legend_plot
PD2182_plots[["legend"]] <- legend_plot
PD2182sequel_plots[["legend"]] <- legend_plot
k_normBy_ab_plots[["legend"]] <- legend_plot
l_normBy_cd_plots[["legend"]] <- legend_plot

#motifs1 <- c("ACATMTGG", "ACGCRTG", "AGGCTT_4", "AGGY_7", "ATCAGCTG", "ATWTTB", "CGRAACCCG", "CTGDAR", "CTGRAA", "GAGD", "GCGC", "GGN_10", "GGN_4", "GGTCTCGC", "GKCGGY", "GTAS", "HKGAGCA", "HRGTA", "RGTAB", "RTAB", "TACTGTAG", "TATAA", "TTTT", "YGGWR", "legend")
#motifs2 <- c("tandem_repeat_s4_01down", "tandem_repeat_s4_01up", "tandem_repeat_s4_02down", "tandem_repeat_s4_02up", "tandem_repeat_s4_03down", "tandem_repeat_s4_03up", "tandem_repeat_s4_04down", "tandem_repeat_s4_04up", "tandem_repeat_s4_05down", "tandem_repeat_s4_05up", "legend")
motifs_high_ipd <- c("ACGCRTG", "ATCAGCTG", "GGN_10", "GGN_4", "HKGAGCA", "HRGTA", "RGTAB", "legend")
motifs_low_ipd <- c("ACATMTGG", "AGGCTT_4", "AGGY_7", "CTGDAR", "RTAB", "TACTGTAG", "legend")

pdf_width <- 10
pdf_height <- 8
pdf("arrange_plot_motif_kinetics.stderr.wide_y.ab.pdf", width = pdf_width, height = pdf_height)
plot_grid(plotlist = ab_plots[motifs_high_ipd], align = "none", ncol = 4)
plot_grid(plotlist = ab_plots[motifs_low_ipd], align = "none", ncol = 4)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.cd.pdf", width = pdf_width, height = pdf_height)
plot_grid(plotlist = cd_plots[motifs_high_ipd], align = "none", ncol = 4)
plot_grid(plotlist = cd_plots[motifs_low_ipd], align = "none", ncol = 4)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.PD2182.pdf", width = pdf_width, height = pdf_height)
plot_grid(plotlist = PD2182_plots[motifs_high_ipd], align = "none", ncol = 4)
plot_grid(plotlist = PD2182_plots[motifs_low_ipd], align = "none", ncol = 4)
invisible(dev.off())
pdf("arrange_plot_motif_kinetics.stderr.wide_y.PD2182sequel.pdf", width = pdf_width, height = pdf_height)
plot_grid(plotlist = PD2182sequel_plots[motifs_high_ipd], align = "none", ncol = 4)
plot_grid(plotlist = PD2182sequel_plots[motifs_low_ipd], align = "none", ncol = 4)
invisible(dev.off())
