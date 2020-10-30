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

plot_motif_kinetics <- function(kinetics, title_text, ylab_text, global_mean, global_var, global_size){
    if (kinetics[, .N] == 0) {
        g1 <- ggplot(NULL) + ggtitle(title_text) + geom_text(aes(x = 0, y = 0, label = "No Data"), size = 12) +
            xlab("") + ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
        return(g1)
    }
    kinetics_summary <- kinetics[value > 0 & is.finite(value)][, .(position, strand, label, log2value = log2(value))][, .(mean = mean(log2value), var = var(log2value), .N), by = .(position, strand, label)]
    #kinetics_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    kinetics_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics_summary$strand <- factor(kinetics_summary$strand, levels = c("+", "-"))
    kinetics_summary$region <- factor(kinetics_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    kinetics_summary[,"position" := .(position - 20)]
    g1 <- ggplot(kinetics_summary[region == "Motif"], aes(position, mean)) + geom_point(alpha = 0.5) + geom_errorbar(aes(ymin = mean - sqrt(var / N), ymax = mean + sqrt(var / N))) +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank()) +
        ggtitle(title_text) + ylab(ylab_text) +
        geom_hline(yintercept = global_mean, linetype = "dashed") +
        scale_x_continuous(breaks = get_integer_breaks)
        #geom_text(aes(y = Inf, label = format(pvalue, digits = 3, scientific = T)), show.legend = F, angle = 90, hjust = 1.2, size = 2, color = "black")
        #scale_x_continuous(breaks = kinetics_summary$position, labels = kinetics_summary$label)
        #scale_x_continuous(breaks = 1:kinetics_summary[,max(position)], labels = kinetics_summary$label)
    return(g1)
}

motifs_high_ipd <- c("AAAABT", "AGAGAGTA", "AGCTATAT", "AGGCAGGC", "ATGGGAYA", "ATTGTTAC",
"CAAACTAC", "CAGYTG", "CCAATCAG", "CRACGAS",
"DCGAGACC", "DGCTTC", "DHNDATCAGCTGATSD", "DNKAMCAGSTGKYMNN",
"GAAGGATC", "GATATRGY", "GCCVCGCC", "GCGACCTA", "GCGCGCGC", "GGHGGY", "GTAGATCA", "GTATCGTA",
"HSDAACAKCTGTCHHW", "KVGACAGCTGTCRDVV", "MRWWNCACCTGCTNDNN",
"NBMAACATMTGGTNKM", "NDNGRCACGTGTCNNND", "NDRHYACGCRTGKYHND", "NMDAACAGCTGTYNHN", "NNAACABCTGTCVNNA", "NNDACCACGTGGTYND",
"RKRCGCGTGTCCCWN", "RNGRCGCGTGTCCYND", "RNMVCASCTGBGBCMT",
"TGACGTCA", "TGGTGSA", "TNDRSCACGTGSBANN",
"VRKRMCACCTGCTNVHH", "WNAHGGCACGTGYCNDA", "WNDATCACGTGAYHWDN",
"YBKRVARWAAADWG", "YGGAR")

motifs_low_ipd <- c("AATMAATA", "ADNRATWAGWANWH", "AGACGCAG", "AWBWWTWAAAANWA",
"CGGYTTGA", "CTGDAR", "CTKCAA", "CWHAAAAAAWAAAA",
"GCGCGTCA", "GTGTGTGY",
"HSDAACAKCTGTCHHW", "KVGACAGCTGTCRDVV",
"NBMAACATMTGGTNKM", "NBNAATATAAWVHT", "NMDAACAGCTGTYNHN", "NNAACABCTGTCVNNA",
"RNMVCASCTGBGBCMT",
"TACCCCKA", "TACCTTGA", "TGACGTCA",
"WAYYAWATWTRKVV", "WHWAATAAAACNVW", "YBKRVARWAAADWG")

motifs_high_ipd <- c(c("GGN_4", "RGTA"), motifs_high_ipd)
motifs_low_ipd <- c("TACTGTAG", motifs_low_ipd)

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

# These constants were derived from deep regions (coverage >= 25).
all_ab_deep_log2_modelPrediction_mean <- -0.217030377571732
all_ab_deep_log2_modelPrediction_var <- 0.400738
all_ab_deep_log2_modelPrediction_size <- 26968723
all_cd_deep_log2_modelPrediction_mean <- -0.218477938679962
all_cd_deep_log2_modelPrediction_var <- 0.401117
all_cd_deep_log2_modelPrediction_size <- 45958624
all_PD2182_deep_log2_modelPrediction_mean <- -0.13304048715978
all_PD2182_deep_log2_modelPrediction_var <- 0.254342
all_PD2182_deep_log2_modelPrediction_size <- 1356981
all_PD2182sequel_deep_log2_modelPrediction_mean <- -0.236980241840712
all_PD2182sequel_deep_log2_modelPrediction_var <- 0.406119
all_PD2182sequel_deep_log2_modelPrediction_size <- 166271875

ab_plots <- list()
cd_plots <- list()
PD2182_plots <- list()
PD2182sequel_plots <- list()
k_normBy_ab_plots <- list()
l_normBy_cd_plots <- list()
legend_plot <- NULL
ylim_ipd <- c(-2, 3.5)
ylab_expr <- bquote(log[2]~IPD~estimate)

for (motif_dir in motif_dirs) {
    ab_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.ab.c_elegans.csv"))
    ab_plots[[motif_dir]] <- plot_motif_kinetics(ab_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_ab_deep_log2_modelPrediction_mean, all_ab_deep_log2_modelPrediction_var, all_ab_deep_log2_modelPrediction_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    cd_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.cd.c_elegans.csv"))
    cd_plots[[motif_dir]] <- plot_motif_kinetics(cd_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_cd_deep_log2_modelPrediction_mean, all_cd_deep_log2_modelPrediction_var, all_cd_deep_log2_modelPrediction_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    PD2182_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182.c_elegans.csv"))
    PD2182_plots[[motif_dir]] <- plot_motif_kinetics(PD2182_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_PD2182_deep_log2_modelPrediction_mean, all_PD2182_deep_log2_modelPrediction_var, all_PD2182_deep_log2_modelPrediction_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    PD2182sequel_modelPrediction <- fread(paste0(motif_dir, "/motif_modelPrediction.PD2182sequel.c_elegans.csv"))
    PD2182sequel_plots[[motif_dir]] <- plot_motif_kinetics(PD2182sequel_modelPrediction, lookup_motif_title(motif_dir), ylab_expr, all_PD2182sequel_deep_log2_modelPrediction_mean, all_PD2182sequel_deep_log2_modelPrediction_var, all_PD2182sequel_deep_log2_modelPrediction_size) + theme(legend.position = "none") + coord_cartesian(ylim = ylim_ipd)
    if (is.null(legend_plot)) {
        legend_plot <- get_legend(plot_motif_kinetics(ab_modelPrediction, motif_dir, ylab_expr, all_ab_deep_log2_modelPrediction_mean, all_ab_deep_log2_modelPrediction_var, all_ab_deep_log2_modelPrediction_size))
    }
}

ab_plots[["legend"]] <- legend_plot
cd_plots[["legend"]] <- legend_plot
PD2182_plots[["legend"]] <- legend_plot
PD2182sequel_plots[["legend"]] <- legend_plot
k_normBy_ab_plots[["legend"]] <- legend_plot
l_normBy_cd_plots[["legend"]] <- legend_plot

pdf_width <- 14
pdf_height_high <- 9 * 2.2
pdf_height_low <- 5 * 2.2
ncol1 <- 5
pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.ab.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = ab_plots[motifs_high_ipd], align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.cd.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = cd_plots[motifs_high_ipd], align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.PD2182.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = PD2182_plots[motifs_high_ipd], align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.PD2182sequel.high.pdf", width = pdf_width, height = pdf_height_high)
plot_grid(plotlist = PD2182sequel_plots[motifs_high_ipd], align = "none", ncol = ncol1)
invisible(dev.off())

pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.ab.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = ab_plots[motifs_low_ipd], align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.cd.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = cd_plots[motifs_low_ipd], align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.PD2182.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = PD2182_plots[motifs_low_ipd], align = "none", ncol = ncol1)
invisible(dev.off())
pdf("arrange_plot_motif_modelPrediction.stderr.wide_y.motif_region.one_row.v3.PD2182sequel.low.pdf", width = pdf_width, height = pdf_height_low)
plot_grid(plotlist = PD2182sequel_plots[motifs_low_ipd], align = "none", ncol = ncol1)
invisible(dev.off())
