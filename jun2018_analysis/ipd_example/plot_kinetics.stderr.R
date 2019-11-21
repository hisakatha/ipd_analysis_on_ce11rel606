library(data.table)
library(ggplot2)

welch.t.pvalue <- function(a_mean, a_variance, a_size, b_mean, b_variance, b_size) {
    a_sq <- (a_variance ^ 2) / a_size
    b_sq <- (b_variance ^ 2) / b_size
    t_stat <- (a_mean - b_mean) / sqrt(a_sq + b_sq)
    dof <- ((a_sq + b_sq) ^ 2) / ((a_sq ^ 2) / (a_size - 1) + (b_sq ^ 2) / (b_size - 1))
    # Two-sided p-value
    pt(-abs(t_stat), dof) * 2
}

plot_motif_kinetics <- function(kinetics, title_text, ylab_text, global_mean, global_var, global_size){
    kinetics_summary <- kinetics[value > 0 & is.finite(value)][, .(position, strand, label, log2value = log2(value))][, .(mean = mean(log2value), var = var(log2value), .N), by = .(position, strand, label)]
    #kinetics_summary[, "label" := list(substring())]
    #kinetics_summary[, "pvalue" := list(welch.t.pvalue(mean, var, N, global_mean, global_var, global_size))]
    #kinetics_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
    kinetics_summary$strand <- factor(kinetics_summary$strand, levels = c("+", "-"))
    #kinetics_summary$region <- factor(kinetics_summary$region, levels = c("Upstream", "Motif", "Downstream", "Unknown"))
    g1 <- ggplot(kinetics_summary, aes(position, mean)) + geom_point() +
        facet_grid(strand ~ ., labeller = label_both) +
        theme(panel.grid = element_blank(), axis.ticks.x = element_blank()) +
        coord_cartesian(ylim = c(-4, 4)) +
        ggtitle(title_text) + ylab(ylab_text) + xlab("Sequence") +
        geom_hline(yintercept = global_mean, linetype = "dashed") +
        scale_x_continuous(breaks = 1:kinetics_summary[,max(position)], labels = unlist(strsplit(dna_seq, "")))
        #geom_text(aes(y = Inf, label = format(pvalue, digits = 3, scientific = T)), show.legend = F, angle = 90, hjust = 1.2, size = 2, color = "black")
        #scale_x_continuous(breaks = kinetics_summary$position, labels = kinetics_summary$label)
        #scale_x_continuous(breaks = 1:kinetics_summary[,max(position)], labels = kinetics_summary$label)
    print(g1)
}

args <- commandArgs(trailingOnly = TRUE)
input_prefix <- args[1]
motif <- "IPD example"
dna_seq <- args[2]

ab_ipd <- fread(paste0(input_prefix, ".ab.csv"))
cd_ipd <- fread(paste0(input_prefix, ".cd.csv"))
k_ipd <- fread(paste0(input_prefix, ".k.csv"))
l_ipd <- fread(paste0(input_prefix, ".l.csv"))

# These constants were derived from deep regions (coverage >= 25).
ab_ipd_mean <- 1.001282
ab_ipd_var <- 0.5264925
ab_ipd_size <- 27095150
cd_ipd_mean <- 1.006526
cd_ipd_var <- 0.4974373
cd_ipd_size <- 46121818
k_normBy_ab_ipdratio_mean <- 1.073165
k_normBy_ab_ipdratio_var <- 0.2105606
k_normBy_ab_ipdratio_size <- 19539347
l_normBy_cd_ipdratio_mean <- 1.072319
l_normBy_cd_ipdratio_var <- 0.1875473
l_normBy_cd_ipdratio_size <- 35514409
ab_log2ipd_mean <- -0.1889734
ab_log2ipd_var <- 0.4976316
ab_log2ipd_size <- 27095150
cd_log2ipd_mean <- -0.1902569
cd_log2ipd_var <- 0.5337555
cd_log2ipd_size <- 46121818
k_log2ipd_mean <- -0.204882933673383
k_log2ipd_var <- 0.778915115928296
k_log2ipd_size <- 134105762
l_log2ipd_mean <- -0.201514159692603
l_log2ipd_var <- 0.769850324083848
l_log2ipd_size <- 140695865
k_normBy_ab_log2ipdratio_mean <- -0.008352887
k_normBy_ab_log2ipdratio_var <- 0.3137371
k_normBy_ab_log2ipdratio_size <- 19539347
l_normBy_cd_log2ipdratio_mean <- -0.001731761
l_normBy_cd_log2ipdratio_var <- 0.2929514
l_normBy_cd_log2ipdratio_size <- 35514409

pdf(paste0(input_prefix, ".pdf"), width = 3, height = 2)
plot_motif_kinetics(ab_ipd, paste0(motif, ": sample AB"), bquote(log[2]~IPD), ab_log2ipd_mean, ab_log2ipd_var, ab_log2ipd_size)
plot_motif_kinetics(cd_ipd, paste0(motif, ": sample CD"), bquote(log[2]~IPD), cd_log2ipd_mean, cd_log2ipd_var, cd_log2ipd_size)
plot_motif_kinetics(k_ipd, paste0(motif, ": sample K"), bquote(log[2]~IPD), k_log2ipd_mean, k_log2ipd_var, k_log2ipd_size)
plot_motif_kinetics(l_ipd, paste0(motif, ": sample L"), bquote(log[2]~IPD), l_log2ipd_mean, l_log2ipd_var, l_log2ipd_size)
#plot_motif_kinetics(k_normBy_ab_ipdratio, paste0(motif, ": log2 IPD ratio; the mean +/- 2 times the standard error of the mean; sample K normalized by AB"), "log2 IPD ratio (the mean +/- 2 times the standard error of the mean)", k_normBy_ab_log2ipdratio_mean, k_normBy_ab_log2ipdratio_var, k_normBy_ab_log2ipdratio_size)
#plot_motif_kinetics(l_normBy_cd_ipdratio, paste0(motif, ": log2 IPD ratio; the mean +/- 2 times the standard error of the mean; sample L normalized by CD"), "log2 IPD ratio (the mean +/- 2 times the standard error of the mean)", l_normBy_cd_log2ipdratio_mean, l_normBy_cd_log2ipdratio_var, l_normBy_cd_log2ipdratio_size)
invisible(dev.off())
