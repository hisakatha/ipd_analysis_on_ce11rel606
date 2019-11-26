# 1-origin genomic positions
#pcr_targets_chr <- c("chrIV", "chrV", "chrIV", "chrIV", "chrV", "chrIII");
#pcr_targets_start <- c(14085970, 3412677, 3178656, 3575682, 19370472, 11298792);
#pcr_targets_end <- c(14087183, 3413924, 3180046, 3576932, 19371630, 11300083);

# 1-origin genomic positions
whole_targets_chr <- c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM", "E._coli_REL606");
whole_targets_start <- c(1, 1, 1, 1, 1, 1, 1, 1);
whole_targets_end <- c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942, 13794, 4629812);

# See https://github.com/mannau/h5/issues/52
library(methods);
library(h5);

# Sample AB
ab_name <- "AB_WGA_VC2010_OP50"
ab_data <- h5file(name = "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_ab_pcr_vc2010_op50/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5", mode = "r");
# Sample CD
cd_name <- "CD_WGA_VC2010"
cd_data <- h5file(name = "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_cd_pcr_vc2010/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5", mode = "r");
# Sample K
k_name <- "K_native_VC2010_OP50"
k_data <- h5file(name = "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_k_nopcr_vc2010_op50/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5", mode = "r");
# Sample L
l_name <- "L_native_VC2010"
l_data <- h5file(name = "/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/jun2018_l_nopcr_vc2010/output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5", mode = "r");

library(data.table);
library(ggplot2);
library(ggExtra);
library(grid);
#library(gridExtra);
library(cowplot);

plots_ipd_per_base <- function(data, name1, name2){
    if (nrow(data) == 0) {
        return(vector("list", 16));
    }
    p <- list();
    lim1 <- range(data$LOG2IPD1) + c(-0.5, 0.5);
    lim2 <- range(data$LOG2IPD2) + c(-0.5, 0.5);
    bases <- c("A", "C", "G", "T");
    first_chr <- bases[1];
    for (chr in bases) {
        data_base <- data[BASE == chr];
        m1 <- mean(data_base$LOG2IPDR1);
        pval1 <- wilcox.test(data_base$LOG2IPDR1, mu = 0)$p.value;
        p_ipdr1_hist <- ggplot(data_base, aes(LOG2IPDR1)) + geom_histogram(binwidth = 0.1) + xlab("log_2 IPD ratio") + ggtitle(paste0("log_2 IPD ratios of ", chr, " in ", name1), subtitle = sprintf("mean=%.3G,p=%.3G", m1, pval1));
        m2 <- mean(data_base$LOG2IPDR2);
        pval2 <- wilcox.test(data_base$LOG2IPDR2, mu = 0)$p.value;
        p_ipdr2_hist <- ggplot(data_base, aes(LOG2IPDR2)) + geom_histogram(binwidth = 0.1) + xlab("log_2 IPD ratio") + ggtitle(paste0("log_2 IPD ratios of ", chr, " in ", name2), subtitle = sprintf("mean=%.3G,p=%.3G", m2, pval2));
        if (nrow(data_base) > 0) {
            rho <- if (nrow(data_base) >= 3) {cor.test(data_base$LOG2IPD1, data_base$LOG2IPD2, method = "spearman")} else {list(p.value = NA, estimate = NA)};
            mse <- mean((data_base$LOG2IPD1 - data_base$LOG2IPD2)^2);
            p_hex <- ggplot(data_base, aes(LOG2IPD1, LOG2IPD2)) + geom_hex(binwidth = 0.2) + coord_equal() + geom_abline(slope = 1, intercept = 0, lty = 2, alpha = 0.3) + scale_x_continuous(limits = lim1) + scale_y_continuous(limits = lim2) + xlab(name1) + ylab(name2) + ggtitle(paste0(chr, "'s log_2 IPD"), subtitle = paste0("N=", nrow(data_base), ", rho=", sprintf("%.3f", rho$estimate[[1]]), ", p=", sprintf("%.3f;\nmse=%.3f", rho$p.value, mse)));
            f1 <- ecdf(data_base$LOG2IPD2 - data_base$LOG2IPD1);
            pval <- wilcox.test(data_base$LOG2IPD1, data_base$LOG2IPD2, paired = TRUE)$p.value;
            p_ecdf <- ggplot(data_base, aes(LOG2IPD2 - LOG2IPD1)) + stat_ecdf(n = 1000) + xlab(paste0("log_2 IPD (", name2, " / ", name1, ")")) + ylab("F(log_2 IPD)") + ggtitle("ECDF of log_2 IPD", subtitle = sprintf("F(-1) = %.3G, 1 - F(1) = %.3G", f1(-1), 1 - f1(1))) + theme_bw() + geom_vline(xintercept = 0, lty = 2) + annotate("text", x=0, y=0, hjust=0, vjust=0, label = paste0(" F(0) = ", sprintf("%0.3f", f1(0)), "\n p = ", sprintf("%0.3G", pval)));
            p <- c(p, list(p_ipdr1_hist, p_ipdr2_hist, p_hex, p_ecdf));
        } else {
            p <- c(p, list(p_ipdr1_hist, p_ipdr2_hist, NULL, NULL));
        }
    }
    return(p);
}

# Note: data2 is normalized by data1
compare_between_two <- function(data1, name1, data2, name2, targets_chr, targets_start, targets_end, coverage_threshold, legend_offset){
    all_data <- data.table(LOG2IPD1=numeric(), LOG2IPD2=numeric(), LOG2IPDR1=numeric(), LOG2IPDR2=numeric(), COV1=numeric(), COV2=numeric(), BASE=character());
    sp <- list();
    for (i in 1:length(targets_chr)) {
        ipdr <- paste0("/", targets_chr[i], "/ipdRatio");
        tmean <- paste0("/", targets_chr[i], "/tMean"); # Raw IPD values
        coverage <- paste0("/", targets_chr[i], "/coverage");
        base <- paste0("/", targets_chr[i], "/base");
        start <- 2 * targets_start[i] - 1;
        end <- 2 * targets_end[i];
        title <- paste0(targets_chr[i], ":", targets_start[i], "-", targets_end[i]);
        #title <- targets_chr[i];
        data <- data.table(LOG2IPD1=log2(data1[tmean][start:end]), LOG2IPD2=log2(data2[tmean][start:end]), LOG2IPDR1=log2(data1[ipdr][start:end]), LOG2IPDR2=log2(data2[ipdr][start:end]), COV1=data1[coverage][start:end], COV2=data2[coverage][start:end], BASE=data1[base][start:end]);
        #data <- data.frame(IPD1=data1[tmean][], IPD2=data2[tmean][], IPDR1=data1[ipdr][], IPDR2=data2[ipdr][], COV1=data1[coverage][], COV2=data2[coverage][], BASE=data1[base][]);
        # You can set a per-target threshold like: thresh <- f(coverage_threshold, data1, data2)
        thresh <- coverage_threshold;
        #mask <- data$COV1 >= thresh & data$COV2 >= thresh;
        masked_data <- data[COV1 >= thresh & COV2 >= thresh];
        if (thresh == coverage_threshold) {
            all_data <- merge(all_data, masked_data, all = TRUE);
        }
        # geom_point is heavy for genome-wide data
        #p <- ggplot(data, aes(COV2, IPD2 / IPD1)) + geom_point() + ggtitle(paste0("IPD (", name2, " / ", name1, ")"), subtitle = title) + ylab(paste0("IPD (", name2, " / ", name1, ")")) + xlab(paste0(name2, " per-strand coverage")) + geom_hline(yintercept = 1, lty = 2);
        data_with_ipd <- data[is.finite(LOG2IPD1) & is.finite(LOG2IPD2)]
        if (data_with_ipd[,.N] > 0) { 
            p <- ggplot(data_with_ipd, aes(pmin(COV1, COV2), LOG2IPD2 - LOG2IPD1)) + geom_hex(binwidth=c(5,0.2)) + ggtitle(paste0("log_2 IPD (", name2, " / ", name1, ")"), subtitle = title) + ylab(paste0("log_2 IPD (", name2, " / ", name1, ")")) + xlab(paste0("The lower per-strand coverage between ", name1, " and ", name2)) + geom_hline(yintercept = 0, lty = 2) + theme(legend.position = "left");
            # FIXME: ggMarginal does not work with geom_hex since ggExtra version 0.9.
            # NOTE: the result of ggMarginal for geom_hex seems wrong (compared to geom_point)
            sp2 <- ggMarginal(p, type = "histogram");
        } else { 
            sp2 <- NULL;
        }
        p <- plots_ipd_per_base(data = masked_data, name1, name2);
        stopifnot(length(p) == 16);
        #grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]], top = textGrob(paste0(title, "(Coverage > ", thresh, ")"), gp = gpar(fontsize = 20)));
        sp1 <- plot_grid(p[[1]], p[[2]], p[[5]], p[[6]], p[[9]], p[[10]], p[[13]], p[[14]], ncol = 2, align = "hv", axis = "lb");
        sp3 <- plot_grid(p[[3]], p[[4]], p[[7]], p[[8]], p[[11]], p[[12]], p[[15]], p[[16]], ncol = 2, align = "hv", axis = "lb");
        p1 <- ggplot(data, aes(COV1)) + geom_histogram(binwidth = 1) + ggtitle(paste0(name1, ": ", title)) + xlab("Per-strand coverage depth");
        p2 <- ggplot(data, aes(COV2)) + geom_histogram(binwidth = 1) + ggtitle(paste0(name2, ": ", title)) + xlab("Per-strand coverage depth");
        sp4 <- plot_grid(p1,p2, ncol=1);
        sp <- c(sp, list(plot_grid(sp4, sp1, sp2, sp3, ncol=2, labels=paste0(i + legend_offset, c("A", "B", "C", "D")))));
    }

    ## Analysis using all the targets
    p <- plots_ipd_per_base(data = all_data, name1, name2);
    p_all_hist <- plot_grid(p[[1]], p[[2]], p[[5]], p[[6]], p[[9]], p[[10]], p[[13]], p[[14]], ncol = 2, align = "hv", axis = "lb", labels = paste0(length(targets_chr) + 1 + legend_offset));
    p_all_ipd <- plot_grid(p[[3]], p[[4]], p[[7]], p[[8]], p[[11]], p[[12]], p[[15]], p[[16]], ncol = 2, align = "hv", axis = "lb", labels = paste0(length(targets_chr) + 2 + legend_offset));
    return(list(sp, p_all_hist, p_all_ipd))
}

thresh <- 25;
pdf_width <- 18;
pdf_height <- 20;
pdf(file = NULL, width = pdf_width, height = pdf_height, onefile = TRUE);
#pdf(sprintf("ipdRatio.correlation.compact.all.cov%d.pdf", thresh), width = 18, height = 20, onefile = TRUE);

### ab vs. cd
result1 <- compare_between_two(ab_data, ab_name, cd_data, cd_name, whole_targets_chr, whole_targets_start, whole_targets_end, thresh, 0);

### ab vs. k
result2 <- compare_between_two(ab_data, ab_name, k_data, k_name, whole_targets_chr, whole_targets_start, whole_targets_end, thresh, length(whole_targets_chr) + 2);

### cd vs. l
result3 <- compare_between_two(cd_data, cd_name, l_data, l_name, whole_targets_chr, whole_targets_start, whole_targets_end, thresh, (length(whole_targets_chr) + 2) * 2);

invisible(dev.off());

#pdf(sprintf("ipdRatio.correlation.compact.all.cov%d_%d_%d.v2.pdf", thresh1, thresh2, thresh3), width = 14, height = 18, onefile = TRUE);
pdf(sprintf("ipdRatio.correlation.compact.ab_cd_k_l.cov%d.pdf", thresh), width = pdf_width, height = pdf_height, onefile = TRUE);
for (sp_tmp in result1[[1]]) { print(sp_tmp); }
print(result1[[2]]);
print(result1[[3]]);
for (sp_tmp in result2[[1]]) { print(sp_tmp); }
print(result2[[2]]);
print(result2[[3]]);
for (sp_tmp in result3[[1]]) { print(sp_tmp); }
print(result3[[2]]);
print(result3[[3]]);
invisible(dev.off());

