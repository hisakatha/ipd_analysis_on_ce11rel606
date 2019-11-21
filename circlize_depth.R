#args <- commandArgs(trailingOnly = TRUE)

width <- 1000
forward_depth_file <- "mapped.alignmentset.merged.bam.forward_depth"
reverse_depth_file <- "mapped.alignmentset.merged.bam.reverse_depth"
depth_lim <- 200

#TODO: add tests for files

# forward_depth, forward_depth_summary, reverse_depth, and reverse_depth_summary will be set.
source("/glusterfs/hisakatha/methylation/summarize_coverage_depth.R")

library(data.table)
ref_info <- fread("/glusterfs/hisakatha/ce11rel606/ce11rel606.fa.fai")[, .(chr = V1, start = 1, end = V2)]
# Enlarge chrM
#modified_sector_names <- c("chrI", "chrII", "chrIII", "chrIV", "chrM (100x)", "chrV", "chrX", "E. coli REL606")
setkey(forward_depth_summary, chr)
setkey(reverse_depth_summary, chr)
modified_sector_names <- paste0(c("chrI", "chrII", "chrIII", "chrIV", "chrM (100x)", "chrV", "chrX", "E. coli REL606"), sprintf(" (mean depth:f=%.3g, r=%.3g)", forward_depth_summary[ref_info$chr, mean_depth], reverse_depth_summary[ref_info$chr, mean_depth]))
modified_sector_width <- ref_info[,end] * c(1,1,1,1,100,1,1,1)

source("/glusterfs/hisakatha/methylation/window_mean_depth.R")

#window_mean_forward_depth <- get_window_mean_depth(forward_depth_file, width)
#window_mean_reverse_depth <- get_window_mean_depth(reverse_depth_file, width)
window_mean_forward_depth <- get_window_mean_depth(depth_data = forward_depth, width = width)
window_mean_reverse_depth <- get_window_mean_depth(depth_data = reverse_depth, width = width)

library(circlize)

pdf("mapped.alignmentset.merged.bam.circlize.pdf", width = 20, height = 20)

circos.par(cell.padding = c(0.02, 0, 0.02, 0))

circos.genomicInitialize(ref_info, sector.names = modified_sector_names, sector.width = modified_sector_width)

#circos.genomicTrack(window_mean_forward_depth, ylim = c(0, coverage_lim), panel.fun = function(region, value, ...){ circos.genomicLines(region, value, type = "segment", ...); circos.yaxis() })
#circos.genomicTrack(window_mean_reverse_depth, ylim = c(0, coverage_lim), panel.fun = function(region, value, ...){ circos.genomicLines(region, value, type = "segment", ...); circos.yaxis() })

#separated_forward_depth <- split(window_mean_forward_depth, window_mean_forward_depth$window_depth > depth_lim)
#separated_reverse_depth <- split(window_mean_reverse_depth, window_mean_reverse_depth$window_depth > depth_lim)
exceed_flags <- window_mean_forward_depth[,window_depth > depth_lim]
window_mean_forward_depth[,c("window_depth", "color") := list(ifelse(exceed_flags, depth_lim, window_depth), ifelse(exceed_flags, "red", "black"))]
exceed_flags <- window_mean_reverse_depth[,window_depth > depth_lim]
window_mean_reverse_depth[,c("window_depth", "color") := list(ifelse(exceed_flags, depth_lim, window_depth), ifelse(exceed_flags, "red", "black"))]
circos.genomicTrack(window_mean_forward_depth, ylim = c(0, depth_lim), panel.fun = function(region, value, ...){ circos.genomicLines(region, value, type = "segment", col = value[[2]], ...); circos.yaxis() })
circos.genomicTrack(window_mean_reverse_depth, ylim = c(0, depth_lim), panel.fun = function(region, value, ...){ circos.genomicLines(region, value, type = "segment", col = value[[2]], ...); circos.yaxis() })

title(sprintf("\n%s (mean depth:f=%.3g, r=%.3g)\n(outer) forward depth, reverse depth (inner)", basename(getwd()), forward_depth_summary["ALL", mean_depth], reverse_depth_summary["ALL", mean_depth]), cex = 2)

invisible(dev.off())
