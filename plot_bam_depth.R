library(data.table)
args <- commandArgs(trailingOnly=TRUE)
depths <- c()
for(arg in args){
    depth_data <- fread(file=arg, sep="\t")
    depths <- c(depths, depth_data$V3)
}
library(ggplot2)
gp <- ggplot(data.frame(DEPTH=depths), aes(DEPTH)) + xlab("Coverage depth per strand")
gp1 <- gp + geom_histogram(binwidth = 1) + coord_cartesian(xlim=c(0,50))
gp2 <- gp + geom_histogram(binwidth = 5)
pdf("bam_depth.pdf")
print(gp1)
print(gp2)
invisible(dev.off())

