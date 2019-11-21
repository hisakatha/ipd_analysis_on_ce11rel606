library(data.table)
depth_cols <- c("refName", "position", "coverage")
forward_depth <- fread("mapped.alignmentset.merged.bam.forward_depth", header = FALSE, col.names = depth_cols)
reverse_depth <- fread("mapped.alignmentset.merged.bam.reverse_depth", header = FALSE, col.names = depth_cols)

ecol_chr <- "E._coli_REL606"

cele_depth <- c(forward_depth[refName != ecol_chr, coverage], reverse_depth[refName != ecol_chr, coverage])
ecol_depth <- c(forward_depth[refName == ecol_chr, coverage], reverse_depth[refName == ecol_chr, coverage])
cat(sprintf("Mean depth: C. elegans: %g\n", mean(cele_depth)))
cat(sprintf("Mean depth: E. coli: %g\n", mean(ecol_depth)))
