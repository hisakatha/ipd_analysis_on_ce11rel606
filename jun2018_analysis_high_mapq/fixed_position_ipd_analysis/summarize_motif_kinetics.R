library(data.table)
args <- commandArgs(trailingOnly = TRUE)
# A path to a csv file (position, strand, label, value, src)
kinetics_path <- args[1]
out_path <- args[2]
position_offset <- 20

kinetics <- fread(kinetics_path)
kinetics_occ <- ifelse(kinetics[, .N] == 0, 0, kinetics[, length(unique(src))])
kinetics_summary <- kinetics[value > 0 & is.finite(value)][, .(mean = mean(value), var = var(value), .N, motif_occ = kinetics_occ), by = .(position, strand, label)]
kinetics_summary[, "region" := list(ifelse(substr(label, 1, 1) == "s", "Upstream", ifelse(substr(label, 1, 1) == "m", "Motif", ifelse(substr(label, 1, 1) == "e", "Downstream", "Unknown"))))]
kinetics_summary[, "position" := .(position - position_offset)]
fwrite(kinetics_summary, file = out_path)
