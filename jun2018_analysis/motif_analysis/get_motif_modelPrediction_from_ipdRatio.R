library(data.table)
args <- commandArgs(trailingOnly = TRUE)

# Raw IPDs
ipds <- fread(args[1], header = TRUE)

# IPD ratios based on PacBio IPD estimates
ipdratios <- fread(args[2], header = TRUE)

output_path <- args[3]

if(ipds[,.N] == 0){ fwrite(data.table(), file = output_path); quit("no", 0) }

stopifnot(all(ipds[, position] == ipdratios[, position]))
stopifnot(all(ipds[, strand] == ipdratios[, strand]))
stopifnot(all((ipds[, value] > 0) == (ipdratios[, value] > 0)))
stopifnot(all(ipds[, label] == ipdratios[, label]))
stopifnot(all(ipds[, src] == ipdratios[, src]))

ipd_estimates <- data.table(position = ipds[, position], strand = ipds[, strand],
                            value = ifelse(ipdratios[, value] > 0, ipds[, value] / ipdratios[, value], 0),
                            label = ipds[, label], src = ipds[, src])
fwrite(ipd_estimates, file = output_path)
