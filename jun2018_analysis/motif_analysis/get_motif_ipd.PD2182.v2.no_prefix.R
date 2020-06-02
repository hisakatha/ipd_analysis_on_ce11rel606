library(data.table)
args <- commandArgs(trailingOnly = TRUE)
motif <- args[1]
motif_length <- nchar(motif)

source("../check_chromosome_order.R")

get_motif_ipd <- function(data, occ, out, ref_skips){
    # data: kinetics data in data.table
    # occ: motif occurrences, 0-based start position
    # out: output file name
    if (data[, .N] == 0 || occ[, .N] == 0) { fwrite(data.table(), file = out); return() }
    check_chromosome_order(data, reference)
    motif_margin <- 20
    # start_margin should be equal to end_margin to handle minus strands
    start_margin <- motif_margin
    end_margin <- motif_margin
    collect_length <- start_margin + motif_length + end_margin
    col_names_start <- paste0("s", rep(1:start_margin, each=2), c("p", "n"))
    col_names_motif <- paste0("m", rep(1:motif_length, each=2), c("p", "n"))
    col_names_end <- paste0("e", rep(1:end_margin, each=2), c("p", "n"))
    col_names <- c(col_names_start, col_names_motif, col_names_end)
    occ_n <- nrow(occ)
    current_table_list <- list()
    for (i in 1:occ_n){
        occ_chr <- occ[i, chr]
        occ_start <- occ[i, start]
        index_start <- occ_start - start_margin
        collect_start <- max(0, index_start)
        occ_end <- occ_start + motif_length - 1
        chr_end <- length(reference[[occ_chr]]) - 1
        index_end <- occ_end + end_margin
        collect_end <- min(chr_end, index_end)
        out_values <- c(rep(0, (collect_start - index_start) * 2), data[ref_skips[occ_chr] + (collect_start * 2 + 1):(collect_end * 2 + 2), tMean], rep(0, (index_end - collect_end) * 2))
        occ_strand <- occ[i, strand]
        if (occ_strand == "+"){
            # no-op
        } else if (occ_strand == "-"){
            # start_margin should be equal to end_margin to handle minus strands in this way
            out_values <- rev(out_values)
        } else {
            print(paste("Unknown strand in merged_occ at line", i))
            print(occ[i])
            stop()
        }
        current_table_list <- c(current_table_list, list(data.table(position = rep(1:collect_length, each = 2),
                                                           strand = rep(c("+","-"), collect_length),
                                                           value = out_values,
                                                           label = col_names,
                                                           src = i)))
    }
    out_data_table <- rbindlist(current_table_list)
    fwrite(out_data_table, file = out)
}

source("../../load_and_save_PD2182_data.saved.fst.PD2182.R", chdir = TRUE)

#bed_col <- c("chr", "start", "end", "name", "score", "strand")
merged_occ_col <- c("chr", "start", "strand")

PD2182_occ <- fread(paste0("motif_occ.deep_kinetics_region.PD2182.bed.merged.sorted.slop20.bed.fa.bed.merged_occ"), sep = " ", col.names = merged_occ_col)

reference_skip <- c(0, cumsum(width(reference)) * 2)
names(reference_skip) <- names(reference)

get_motif_ipd(PD2182_data, PD2182_occ, "motif_ipd.PD2182.csv", reference_skip)
