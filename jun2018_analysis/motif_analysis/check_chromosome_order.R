library(Biostrings)
library(data.table)

check_chromosome_order <- function(data, reference) {
    # data: data.table data with column 'refName'
    # reference: Biostrings fasta
    ref_order <- names(reference)
    data_order <- data[, rle(refName)$values]
    if (!identical(ref_order, data_order)) {
        stop("The order of chromosomes in data is not consistent with that in the reference\n")
    }
}
