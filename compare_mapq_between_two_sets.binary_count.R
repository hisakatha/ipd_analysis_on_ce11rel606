library(data.table)

binary_count <- function(data) {
    data[, .(count = .N), keyby = .("top_mapq_celegans_over_127" = top_mapq_celegans > 127, "top_mapq_ecoli_over_127" = top_mapq_ecoli > 127)]
}

ratio_in_celegans <- function(data) {
    both_high <- data[top_mapq_celegans_over_127 == TRUE & top_mapq_ecoli_over_127 == TRUE, count]
    stopifnot(length(both_high) == 1)
    celegans_high <- data[top_mapq_celegans_over_127 == TRUE, sum(count)]
    stopifnot(length(celegans_high) == 1)
    data.table(celegans_high = celegans_high, both_high = both_high, ratio = both_high / celegans_high)
}

path <- "mapped.alignmentset.merged.sorted_by_name.bam.compare_celegans_ecoli.csv"
d1 <- fread(path)
d1_count <- binary_count(d1)
fwrite(d1_count, file = "compare_mapq_between_two_sets.binary_count.full.csv")
fwrite(ratio_in_celegans(d1_count))

d1_subset <- d1[count_celegans > 0 & count_ecoli > 0]
d1_subset_count <- binary_count(d1_subset)
fwrite(d1_subset_count, file = "compare_mapq_between_two_sets.binary_count.aligned_to_both.csv")
