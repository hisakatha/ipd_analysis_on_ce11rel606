library(data.table)
library(ggplot2)
source("/glusterfs/hisakatha/methylation/region.t.test.R")
source("/glusterfs/hisakatha/methylation/window_methylation_density.R")

source("load_jun2018_data.saved.fst.R")

#bed_col <- c("chr", "start", "end", "name", "score", "strand")
#gff_col <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
#intron_outside_exon <- fread("/glusterfs/hisakatha/wormbase/c_elegans.PRJNA13758.WS267.annotations.gff3.gz.ucsc_chr.intron_outside_exon", sep = "\t", col.names = bed_col)

coverage_thres <- 25
ipdRatio_thres <- 2.0
score_thres <- 20
running_width <- 2000
target_base <- "A"

ab_modified <- simple_data_modified(ab_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
cd_modified <- simple_data_modified(cd_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
k_modified <- simple_data_modified(k_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
l_modified <- simple_data_modified(l_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
abcd_modified <- simple_data_modified(abcd_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
kl_modified <- simple_data_modified(kl_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
k_normBy_ab_modified <- normalized_data_modified(k_normBy_ab_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
l_normBy_cd_modified <- normalized_data_modified(l_normBy_cd_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
kl_normBy_abcd_modified <- normalized_data_modified(kl_normBy_abcd_data, coverage_thres, ipdRatio_thres, score_thres, target_base)
ab_callable <- simple_data_callable(ab_data, coverage_thres, target_base)
cd_callable <- simple_data_callable(cd_data, coverage_thres, target_base)
k_callable <- simple_data_callable(k_data, coverage_thres, target_base)
l_callable <- simple_data_callable(l_data, coverage_thres, target_base)
abcd_callable <- simple_data_callable(abcd_data, coverage_thres, target_base)
kl_callable <- simple_data_callable(kl_data, coverage_thres, target_base)
k_normBy_ab_callable <- normalized_data_callable(k_normBy_ab_data, coverage_thres, target_base)
l_normBy_cd_callable <- normalized_data_callable(l_normBy_cd_data, coverage_thres, target_base)
kl_normBy_abcd_callable <- normalized_data_callable(kl_normBy_abcd_data, coverage_thres, target_base)
ab_data <- running_methylation_density(ab_data, ab_modified, ab_callable, running_width)
cd_data <- running_methylation_density(cd_data, cd_modified, cd_callable, running_width)
k_data <- running_methylation_density(k_data, k_modified, k_callable, running_width)
l_data <- running_methylation_density(l_data, l_modified, l_callable, running_width)
abcd_data <- running_methylation_density(abcd_data, abcd_modified, abcd_callable, running_width)
kl_data <- running_methylation_density(kl_data, kl_modified, kl_callable, running_width)
k_normBy_ab_data <- running_methylation_density(k_normBy_ab_data, k_normBy_ab_modified, k_normBy_ab_callable, running_width)
l_normBy_cd_data <- running_methylation_density(l_normBy_cd_data, l_normBy_cd_modified, l_normBy_cd_callable, running_width)
kl_normBy_abcd_data <- running_methylation_density(kl_normBy_abcd_data, kl_normBy_abcd_modified, kl_normBy_abcd_callable, running_width)

ab_data_valid <- ab_data[coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
cd_data_valid <- cd_data[coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
k_data_valid <- k_data[coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
l_data_valid <- l_data[coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
abcd_data_valid <- abcd_data[coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
kl_data_valid <- kl_data[coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
k_normBy_ab_data_valid <- k_normBy_ab_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
l_normBy_cd_data_valid <- l_normBy_cd_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]
kl_normBy_abcd_data_valid <- kl_normBy_abcd_data[native_coverage >= coverage_thres & control_coverage >= coverage_thres & base == target_base & is.finite(modified_with_callable)]

union <- cbind(a = k_normBy_ab_data, b = l_normBy_cd_data)
union_valid <- union[a.native_coverage >= coverage_thres & a.control_coverage >= coverage_thres & a.base == target_base & is.finite(a.modified_with_callable) &
    b.native_coverage >= coverage_thres & b.control_coverage >= coverage_thres & b.base == target_base & is.finite(b.modified_with_callable)]
union_valid[, c("refName", "tpl", "strand", "base", "modified_with_callable") := .(a.refName, a.tpl, a.strand, a.base, a.modified_with_callable * b.modified_with_callable)]

is_within_region_without_strand_gff <- function(data, chr, gff_starts, gff_ends){
    bed_starts <- gff_starts - 1
    bed_ends <- gff_ends
    is_within_region_without_strand(data, data.table(chr, start = bed_starts, end = bed_ends))
}

database <- fread("get_enriched_gene.feature_database")
summary_data <- fread("calc_enrichment.both.csv")

get_enriched_gene_OLD <- function(data, data_name, feature_name, database, summary_data){
    # "data" should be a valid data, that is, consist of callable bases.
    all_callable <- data[, .N]
    all_modified <- data[modified_with_callable == 1, .N]
    all_unmodified <- all_callable - all_modified
    #database <- database[
    #         feature == feature_name, .(sample = data_name, feature = feature[1], flags = list(is_within_region_without_strand_gff(data, chr, start, end))), by = .(gene)][
    #         , .(sample, feature, gene, callable = data[flags[[1]], .N], modified = data[flags[[1]] & (modified_with_callable == 1), .N])][
    #         , c("fraction") := .(modified / callable)
    #        ]
    database <- database[
             feature == feature_name, .(sample = data_name, feature = feature[1], flags = list(is_within_region_without_strand_gff(data, chr, start, end))), by = .(gene)]
    print(database)
    print(length(database[1, flags]))
    database <- database[, .(sample, feature, gene, callable = data[flags[[1]], .N], modified = data[flags[[1]] & (modified_with_callable == 1), .N])][
             , c("fraction") := .(modified / callable)
            ]
    print(database)
    all_feature_callable <- database[, sum(callable)]
    all_feature_modified <- database[, sum(modified)]
    all_feature_unmodified <- all_feature_callable - all_feature_modified
    merged_feature_callable <- summary_data[sample == data_name & region == feature_name, callableNum]
    merged_feature_modified <- summary_data[sample == data_name & region == feature_name, modifiedNum]
    merged_feature_unmodified <- merged_feature_callable - merged_feature_modified
    database[, c("pvalue_using_all", "pvalue_using_features", "pvalue_using_merged_features") := .(
                    fisher.test(rbind(c(all_unmodified, all_modified), c(callable - modified, modified)))$p.value,
                    fisher.test(rbind(c(all_feature_unmodified, all_feature_modified), c(callable - modified, modified)))$p.value,
                    fisher.test(rbind(c(merged_feature_unmodified, merged_feature_modified), c(callable - modified, modified)))$p.value
                                                                                                   ), by = gene]
    database[order(-fraction)]
}

get_enriched_gene <- function(data, data_name, feature_name, database, summary_data){
    # "data" should be a valid data, that is, consist of callable bases.
    cat(sprintf("Start get_enriched_gene() for data: %s; feature: %s\n", data_name, feature_name))
    all_callable <- data[, .N]
    all_modified <- data[modified_with_callable == 1, .N]
    cat(sprintf("all_callable: %d; all_modified: %d\n", all_callable, all_modified))
    all_unmodified <- all_callable - all_modified
    database <- database[feature == feature_name]
    genes <- database[, .N, by = gene][, gene]
    result <- NULL
    for(g in genes){
        g_info <- database[gene == g]
        flags = is_within_region_without_strand_gff(data, g_info[,chr], g_info[,start], g_info[,end])
        if(length(flags) != data[, .N]){ cat(sprintf("gene: %s; length(flags): %d; data[, .N]: %d\n", g, length(flags), data[, .N])) }
        g_result <- data.table(data = data_name, feature = feature_name, gene = g, callable = data[flags, .N], modified = data[flags & (modified_with_callable == 1), .N])[
             , c("fraction") := .(modified / callable)
            ]
        result <- rbind(result, g_result)
    }
    all_feature_callable <- result[, sum(callable)]
    all_feature_modified <- result[, sum(modified)]
    cat(sprintf("all_feature_callable: %d; all_feature_modified: %d\n", all_feature_callable, all_feature_modified))
    all_feature_unmodified <- all_feature_callable - all_feature_modified
    merged_feature_callable <- summary_data[sample == data_name & region == feature_name, callableNum]
    merged_feature_modified <- summary_data[sample == data_name & region == feature_name, modifiedNum]
    cat(sprintf("merged_feature_callable: %d; merged_feature_modified: %d\n", merged_feature_callable, merged_feature_modified))
    merged_feature_unmodified <- merged_feature_callable - merged_feature_modified
    result[, c("pvalue_using_all", "pvalue_using_features", "pvalue_using_merged_features") := .(
                    fisher.test(rbind(c(all_unmodified, all_modified), c(callable - modified, modified)))$p.value,
                    fisher.test(rbind(c(all_feature_unmodified, all_feature_modified), c(callable - modified, modified)))$p.value,
                    fisher.test(rbind(c(merged_feature_unmodified, merged_feature_modified), c(callable - modified, modified)))$p.value
                                                                                                   ), by = gene]
    result[order(-fraction)]
}

summarize.enriched.gene <- function(data, data_name, database, summary_data){
    feature_list <- c("intron", "exon")
    d_all <- NULL
    for(f in feature_list){
        d <- get_enriched_gene(data, data_name, f, database, summary_data)
        d_all <- rbind(d_all, d)
    }
    return(d_all)
}

d_all <- NULL

d <- summarize.enriched.gene(k_normBy_ab_data_valid, "k_normBy_ab", database, summary_data)
d_all <- rbind(d_all, d)

d <- summarize.enriched.gene(l_normBy_cd_data_valid, "l_normBy_cd", database, summary_data)
#d <- summarize.enriched.gene(l_normBy_cd_data_valid[100001:300000], "l_normBy_cd", database, summary_data)
d_all <- rbind(d_all, d)

d <- summarize.enriched.gene(kl_normBy_abcd_data_valid, "kl_normBy_abcd", database, summary_data)
d_all <- rbind(d_all, d)

d <- summarize.enriched.gene(union_valid, "intersection of callable k_normBy_ab and l_normBy_cd", database, summary_data)
d_all <- rbind(d_all, d)

fwrite(d_all, file = "get_enriched_gene.both.csv")

