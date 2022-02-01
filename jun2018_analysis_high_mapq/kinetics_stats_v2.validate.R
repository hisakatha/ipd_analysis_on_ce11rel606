library(data.table)
ipds_by_r <- fread("kinetics_stats.per_chr.csv")
ipds_by_c <- fread("kinetics_stats_v2.ab_deep.csv")
merged <- merge(ipds_by_c[,.(kmer_string,chromosome,ipd_sum,count_c=count,mean_c=ipd_sum/count)],
                ipds_by_r[base %in% c("A","C","G","T") & sample=="ab_deep" & type=="IPD", .(kmer_string=base,chromosome=refName,mean_r=mean,count_r=N)],
                by=c("kmer_string","chromosome"))
print(merged)
all(merged[, abs(mean_c - mean_r) < 0.00001])
all(merged[, count_c == count_r])
