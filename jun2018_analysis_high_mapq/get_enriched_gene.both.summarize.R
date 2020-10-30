library(data.table)

enrich <- NULL
enrich_tmp <- fread("get_enriched_gene.both.csv")
enrich <- rbind(enrich, enrich_tmp)
enrich_tmp <- fread("get_enriched_gene.both.tts_down500.csv")
enrich <- rbind(enrich, enrich_tmp)

# You may want the number of features
feature_database <- fread("get_enriched_gene.feature_database")

enrich[, "BH_using_all" := .(p.adjust(pvalue_using_all, method = "BH"))]

#filtered <- enrich[(data == "k_normBy_ab" | data == "l_normBy_cd") & pvalue_using_all <= 10^-5, .SD, by = .(gene)]
# Extract genes reproducibly enriched with DNA modification
filtered <- enrich[(data == "k_normBy_ab" | data == "l_normBy_cd") & BH_using_all <= 0.05, c(.SD, .N), by = .(gene, feature)][N > 1]
fwrite(x = filtered, file = "get_enriched_gene.both.summarize.csv")

