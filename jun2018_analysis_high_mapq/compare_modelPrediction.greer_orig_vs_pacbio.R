library(data.table)
library(tidyverse)
# load greer_orig_data
source("load_and_save_greer_data.saved.fst.greer_orig.R")
# load greer_orig_pacbio_data
source("load_and_save_greer_data.saved.fst.greer_orig_pacbio.R")
# load greer_pacbio_data
source("load_and_save_greer_data.saved.fst.greer_pacbio.R")
# load greer_pacbio_orig_data
source("load_and_save_greer_data.saved.fst.greer_pacbio_orig.R")

min_cov_thres <- 3

flags_p6 <- greer_orig_pacbio_data[, coverage >= min_cov_thres] &
    greer_pacbio_data[, coverage >= min_cov_thres] &
    greer_pacbio_orig_data[, coverage >= min_cov_thres]
greer_orig_pacbio_modelPrediction <- greer_orig_pacbio_data[flags_p6, modelPrediction]
flag_eq1 <- all(greer_orig_pacbio_modelPrediction == greer_pacbio_data[flags_p6, modelPrediction])
cat(sprintf("orig_pacbio == pacbio: %s\n", ifelse(flag_eq1, "True", "False")))
flag_eq2 <- all(greer_orig_pacbio_modelPrediction == greer_pacbio_orig_data[flags_p6, modelPrediction])
cat(sprintf("orig_pacbio == pacbio_orig: %s\n", ifelse(flag_eq2, "True", "False")))

flags_common <- greer_orig_data[, coverage >= min_cov_thres] & greer_pacbio_data[, coverage >= min_cov_thres]
comparison_data <- greer_orig_data[flags_common, .(refName, tpl, log2_modelPrediction_orig = log2(modelPrediction),
                                                   log2_modelPrediction_pacbio = greer_pacbio_data[flags_common, log2(modelPrediction)])]

g1 <- ggplot(comparison_data) + geom_hex(aes(log2_modelPrediction_pacbio, log2_modelPrediction_orig), binwidth = 0.1) +
    xlab("log2 of predicted IPD in the public PacBio data (binding kit P6)") +
    ylab("log2 of predicted IPD in the Greer's original data (binding kit P4)") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")

pdf("compare_modelPrediction.greer_orig_vs_pacbio.pdf", height = 6, width = 6)
print(g1)
print(g1 + scale_fill_continuous(trans = "log10"))
invisible(dev.off())
