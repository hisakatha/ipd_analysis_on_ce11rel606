SHELL := bash
COLLECT_BIN := /glusterfs/hisakatha/methylation/svr_ipd/collect_ipd_per_kmer/collect_ipd
k := 2
l := 20
COV_THRES := 25

.DELETE_ON_ERROR:
.SECONDEXPANSION:

target =
all: $$(target)

.PHONY: all

input_file := ipd_summary_high_mapq.h5
output_prefix := dimer_ipd_stats

output_ab := $(output_prefix).ab.csv
$(output_ab): ../../jun2018_ab_pcr_vc2010_op50_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<
output_cd := $(output_prefix).cd.csv
$(output_cd): ../../jun2018_cd_pcr_vc2010_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<
output_k := $(output_prefix).k.csv
$(output_k): ../../jun2018_k_nopcr_vc2010_op50_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<
output_l := $(output_prefix).l.csv
$(output_l): ../../jun2018_l_nopcr_vc2010_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<
output_abcd := $(output_prefix).abcd.csv
$(output_abcd): ../../jun2018_abcd_pcr_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<
output_kl := $(output_prefix).kl.csv
$(output_kl): ../../jun2018_kl_nopcr_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<
output_PD2182 := $(output_prefix).PD2182.csv
$(output_PD2182): ../../PD2182-20150119_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<
output_PD2182sequel := $(output_prefix).PD2182sequel.csv
$(output_PD2182sequel): ../../PD2182-sequel_no_chunk/$(input_file)
	$(COLLECT_BIN) -k $(k) -l $(l) -t $(COV_THRES) -o $@ $<

target += $(output_ab) $(output_cd) $(output_k) $(output_l) $(output_PD2182sequel) $(output_abcd) $(output_kl)
target += $(output_PD2182)
