SHELL := bash

PATTERN := $(shell cat PATTERN)
ifndef PATTERN
$(error PATTERN is undefined)
endif

occ_prefix := motif_occ

motif_ipds := motif_ipd.ab.csv motif_ipd.cd.csv
motif_wga_ipdratios := motif_ipdratio.ab.csv motif_ipdratio.cd.csv
motif_hawaiian_ipdratios := motif_ipdratio.PD2182.csv motif_ipdratio.PD2182sequel.csv
motif_ipdratios := motif_ipdratio.k_normBy_ab.csv motif_ipdratio.l_normBy_cd.csv
motif_ipdratios2 := motif_ipdratio.k.csv motif_ipdratio.l.csv
motif_kinetics_csv := $(motif_ipds) motif_ipd.PD2182.csv motif_ipd.PD2182sequel.csv $(motif_ipdratios) motif_ipd.k.csv motif_ipd.l.csv $(motif_ipdratios2)
motif_kinetics_csv += motif_ipd.abcd.csv motif_ipdratio.abcd.csv motif_ipd.kl.csv motif_ipdratio.kl.csv motif_ipdratio.kl_normBy_abcd.csv
motif_kinetics_csv_celegans := $(motif_kinetics_csv:.csv=.c_elegans.csv)
motif_kinetics_csv_celegans += $(motif_wga_ipdratios:.csv=.c_elegans.csv)
motif_kinetics_csv_celegans += $(motif_hawaiian_ipdratios:.csv=.c_elegans.csv)
motif_kinetics_csv_ecoli := $(motif_kinetics_csv_celegans:.c_elegans.csv=.e_coli.csv)
motif_modelPrediction := motif_modelPrediction.ab.csv motif_modelPrediction.cd.csv motif_modelPrediction.PD2182.csv motif_modelPrediction.PD2182sequel.csv
motif_modelPrediction += motif_modelPrediction.k.csv motif_modelPrediction.l.csv
motif_modelPrediction += motif_modelPrediction.abcd.csv motif_modelPrediction.kl.csv
motif_modelPrediction_celegans := $(motif_modelPrediction:.csv=.c_elegans.csv)
motif_modelPrediction_ecoli := $(motif_modelPrediction:.csv=.e_coli.csv)
motif_summary_input := $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli) $(motif_modelPrediction_celegans) $(motif_modelPrediction_ecoli)
motif_summary_target := $(motif_summary_input:.csv=.summary.csv)
motif_kinetics_pdf := plot_motif_kinetics.pdf
motif_kinetics_stderr_pdf := plot_motif_kinetics.stderr.pdf
motif_kinetics_stderr_widey_pdf := plot_motif_kinetics.stderr.wide_y.pdf
cor_pdf := plot_kinetics_correlation.pdf
cor_pdf2 := plot_ipdratio_correlation.pdf
cor_pdf_celegans := plot_kinetics_correlation.c_elegans.pdf
#cor_pdf2_celegans := plot_ipdratio_correlation.c_elegans.pdf
cor_pdf_celegans_ecoli := plot_kinetics_correlation.celegans_ecoli.pdf
cor_pdf_celegans_ecoli_subset1 := plot_kinetics_correlation.celegans_ecoli_subset1.pdf
cor_pdf_celegans_ecoli_subset1_positive_strand := plot_kinetics_correlation.celegans_ecoli_subset1_positive_strand.pdf

.SECONDEXPANSION:
target := $(motif_kinetics_pdf) $(motif_kinetics_stderr_pdf) $(motif_kinetics_stderr_widey_pdf) $(cor_pdf) $(motif_wga_ipdratios) $(motif_hawaiian_ipdratios) $(cor_pdf2)
target += $(motif_kinetics_csv_celegans) $(cor_pdf_celegans)
target += $(motif_modelPrediction) $(motif_modelPrediction_celegans)
target += $(motif_kinetics_csv_ecoli) $(motif_modelPrediction_ecoli) $(cor_pdf_celegans_ecoli) $(cor_pdf_celegans_ecoli_subset1) $(cor_pdf_celegans_ecoli_subset1_positive_strand)
target += $(motif_summary_target)
all: $$(target)

.DELETE_ON_ERROR:

# NOTE: It is important to use non-merged files after slopping/extension (*.slop20.bed.fa),
# rather than merged files (*.slop20.merged.bed.fa) for locating motifs near feature loci.
# I want to exclude occurrences of motifs at the edges of extended regions:
# --xxxxxXxxxxx------------
# ------------xxxxxXxxxxx--
# ---YYY-----YYY-----------
# X: feature locus
# x: extended region
# YYY: target motif
# The first motif should be captured, but the second should be ignored.

# TODO?: replace deep_region_cov25, which is not high mapq, into deep_kinetics_region
# pattern occurrence for background
occ_background := $(occ_prefix).ab.bed $(occ_prefix).cd.bed $(occ_prefix).ce11rel606.bed $(occ_prefix).k_normBy_ab.bed $(occ_prefix).l_normBy_cd.bed $(occ_prefix).PD2182.bed $(occ_prefix).PD2182sequel.bed
occ_background += $(occ_prefix).abcd.bed $(occ_prefix).kl.bed $(occ_prefix).kl_normBy_abcd.bed
$(occ_prefix).ce11rel606.bed: /glusterfs/hisakatha/ce11rel606/ce11rel606.fa
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
DIR := /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606
FILE := mapped.alignmentset.merged.high_mapq.bam.cov25.slop20.bed.fa
$(occ_prefix).ab.bed: $(DIR)/jun2018_ab_pcr_vc2010_op50_no_chunk/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).cd.bed: $(DIR)/jun2018_cd_pcr_vc2010_no_chunk/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).PD2182.bed: $(DIR)/PD2182-20150119_no_chunk/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).PD2182sequel.bed: $(DIR)/PD2182-sequel_no_chunk/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).k_normBy_ab.bed: $(DIR)/jun2018_analysis_high_mapq/deep_region_cov25.k_normBy_ab.slop20.bed.fa
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).l_normBy_cd.bed: $(DIR)/jun2018_analysis_high_mapq/deep_region_cov25.l_normBy_cd.slop20.bed.fa
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).abcd.bed: $(DIR)/jun2018_abcd_pcr_no_chunk/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).kl.bed: $(DIR)/jun2018_kl_nopcr_no_chunk/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(occ_prefix).kl_normBy_abcd.bed: $(DIR)/jun2018_analysis_high_mapq/deep_region_cov25.kl_normBy_abcd.slop20.bed.fa
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

deep_kinetics_beds := $(wildcard ../../deep_kinetics_region.*.slop20.bed.fa)
occ_deep_kinetics := $(deep_kinetics_beds:../../%=$(occ_prefix).%.bed)
$(occ_deep_kinetics): $(occ_prefix).%.bed: ../../%
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

# pattern occurrence for features (such as DNA modification and high IPD)
occ_feature_fa := $(wildcard ../../*.slop20.fullLength.gff.fa)
occ_feature := $(occ_feature_fa:../../%=$(occ_prefix).%.bed)
$(occ_feature): $(occ_prefix).%.bed: ../../%
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

occ := $(occ_background) $(occ_feature) $(occ_deep_kinetics)
merged_occ := $(occ:=.merged_occ)
$(merged_occ): %.merged_occ: %
	cat $< | ../merge_occurrence.sh > $@ && touch $@

target += $(merged_occ)

# These lines require a lot of RAM
$(motif_ipds): get_motif_ipd.R.log
	touch $@
get_motif_ipd.R.log: $(occ_prefix).deep_kinetics_region.ab.bed.merged.sorted.slop20.bed.fa.bed.merged_occ $(occ_prefix).deep_kinetics_region.cd.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.v2.R $(PATTERN) &> $@ && touch $@
$(motif_wga_ipdratios): get_motif_wga_ipdratio.R.log
	touch $@
get_motif_wga_ipdratio.R.log: $(occ_prefix).deep_kinetics_region.ab.bed.merged.sorted.slop20.bed.fa.bed.merged_occ $(occ_prefix).deep_kinetics_region.cd.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_wga_ipdratio.R $(PATTERN) &> $@ && touch $@
motif_ipd.PD2182.csv: $(occ_prefix).deep_kinetics_region.PD2182.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.PD2182.v2.R $(PATTERN) &> get_motif_ipd.PD2182.R.log && touch $@
motif_ipd.PD2182sequel.csv: $(occ_prefix).deep_kinetics_region.PD2182sequel.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.PD2182sequel.v2.R $(PATTERN) &> get_motif_ipd.PD2182sequel.R.log && touch $@
$(motif_ipdratios): get_motif_ipdratio.R.log
	touch $@
get_motif_ipdratio.R.log: $(occ_prefix).deep_kinetics_region.k_normBy_ab.bed.slop20.bed.fa.bed.merged_occ $(occ_prefix).deep_kinetics_region.l_normBy_cd.bed.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.v2.R $(PATTERN) &> $@ && touch $@

motif_ipdratio.PD2182.csv: $(occ_prefix).deep_kinetics_region.PD2182.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.PD2182.R $(PATTERN) &> get_motif_ipdratio.PD2182.R.log && touch $@
motif_ipdratio.PD2182sequel.csv: $(occ_prefix).deep_kinetics_region.PD2182sequel.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.PD2182sequel.R $(PATTERN) &> get_motif_ipdratio.PD2182sequel.R.log && touch $@

motif_ipdratio.k.csv: $(occ_prefix).deep_kinetics_region.k.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.k.R $(PATTERN) &> get_motif_ipdratio.k.R.log && touch $@
motif_ipdratio.l.csv: $(occ_prefix).deep_kinetics_region.l.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.l.R $(PATTERN) &> get_motif_ipdratio.l.R.log && touch $@

motif_ipd.k.csv: $(occ_prefix).deep_kinetics_region.k.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.k.v2.R $(PATTERN) &> get_motif_ipd.k.R.log && touch $@
motif_ipd.l.csv: $(occ_prefix).deep_kinetics_region.l.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.l.v2.R $(PATTERN) &> get_motif_ipd.l.R.log && touch $@

motif_ipd.abcd.csv: $(occ_prefix).deep_kinetics_region.abcd.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.abcd.v2.R $(PATTERN) &> get_motif_ipd.abcd.R.log && touch $@
motif_ipd.kl.csv: $(occ_prefix).deep_kinetics_region.kl.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipd.kl.v2.R $(PATTERN) &> get_motif_ipd.kl.R.log && touch $@
motif_ipdratio.abcd.csv: $(occ_prefix).deep_kinetics_region.abcd.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.abcd.R $(PATTERN) &> get_motif_ipdratio.abcd.R.log && touch $@
motif_ipdratio.kl.csv: $(occ_prefix).deep_kinetics_region.kl.bed.merged.sorted.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.kl.R $(PATTERN) &> get_motif_ipdratio.kl.R.log && touch $@
motif_ipdratio.kl_normBy_abcd.csv: $(occ_prefix).deep_kinetics_region.kl_normBy_abcd.bed.slop20.bed.fa.bed.merged_occ
	Rscript ../get_motif_ipdratio.kl_normBy_abcd.R $(PATTERN) &> get_motif_ipdratio.kl_normBy_abcd.R.log && touch $@

$(motif_kinetics_csv_celegans): %.c_elegans.csv: %.csv
	sample=$$(echo $< | sed -E "s/motif[_0-9a-zA-Z]+\.([_0-9a-zA-Z]+)\.csv/\1/"); \
	Rscript ../extract_celegans_data.R $< $(occ_prefix).deep_kinetics_region.$${sample}.bed*.slop20.bed.fa.bed.merged_occ $@

$(motif_kinetics_csv_ecoli): %.e_coli.csv: %.csv
	sample=$$(echo $< | sed -E "s/motif[_0-9a-zA-Z]+\.([_0-9a-zA-Z]+)\.csv/\1/"); \
	Rscript ../extract_ecoli_data.R $< $(occ_prefix).deep_kinetics_region.$${sample}.bed*.slop20.bed.fa.bed.merged_occ $@

$(motif_summary_target): %.summary.csv: %.csv
	Rscript ../summarize_motif_kinetics.R $< $@

$(motif_kinetics_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.R $(PATTERN) && touch $@

$(motif_kinetics_stderr_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.stderr.R $(PATTERN) && touch $@

$(motif_kinetics_stderr_widey_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.stderr.wide_y.R $(PATTERN) && touch $@

$(cor_pdf): $(motif_kinetics_csv)
	Rscript ../plot_kinetics_correlation.R

$(cor_pdf_celegans): $(motif_kinetics_csv_celegans)
	Rscript ../plot_kinetics_correlation.c_elegans.R

$(cor_pdf2): $(motif_wga_ipdratios) $(motif_hawaiian_ipdratios)
	Rscript ../plot_ipdratio_correlation.R

$(cor_pdf_celegans_ecoli): $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli)
	Rscript ../plot_kinetics_correlation.celegans_ecoli.R

$(cor_pdf_celegans_ecoli_subset1): $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli)
	Rscript ../plot_kinetics_correlation.celegans_ecoli_subset1.R

$(cor_pdf_celegans_ecoli_subset1_positive_strand): $(motif_kinetics_csv_celegans) $(motif_kinetics_csv_ecoli)
	Rscript ../plot_kinetics_correlation.celegans_ecoli_subset1_positive_strand.R

$(motif_modelPrediction): motif_modelPrediction.%.csv: motif_ipd.%.csv motif_ipdratio.%.csv
	Rscript ../get_motif_modelPrediction_from_ipdRatio.R $< $(word 2,$^) $@ && touch $@
$(motif_modelPrediction_celegans): motif_modelPrediction.%.c_elegans.csv: motif_ipd.%.c_elegans.csv motif_ipdratio.%.c_elegans.csv
	Rscript ../get_motif_modelPrediction_from_ipdRatio.R $< $(word 2,$^) $@ && touch $@
$(motif_modelPrediction_ecoli): motif_modelPrediction.%.e_coli.csv: motif_ipd.%.e_coli.csv motif_ipdratio.%.e_coli.csv
	Rscript ../get_motif_modelPrediction_from_ipdRatio.R $< $(word 2,$^) $@ && touch $@
