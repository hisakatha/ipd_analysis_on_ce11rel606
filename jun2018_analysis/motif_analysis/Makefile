SHELL := bash

#PATTERN := HRGTA
PATTERN := $(shell cat PATTERN)
ifndef PATTERN
$(error PATTERN is undefined)
endif

motif_ipds := motif_ipd.ab.csv motif_ipd.cd.csv
motif_ipdratios := motif_ipdratio.k_normBy_ab.csv motif_ipdratio.l_normBy_cd.csv
motif_kinetics_csv := $(motif_ipds) motif_ipd.PD2182.csv motif_ipd.PD2182sequel.csv $(motif_ipdratios)
motif_kinetics_pdf := plot_motif_kinetics.pdf
motif_kinetics_stderr_pdf := plot_motif_kinetics.stderr.pdf
motif_kinetics_stderr_widey_pdf := plot_motif_kinetics.stderr.wide_y.pdf
cor_pdf := plot_kinetics_correlation.pdf

.SECONDEXPANSION:
target := $(motif_kinetics_pdf) $(motif_kinetics_stderr_pdf) $(motif_kinetics_stderr_widey_pdf) $(cor_pdf)
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

# pattern occurrence for background
occ_background := $(PATTERN).ab.bed $(PATTERN).cd.bed $(PATTERN).ce11rel606.bed $(PATTERN).k_normBy_ab.bed $(PATTERN).l_normBy_cd.bed $(PATTERN).PD2182.bed $(PATTERN).PD2182sequel.bed
$(PATTERN).ce11rel606.bed: /glusterfs/hisakatha/ce11rel606/ce11rel606.fa
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
DIR := /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606
FILE := mapped.alignmentset.merged.bam.cov25.slop20.bed.fa
$(PATTERN).ab.bed: $(DIR)/jun2018_ab_pcr_vc2010_op50/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(PATTERN).cd.bed: $(DIR)/jun2018_cd_pcr_vc2010/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(PATTERN).PD2182.bed: $(DIR)/PD2182-20150119/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(PATTERN).PD2182sequel.bed: $(DIR)/PD2182-sequel/$(FILE)
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(PATTERN).k_normBy_ab.bed: $(DIR)/jun2018_analysis/deep_region_cov25.k_normBy_ab.slop20.bed.fa
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@
$(PATTERN).l_normBy_cd.bed: $(DIR)/jun2018_analysis/deep_region_cov25.l_normBy_cd.slop20.bed.fa
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

deep_kinetics_beds := $(wildcard ../../deep_kinetics_region.*.slop20.bed.fa)
occ_deep_kinetics := $(deep_kinetics_beds:../../%=$(PATTERN).%.bed)
$(occ_deep_kinetics): $(PATTERN).%.bed: ../../%
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

# pattern occurrence for features (such as DNA modification and high IPD)
occ_feature_fa := $(wildcard ../../*.slop20.fullLength.gff.fa)
occ_feature := $(occ_feature_fa:../../%=$(PATTERN).%.bed)
$(occ_feature): $(PATTERN).%.bed: ../../%
	seqkit locate --bed -i -d -p $(PATTERN) $< > $@ && touch $@

occ := $(occ_background) $(occ_feature) $(occ_deep_kinetics)
merged_occ := $(occ:=.merged_occ)
$(merged_occ): %.merged_occ: %
	cat $< | ../merge_occurrence.sh > $@ && touch $@

target += $(merged_occ)

# This line requires a lot of RAM
$(motif_ipds): get_motif_ipd.R.log
	touch $@
get_motif_ipd.R.log: $(PATTERN).ab.bed.merged_occ $(PATTERN).cd.bed.merged_occ
	Rscript ../get_motif_ipd.v2.R $(PATTERN) &> $@ && touch $@
motif_ipd.PD2182.csv: $(PATTERN).PD2182.bed.merged_occ
	Rscript ../get_motif_ipd.PD2182.v2.R $(PATTERN) &> get_motif_ipd.PD2182.R.log && touch $@
motif_ipd.PD2182sequel.csv: $(PATTERN).PD2182sequel.bed.merged_occ
	Rscript ../get_motif_ipd.PD2182sequel.v2.R $(PATTERN) &> get_motif_ipd.PD2182sequel.R.log && touch $@

# This line requires a lot of RAM
$(motif_ipdratios): get_motif_ipdratio.R.log
	touch $@
get_motif_ipdratio.R.log: $(PATTERN).k_normBy_ab.bed.merged_occ $(PATTERN).l_normBy_cd.bed.merged_occ
	Rscript ../get_motif_ipdratio.v2.R $(PATTERN) &> $@ && touch $@

$(motif_kinetics_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.R $(PATTERN) && touch $@

$(motif_kinetics_stderr_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.stderr.R $(PATTERN) && touch $@

$(motif_kinetics_stderr_widey_pdf): $(motif_kinetics_csv)
	Rscript ../plot_motif_kinetics.stderr.wide_y.R $(PATTERN) && touch $@

$(cor_pdf): $(motif_kinetics_csv)
	Rscript ../plot_kinetics_correlation.R
