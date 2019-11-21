SHELL := /bin/bash

REGION_SH := get_enriched_gene.name_to_region.sh
REGION_GFF := get_enriched_gene.both.summarize.gene_region.gff
IGV_SCRIPT := get_enriched_gene.both.summarize.gene_region.igv_script
SNAPSHOT_DIR := igv_snapshot

ENRICH_BOTH_CSV := get_enriched_gene.both.summarize.csv
ENRICH_BOTH_NAMES := get_enriched_gene.both.summarize.gene_names

FA := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa
FAI := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa.fai

#EXTREME_IPD := $(wildcard extreme_ipd.*.gff)
EXTREME_IPD := $(shell ls extreme_ipd.*.gff | grep -E 'extreme_ipd\.[^\.]+\.[^\.]+(\.coverage[0-9]+)?.gff')
EXTREME_IPD_SORTED_GFF := $(EXTREME_IPD:.gff=.sorted.gff)
SLOP := 20
EXTREME_IPD_SLOP := $(EXTREME_IPD:.gff=.slop$(SLOP).gff)
EXTREME_IPD_FA := $(EXTREME_IPD_SLOP:=.fa)
EXTREME_IPD_FULL := $(EXTREME_IPD_SLOP:.gff=.fullLength.gff)
EXTREME_IPD_FULL_FA := $(EXTREME_IPD_FULL:=.fa)
EXTREME_IPD_FULL_MERGED := $(EXTREME_IPD_FULL:=.merged)
EXTREME_IPD_FULL_MERGED_FA := $(EXTREME_IPD_FULL_MERGED:=.fa)
EXTREME_IPD_CONTEXT := $(EXTREME_IPD:.gff=.withContext.gff)
#EXTREME_IPD_SORTED := $(EXTREME_IPD_CONTEXT:.gff=.sorted.gff)
#EXTREME_IPD_GZ := $(EXTREME_IPD_SORTED:=.gz)
#EXTREME_IPD_TBI := $(EXTREME_IPD_GZ:=.tbi)
EXTREME_IPD_SORTED_GZ := $(EXTREME_IPD_CONTEXT:.gff=.sorted.gff.gz)
EXTREME_IPD_TBI := $(EXTREME_IPD_SORTED_GZ:=.tbi)

.DELETE_ON_ERROR:
.SECONDEXPANSION:
target = $(TTS_DOWN_MERGED) $(IGV_SCRIPT) $(ENRICH_BOTH_NAMES) $(EXTREME_IPD_TBI) $(EXTREME_IPD_FULL_FA) $(EXTREME_IPD_FULL_MERGED_FA) $(EXTREME_IPD_SORTED_GFF) $(MOD_TBI) $(MOD_FULL_FA) $(MOD_FULL_MERGED_FA) $(DEEP_REGION_SLOP_FA) $(DEEP_REGION_SLOP_MERGED_FA) $(MOD_FULL_CONTEXT)
all: $$(target)
.PHONY: all

$(REGION_GFF): $(REGION_SH) $(ENRICH_BOTH_CSV)
	./$< > $@

$(IGV_SCRIPT): $(REGION_GFF)
	bedtools igv -slop 500 -img svg -i $< > $@
	echo "Call igv_script at $(SNAPSHOT_DIR) in a IGV session to take snapshots"

$(ENRICH_BOTH_CSV): get_enriched_gene.both.csv get_enriched_gene.both.tts_down500.csv
	Rscript get_enriched_gene.both.summarize.R

$(ENRICH_BOTH_NAMES): $(ENRICH_BOTH_CSV)
	cat $< | tail -n +2 | cut -f1 -d, | uniq > $@

$(EXTREME_IPD_SORTED_GFF): %.sorted.gff: %.gff
	bedtools sort -header -faidx $(FAI) -i $< > $@

$(EXTREME_IPD_SLOP): %.slop$(SLOP).gff: %.gff
	bedtools slop -header -i $< -g $(FAI) -b $(SLOP) > $@

$(EXTREME_IPD_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@

$(EXTREME_IPD_FULL): %.fullLength.gff: %.gff
	cat $< | grep -v '^#' | awk -v EXPECTED=$$(($(SLOP) * 2 + 1)) '{if($$5 - $$4 + 1 == EXPECTED){print}}' > $@

$(EXTREME_IPD_FULL_FA): %.fullLength.gff.fa: %.fullLength.gff %.gff.fa
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@
	@cmp $@ $(word 2,$^) || echo "WARNING: There are partial-length regions in $<. You have to use $@ for MEME." >&2

$(EXTREME_IPD_FULL_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 3,6,7 -o distinct,mean,distinct > $@

$(EXTREME_IPD_FULL_MERGED_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@

# bedtools getfasta (v2.27.0) returns one-line sequence per region,
# so you can merge a gff and a fasta as below.
$(EXTREME_IPD_CONTEXT): %.withContext.gff: %.gff %.slop$(SLOP).gff.fa
	$(eval p1 = $<)
	$(eval p2 = $(word 2, $^))
	(grep '^#' $(p1); paste -d "\0" <(grep -v '^#' $(p1) | sed -e "s/$$/;context=/") <(grep -v '>' $(p2))) > $@

#$(EXTREME_IPD_SORTED): %.sorted.gff: %.gff
#	(grep ^"#" $<; grep -v ^"#" $< | grep -v "^$$" | grep "\t" | sort -k1,1 -k4,4n) > $@

#$(EXTREME_IPD_GZ): %.gz: %
#	bgzip $<

$(EXTREME_IPD_SORTED_GZ): %.sorted.gff.gz: %.gff
	(grep ^"#" $<; grep -v ^"#" $< | grep -v "^$$" | grep "\t" | sort -k1,1 -k4,4n) | bgzip -c > $@

$(EXTREME_IPD_TBI): %.tbi: %
	tabix -p gff $<

# Rules for modifications.*.gff
MOD := $(shell ls modifications.*.gff | grep -E 'modifications\.[^\.]+\.[^\.]+(\.[ACGT])?\.gff')
MSLOP := 20

MOD_SLOP := $(MOD:.gff=.slop$(MSLOP).gff)
$(MOD_SLOP): %.slop$(MSLOP).gff: %.gff
	bedtools slop -header -i $< -g $(FAI) -b $(MSLOP) > $@

MOD_FA := $(MOD_SLOP:=.fa)
$(MOD_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@

MOD_FULL := $(MOD_SLOP:.gff=.fullLength.gff)
$(MOD_FULL): %.fullLength.gff: %.gff
	cat $< | grep -v '^#' | awk -v EXPECTED=$$(($(MSLOP) * 2 + 1)) '{if($$5 - $$4 + 1 == EXPECTED){print}}' > $@

MOD_FULL_FA := $(MOD_FULL:=.fa)
$(MOD_FULL_FA): %.fullLength.gff.fa: %.fullLength.gff %.gff.fa
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@
	@cmp $@ $(word 2,$^) || echo "WARNING: There are partial-length regions in $<. You have to use $@ for MEME." >&2

MOD_FULL_MERGED := $(MOD_FULL:=.merged)
$(MOD_FULL_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 3,6,7 -o distinct,mean,distinct > $@

MOD_FULL_MERGED_FA := $(MOD_FULL_MERGED:=.fa)
$(MOD_FULL_MERGED_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@

MOD_CONTEXT := $(MOD:.gff=.withContext.gff)
$(MOD_CONTEXT): %.withContext.gff: %.gff %.slop$(MSLOP).gff.fa
	$(eval p1 = $<)
	$(eval p2 = $(word 2, $^))
	(grep '^#' $(p1); paste -d "\0" <(grep -v '^#' $(p1) | sed -e "s/$$/;context=/") <(grep -v '>' $(p2))) > $@

MOD_FULL_CONTEXT := $(MOD_FULL:.gff=.withContext.gff)
$(MOD_FULL_CONTEXT): %.withContext.gff: %.gff %.gff.fa
	$(eval p1 = $<)
	$(eval p2 = $(word 2, $^))
	(grep '^#' $(p1); paste -d "\0" <(grep -v '^#' $(p1) | sed -e "s/$$/;context=/") <(grep -v '>' $(p2))) > $@

MOD_SORTED_GZ := $(MOD_CONTEXT:.gff=.sorted.gff.gz)
$(MOD_SORTED_GZ): %.sorted.gff.gz: %.gff
	(grep ^"#" $<; grep -v ^"#" $< | grep -v "^$$" | grep "\t" | sort -k1,1 -k4,4n) | bgzip -c > $@

MOD_TBI := $(MOD_SORTED_GZ:=.tbi)
$(MOD_TBI): %.tbi: %
	tabix -p gff $<

## Deep regions
DEEP_PREFIX := /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606
DEEP_BED_NAME := mapped.alignmentset.merged.bam.cov25
DSLOP := $(MSLOP)
deep_k_ab := deep_region_cov25.k_normBy_ab.bed
deep_l_cd := deep_region_cov25.l_normBy_cd.bed
DEEP_REGION := $(deep_k_ab) $(deep_l_cd)
$(deep_k_ab): $(DEEP_PREFIX)/jun2018_ab_pcr_vc2010_op50/$(DEEP_BED_NAME) $(DEEP_PREFIX)/jun2018_k_nopcr_vc2010_op50/$(DEEP_BED_NAME)
	bedtools intersect -a $< -b $(word 2,$^) -s -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@

$(deep_l_cd): $(DEEP_PREFIX)/jun2018_cd_pcr_vc2010/$(DEEP_BED_NAME) $(DEEP_PREFIX)/jun2018_l_nopcr_vc2010/$(DEEP_BED_NAME)
	bedtools intersect -a $< -b $(word 2,$^) -s -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@

# Note that fasta headers chr:start-end are BED-format, that is, 0-based inclusive start and 0-based exclusive end.
DEEP_REGION_FA := $(DEEP_REGION:=.fa)
$(DEEP_REGION_FA): %.fa: %
	bedtools getfasta -fi $(FA) -bed $< > $@

# These files were found to be useless for motif counting
#DEEP_REGION_MERGED := $(DEEP_REGION:.bed=.merged.bed)
#DEEP_REGION_MERGED_FA := $(DEEP_REGION_MERGED:=.fa)
#DEEP_REGION_OVERLAP_MERGED := $(DEEP_REGION:.bed=.overlap41_merged.bed)
#DEEP_REGION_OVERLAP_MERGED_FA := $(DEEP_REGION_OVERLAP_MERGED:=.fa)
#$(DEEP_REGION_MERGED): %.merged.bed: %.bed
#	bedtools merge -i $< -c 4,5,6 -o first,min,distinct | sed -E 's/[-\+],[-\+]/./' > $@
#$(DEEP_REGION_MERGED_FA): %.fa: %
#	bedtools getfasta -fi $(FA) -bed $< > $@
# Specify 41 bp overlap
#$(DEEP_REGION_OVERLAP_MERGED): %.overlap41_merged.bed: %.bed
#	bedtools merge -i $< -c 4,5,6 -o first,min,distinct -d -41 | sed -E 's/[-\+],[-\+]/./' > $@
#$(DEEP_REGION_OVERLAP_MERGED_FA): %.fa: %
#	bedtools getfasta -fi $(FA) -bed $< > $@

DEEP_REGION_SLOP := $(DEEP_REGION:.bed=.slop$(DSLOP).bed)
$(DEEP_REGION_SLOP): %.slop$(DSLOP).bed: %.bed
	bedtools slop -i $< -g $(FAI) -b $(DSLOP) > $@
DEEP_REGION_SLOP_FA := $(DEEP_REGION_SLOP:=.fa)
$(DEEP_REGION_SLOP_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) > $@

# For counting the size of the regions
DEEP_REGION_SLOP_MERGED := $(DEEP_REGION_SLOP:=.merged)
$(DEEP_REGION_SLOP_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > $@
DEEP_REGION_SLOP_MERGED_FA := $(DEEP_REGION_SLOP_MERGED:=.fa)
$(DEEP_REGION_SLOP_MERGED_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@

## Shuffled regions in deep regions
SHUF_BIN := /home/hisakatha/bedtools_util/shuffle_with_strand.sh
MOD_SHUF := $(MOD:.gff=.shuffled.gff)
MOD_SHUF_in_k_normBy_ab := $(filter modifications.k_normBy_ab.%,$(MOD_SHUF))
$(MOD_SHUF_in_k_normBy_ab): %.shuffled.gff: %.gff $(deep_k_ab)
	$(SHUF_BIN) $(FAI) $(word 2,$^) $< | bedtools sort -faidx $(FAI) > $@

MOD_SHUF_in_l_normBy_cd := $(filter modifications.l_normBy_cd.%,$(MOD_SHUF))
#MOD_SHUF_in_l_normBy_cd := $(shell echo $(MOD_SHUF) | xargs -n1 | grep l_normBy_cd)
$(MOD_SHUF_in_l_normBy_cd): %.shuffled.gff: %.gff $(deep_l_cd)
	$(SHUF_BIN) $(FAI) $(word 2,$^) $< | bedtools sort -faidx $(FAI) > $@

MOD_SHUF_target := $(MOD_SHUF_in_k_normBy_ab) $(MOD_SHUF_in_l_normBy_cd)

MOD_SHUF_SLOP := $(MOD_SHUF_target:.gff=.slop$(MSLOP).gff)
$(MOD_SHUF_SLOP): %.slop$(MSLOP).gff: %.gff
	bedtools slop -header -i $< -g $(FAI) -b $(MSLOP) > $@

MOD_SHUF_SLOP_FULL := $(MOD_SHUF_SLOP:.gff=.fullLength.gff)
$(MOD_SHUF_SLOP_FULL): %.fullLength.gff: %.gff
	cat $< | grep -v '^#' | awk -v EXPECTED=$$(($(MSLOP) * 2 + 1)) '{if($$5 - $$4 + 1 == EXPECTED){print}}' > $@
MOD_SHUF_SLOP_FULL_FA := $(MOD_SHUF_SLOP_FULL:=.fa)
$(MOD_SHUF_SLOP_FULL_FA): %.fa: %
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@
target += $(MOD_SHUF_SLOP_FULL_FA)

# For counting the size of the regions
MOD_SHUF_SLOP_FULL_MERGED := $(MOD_SHUF_SLOP_FULL:=.merged)
$(MOD_SHUF_SLOP_FULL_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 3,6,7 -o distinct,mean,distinct > $@
MOD_SHUF_SLOP_FULL_MERGED_FA := $(MOD_SHUF_SLOP_FULL_MERGED:=.fa)
$(MOD_SHUF_SLOP_FULL_MERGED_FA): %.fa: %
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@
target += $(MOD_SHUF_SLOP_FULL_MERGED_FA)

## Deep regions in kinetics data
DEEP_KINETICS := deep_kinetics_region.ab.bed deep_kinetics_region.cd.bed deep_kinetics_region.k.bed deep_kinetics_region.l.bed \
	deep_kinetics_region.abcd.bed deep_kinetics_region.kl.bed \
	deep_kinetics_region.PD2182.bed deep_kinetics_region.PD2182sequel.bed
$(DEEP_KINETICS): convert_kinetics_deep_region_into_bed.R.log
convert_kinetics_deep_region_into_bed.R.log:
	Rscript convert_kinetics_deep_region_into_bed.R >& $@ && touch $@
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_MERGED_SORTED := $(DEEP_KINETICS:=.merged.sorted)
$(DEEP_KINETICS_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@
DEEP_KINETICS2 := deep_kinetics_region.k_normBy_ab.bed deep_kinetics_region.l_normBy_cd.bed deep_kinetics_region.kl_normBy_abcd.bed
deep_kinetics_region.k_normBy_ab.bed: deep_kinetics_region.k.bed.merged.sorted deep_kinetics_region.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@
deep_kinetics_region.l_normBy_cd.bed: deep_kinetics_region.l.bed.merged.sorted deep_kinetics_region.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@
deep_kinetics_region.kl_normBy_abcd.bed: deep_kinetics_region.kl.bed.merged.sorted deep_kinetics_region.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@
target += $(DEEP_KINETICS_MERGED_SORTED) $(DEEP_KINETICS2)

## Debug
print-%  : ; @echo $* = $($*)

