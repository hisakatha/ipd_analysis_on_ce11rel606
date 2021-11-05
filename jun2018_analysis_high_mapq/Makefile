SHELL := /bin/bash

REGION_SH := get_enriched_gene.name_to_region.sh
REGION_GFF := get_enriched_gene.both.summarize.gene_region.gff
IGV_SCRIPT := get_enriched_gene.both.summarize.gene_region.igv_script
SNAPSHOT_DIR := igv_snapshot

ENRICH_BOTH_CSV := get_enriched_gene.both.summarize.csv
ENRICH_BOTH_NAMES := get_enriched_gene.both.summarize.gene_names

### Variables you have to set
FA := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa
FAI := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa.fai

#EXTREME_IPD := $(wildcard extreme_ipd.*.gff)
EXTREME_IPD := $(shell ls extreme_ipd.*.gff | grep -E 'extreme_ipd\.[^\.]+\.[^\.]+(\.coverage[0-9]+)?.gff')
EXTREME_IPD_SORTED_GFF := $(EXTREME_IPD:.gff=.sorted.gff)
EXTREME_IPD_SORTED_Celegans_GFF := $(EXTREME_IPD_SORTED_GFF:.gff=.c_elegans.gff)
EXTREME_IPD_SORTED_Ecoli_GFF := $(EXTREME_IPD_SORTED_GFF:.gff=.e_coli.gff)
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

EXTREME_IPD_FULL_Celegans := $(EXTREME_IPD_FULL:.gff=.c_elegans.gff)
EXTREME_IPD_FULL_Ecoli := $(EXTREME_IPD_FULL:.gff=.e_coli.gff)
EXTREME_IPD_FULL_Celegans_FA := $(EXTREME_IPD_FULL_Celegans:=.fa)
EXTREME_IPD_FULL_Ecoli_FA := $(EXTREME_IPD_FULL_Ecoli:=.fa)

.DELETE_ON_ERROR:
.SECONDEXPANSION:
target = $(EXTREME_IPD_TBI) $(EXTREME_IPD_FULL_FA) $(EXTREME_IPD_FULL_MERGED_FA) $(EXTREME_IPD_SORTED_GFF) $(MOD_TBI) $(MOD_FULL_FA) $(MOD_FULL_MERGED_FA) $(DEEP_REGION_SLOP_FA) $(DEEP_REGION_SLOP_MERGED_FA) $(MOD_FULL_CONTEXT) $(EXTREME_IPD_FULL_Celegans_FA) $(EXTREME_IPD_FULL_Ecoli_FA) $(EXTREME_IPD_SORTED_Celegans_GFF) $(EXTREME_IPD_SORTED_Ecoli_GFF)
#target += $(ENRICH_BOTH_NAMES) $(IGV_SCRIPT)
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

$(EXTREME_IPD_SORTED_Celegans_GFF): %.c_elegans.gff: %.gff
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@

$(EXTREME_IPD_SORTED_Ecoli_GFF): %.e_coli.gff: %.gff
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@

$(EXTREME_IPD_SLOP): %.slop$(SLOP).gff: %.gff
	bedtools slop -header -i $< -g $(FAI) -b $(SLOP) > $@

$(EXTREME_IPD_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@

$(EXTREME_IPD_FULL): %.fullLength.gff: %.gff
	cat $< | grep -v '^#' | awk -v EXPECTED=$$(($(SLOP) * 2 + 1)) '{if($$5 - $$4 + 1 == EXPECTED){print}}' > $@

$(EXTREME_IPD_FULL_FA): %.fullLength.gff.fa: %.fullLength.gff %.gff.fa
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@
	@cmp $@ $(word 2,$^) || echo "WARNING: There are partial-length regions in $<. You have to use $@ for MEME." >&2

EXTREME_IPD_FULL_SORTED_BED := $(EXTREME_IPD_FULL:=.sorted.bed)
$(EXTREME_IPD_FULL_SORTED_BED): %.sorted.bed: %
	bedtools sort -i $< -faidx $(FAI) | awk -v OFS="\t" '{print $$1,$$4-1,$$5,$$3,$$6,$$7}' > $@
target += $(EXTREME_IPD_FULL_SORTED_BED)

$(EXTREME_IPD_FULL_Celegans): %.c_elegans.gff: %.gff
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@

$(EXTREME_IPD_FULL_Celegans_FA): %.fa: %
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@

$(EXTREME_IPD_FULL_Ecoli): %.e_coli.gff: %.gff
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@

$(EXTREME_IPD_FULL_Ecoli_FA): %.fa: %
	cat $< | bedtools getfasta -bed - -fi $(FA) -s > $@

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

#MOD_Celegans := $(MOD:.gff=.c_elegans.gff)
MOD_Celegans := modifications.k_normBy_ab.score20_ipdratio2_coverage25.c_elegans.gff \
	modifications.l_normBy_cd.score20_ipdratio2_coverage25.c_elegans.gff
MOD_Celegans += modifications.kl_normBy_abcd.score20_ipdratio2_coverage25.c_elegans.gff
$(MOD_Celegans): %.c_elegans.gff: %.gff
	(cat $< | grep "^#"; cat $< | grep -v "^#" | grep -v "E._coli" || [[ $$? == 1 ]]) > $@
target += $(MOD_Celegans)

MOD_SLOP := $(MOD:.gff=.slop$(MSLOP).gff)
$(MOD_SLOP): %.slop$(MSLOP).gff: %.gff
	bedtools slop -header -i $< -g $(FAI) -b $(MSLOP) > $@

MOD_FA := $(MOD_SLOP:=.fa)
$(MOD_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@

MOD_FULL := $(MOD_SLOP:.gff=.fullLength.gff)
$(MOD_FULL): %.fullLength.gff: %.gff
	cat $< | grep -v '^#' | awk -v EXPECTED=$$(($(MSLOP) * 2 + 1)) '{if($$5 - $$4 + 1 == EXPECTED){print}}' > $@

MOD_FULL_SORTED_BED := $(MOD_FULL:=.sorted.bed)
$(MOD_FULL_SORTED_BED): %.sorted.bed: %
	bedtools sort -i $< -faidx $(FAI) | awk -v OFS="\t" '{print $$1,$$4-1,$$5,$$3,$$6,$$7}' > $@
target += $(MOD_FULL_SORTED_BED)

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
DEEP_PREFIX := ..
DEEP_BED_NAME := mapped.alignmentset.merged.high_mapq.bam.cov25
DSLOP := $(MSLOP)
deep_k_ab := deep_region_cov25.k_normBy_ab.bed
deep_l_cd := deep_region_cov25.l_normBy_cd.bed
deep_kl_abcd := deep_region_cov25.kl_normBy_abcd.bed
DEEP_REGION := $(deep_k_ab) $(deep_l_cd) $(deep_kl_abcd)
$(deep_k_ab): $(DEEP_PREFIX)/jun2018_ab_pcr_vc2010_op50_no_chunk/$(DEEP_BED_NAME) $(DEEP_PREFIX)/jun2018_k_nopcr_vc2010_op50_no_chunk/$(DEEP_BED_NAME)
	bedtools intersect -a $< -b $(word 2,$^) -s -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@

$(deep_l_cd): $(DEEP_PREFIX)/jun2018_cd_pcr_vc2010_no_chunk/$(DEEP_BED_NAME) $(DEEP_PREFIX)/jun2018_l_nopcr_vc2010_no_chunk/$(DEEP_BED_NAME)
	bedtools intersect -a $< -b $(word 2,$^) -s -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@

$(deep_kl_abcd): $(DEEP_PREFIX)/jun2018_abcd_pcr_no_chunk/$(DEEP_BED_NAME) $(DEEP_PREFIX)/jun2018_kl_nopcr_no_chunk/$(DEEP_BED_NAME)
	bedtools intersect -a $< -b $(word 2,$^) -s -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@

# Note that fasta headers chr:start-end are BED-format, that is, 0-based inclusive start and 0-based exclusive end.
DEEP_REGION_FA := $(DEEP_REGION:=.fa)
$(DEEP_REGION_FA): %.fa: %
	bedtools getfasta -fi $(FA) -bed $< > $@ && touch $@

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
	bedtools slop -i $< -g $(FAI) -b $(DSLOP) > $@ && touch $@
DEEP_REGION_SLOP_FA := $(DEEP_REGION_SLOP:=.fa)
$(DEEP_REGION_SLOP_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) > $@ && touch $@

# For counting the size of the regions
DEEP_REGION_SLOP_MERGED := $(DEEP_REGION_SLOP:=.merged)
$(DEEP_REGION_SLOP_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > $@ && touch $@
DEEP_REGION_SLOP_MERGED_FA := $(DEEP_REGION_SLOP_MERGED:=.fa)
$(DEEP_REGION_SLOP_MERGED_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) -s > $@ && touch $@

## Shuffled regions in deep regions
#SHUF_BIN := /home/hisakatha/bedtools_util/shuffle_with_strand.sh
SHUF_BIN := shuffle_with_strand.sh
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

MOD_SHUF_SLOP_FULL_SORTED_BED := $(MOD_SHUF_SLOP_FULL:=.sorted.bed)
$(MOD_SHUF_SLOP_FULL_SORTED_BED): %.sorted.bed: %
	bedtools sort -i $< -faidx $(FAI) | awk -v OFS="\t" '{print $$1,$$4-1,$$5,$$3,$$6,$$7}' > $@
target += $(MOD_SHUF_SLOP_FULL_SORTED_BED)

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
	deep_kinetics_region.PD2182.bed deep_kinetics_region.PD2182sequel.bed
DEEP_KINETICS += deep_kinetics_region.abcd.bed deep_kinetics_region.kl.bed
$(DEEP_KINETICS)&:
	Rscript convert_kinetics_deep_region_into_bed.R
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_MERGED_SORTED := $(DEEP_KINETICS:=.merged.sorted)
$(DEEP_KINETICS_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS2 := deep_kinetics_region.k_normBy_ab.bed deep_kinetics_region.l_normBy_cd.bed
DEEP_KINETICS2 += deep_kinetics_region.kl_normBy_abcd.bed
deep_kinetics_region.k_normBy_ab.bed: deep_kinetics_region.k.bed.merged.sorted deep_kinetics_region.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region.l_normBy_cd.bed: deep_kinetics_region.l.bed.merged.sorted deep_kinetics_region.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region.kl_normBy_abcd.bed: deep_kinetics_region.kl.bed.merged.sorted deep_kinetics_region.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_BASIC := $(DEEP_KINETICS_MERGED_SORTED) $(DEEP_KINETICS2)

DEEP_KINETICS_Celegans := $(DEEP_KINETICS_BASIC:=.c_elegans.bed)
$(DEEP_KINETICS_Celegans): %.c_elegans.bed: %
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
DEEP_KINETICS_Ecoli := $(DEEP_KINETICS_BASIC:=.e_coli.bed)
$(DEEP_KINETICS_Ecoli): %.e_coli.bed: %
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
target += $(DEEP_KINETICS_Celegans) $(DEEP_KINETICS_Ecoli)

## Deep kinetics region (threshold = 100)
DEEP_KINETICS_COV100 := deep_kinetics_region_cov100.ab.bed deep_kinetics_region_cov100.cd.bed deep_kinetics_region_cov100.k.bed deep_kinetics_region_cov100.l.bed \
	deep_kinetics_region_cov100.PD2182.bed deep_kinetics_region_cov100.PD2182sequel.bed
DEEP_KINETICS_COV100 += deep_kinetics_region_cov100.abcd.bed deep_kinetics_region_cov100.kl.bed
$(DEEP_KINETICS_COV100)&:
	Rscript convert_kinetics_deep_region_into_bed.cov100.R
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_COV100_MERGED_SORTED := $(DEEP_KINETICS_COV100:=.merged.sorted)
$(DEEP_KINETICS_COV100_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS2_COV100 := deep_kinetics_region_cov100.k_normBy_ab.bed deep_kinetics_region_cov100.l_normBy_cd.bed
DEEP_KINETICS2_COV100 += deep_kinetics_region_cov100.kl_normBy_abcd.bed
deep_kinetics_region_cov100.k_normBy_ab.bed: deep_kinetics_region_cov100.k.bed.merged.sorted deep_kinetics_region_cov100.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov100.l_normBy_cd.bed: deep_kinetics_region_cov100.l.bed.merged.sorted deep_kinetics_region_cov100.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov100.kl_normBy_abcd.bed: deep_kinetics_region_cov100.kl.bed.merged.sorted deep_kinetics_region_cov100.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_COV100_BASIC := $(DEEP_KINETICS_COV100_MERGED_SORTED) $(DEEP_KINETICS2_COV100)

DEEP_KINETICS_COV100_Celegans := $(DEEP_KINETICS_COV100_BASIC:=.c_elegans.bed)
$(DEEP_KINETICS_COV100_Celegans): %.c_elegans.bed: %
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
DEEP_KINETICS_COV100_Ecoli := $(DEEP_KINETICS_COV100_BASIC:=.e_coli.bed)
$(DEEP_KINETICS_COV100_Ecoli): %.e_coli.bed: %
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
target += $(DEEP_KINETICS_COV100_Celegans) $(DEEP_KINETICS_COV100_Ecoli)

## Deep kinetics region (threshold = 200)
DEEP_KINETICS_COV200 := deep_kinetics_region_cov200.ab.bed deep_kinetics_region_cov200.cd.bed deep_kinetics_region_cov200.k.bed deep_kinetics_region_cov200.l.bed \
	deep_kinetics_region_cov200.PD2182.bed deep_kinetics_region_cov200.PD2182sequel.bed
DEEP_KINETICS_COV200 += deep_kinetics_region_cov200.abcd.bed deep_kinetics_region_cov200.kl.bed
$(DEEP_KINETICS_COV200)&:
	Rscript convert_kinetics_deep_region_into_bed.cov200.R
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_COV200_MERGED_SORTED := $(DEEP_KINETICS_COV200:=.merged.sorted)
$(DEEP_KINETICS_COV200_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS2_COV200 := deep_kinetics_region_cov200.k_normBy_ab.bed deep_kinetics_region_cov200.l_normBy_cd.bed
DEEP_KINETICS2_COV200 += deep_kinetics_region_cov200.kl_normBy_abcd.bed
deep_kinetics_region_cov200.k_normBy_ab.bed: deep_kinetics_region_cov200.k.bed.merged.sorted deep_kinetics_region_cov200.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov200.l_normBy_cd.bed: deep_kinetics_region_cov200.l.bed.merged.sorted deep_kinetics_region_cov200.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov200.kl_normBy_abcd.bed: deep_kinetics_region_cov200.kl.bed.merged.sorted deep_kinetics_region_cov200.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_COV200_BASIC := $(DEEP_KINETICS_COV200_MERGED_SORTED) $(DEEP_KINETICS2_COV200)

DEEP_KINETICS_COV200_Celegans := $(DEEP_KINETICS_COV200_BASIC:=.c_elegans.bed)
$(DEEP_KINETICS_COV200_Celegans): %.c_elegans.bed: %
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
DEEP_KINETICS_COV200_Ecoli := $(DEEP_KINETICS_COV200_BASIC:=.e_coli.bed)
$(DEEP_KINETICS_COV200_Ecoli): %.e_coli.bed: %
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
target += $(DEEP_KINETICS_COV200_Celegans) $(DEEP_KINETICS_COV200_Ecoli)

## Deep kinetics region (threshold = 300)
DEEP_KINETICS_COV300 := deep_kinetics_region_cov300.ab.bed deep_kinetics_region_cov300.cd.bed deep_kinetics_region_cov300.k.bed deep_kinetics_region_cov300.l.bed \
	deep_kinetics_region_cov300.PD2182.bed deep_kinetics_region_cov300.PD2182sequel.bed
DEEP_KINETICS_COV300 += deep_kinetics_region_cov300.abcd.bed deep_kinetics_region_cov300.kl.bed
$(DEEP_KINETICS_COV300)&:
	Rscript convert_kinetics_deep_region_into_bed.cov300.R
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_COV300_MERGED_SORTED := $(DEEP_KINETICS_COV300:=.merged.sorted)
$(DEEP_KINETICS_COV300_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS2_COV300 := deep_kinetics_region_cov300.k_normBy_ab.bed deep_kinetics_region_cov300.l_normBy_cd.bed
DEEP_KINETICS2_COV300 += deep_kinetics_region_cov300.kl_normBy_abcd.bed
deep_kinetics_region_cov300.k_normBy_ab.bed: deep_kinetics_region_cov300.k.bed.merged.sorted deep_kinetics_region_cov300.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov300.l_normBy_cd.bed: deep_kinetics_region_cov300.l.bed.merged.sorted deep_kinetics_region_cov300.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov300.kl_normBy_abcd.bed: deep_kinetics_region_cov300.kl.bed.merged.sorted deep_kinetics_region_cov300.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_COV300_BASIC := $(DEEP_KINETICS_COV300_MERGED_SORTED) $(DEEP_KINETICS2_COV300)

DEEP_KINETICS_COV300_Celegans := $(DEEP_KINETICS_COV300_BASIC:=.c_elegans.bed)
$(DEEP_KINETICS_COV300_Celegans): %.c_elegans.bed: %
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
DEEP_KINETICS_COV300_Ecoli := $(DEEP_KINETICS_COV300_BASIC:=.e_coli.bed)
$(DEEP_KINETICS_COV300_Ecoli): %.e_coli.bed: %
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
target += $(DEEP_KINETICS_COV300_Celegans) $(DEEP_KINETICS_COV300_Ecoli)

## Deep kinetics region (threshold = 400)
DEEP_KINETICS_COV400 := deep_kinetics_region_cov400.ab.bed deep_kinetics_region_cov400.cd.bed deep_kinetics_region_cov400.k.bed deep_kinetics_region_cov400.l.bed \
	deep_kinetics_region_cov400.PD2182.bed deep_kinetics_region_cov400.PD2182sequel.bed
DEEP_KINETICS_COV400 += deep_kinetics_region_cov400.abcd.bed deep_kinetics_region_cov400.kl.bed
$(DEEP_KINETICS_COV400)&:
	Rscript convert_kinetics_deep_region_into_bed.cov400.R
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_COV400_MERGED_SORTED := $(DEEP_KINETICS_COV400:=.merged.sorted)
$(DEEP_KINETICS_COV400_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS2_COV400 := deep_kinetics_region_cov400.k_normBy_ab.bed deep_kinetics_region_cov400.l_normBy_cd.bed
DEEP_KINETICS2_COV400 += deep_kinetics_region_cov400.kl_normBy_abcd.bed
deep_kinetics_region_cov400.k_normBy_ab.bed: deep_kinetics_region_cov400.k.bed.merged.sorted deep_kinetics_region_cov400.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov400.l_normBy_cd.bed: deep_kinetics_region_cov400.l.bed.merged.sorted deep_kinetics_region_cov400.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov400.kl_normBy_abcd.bed: deep_kinetics_region_cov400.kl.bed.merged.sorted deep_kinetics_region_cov400.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_COV400_BASIC := $(DEEP_KINETICS_COV400_MERGED_SORTED) $(DEEP_KINETICS2_COV400)

DEEP_KINETICS_COV400_Celegans := $(DEEP_KINETICS_COV400_BASIC:=.c_elegans.bed)
$(DEEP_KINETICS_COV400_Celegans): %.c_elegans.bed: %
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
DEEP_KINETICS_COV400_Ecoli := $(DEEP_KINETICS_COV400_BASIC:=.e_coli.bed)
$(DEEP_KINETICS_COV400_Ecoli): %.e_coli.bed: %
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
target += $(DEEP_KINETICS_COV400_Celegans) $(DEEP_KINETICS_COV400_Ecoli)

## Deep kinetics region (threshold = 500)
DEEP_KINETICS_COV500 := deep_kinetics_region_cov500.ab.bed deep_kinetics_region_cov500.cd.bed deep_kinetics_region_cov500.k.bed deep_kinetics_region_cov500.l.bed \
	deep_kinetics_region_cov500.PD2182.bed deep_kinetics_region_cov500.PD2182sequel.bed
DEEP_KINETICS_COV500 += deep_kinetics_region_cov500.abcd.bed deep_kinetics_region_cov500.kl.bed
$(DEEP_KINETICS_COV500)&:
	Rscript convert_kinetics_deep_region_into_bed.cov500.R
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_COV500_MERGED_SORTED := $(DEEP_KINETICS_COV500:=.merged.sorted)
$(DEEP_KINETICS_COV500_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS2_COV500 := deep_kinetics_region_cov500.k_normBy_ab.bed deep_kinetics_region_cov500.l_normBy_cd.bed
DEEP_KINETICS2_COV500 += deep_kinetics_region_cov500.kl_normBy_abcd.bed
deep_kinetics_region_cov500.k_normBy_ab.bed: deep_kinetics_region_cov500.k.bed.merged.sorted deep_kinetics_region_cov500.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov500.l_normBy_cd.bed: deep_kinetics_region_cov500.l.bed.merged.sorted deep_kinetics_region_cov500.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov500.kl_normBy_abcd.bed: deep_kinetics_region_cov500.kl.bed.merged.sorted deep_kinetics_region_cov500.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_COV500_BASIC := $(DEEP_KINETICS_COV500_MERGED_SORTED) $(DEEP_KINETICS2_COV500)

DEEP_KINETICS_COV500_Celegans := $(DEEP_KINETICS_COV500_BASIC:=.c_elegans.bed)
$(DEEP_KINETICS_COV500_Celegans): %.c_elegans.bed: %
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
DEEP_KINETICS_COV500_Ecoli := $(DEEP_KINETICS_COV500_BASIC:=.e_coli.bed)
$(DEEP_KINETICS_COV500_Ecoli): %.e_coli.bed: %
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
target += $(DEEP_KINETICS_COV500_Celegans) $(DEEP_KINETICS_COV500_Ecoli)

## Deep kinetics region (threshold = 600)
DEEP_KINETICS_COV600 := deep_kinetics_region_cov600.ab.bed deep_kinetics_region_cov600.cd.bed deep_kinetics_region_cov600.k.bed deep_kinetics_region_cov600.l.bed \
	deep_kinetics_region_cov600.PD2182.bed deep_kinetics_region_cov600.PD2182sequel.bed
DEEP_KINETICS_COV600 += deep_kinetics_region_cov600.abcd.bed deep_kinetics_region_cov600.kl.bed
$(DEEP_KINETICS_COV600)&:
	Rscript convert_kinetics_deep_region_into_bed.cov600.R
# DEEP_KINETICS are supposed to be sorted
DEEP_KINETICS_COV600_MERGED_SORTED := $(DEEP_KINETICS_COV600:=.merged.sorted)
$(DEEP_KINETICS_COV600_MERGED_SORTED): %.merged.sorted: %
	bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i $< | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS2_COV600 := deep_kinetics_region_cov600.k_normBy_ab.bed deep_kinetics_region_cov600.l_normBy_cd.bed
DEEP_KINETICS2_COV600 += deep_kinetics_region_cov600.kl_normBy_abcd.bed
deep_kinetics_region_cov600.k_normBy_ab.bed: deep_kinetics_region_cov600.k.bed.merged.sorted deep_kinetics_region_cov600.ab.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov600.l_normBy_cd.bed: deep_kinetics_region_cov600.l.bed.merged.sorted deep_kinetics_region_cov600.cd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
deep_kinetics_region_cov600.kl_normBy_abcd.bed: deep_kinetics_region_cov600.kl.bed.merged.sorted deep_kinetics_region_cov600.abcd.bed.merged.sorted
	bedtools intersect -s -a $< -b $(word 2,$^) -sorted -g $(FAI) | bedtools sort -faidx $(FAI) > $@ && touch $@
DEEP_KINETICS_COV600_BASIC := $(DEEP_KINETICS_COV600_MERGED_SORTED) $(DEEP_KINETICS2_COV600)

DEEP_KINETICS_COV600_Celegans := $(DEEP_KINETICS_COV600_BASIC:=.c_elegans.bed)
$(DEEP_KINETICS_COV600_Celegans): %.c_elegans.bed: %
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
DEEP_KINETICS_COV600_Ecoli := $(DEEP_KINETICS_COV600_BASIC:=.e_coli.bed)
$(DEEP_KINETICS_COV600_Ecoli): %.e_coli.bed: %
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) > $@ && touch $@
target += $(DEEP_KINETICS_COV600_Celegans) $(DEEP_KINETICS_COV600_Ecoli)


# Intentionally, these are not merged for motif analysis
DEEP_KINETICS_SLOP := $(DEEP_KINETICS_BASIC:=.slop$(DSLOP).bed) $(DEEP_KINETICS_COV100_BASIC:=.slop$(DSLOP).bed) $(DEEP_KINETICS_COV200_BASIC:=.slop$(DSLOP).bed)
DEEP_KINETICS_SLOP += $(DEEP_KINETICS_COV300_BASIC:=.slop$(DSLOP).bed) $(DEEP_KINETICS_COV400_BASIC:=.slop$(DSLOP).bed) $(DEEP_KINETICS_COV500_BASIC:=.slop$(DSLOP).bed) $(DEEP_KINETICS_COV600_BASIC:=.slop$(DSLOP).bed)
$(DEEP_KINETICS_SLOP): %.slop$(DSLOP).bed: %
	bedtools slop -i $< -g $(FAI) -b $(DSLOP) > $@ && touch $@
DEEP_KINETICS_SLOP_FA := $(DEEP_KINETICS_SLOP:=.fa)
$(DEEP_KINETICS_SLOP_FA): %.fa: %
	bedtools getfasta -bed $< -fi $(FA) > $@ && touch $@
DEEP_KINETICS_SLOP_Celegans_FA := $(DEEP_KINETICS_SLOP:.bed=.c_elegans.bed.fa)
$(DEEP_KINETICS_SLOP_Celegans_FA): %.c_elegans.bed.fa: %.bed
	(cat $< | grep -v "E._coli" || [[ $$? == 1 ]]) | bedtools getfasta -bed - -fi $(FA) > $@ && touch $@
DEEP_KINETICS_SLOP_Ecoli_FA := $(DEEP_KINETICS_SLOP:.bed=.e_coli.bed.fa)
$(DEEP_KINETICS_SLOP_Ecoli_FA): %.e_coli.bed.fa: %.bed
	(cat $< | grep "E._coli" || [[ $$? == 1 ]]) | bedtools getfasta -bed - -fi $(FA) > $@ && touch $@
# For counting the size of the regions
DEEP_KINETICS_SLOP_MERGED := $(DEEP_KINETICS_SLOP:=.merged)
$(DEEP_KINETICS_SLOP_MERGED): %.merged: %
	bedtools sort -i $< -faidx $(FAI) | bedtools merge -s -c 4,5,6 -o distinct,mean,distinct > $@ && touch $@
target += $(DEEP_KINETICS_BASIC) $(DEEP_KINETICS_SLOP_FA) $(DEEP_KINETICS_SLOP_MERGED) $(DEEP_KINETICS_SLOP_Celegans_FA) $(DEEP_KINETICS_SLOP_Ecoli_FA)

## Debug
print-%  : ; @echo $* = $($*)

