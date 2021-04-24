SHELL := bash
# BASE:=$(basename $1 .bam)
BASE := mapped.alignmentset
INPUTS := $(wildcard output/tasks/pbalign.tasks.pbalign-*/$(BASE).bam)
COV_THRES ?= 25
MAPQ_THRES ?= 128

### Variables you have to set
# Relative paths should be based on the sample directories

REF := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa
FAI := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa.fai
pbindex := /bio/package/pacbio/smrtlink/smrtcmds/bin/pbindex
IGVT := /home/hisakatha/IGVTools/igvtools

# `cargo build` in a clone of https://github.com/hisakatha/compare-mapq-between-two-sets
comp_bin := ~/compare_mapq_between_two_sets/target/release/compare_mapq_between_two_sets

# `cargo build` in a clone of https://github.com/hisakatha/extract-bam-unique-alignment
extract_bin := ~/extract_bam_unique_alignment/target/release/extract_bam_unique_alignment

# UCSC Utilities (https://genome.ucsc.edu/util.html)
bedg2bw := /home/hisakatha/ucsc/bedGraphToBigWig

# samtools doesn't output coverage 0; bedtools does.
BEDT := /home/hisakatha/bedtools2/bin/bedtools

###

PLOT_SCRIPT := ../plot_bam_depth.R
MERGED := $(BASE).merged.bam
BAI := $(MERGED).bai
MERGED_MD := $(BASE).merged.md.bam
BAI_MD := $(MERGED_MD).bai
FDEPTH := $(MERGED).forward_depth
RDEPTH := $(MERGED).reverse_depth
DEPTH_PDF := bam_depth.pdf
BEDGRAPH := $(MERGED).coverage.bedgraph
BEDGRAPH_SORTED := $(MERGED).coverage.bedgraph.sorted
FWD_BEDGRAPH := $(MERGED).forward_depth.bedgraph
REV_BEDGRAPH := $(MERGED).reverse_depth.bedgraph
NAME_SORTED := $(MERGED:.bam=.sorted_by_name.bam)
HQ_UNIQUE := $(NAME_SORTED:.bam=.hq_unique.bam)
HQ_UNIQUE_SORTED := $(HQ_UNIQUE:.bam=.sorted.bam)
HQ_UNIQUE_SORTED_BAI := $(HQ_UNIQUE_SORTED:=.bai)
HQ_UNIQUE_SORTED_PBI := $(HQ_UNIQUE_SORTED:=.pbi)
HQ_UNIQUE_SORTED_FDEPTH := $(HQ_UNIQUE_SORTED).forward_depth
HQ_UNIQUE_SORTED_RDEPTH := $(HQ_UNIQUE_SORTED).reverse_depth
MAPQ_COMP_CSV := $(NAME_SORTED:=.compare_celegans_ecoli.csv)
BOTH_HIGH_MAPQ_READ := $(MAPQ_COMP_CSV:.csv=.both_high_mapq_read.csv)

FWD_DEEP_REGION := $(FWD_BEDGRAPH).cov$(COV_THRES)
REV_DEEP_REGION := $(REV_BEDGRAPH).cov$(COV_THRES)
DEEP_REGION := $(MERGED).cov$(COV_THRES)
DEEP_REGION_FA := $(DEEP_REGION).fa
SLOP := 20
DEEP_REGION_SLOP := $(DEEP_REGION).slop$(SLOP).bed
DEEP_REGION_SLOP_FA := $(DEEP_REGION_SLOP).fa
DEEP_REGION_SLOP_MERGED := $(DEEP_REGION_SLOP).merged
DEEP_REGION_SLOP_MERGED_FA := $(DEEP_REGION_SLOP_MERGED).fa

HQ := $(MERGED:.bam=.high_mapq.bam)
HQ_BEDGRAPH := $(HQ).coverage.bedgraph
HQ_BEDGRAPH_SORTED := $(HQ).coverage.bedgraph.sorted
HQ_FDEPTH := $(HQ).forward_depth
HQ_RDEPTH := $(HQ).reverse_depth
HQ_FWD_BEDGRAPH := $(HQ_FDEPTH).bedgraph
HQ_REV_BEDGRAPH := $(HQ_RDEPTH).bedgraph
HQ_FWD_DEEP_REGION := $(HQ_FWD_BEDGRAPH).cov$(COV_THRES)
HQ_REV_DEEP_REGION := $(HQ_REV_BEDGRAPH).cov$(COV_THRES)
HQ_DEEP_REGION := $(HQ).cov$(COV_THRES)
HQ_DEEP_REGION_SLOP := $(HQ_DEEP_REGION).slop$(SLOP).bed
HQ_DEEP_REGION_SLOP_FA := $(HQ_DEEP_REGION_SLOP).fa

# These files are useless
#DEEP_REGION_MERGED := $(DEEP_REGION).merged.bed
#DEEP_REGION_MERGED_FA := $(DEEP_REGION_MERGED).fa
#DEEP_REGION_OVERLAP_MERGED := $(DEEP_REGION).overlap41_merged.bed
#DEEP_REGION_OVERLAP_MERGED_FA := $(DEEP_REGION_OVERLAP_MERGED).fa

targets_DEEP := $(FWD_DEEP_REGION) $(REV_DEEP_REGION) $(DEEP_REGION) $(DEEP_REGION_FA) $(DEEP_REGION_SLOP_FA) $(DEEP_REGION_SLOP_MERGED_FA)

BIGWIG := $(BEDGRAPH).ce11.bigWig
BIGWIG_ALL := $(BEDGRAPH_SORTED).bigWig
HQ_BIGWIG_ALL := $(HQ_BEDGRAPH_SORTED).bigWig
TDF := $(MERGED).tdf
CIRCLIZE := $(MERGED).circlize.pdf

DIMER_IPD := dimer_ipd_stats.csv
TETRAMER_IPD := tetramer_ipd_stats.csv

ALIGNMENT_LEN := $(MERGED).alignment_length
MERGED_STATS := $(MERGED).stats
HQ_ALIGNMENT_LEN := $(HQ).alignment_length
HQ_MERGED_STATS := $(HQ).stats

.DELETE_ON_ERROR:
.SECONDEXPANSION:

targets := $(MERGED) $(BAI) $(BAI_MD) $(FDEPTH) $(RDEPTH) $(DEPTH_PDF) $(BEDGRAPH) $(FWD_BEDGRAPH) $(REV_BEDGRAPH) $(BIGWIG) $(TDF) $(CIRCLIZE) $(targets_DEEP) $(DIMER_IPD) $(BIGWIG_ALL) $(HQ_BIGWIG_ALL)
targets += $(ALIGNMENT_LEN) $(MERGED_STATS)
targets += $(HQ_ALIGNMENT_LEN) $(HQ_MERGED_STATS) $(HQ_FDEPTH) $(HQ_RDEPTH) $(HQ_DEEP_REGION_SLOP_FA)
targets += $(TETRAMER_IPD)
all: $$(targets)

mapq_comp_csv: $(MAPQ_COMP_CSV) $(BOTH_HIGH_MAPQ_READ)
hq_unique: mapq_comp_csv $(HQ_UNIQUE_SORTED_FDEPTH) $(HQ_UNIQUE_SORTED_RDEPTH) $(HQ_UNIQUE_SORTED_BAI) $(HQ_UNIQUE_SORTED_PBI)

.PHONY: all mapq_comp_csv hq_unique

$(MERGED): $(INPUTS)
	if [[ $(words $^) == 1 ]]; then ln -s $< $@; else samtools merge -@ 8 -c $@ $^; fi

$(BAI): $(MERGED)
	samtools index -@ 4 $<

$(MERGED_MD): $(MERGED)
	samtools calmd -@ 4 -b $< $(REF) > $@

$(BAI_MD): $(MERGED_MD)
	samtools index -@ 4 $<

$(NAME_SORTED): $(MERGED)
	samtools sort -n -@ 8 $< > $@

$(MAPQ_COMP_CSV): $(NAME_SORTED)
	$(comp_bin) $< celegans chrI,chrII,chrIII,chrIV,chrV,chrX,chrM ecoli E._coli_REL606 > $@ 2> $(@:.csv=.log)

$(BOTH_HIGH_MAPQ_READ): $(MAPQ_COMP_CSV)
	Rscript ../extract_both_high_mapq.R $< $@

$(HQ_UNIQUE): $(NAME_SORTED)
	$(extract_bin) $< $(MAPQ_THRES) $@

$(HQ_UNIQUE_SORTED): $(HQ_UNIQUE)
	samtools sort -@ 4 $< > $@

$(HQ_UNIQUE_SORTED_BAI): $(HQ_UNIQUE_SORTED)
	samtools index $<

$(HQ_UNIQUE_SORTED_PBI): $(HQ_UNIQUE_SORTED)
	$(pbindex) $<

aln_len_bin := ../get_alignment_length.sh
$(ALIGNMENT_LEN): $(MERGED)
	samtools view -@ 4 $< | $(aln_len_bin) > $@

$(MERGED_STATS): $(MERGED)
	samtools stats -@ 4 $< > $@

# SAMT = samtools
## $SAMT view -b -f 0x10 $BASE.merged.bam > $BASE.merged.reverse.bam
## $SAMT view -b -F 0x10 $BASE.merged.bam > $BASE.merged.forward.bam
## $SAMT depth $BASE.merged.reverse.bam > $BASE.merged.reverse.bam.depth
## $SAMT depth $BASE.merged.forward.bam > $BASE.merged.forward.bam.depth

# 1-based coverage info
$(FDEPTH): $(MERGED)
	$(BEDT) genomecov -ibam $< -d -strand + > $@
$(RDEPTH): $(MERGED)
	$(BEDT) genomecov -ibam $< -d -strand - > $@

$(HQ_FDEPTH): $(HQ)
	$(BEDT) genomecov -ibam $< -d -strand + > $@
$(HQ_RDEPTH): $(HQ)
	$(BEDT) genomecov -ibam $< -d -strand - > $@

$(HQ_UNIQUE_SORTED_FDEPTH): $(HQ_UNIQUE_SORTED)
	$(BEDT) genomecov -ibam $< -d -strand + > $@
$(HQ_UNIQUE_SORTED_RDEPTH): $(HQ_UNIQUE_SORTED)
	$(BEDT) genomecov -ibam $< -d -strand - > $@

$(DEPTH_PDF): $(FDEPTH) $(RDEPTH)
	Rscript ${PLOT_SCRIPT} $^

# coverage info in bedgraph format (0-based)
$(BEDGRAPH): $(MERGED)
	$(BEDT) genomecov -ibam $^ -bga > $@
$(FWD_BEDGRAPH): $(MERGED)
	$(BEDT) genomecov -ibam $< -bga -strand + > $@
$(REV_BEDGRAPH): $(MERGED)
	$(BEDT) genomecov -ibam $< -bga -strand - > $@
$(HQ_FWD_BEDGRAPH): $(HQ)
	$(BEDT) genomecov -ibam $< -bga -strand + > $@
$(HQ_REV_BEDGRAPH): $(HQ)
	$(BEDT) genomecov -ibam $< -bga -strand - > $@

$(FWD_DEEP_REGION): $(FWD_BEDGRAPH)
	cat $< | awk -v OFS="\t" '$$4 >= $(COV_THRES){print $$1,$$2,$$3,"deep_region",$$4,"+"}' | $(BEDT) merge -s -c 4,5,6 -o distinct,mean,distinct > $@
$(REV_DEEP_REGION): $(REV_BEDGRAPH)
	cat $< | awk -v OFS="\t" '$$4 >= $(COV_THRES){print $$1,$$2,$$3,"deep_region",$$4,"-"}' | $(BEDT) merge -s -c 4,5,6 -o distinct,mean,distinct > $@
$(HQ_FWD_DEEP_REGION): $(HQ_FWD_BEDGRAPH)
	cat $< | awk -v OFS="\t" '$$4 >= $(COV_THRES){print $$1,$$2,$$3,"deep_region",$$4,"+"}' | $(BEDT) merge -s -c 4,5,6 -o distinct,mean,distinct > $@
$(HQ_REV_DEEP_REGION): $(HQ_REV_BEDGRAPH)
	cat $< | awk -v OFS="\t" '$$4 >= $(COV_THRES){print $$1,$$2,$$3,"deep_region",$$4,"-"}' | $(BEDT) merge -s -c 4,5,6 -o distinct,mean,distinct > $@

$(BIGWIG): $(BEDGRAPH)
	cat $< | grep -v "E._coli_REL606" > $<.ce11 && \
	$(bedg2bw) $<.ce11 $(FAI) $@ && \
	rm $<.ce11

$(DEEP_REGION): $(FWD_DEEP_REGION) $(REV_DEEP_REGION)
	cat $^ | $(BEDT) sort -faidx $(FAI) > $@
$(HQ_DEEP_REGION): $(HQ_FWD_DEEP_REGION) $(HQ_REV_DEEP_REGION)
	cat $^ | $(BEDT) sort -faidx $(FAI) > $@
$(DEEP_REGION_FA): $(DEEP_REGION)
	$(BEDT) getfasta -fi $(REF) -bed $< > $@

$(DEEP_REGION_SLOP): $(DEEP_REGION)
	cat $^ | $(BEDT) slop -b $(SLOP) -g $(FAI) > $@
$(DEEP_REGION_SLOP_FA): $(DEEP_REGION_SLOP)
	$(BEDT) getfasta -fi $(REF) -bed $< > $@
$(HQ_DEEP_REGION_SLOP): $(HQ_DEEP_REGION)
	cat $^ | $(BEDT) slop -b $(SLOP) -g $(FAI) > $@
$(HQ_DEEP_REGION_SLOP_FA): $(HQ_DEEP_REGION_SLOP)
	$(BEDT) getfasta -fi $(REF) -bed $< > $@

## DEEP_REGION_MERGED is a merged regions regardless of strands
# Note that DEEP_REGION is sorted, and therefore, DEEP_REGION_MERGED is sorted.
#$(DEEP_REGION_MERGED): $(DEEP_REGION)
#	$(BEDT) merge -i $< -c 4,5,6 -o first,min,distinct | sed -E 's/[-\+],[-\+]/./' > $@
#$(DEEP_REGION_MERGED_FA): $(DEEP_REGION_MERGED)
#	$(BEDT) getfasta -fi $(REF) -bed $< > $@

# Specify 41 bp overlap
#$(DEEP_REGION_OVERLAP_MERGED): $(DEEP_REGION)
#	$(BEDT) merge -i $< -c 4,5,6 -o first,min,distinct -d -41 | sed -E 's/[-\+],[-\+]/./' > $@
#$(DEEP_REGION_OVERLAP_MERGED_FA): $(DEEP_REGION_OVERLAP_MERGED)
#	$(BEDT) getfasta -fi $(REF) -bed $< > $@
$(DEEP_REGION_SLOP_MERGED): $(DEEP_REGION_SLOP)
	$(BEDT) merge -i $< -c 4,5,6 -o first,min,distinct -s | sed -E 's/[-\+],[-\+]/./' > $@
$(DEEP_REGION_SLOP_MERGED_FA): $(DEEP_REGION_SLOP_MERGED)
	$(BEDT) getfasta -fi $(REF) -bed $< -s > $@

$(BEDGRAPH_SORTED): $(BEDGRAPH)
	$(BEDT) sort -i $< > $@

$(BIGWIG_ALL): $(BEDGRAPH_SORTED)
	$(bedg2bw) $< $(FAI) $@

$(HQ): $(MERGED)
	samtools view -q 128 -b $< > $@

$(HQ_BEDGRAPH): $(HQ)
	$(BEDT) genomecov -ibam $< -bga > $@

$(HQ_BEDGRAPH_SORTED): $(HQ_BEDGRAPH)
	$(BEDT) sort -i $< > $@

$(HQ_BIGWIG_ALL): $(HQ_BEDGRAPH_SORTED)
	$(bedg2bw) $< $(FAI) $@

$(HQ_ALIGNMENT_LEN): $(HQ)
	samtools view -@ 4 $< | $(aln_len_bin) > $@

$(HQ_MERGED_STATS): $(HQ)
	samtools stats -@ 4 $< > $@

$(TDF): $(MERGED)
	$(IGVT) count -z 9 -w 1 $< $@ $(FAI)

CIRCLIZE_SCRIPT := /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/circlize_depth.R
$(CIRCLIZE): $(FDEPTH) $(RDEPTH)
	Rscript $(CIRCLIZE_SCRIPT)

gathered_ipd := output/tasks/kinetics_tools.tasks.gather_kinetics_h5-1/file.h5
scattered_ipds := $(wildcard output/tasks/kinetics_tools.tasks.ipd_summary-*/basemods.h5)
$(gathered_ipd): $(scattered_ipds)
	if [[ $(words $(scattered_ipds)) == 1 ]]; then mkdir -p $$(dirname $@) && ln -s ../kinetics_tools.tasks.ipd_summary-0/basemods.h5 $@; else echo "ERROR: $@ should have been made by pbsmrtpipe" && exit 1; fi

COLLECT_BIN := /glusterfs/hisakatha/methylation/svr_ipd/collect_ipd_per_kmer/collect_ipd
$(DIMER_IPD): $(gathered_ipd)
	$(COLLECT_BIN) -k 2 -l 20 -t 25 -o $@ $<
$(TETRAMER_IPD): $(gathered_ipd)
	$(COLLECT_BIN) -k 4 -l 20 -t 25 -o $@ $<
