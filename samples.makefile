# BASE:=$(basename $1 .bam)
BASE := mapped.alignmentset
INPUTS := $(wildcard output/tasks/pbalign.tasks.pbalign-*/$(BASE).bam)
COV_THRES ?= 25

MERGED := $(BASE).merged.bam
BAI := $(MERGED).bai
FDEPTH := $(MERGED).forward_depth
RDEPTH := $(MERGED).reverse_depth
DEPTH_PDF := bam_depth.pdf
BEDGRAPH := $(MERGED).coverage.bedgraph
FWD_BEDGRAPH := $(MERGED).forward_depth.bedgraph
REV_BEDGRAPH := $(MERGED).reverse_depth.bedgraph

FWD_DEEP_REGION := $(FWD_BEDGRAPH).cov$(COV_THRES)
REV_DEEP_REGION := $(REV_BEDGRAPH).cov$(COV_THRES)
DEEP_REGION := $(MERGED).cov$(COV_THRES)
DEEP_REGION_FA := $(DEEP_REGION).fa
SLOP := 20
DEEP_REGION_SLOP := $(DEEP_REGION).slop$(SLOP).bed
DEEP_REGION_SLOP_FA := $(DEEP_REGION_SLOP).fa
DEEP_REGION_SLOP_MERGED := $(DEEP_REGION_SLOP).merged
DEEP_REGION_SLOP_MERGED_FA := $(DEEP_REGION_SLOP_MERGED).fa

# These files are useless
#DEEP_REGION_MERGED := $(DEEP_REGION).merged.bed
#DEEP_REGION_MERGED_FA := $(DEEP_REGION_MERGED).fa
#DEEP_REGION_OVERLAP_MERGED := $(DEEP_REGION).overlap41_merged.bed
#DEEP_REGION_OVERLAP_MERGED_FA := $(DEEP_REGION_OVERLAP_MERGED).fa

targets_DEEP := $(FWD_DEEP_REGION) $(REV_DEEP_REGION) $(DEEP_REGION) $(DEEP_REGION_FA) $(DEEP_REGION_SLOP_FA) $(DEEP_REGION_SLOP_MERGED_FA)

BIGWIG := $(BEDGRAPH).ce11.bigWig
BIGWIG_ALL := $(BEDGRAPH).bigWig
TDF := $(MERGED).tdf
CIRCLIZE := $(MERGED).circlize.pdf

.DELETE_ON_ERROR:

targets := $(MERGED) $(BAI) $(FDEPTH) $(RDEPTH) $(DEPTH_PDF) $(BEDGRAPH) $(FWD_BEDGRAPH) $(REV_BEDGRAPH) $(BIGWIG) $(TDF) $(CIRCLIZE) $(targets_DEEP)
all: $(targets)

$(MERGED): $(INPUTS)
	samtools merge -@ 8 $@ $?

$(BAI): $(MERGED)
	samtools index -@ 4 $?

# SAMT = samtools
## $SAMT view -b -f 0x10 $BASE.merged.bam > $BASE.merged.reverse.bam
## $SAMT view -b -F 0x10 $BASE.merged.bam > $BASE.merged.forward.bam
## $SAMT depth $BASE.merged.reverse.bam > $BASE.merged.reverse.bam.depth
## $SAMT depth $BASE.merged.forward.bam > $BASE.merged.forward.bam.depth

# samtools doesn't output coverage 0; bedtools does.
BEDT := /home/hisakatha/bedtools2/bin/bedtools

# 1-based coverage info
$(FDEPTH): $(MERGED)
	$(BEDT) genomecov -ibam $? -d -strand + > $@
$(RDEPTH): $(MERGED)
	$(BEDT) genomecov -ibam $? -d -strand - > $@

PLOT_SCRIPT := /home/hisakatha/glusterfs/methylation/svr_ipd/plot_bam_depth.R

$(DEPTH_PDF): $(FDEPTH) $(RDEPTH)
	Rscript ${PLOT_SCRIPT} $?

# coverage info in bedgraph format (0-based)
$(BEDGRAPH): $(MERGED)
	$(BEDT) genomecov -ibam $? -bga > $@
$(FWD_BEDGRAPH): $(MERGED)
	$(BEDT) genomecov -ibam $< -bga -strand + > $@
$(REV_BEDGRAPH): $(MERGED)
	$(BEDT) genomecov -ibam $< -bga -strand - > $@

$(FWD_DEEP_REGION): $(FWD_BEDGRAPH)
	cat $< | awk -v OFS="\t" '$$4 >= $(COV_THRES){print $$1,$$2,$$3,"deep_region",$$4,"+"}' | $(BEDT) merge -s -c 4,5,6 -o distinct,mean,distinct > $@
$(REV_DEEP_REGION): $(REV_BEDGRAPH)
	cat $< | awk -v OFS="\t" '$$4 >= $(COV_THRES){print $$1,$$2,$$3,"deep_region",$$4,"-"}' | $(BEDT) merge -s -c 4,5,6 -o distinct,mean,distinct > $@

REF := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa
FAI := /glusterfs/hisakatha/ce11rel606/ce11rel606.fa.fai
$(BIGWIG): $(BEDGRAPH)
	cat $< | grep -v "E._coli_REL606" > $<.ce11 && \
	/home/hisakatha/ucsc/bedGraphToBigWig $<.ce11 $(FAI) $@ && \
	rm $<.ce11

$(DEEP_REGION): $(FWD_DEEP_REGION) $(REV_DEEP_REGION)
	cat $^ | $(BEDT) sort -faidx $(FAI) > $@
$(DEEP_REGION_FA): $(DEEP_REGION)
	$(BEDT) getfasta -fi $(REF) -bed $< > $@

$(DEEP_REGION_SLOP): $(DEEP_REGION)
	cat $^ | $(BEDT) slop -b $(SLOP) -g $(FAI) > $@
$(DEEP_REGION_SLOP_FA): $(DEEP_REGION_SLOP)
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

# Need to sort BEDGRAPH
$(BIGWIG_ALL): $(BEDGRAPH)
	/home/hisakatha/ucsc/bedGraphToBigWig $< $(FAI) $@

IGVT := /home/hisakatha/IGVTools/igvtools
$(TDF): $(MERGED)
	$(IGVT) count -z 9 -w 1 $< $@ $(FAI)

CIRCLIZE_SCRIPT := /glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606/circlize_depth.R
$(CIRCLIZE): $(FDEPTH) $(RDEPTH)
	Rscript $(CIRCLIZE_SCRIPT)

