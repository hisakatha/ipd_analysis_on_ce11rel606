SHELL := bash
.DELETE_ON_ERROR:
.SECONGEXPANSION:

motif_dirs := ACATMTGG ACGCRTG AGGCTT_4 AGGY_7 ATCAGCTG ATWTTB CGRAACCCG CTGDAR CTGRAA GAGD GATC GCGC GGN_10 GGN_4 GGTCTCGC GKCGGY GTAS HKGAGCA HRGTA RGTA RGTAB RTAB TACTGTAG TATAA TTTT YGGWR
tandem_dirs := tandem_repeat_s4_01down tandem_repeat_s4_01up tandem_repeat_s4_02down tandem_repeat_s4_02up tandem_repeat_s4_03down tandem_repeat_s4_03up \
	tandem_repeat_s4_04down tandem_repeat_s4_04up tandem_repeat_s4_05down tandem_repeat_s4_05up
all_dirs := $(motif_dirs) $(tandem_dirs)

all: $(all_dirs)
motif: $(motif_dirs)
tandem: $(tandem_dirs)

.PHONY: all motif tandem $(all_dirs)

$(motif_dirs):
	ln -sf ../Makefile.motif $@/Makefile
	$(MAKE) -C $@ occ_prefix=$$(cat $@/PATTERN)

$(tandem_dirs):
	ln -sf ../Makefile.motif $@/Makefile
	$(MAKE) -C $@
