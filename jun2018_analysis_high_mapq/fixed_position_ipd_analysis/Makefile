SHELL := bash
.DELETE_ON_ERROR:
.SECONGEXPANSION:

all:

.PHONY: all

fixed_m6A_dirs := greer_pacbio_m6A_cov50 greer_pacbio_orig_m6A_cov50
fixed_m6A: $(fixed_m6A_dirs)
.PHONY: fixed_m6A $(fixed_m6A_dirs)
#greer_pacbio_m6A_cov50: ../../greer_pacbio_no_chunk/ipd_summary.m6A_cov50.gff
#greer_pacbio_orig_m6A_cov50: ../../greer_pacbio_orig_no_chunk/ipd_summary.m6A_cov50.gff
$(fixed_m6A_dirs):
	$(MAKE) -C $@
