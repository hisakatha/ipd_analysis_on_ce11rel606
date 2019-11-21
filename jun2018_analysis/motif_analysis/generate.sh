#!/usr/bin/env bash
REF=/glusterfs/hisakatha/ce11rel606/ce11rel606.fa
REF_PREFIX=$(basename $REF .fa)

#PATTERN=CTGDAR
PATTERN=$1

seqkit locate --bed -i -d -p $PATTERN $REF > $PATTERN.${REF_PREFIX}.bed

DIR=/glusterfs/hisakatha/methylation/smrtpipe/vs_ce11rel606
FILE=mapped.alignmentset.merged.bam.cov25.slop20.bed.fa

seqkit locate --bed -i -d -p $PATTERN $DIR/jun2018_ab_pcr_vc2010_op50/$FILE > $PATTERN.ab.bed
seqkit locate --bed -i -d -p $PATTERN $DIR/jun2018_cd_pcr_vc2010/$FILE > $PATTERN.cd.bed
seqkit locate --bed -i -d -p $PATTERN $DIR/jun2018_analysis/deep_region_cov25.k_normBy_ab.slop20.bed.fa > $PATTERN.k_normBy_ab.bed
seqkit locate --bed -i -d -p $PATTERN $DIR/jun2018_analysis/deep_region_cov25.l_normBy_cd.slop20.bed.fa > $PATTERN.l_normBy_cd.bed

