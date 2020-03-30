#!/usr/bin/env bash
OUTDIR=ame.high_ipd.cd.with_low_motifs.out
CONTROL=../../deep_kinetics_region.cd.bed.merged.sorted.slop20.c_elegans.bed.fa
INPUT=../../extreme_ipd.high.cd.coverage25.slop20.fullLength.c_elegans.gff.fa
MOTIF=motifs_low_ipd.meme
~/memesuite_5_1_1_install/bin/ame -oc $OUTDIR --control $CONTROL $INPUT $MOTIF
